//! BayesTraits variable-rates output (`*.VarRates.txt`): per-branch rate scalars.
//!
//! File layout: the taxa (count, then index and name per line), the nodes
//! (count, then per line: node id, original branch length (-1 for the root),
//! number of taxa, taxon indices), then the MCMC log, whose header contains
//! repeating groups "Node ID, Scaler, Creat It, Node / Branch".
//!
//! In each iteration a *Node* scalar multiplies the branch leading to that node
//! and every branch below it; a *Branch* scalar multiplies that branch only.
//! Scalars that overlap multiply. Each branch's mean scalar is the average of
//! its combined scalar over iterations; "posterior sample" is the percentage
//! of iterations in which that combined scalar is not 1.

use crate::tree::{Attr, NodeId, Tree};
use anyhow::{bail, Context, Result};
use std::collections::{BTreeSet, HashMap};

#[derive(Clone)]
pub struct VrNode {
    pub id: usize,
    /// Original branch length (-1 for the root).
    pub length: f64,
    pub taxa: BTreeSet<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// Rate scalar on a node's clade, including its stem branch.
    Node,
    /// Rate scalar on one branch.
    Branch,
}

#[derive(Clone)]
pub struct Scalar {
    pub node: usize,
    pub value: f64,
    pub kind: Kind,
}

pub struct VarRates {
    pub taxa: Vec<String>,
    pub nodes: Vec<VrNode>,
    pub iterations: Vec<Vec<Scalar>>,
    /// Log-likelihood of each iteration ("Lh" column; NaN if absent).
    pub lh: Vec<f64>,
}

/// Summary of one branch (the branch leading to node `id`).
pub struct BranchRate {
    pub id: usize,
    pub n_taxa: usize,
    pub original_length: f64,
    pub is_root: bool,
    /// % of iterations in which the branch's combined scalar is not 1.
    pub pct_scaled: f64,
    pub mean: f64,
    pub median: f64,
    /// 95% HPD interval of the branch's scalar over all iterations (0 for the root).
    pub hpd: (f64, f64),
    /// Posterior probability (proportion of iterations) of a Node scalar placed at
    /// this node (rescaling its clade).
    pub node_pp: f64,
    /// Posterior probability of a Branch scalar placed on this branch.
    pub branch_pp: f64,
    /// Mean and median value of the node scalar placed at this node, over the
    /// iterations with one (1, no rescaling, if none).
    pub node_mean: f64,
    pub node_median: f64,
    /// 95% HPD interval of the node scalar, over the iterations with one (1 if none).
    pub node_hpd: (f64, f64),
    /// The matching node in the Canopy tree, if found.
    pub tree_node: Option<NodeId>,
}

pub fn parse(text: &str) -> Result<VarRates> {
    let mut lines = text.lines().map(str::trim_end).filter(|l| !l.trim().is_empty());
    let mut next = |what: &str| lines.next().with_context(|| format!("file ended while reading {}", what));
    let ntaxa: usize = next("the taxon count")?.trim().parse().context("first line should be the number of taxa")?;
    let mut taxa = vec![String::new(); ntaxa];
    for _ in 0..ntaxa {
        let l = next("taxa")?;
        let mut f = l.split('\t');
        let i: usize = f.next().unwrap_or("").trim().parse().context("taxon line should start with its index")?;
        let name = f.next().unwrap_or("").trim().to_string();
        if i >= ntaxa {
            bail!("taxon index {} out of range", i);
        }
        taxa[i] = name;
    }
    let nnodes: usize = next("the node count")?.trim().parse().context("expected the number of nodes after the taxa")?;
    let mut nodes = Vec::with_capacity(nnodes);
    for _ in 0..nnodes {
        let l = next("nodes")?;
        let f: Vec<&str> = l.split('\t').map(str::trim).filter(|s| !s.is_empty()).collect();
        if f.len() < 3 {
            bail!("malformed node line: {}", l);
        }
        let id: usize = f[0].parse().context("node id")?;
        let length: f64 = f[1].parse().context("node branch length")?;
        let count: usize = f[2].parse().context("node taxon count")?;
        let set: BTreeSet<usize> = f[3..].iter().map(|s| s.parse::<usize>()).collect::<std::result::Result<_, _>>().context("node taxa")?;
        if set.len() != count {
            bail!("node {} lists {} taxa but says {}", id, set.len(), count);
        }
        nodes.push(VrNode { id, length, taxa: set });
    }
    let header = next("the MCMC header")?;
    let cols: Vec<&str> = header.split('\t').map(str::trim).collect();
    let Some(first) = cols.iter().position(|c| c.eq_ignore_ascii_case("Node ID")) else { bail!("no \"Node ID\" column in the MCMC header") };
    let mut iterations = Vec::new();
    let lh_col = cols.iter().position(|c| c.eq_ignore_ascii_case("Lh"));
    let mut lh = Vec::new();
    for l in lines {
        let f: Vec<&str> = l.split('\t').map(str::trim).collect();
        lh.push(lh_col.and_then(|c| f.get(c)).and_then(|s| s.parse::<f64>().ok()).unwrap_or(f64::NAN));
        let mut row = Vec::new();
        let mut k = first;
        while k + 3 < f.len() + 1 && k < f.len() && !f[k].is_empty() {
            let node: usize = f[k].parse().with_context(|| format!("bad Node ID '{}'", f[k]))?;
            let value: f64 = f.get(k + 1).unwrap_or(&"").parse().with_context(|| format!("bad scaler for node {}", node))?;
            let kind = match f.get(k + 3).copied().unwrap_or("").to_ascii_lowercase().as_str() {
                "node" => Kind::Node,
                "branch" => Kind::Branch,
                other => bail!("expected Node or Branch, found '{}'", other),
            };
            row.push(Scalar { node, value, kind });
            k += 4;
        }
        iterations.push(row);
    }
    if iterations.is_empty() {
        bail!("no MCMC iterations found");
    }
    Ok(VarRates { taxa, nodes, iterations, lh })
}

impl VarRates {
    /// Several runs of the same analysis pooled into one posterior: their
    /// (post-burn-in) iterations one after another. The runs must have the same
    /// taxa and nodes, i.e. come from the same tree and data.
    pub fn combine(runs: &[VarRates]) -> Result<VarRates> {
        let Some(first) = runs.first() else { bail!("no runs to combine") };
        for (i, r) in runs.iter().enumerate().skip(1) {
            if r.taxa != first.taxa {
                bail!("run {} has different taxa from run 1; runs can only be combined when they analyse the same tree and data", i + 1);
            }
            let same = r.nodes.len() == first.nodes.len() && r.nodes.iter().zip(&first.nodes).all(|(a, b)| a.id == b.id && a.taxa == b.taxa);
            if !same {
                bail!("run {} has different tree nodes from run 1; runs can only be combined when they analyse the same tree", i + 1);
            }
        }
        Ok(VarRates {
            taxa: first.taxa.clone(),
            nodes: first.nodes.clone(),
            iterations: runs.iter().flat_map(|r| r.iterations.iter().cloned()).collect(),
            lh: runs.iter().flat_map(|r| r.lh.iter().copied()).collect(),
        })
    }

    /// Each branch's combined scalar (product of the scalars applying to it) in
    /// every iteration, in node order.
    pub fn branch_traces(&self) -> Vec<Vec<f64>> {
        let by_id: HashMap<usize, &VrNode> = self.nodes.iter().map(|n| (n.id, n)).collect();
        let mut per: Vec<Vec<f64>> = vec![Vec::with_capacity(self.iterations.len()); self.nodes.len()];
        for row in &self.iterations {
            for (j, b) in self.nodes.iter().enumerate() {
                let mut s = 1.0;
                for sc in row {
                    let Some(a) = by_id.get(&sc.node) else { continue };
                    let applies = match sc.kind {
                        Kind::Node => b.taxa.is_subset(&a.taxa),
                        Kind::Branch => sc.node == b.id,
                    };
                    if applies {
                        s *= sc.value;
                    }
                }
                per[j].push(s);
            }
        }
        per
    }

    /// Mean and median combined scalar per branch, and how often it was rescaled.
    pub fn summarize(&self) -> Vec<BranchRate> {
        let per = self.branch_traces();
        self.nodes
            .iter()
            .zip(per)
            .map(|(n, v)| {
                let is_root = n.length < 0.0;
                let k = v.len().max(1) as f64;
                let scaled = v.iter().filter(|&&x| (x - 1.0).abs() > 1e-12).count() as f64;
                let placed = |kind: Kind| self.iterations.iter().filter(|row| row.iter().any(|sc| sc.node == n.id && sc.kind == kind)).count() as f64;

                let median = median_of(&v).unwrap_or(1.0);
                let at_node = node_scalar_values(&self.iterations, n.id, |sc| (sc.node, sc.kind == Kind::Node, sc.value));
                BranchRate {
                    id: n.id,
                    n_taxa: n.taxa.len(),
                    original_length: if is_root { 0.0 } else { n.length },
                    is_root,
                    // The root has no branch: reported as 0, like BayesTraits' post-processor.
                    pct_scaled: if is_root { 0.0 } else { 100.0 * scaled / k },
                    mean: if is_root { 0.0 } else { v.iter().sum::<f64>() / k },
                    median: if is_root { 0.0 } else { median },
                    hpd: if is_root { (0.0, 0.0) } else { crate::consensus::hpd(&v, 0.95).unwrap_or((1.0, 1.0)) },
                    node_pp: placed(Kind::Node) / k,
                    branch_pp: if is_root { 0.0 } else { placed(Kind::Branch) / k },
                    node_mean: mean_or_one(&at_node),
                    node_median: median_of(&at_node).unwrap_or(1.0),
                    node_hpd: crate::consensus::hpd(&at_node, 0.95).unwrap_or((1.0, 1.0)),
                    tree_node: None,
                }
            })
            .collect()
    }

    /// Find each VarRates node in `tree` by its species, and write the summaries
    /// onto the tree as `rate_scalar`, `rate_scalar_median` and `rate_scaled_pct`.
    /// Returns warnings (unmatched nodes, branch-length mismatches).
    pub fn apply(&self, rates: &mut [BranchRate], tree: &mut Tree) -> Vec<String> {
        let mut warnings = Vec::new();
        let tips: HashMap<String, NodeId> = tree.tips().into_iter().map(|t| (norm(tree.label(t)), t)).collect();
        let by_id: HashMap<usize, &VrNode> = self.nodes.iter().map(|n| (n.id, n)).collect();
        for r in rates.iter_mut() {
            let n = by_id[&r.id];
            let ids: Option<Vec<NodeId>> = n.taxa.iter().map(|&i| self.taxa.get(i).and_then(|name| tips.get(&norm(name)).copied())).collect();
            let Some(ids) = ids else {
                warnings.push(format!("node {}: some of its species are not in the tree", r.id));
                continue;
            };
            let Some(m) = tree.mrca(&ids) else { continue };
            if tree.tips_below(m).len() != ids.len() {
                warnings.push(format!("node {}: its species do not form a clade in this tree", r.id));
                continue;
            }
            r.tree_node = Some(m);
            if r.is_root {
                continue;
            }
            if let Some(l) = tree.nodes[m].length {
                if (l - r.original_length).abs() > 1e-4 * l.abs().max(1.0) {
                    warnings.push(format!("node {}: branch length {} in the tree vs {} in the file", r.id, crate::tree::format_num(l, 6), crate::tree::format_num(r.original_length, 6)));
                }
            }
            let a = &mut tree.nodes[m].attrs;
            a.insert("rate_scalar".into(), Attr::Num(r.mean));
            a.insert("rate_scalar_median".into(), Attr::Num(r.median));
            a.insert("rate_scaled_pct".into(), Attr::Num(r.pct_scaled));
            a.insert("rate_node_pp".into(), Attr::Num(r.node_pp));
            a.insert("rate_branch_pp".into(), Attr::Num(r.branch_pp));
            a.insert("rate_node_mean".into(), Attr::Num(r.node_mean));
            a.insert("rate_node_median".into(), Attr::Num(r.node_median));
            a.insert("rate_scalar_hpd".into(), Attr::Range(r.hpd.0, r.hpd.1));
            a.insert("rate_node_hpd".into(), Attr::Range(r.node_hpd.0, r.node_hpd.1));
        }
        warnings
    }
}

/// Value of the node scalar(s) placed at node `id` in each iteration that has
/// one (several in one iteration multiply). `f` gives (node, is a node scalar, value).
pub(crate) fn node_scalar_values<E>(iterations: &[Vec<E>], id: usize, f: impl Fn(&E) -> (usize, bool, f64)) -> Vec<f64> {
    iterations
        .iter()
        .filter_map(|row| {
            let v: Vec<f64> = row.iter().map(&f).filter(|&(n, node, _)| node && n == id).map(|x| x.2).collect();
            (!v.is_empty()).then(|| v.iter().product())
        })
        .collect()
}

/// Mean of `v`, or 1 (no rescaling) if it is empty.
pub(crate) fn mean_or_one(v: &[f64]) -> f64 {
    if v.is_empty() {
        1.0
    } else {
        v.iter().sum::<f64>() / v.len() as f64
    }
}

fn median_of(v: &[f64]) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let m = s.len() / 2;
    Some(if s.len() % 2 == 1 { s[m] } else { (s[m - 1] + s[m]) / 2.0 })
}

fn norm(s: &str) -> String {
    s.trim().replace(' ', "_").to_ascii_lowercase()
}


/// MCMC diagnostics of one or more runs (named `names`): the log-likelihood and
/// each branch's rate scalar per run (the root, never rescaled, left out), plus the
/// comparison of posterior probabilities across runs when there are several.
pub fn diagnostics(runs: &[VarRates], names: &[String]) -> crate::mcmc::Diagnostics {
    use crate::mcmc::Param;
    let traces: Vec<Vec<Vec<f64>>> = runs.iter().map(|r| r.branch_traces()).collect();
    let mut params = vec![Param { name: "Log-likelihood".into(), chains: runs.iter().map(|r| r.lh.clone()).collect() }];
    for (j, n) in runs[0].nodes.iter().enumerate() {
        if n.length < 0.0 {
            continue;
        }
        let what = match n.taxa.iter().next() {
            Some(&t) if n.taxa.len() == 1 => runs[0].taxa.get(t).map(|s| s.replace('_', " ")).unwrap_or_default(),
            _ => format!("{} taxa", n.taxa.len()),
        };
        params.push(Param { name: format!("Rate scalar: node {} ({})", n.id, what), chains: traces.iter().map(|t| t[j].clone()).collect() });
    }
    crate::mcmc::Diagnostics {
        title: String::new(),
        runs: runs.iter().zip(names).map(|(r, nm)| (nm.clone(), r.iterations.len())).collect(),
        params,
        report: (runs.len() > 1).then(|| convergence(runs, names)),
    }
}
/// Node-scalar posterior of node `n`: the proportion (0–1) of the sample with a
/// node scalar placed there, if a VarRates file has been applied.
/// Convergence of several variable-rates runs (named `names`): the
/// log-likelihood and each branch's rate scalar compared across runs (R-hat, ESS),
/// and the posterior probabilities of branch and node scalars in each run.
pub fn convergence(runs: &[VarRates], names: &[String]) -> crate::mcmc::Report {
    let summaries: Vec<Vec<BranchRate>> = runs.iter().map(|r| r.summarize()).collect();
    let traces: Vec<Vec<Vec<f64>>> = runs.iter().map(|r| r.branch_traces()).collect();
    let lh: Vec<Vec<f64>> = runs.iter().map(|r| r.lh.clone()).collect();
    let branches = runs[0]
        .nodes
        .iter()
        .enumerate()
        .map(|(j, n)| {
            let chains: Vec<Vec<f64>> = traces.iter().map(|t| t[j].clone()).collect();
            crate::mcmc::BranchCheck {
                id: n.id,
                n_taxa: n.taxa.len(),
                pp: summaries.iter().map(|s| s[j].branch_pp).collect(),
                node_pp: summaries.iter().map(|s| s[j].node_pp).collect(),
                rhat: crate::mcmc::rhat(&chains),
                ess: crate::mcmc::ess(&chains),
            }
        })
        .collect();
    crate::mcmc::Report {
        runs: runs.iter().zip(names).map(|(r, nm)| (nm.clone(), r.iterations.len(), r.lh.iter().sum::<f64>() / r.lh.len().max(1) as f64)).collect(),
        lh_rhat: crate::mcmc::rhat(&lh),
        lh_ess: crate::mcmc::ess(&lh),
        branches,
        quantity: "rate scalar",
    }
}

pub fn node_posterior(tree: &Tree, n: NodeId) -> Option<f64> {
    let a = &tree.nodes[n].attrs;
    a.get("rate_node_pp").and_then(Attr::as_f64).or_else(|| a.get("fabric_node_prop").and_then(Attr::as_f64))
}

/// The branch-posterior attribute for a BayesTraits branch value: the share of
/// the sample with a branch scalar on the branch (variable rates) or an LS_Beta
/// (Fabric). None for other attributes.
pub fn branch_pp_key(attr: &str) -> Option<&'static str> {
    if attr.starts_with("fabric_beta") {
        Some("fabric_beta_prop")
    } else if attr.starts_with("rate_scal") {
        Some("rate_branch_pp")
    } else {
        None
    }
}

/// Whether branch `n` shows `attr` at posterior `min_pp`: always for other
/// attributes or a zero threshold, else when its branch posterior is at least `min_pp`.
pub fn branch_passes(tree: &Tree, n: NodeId, attr: &str, min_pp: f64) -> bool {
    match branch_pp_key(attr) {
        Some(k) if min_pp > 0.0 => tree.nodes[n].attrs.get(k).and_then(Attr::as_f64).map(|p| p >= min_pp - 1e-12).unwrap_or(false),
        _ => true,
    }
}

/// Copy of `tree` with each branch length multiplied by its `rate_scalar`.
pub fn scaled_tree(tree: &Tree) -> Tree {
    let mut t = tree.compact();
    for n in t.preorder() {
        let s = t.nodes[n].attrs.get("rate_scalar").and_then(Attr::as_f64);
        if let (Some(s), Some(l)) = (s, t.nodes[n].length.as_mut()) {
            *l *= s;
        }
    }
    t.name = Some("mean_rate_scaled".into());
    t
}

/// Column names and values per branch: Node ID, No Taxa, Original BL, Mean
/// Scalar, Median Scalar, Scalar 95% HPD, Branch Scalar PP, Mean / Median Node
/// Scalar, Node Scalar 95% HPD, Node Scalar PP.
pub fn columns(rates: &[BranchRate]) -> (Vec<String>, Vec<Vec<f64>>, Vec<usize>) {
    let header = [
        "Node ID",
        "No Taxa",
        "Original BL",
        "Mean Scalar",
        "Median Scalar",
        "Scalar 95% HPD Lower",
        "Scalar 95% HPD Upper",
        "Branch Scalar PP",
        "Mean Node Scalar",
        "Median Node Scalar",
        "Node Scalar 95% HPD Lower",
        "Node Scalar 95% HPD Upper",
        "Node Scalar PP",
    ];
    let rows = rates
        .iter()
        .map(|r| vec![r.id as f64, r.n_taxa as f64, r.original_length, r.mean, r.median, r.hpd.0, r.hpd.1, r.branch_pp, r.node_mean, r.node_median, r.node_hpd.0, r.node_hpd.1, r.node_pp])
        .collect();
    // Decimals per column: 6, with 9 for the posterior probabilities (as in the Fabric table).
    (header.iter().map(|s| s.to_string()).collect(), rows, vec![0, 0, 6, 6, 6, 6, 6, 9, 6, 6, 6, 6, 9])
}

/// Table in the layout of BayesTraits' post-processor (6 decimals; 9 for the PP columns).
pub fn table(rates: &[BranchRate]) -> (Vec<String>, Vec<Vec<String>>) {
    let (header, rows, digits) = columns(rates);
    let rows = rows.iter().map(|r| r.iter().zip(&digits).map(|(&x, &d)| crate::tree::format_num(x, d)).collect()).collect();
    (header, rows)
}
#[cfg(test)]
mod tests {
    use super::*;

    const FILE: &str = "9\n0\tHomo_sapiens\n1\tPan_troglodytes\n2\tGorilla_gorilla\n3\tPongo_abelii\n4\tHylobates_lar\n5\tMacaca_mulatta\n6\tPapio_anubis\n7\tColobus_guereza\n8\tNasalis_larvatus\n17\n0\t-1.000000\t9\t0\t1\t2\t3\t4\t5\t6\t7\t8\t\n1\t10.195000\t5\t0\t1\t2\t3\t4\t\n2\t3.851600\t4\t0\t1\t2\t3\t\n3\t7.325500\t3\t0\t1\t2\t\n4\t1.993571\t2\t0\t1\t\n5\t6.591829\t1\t0\t\n6\t6.591829\t1\t1\t\n7\t8.585400\t1\t2\t\n8\t15.910500\t1\t3\t\n9\t19.762500\t1\t4\t\n10\t10.898300\t4\t5\t6\t7\t8\t\n11\t8.786356\t2\t5\t6\t\n12\t10.273344\t1\t5\t\n13\t10.273344\t1\t6\t\n14\t7.599467\t2\t7\t8\t\n15\t11.460033\t1\t7\t\n16\t11.460033\t1\t8\t\nIt\tLh\tLh + Prior\tNo Pram\tAlpha\tSigma^2\tAlpha Scale Prior\tNode ID\tScaler\tCreat It\tNode / Branch\t\n10000\t-43.720198\t-58.088139\t2\t2.437032\t42.512694\t1.200000\t2\t4.239166\t9874\tNode\t2\t1.415276\t9902\tBranch\t\n11000\t-42.969946\t-55.619873\t1\t30.747735\t13.708323\t1.200000\t3\t6.351064\t10895\tNode\t\n12000\t-44.512828\t-63.049036\t2\t-15.427550\t20.768136\t1.200000\t2\t26.807782\t11804\tNode\t5\t1.737369\t11992\tBranch\t\n";

    #[test]
    fn runs_combine_and_compare() {
        let a = parse(FILE).unwrap();
        let b = parse(FILE).unwrap();
        assert_eq!(a.lh.len(), 3);
        assert!((a.lh[0] - FILE.lines().nth(29).unwrap().split('\t').nth(1).unwrap().parse::<f64>().unwrap()).abs() < 1e-12);
        let both = VarRates::combine(&[a, b]).unwrap();
        assert_eq!((both.iterations.len(), both.lh.len()), (6, 6));
        // Identical runs: the same posterior probabilities in each, so no spread.
        let runs = [parse(FILE).unwrap(), parse(FILE).unwrap()];
        let r = convergence(&runs, &["a".to_string(), "b".to_string()]);
        assert_eq!(r.runs.len(), 2);
        assert_eq!(r.asd_pp(), (0.0, 0.0));
        assert_eq!(r.table().1.len(), runs[0].nodes.len());
        // Runs from a different tree can't be pooled.
        let other = FILE.replacen("3\t7.325500\t3\t0\t1\t2", "3\t7.325500\t3\t0\t1\t3", 1);
        assert!(VarRates::combine(&[parse(FILE).unwrap(), parse(&other).unwrap()]).is_err());
    }

    #[test]
    fn matches_bayestraits_postprocessor() {
        let vr = parse(FILE).unwrap();
        assert_eq!(vr.iterations.len(), 3);
        let s = vr.summarize();
        // (node id, posterior sample %, mean scalar) from the BayesTraits target.
        let target = [
            (0, 0.0, 0.0),
            (1, 0.0, 1.0),
            (2, 66.666667, 11.269124),
            (3, 100.0, 12.466004),
            (4, 100.0, 12.466004),
            (5, 100.0, 19.05508),
            (6, 100.0, 12.466004),
            (7, 100.0, 12.466004),
            (8, 66.666667, 10.682316),
            (9, 0.0, 1.0),
            (10, 0.0, 1.0),
            (16, 0.0, 1.0),
        ];
        for (id, pct, mean) in target {
            let r = s.iter().find(|r| r.id == id).unwrap();
            assert!((r.pct_scaled - pct).abs() < 1e-5, "node {} pct {}", id, r.pct_scaled);
            assert!((r.mean - mean).abs() < 1e-5, "node {} mean {}", id, r.mean);
        }
        // Where scalars were placed: node 2 (Node) in 2 of 3 iterations, node 3 in 1;
        // Branch scalars on node 2's and node 5's branches once each.
        let pct = |id: usize| s.iter().find(|r| r.id == id).map(|r| (r.node_pp, r.branch_pp)).unwrap();
        let close = |a: f64, b: f64| (a - b).abs() < 1e-5;
        assert!(close(pct(2).0, 2.0 / 3.0) && close(pct(2).1, 1.0 / 3.0));
        assert!(close(pct(3).0, 1.0 / 3.0) && close(pct(3).1, 0.0));
        assert!(close(pct(5).0, 0.0) && close(pct(5).1, 1.0 / 3.0));
    }

    #[test]
    fn applies_to_tree_by_species() {
        let mut t = crate::io::newick::parse_newick(
            "(((((Homo_sapiens:6.5918285714,Pan_troglodytes:6.5918285714):1.9935714286,Gorilla_gorilla:8.5854):7.3255,Pongo_abelii:15.9105):3.8516,Hylobates_lar:19.7625):10.195,((Macaca_mulatta:10.2733444444,Papio_anubis:10.2733444444):8.7863555556,(Colobus_guereza:11.4600333333,Nasalis_larvatus:11.4600333333):7.5994666667):10.8983);",
        )
        .unwrap();
        let vr = parse(FILE).unwrap();
        let mut s = vr.summarize();
        let warnings = vr.apply(&mut s, &mut t);
        assert!(warnings.is_empty(), "{:?}", warnings);
        let homo = t.find_label("Homo_sapiens").unwrap();
        assert!((t.nodes[homo].attrs["rate_scalar"].as_f64().unwrap() - 19.05508).abs() < 1e-5);
        let scaled = scaled_tree(&t);
        let h = scaled.find_label("Homo_sapiens").unwrap();
        assert!((scaled.nodes[h].length.unwrap() - 6.5918285714 * 19.05508).abs() < 1e-3);
    }

    #[test]
    fn labels_show_scalars_and_node_percentages() {
        use crate::scene::{build, ApproxEnv, BuildInput, Prim};
        use crate::style::*;
        let mut t = crate::io::newick::parse_newick(
            "(((((Homo_sapiens:6.59,Pan_troglodytes:6.59):1.99,Gorilla_gorilla:8.59):7.33,Pongo_abelii:15.91):3.85,Hylobates_lar:19.76):10.2,((Macaca_mulatta:10.27,Papio_anubis:10.27):8.79,(Colobus_guereza:11.46,Nasalis_larvatus:11.46):7.6):10.9);",
        )
        .unwrap();
        let vr = parse(FILE).unwrap();
        let mut s = vr.summarize();
        vr.apply(&mut s, &mut t);
        let layers = vec![
            LayerEntry::new(Layer::BranchLengths(BranchLengthStyle { size: 7.0, color: Color::BLACK, digits: 2, attr: Some("rate_scalar".into()), hide_value: Some(1.0), min_pp: 0.0, min_value: None, min_abs: false })),
            LayerEntry::new(Layer::default_node_scalars()),
        ];
        let view = ViewState::default();
        let sc = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        let texts: Vec<String> = sc.prims.iter().filter_map(|p| if let Prim::Text { text, .. } = p { Some(text.clone()) } else { None }).collect();
        // Scaled branches only: 7 rescaled branches (nodes 2-8), none of the unscaled (= 1) ones.
        assert!(texts.contains(&"19.06".to_string()) && texts.contains(&"11.27".to_string()));
        assert!(!texts.contains(&"1".to_string()));
        // Node scalars: circles at node 2 (posterior 0.667) and node 3 (0.333), in
        // different colors, with a legend titled for them.
        let circles = |sc: &crate::scene::Scene| -> Vec<Color> { sc.prims.iter().filter_map(|p| if let Prim::Circle { fill, .. } = p { Some(*fill) } else { None }).collect() };
        let all = circles(&sc);
        assert_eq!(all.len(), 2, "{:?}", all);
        assert_ne!(all[0], all[1]);
        assert!(texts.iter().any(|t| t == "Node scalar posterior"));
        // A minimum posterior filters them: 0.5 keeps node 2 only, 0.2 keeps both.
        let with_threshold = |th: f64| {
            let mut ls = layers.clone();
            if let Layer::NodeScalars(s) = &mut ls[1].layer {
                s.threshold = th;
            }
            let sc = build(&BuildInput { tree: &t, view: &view, layers: &ls, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
            circles(&sc).len()
        };
        assert_eq!(with_threshold(0.5), 1);
        assert_eq!(with_threshold(0.2), 2);
        assert_eq!(with_threshold(0.7), 0);
        // One color, no color key; sized by posterior, the 0.667 node gets the larger circle.
        let mut ls = layers.clone();
        if let Layer::NodeScalars(s) = &mut ls[1].layer {
            s.single_color = Some(Color::rgb(200, 0, 0));
            s.size_by_pp = true;
        }
        let sc1 = build(&BuildInput { tree: &t, view: &view, layers: &ls, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        let dots: Vec<(Color, f32)> = sc1.prims.iter().filter_map(|p| if let Prim::Circle { fill, r, .. } = p { Some((*fill, *r)) } else { None }).collect();
        assert_eq!(dots.len(), 2);
        assert!(dots.iter().all(|d| d.0 == Color::rgb(200, 0, 0)));
        assert!((dots[0].1 - dots[1].1).abs() > 0.5);
        assert!(!sc1.prims.iter().any(|p| matches!(p, Prim::Text { text, .. } if text == "Node scalar posterior")));
        // Branch labels at a minimum branch posterior: only branches with a branch
        // scalar that often (and a rate other than 1) keep their label.
        let labels_at = |th: f64| {
            let mut ls = layers.clone();
            if let Layer::BranchLengths(s) = &mut ls[0].layer {
                s.min_pp = th;
            }
            ls.truncate(1);
            let sc = build(&BuildInput { tree: &t, view: &view, layers: &ls, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
            sc.prims.iter().filter(|p| matches!(p, Prim::Text { .. })).count()
        };
        let expected = |th: f64| {
            t.preorder()
                .into_iter()
                .filter(|&n| t.nodes[n].parent.is_some())
                .filter(|&n| t.nodes[n].attrs.get("rate_scalar").and_then(Attr::as_f64).map(|v| (v - 1.0).abs() > 1e-9).unwrap_or(false))
                .filter(|&n| t.nodes[n].attrs.get("rate_branch_pp").and_then(Attr::as_f64).unwrap_or(0.0) >= th - 1e-12)
                .count()
        };
        for th in [0.0, 0.2, 0.5, 0.9] {
            assert_eq!(labels_at(th), expected(th), "threshold {}", th);
        }
        assert!(expected(0.5) < expected(0.0));
        assert_eq!((branch_pp_key("fabric_beta_median"), branch_pp_key("rate_scalar"), branch_pp_key("rate_node_pp")), (Some("fabric_beta_prop"), Some("rate_branch_pp"), None));
        // A minimum value keeps only labels at least that large.
        let mut ls = layers.clone();
        if let Layer::BranchLengths(s) = &mut ls[0].layer {
            s.min_value = Some(15.0);
        }
        ls.truncate(1);
        let sc = build(&BuildInput { tree: &t, view: &view, layers: &ls, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        let big = t.preorder().into_iter().filter(|&n| t.nodes[n].parent.is_some() && t.nodes[n].attrs.get("rate_scalar").and_then(Attr::as_f64).map(|v| v >= 15.0).unwrap_or(false)).count();
        assert!(big > 0 && big < expected(0.0));
        assert_eq!(sc.prims.iter().filter(|p| matches!(p, Prim::Text { .. })).count(), big);
        // Branch colors: below the threshold a branch keeps the line color.
        let mut ts = TreeStyle { width: 1.0, color: Color::BLACK, color_by: Some("rate_scalar".into()), palette: Palette::Viridis, reverse_palette: false, min_pp: 0.0, min_value: None, min_abs: false, width_by: None, width_min: 0.5 };
        let black = |ts: &TreeStyle| {
            let ls = vec![LayerEntry::new(Layer::Tree(ts.clone()))];
            let sc = build(&BuildInput { tree: &t, view: &view, layers: &ls, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
            sc.prims.iter().filter(|p| matches!(p, Prim::Path { color, .. } if *color == Color::BLACK)).count()
        };
        let none = black(&ts);
        ts.min_pp = 0.5;
        assert!(black(&ts) > none);
    }

    #[test]
    fn table_layout() {
        let (header, rows) = table(&parse(FILE).unwrap().summarize());
        assert_eq!(
            header,
            [
                "Node ID",
                "No Taxa",
                "Original BL",
                "Mean Scalar",
                "Median Scalar",
                "Scalar 95% HPD Lower",
                "Scalar 95% HPD Upper",
                "Branch Scalar PP",
                "Mean Node Scalar",
                "Median Node Scalar",
                "Node Scalar 95% HPD Lower",
                "Node Scalar 95% HPD Upper",
                "Node Scalar PP"
            ]
        );
        // Node 2's scalar is 5.99959, 1 and 26.807782 in the three iterations (95% of
        // three samples is all three); its node scalars are 4.239166 and 26.807782.
        assert_eq!(rows[2].join("\t"), "2\t4\t3.8516\t11.269124\t5.99959\t1\t26.807782\t0.333333333\t15.523474\t15.523474\t4.239166\t26.807782\t0.666666667");
        // Node 3: one node scalar (6.351064); node 5: none, so 1.
        assert!(rows[3].join("\t").ends_with("\t0\t6.351064\t6.351064\t6.351064\t6.351064\t0.333333333"), "{:?}", rows[3]);
        assert!(rows[5].join("\t").ends_with("\t0.333333333\t1\t1\t1\t1\t0"), "{:?}", rows[5]);
    }
}