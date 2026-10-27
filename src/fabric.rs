//! BayesTraits Fabric-model output (`*.VarRates.txt` from a Fabric run).
//!
//! File layout: the taxa (count, then index and name per line), the nodes
//! (count, then per line: node id, original branch length (-1 for the root),
//! number of taxa, taxon indices), then the MCMC log, whose header contains
//! repeating groups "Node ID, Scaler, Creat It, Node / Branch". The last field
//! of each group is "Node" (a node scalar) or "LS_Beta" (a value on a single
//! branch).
//!
//! Per branch: the mean and median LS_Beta over the iterations in which that
//! branch carried one (0 if it never did), the proportion of iterations with an
//! LS_Beta on it, and the proportion of iterations with a node scalar placed at
//! that node.

use crate::tree::{Attr, NodeId, Tree};
use anyhow::{bail, Context, Result};
use std::collections::{BTreeSet, HashMap};

#[derive(Clone)]
pub struct FabNode {
    pub id: usize,
    /// Original branch length (-1 for the root).
    pub length: f64,
    pub taxa: BTreeSet<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Node,
    LsBeta,
}

#[derive(Clone)]
pub struct Entry {
    pub node: usize,
    pub value: f64,
    pub kind: Kind,
}

pub struct Fabric {
    pub taxa: Vec<String>,
    pub nodes: Vec<FabNode>,
    pub iterations: Vec<Vec<Entry>>,
    /// Log-likelihood of each iteration ("Lh" column; NaN if absent).
    pub lh: Vec<f64>,
}

/// Summary of one branch (the branch leading to node `id`).
pub struct FabricBranch {
    pub id: usize,
    pub n_taxa: usize,
    /// The species name when the node is a single tip (empty otherwise).
    pub tip: String,
    pub original_length: f64,
    pub is_root: bool,
    /// Mean and median LS_Beta over the iterations with one (0 if none).
    pub beta_mean: f64,
    pub beta_median: f64,
    /// 95% HPD interval of LS_Beta over the iterations with one (0 if none).
    pub beta_hpd: (f64, f64),
    /// Proportion of iterations with an LS_Beta on this branch.
    pub beta_prop: f64,
    /// Proportion of iterations with a node scalar placed at this node.
    pub node_prop: f64,
    /// Mean and median value of the node scalar placed at this node, over the
    /// iterations with one (1, no rescaling, if none).
    pub node_mean: f64,
    pub node_median: f64,
    /// 95% HPD interval of the node scalar over the iterations with one (1 if none).
    pub node_hpd: (f64, f64),
    /// Mean product of the node scalars applying to this branch (those placed at
    /// this node or an ancestor), 1 in iterations with none; 1 for the root.
    pub node_scalar: f64,
    /// The matching node in the Canopy tree, if found.
    pub tree_node: Option<NodeId>,
}

pub fn parse(text: &str) -> Result<Fabric> {
    let mut lines = text.lines().map(str::trim_end).filter(|l| !l.trim().is_empty());
    let mut next = |what: &str| lines.next().with_context(|| format!("file ended while reading {}", what));
    let ntaxa: usize = next("the taxon count")?.trim().parse().context("first line should be the number of taxa")?;
    let mut taxa = vec![String::new(); ntaxa];
    for _ in 0..ntaxa {
        let l = next("taxa")?;
        let mut f = l.split('\t');
        let i: usize = f.next().unwrap_or("").trim().parse().context("taxon line should start with its index")?;
        if i >= ntaxa {
            bail!("taxon index {} out of range", i);
        }
        taxa[i] = f.next().unwrap_or("").trim().to_string();
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
        nodes.push(FabNode { id, length, taxa: set });
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
        while k < f.len() && !f[k].is_empty() {
            let node: usize = f[k].parse().with_context(|| format!("bad Node ID '{}'", f[k]))?;
            let value: f64 = f.get(k + 1).unwrap_or(&"").parse().with_context(|| format!("bad scaler for node {}", node))?;
            let kind = match f.get(k + 3).copied().unwrap_or("").to_ascii_lowercase().as_str() {
                "node" => Kind::Node,
                "ls_beta" => Kind::LsBeta,
                other => bail!("expected Node or LS_Beta, found '{}' (is this a Fabric run?)", other),
            };
            row.push(Entry { node, value, kind });
            k += 4;
        }
        iterations.push(row);
    }
    if iterations.is_empty() {
        bail!("no MCMC iterations found");
    }
    Ok(Fabric { taxa, nodes, iterations, lh })
}

fn median(v: &[f64]) -> Option<f64> {
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

impl Fabric {
    /// Several runs of the same analysis pooled into one posterior: their
    /// (post-burn-in) iterations one after another. The runs must have the same
    /// taxa and nodes, i.e. come from the same tree and data.
    pub fn combine(runs: &[Fabric]) -> Result<Fabric> {
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
        Ok(Fabric {
            taxa: first.taxa.clone(),
            nodes: first.nodes.clone(),
            iterations: runs.iter().flat_map(|r| r.iterations.iter().cloned()).collect(),
            lh: runs.iter().flat_map(|r| r.lh.iter().copied()).collect(),
        })
    }

    /// Each branch's LS_Beta in every iteration (0 when it has none), in node order.
    pub fn beta_traces(&self) -> Vec<Vec<f64>> {
        self.nodes
            .iter()
            .map(|n| self.iterations.iter().map(|row| row.iter().filter(|e| e.node == n.id && e.kind == Kind::LsBeta).map(|e| e.value).sum()).collect())
            .collect()
    }

    pub fn summarize(&self) -> Vec<FabricBranch> {
        let k = self.iterations.len().max(1) as f64;
        let by_id: HashMap<usize, &FabNode> = self.nodes.iter().map(|n| (n.id, n)).collect();
        // Each branch on its own, over the cores.
        crate::mcmc::par_map(&self.nodes, |n| {
            let is_root = n.length < 0.0;
            // LS_Beta on this branch, one value per iteration that has one.
            let betas: Vec<f64> = self
                .iterations
                .iter()
                .filter_map(|row| {
                    let v: Vec<f64> = row.iter().filter(|e| e.node == n.id && e.kind == Kind::LsBeta).map(|e| e.value).collect();
                    (!v.is_empty()).then(|| v.iter().sum())
                })
                .collect();
            let node_hits = self.iterations.iter().filter(|row| row.iter().any(|e| e.node == n.id && e.kind == Kind::Node)).count();
            let at_node = crate::varrates::node_scalar_values(&self.iterations, n.id, |e| (e.node, e.kind == Kind::Node, e.value));
            // Node scalars rescale the clade, stem branch included; overlapping ones multiply.
            let products: f64 = self
                .iterations
                .iter()
                .map(|row| row.iter().filter(|e| e.kind == Kind::Node && by_id.get(&e.node).map(|a| n.taxa.is_subset(&a.taxa)).unwrap_or(false)).map(|e| e.value).product::<f64>())
                .sum();
            FabricBranch {
                id: n.id,
                n_taxa: n.taxa.len(),
                tip: crate::varrates::tip_name(&self.taxa, &n.taxa),
                original_length: if is_root { 0.0 } else { n.length },
                is_root,
                beta_mean: if betas.is_empty() { 0.0 } else { betas.iter().sum::<f64>() / betas.len() as f64 },
                beta_median: median(&betas).unwrap_or(0.0),
                beta_hpd: crate::consensus::hpd(&betas, 0.95).unwrap_or((0.0, 0.0)),
                beta_prop: betas.len() as f64 / k,
                node_prop: node_hits as f64 / k,
                node_mean: crate::varrates::mean_or_one(&at_node),
                node_median: median(&at_node).unwrap_or(1.0),
                node_hpd: crate::consensus::hpd(&at_node, 0.95).unwrap_or((1.0, 1.0)),
                node_scalar: if is_root { 1.0 } else { products / k },
                tree_node: None,
            }
        })
    }

    /// Find each node in `tree` by its species and store the summaries on it as
    /// `fabric_beta_mean`, `fabric_beta_median`, `fabric_beta_prop` and
    /// `fabric_node_prop`. Returns warnings (unmatched nodes, length mismatches).
    pub fn apply(&self, rows: &mut [FabricBranch], tree: &mut Tree) -> Vec<String> {
        let mut warnings = Vec::new();
        let tips: HashMap<String, NodeId> = tree.tips().into_iter().map(|t| (norm(tree.label(t)), t)).collect();
        let by_id: HashMap<usize, &FabNode> = self.nodes.iter().map(|n| (n.id, n)).collect();
        for r in rows.iter_mut() {
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
            a.insert("fabric_beta_mean".into(), Attr::Num(r.beta_mean));
            a.insert("fabric_beta_median".into(), Attr::Num(r.beta_median));
            a.insert("fabric_beta_prop".into(), Attr::Num(r.beta_prop));
            a.insert("fabric_node_prop".into(), Attr::Num(r.node_prop));
            a.insert("fabric_node_mean".into(), Attr::Num(r.node_mean));
            a.insert("fabric_node_median".into(), Attr::Num(r.node_median));
            a.insert("fabric_beta_hpd".into(), Attr::Range(r.beta_hpd.0, r.beta_hpd.1));
            a.insert("fabric_node_hpd".into(), Attr::Range(r.node_hpd.0, r.node_hpd.1));
            a.insert("fabric_node_scalar".into(), Attr::Num(r.node_scalar));
        }
        warnings
    }
}


/// MCMC diagnostics of one or more runs (named `names`): the log-likelihood and
/// each branch's β per run (the root, never rescaled, left out), plus the
/// comparison of posterior probabilities across runs when there are several.
pub fn diagnostics(runs: &[Fabric], names: &[String]) -> crate::mcmc::Diagnostics {
    use crate::mcmc::Param;
    let traces: Vec<Vec<Vec<f64>>> = runs.iter().map(|r| r.beta_traces()).collect();
    let mut params = vec![Param { name: "Log-likelihood".into(), chains: runs.iter().map(|r| r.lh.clone()).collect() }];
    for (j, n) in runs[0].nodes.iter().enumerate() {
        if n.length < 0.0 {
            continue;
        }
        let what = match n.taxa.iter().next() {
            Some(&t) if n.taxa.len() == 1 => runs[0].taxa.get(t).map(|s| s.replace('_', " ")).unwrap_or_default(),
            _ => format!("{} taxa", n.taxa.len()),
        };
        params.push(Param { name: format!("β: node {} ({})", n.id, what), chains: traces.iter().map(|t| t[j].clone()).collect() });
    }
    crate::mcmc::Diagnostics {
        title: String::new(),
        runs: runs.iter().zip(names).map(|(r, nm)| (nm.clone(), r.iterations.len())).collect(),
        params,
        report: (runs.len() > 1).then(|| convergence(runs, names)),
        cache: Default::default(),
    }
}
/// Copy of `tree` with each branch length multiplied by its mean node scalar
/// (`fabric_node_scalar`); LS_Beta shifts the trait, not the rate, so it is left out.
/// Convergence of several Fabric runs (named `names`): the log-likelihood and
/// each branch's LS_Beta compared across runs (R-hat, ESS), and the posterior
/// probabilities of an LS_Beta and of a node scalar in each run.
pub fn convergence(runs: &[Fabric], names: &[String]) -> crate::mcmc::Report {
    let summaries: Vec<Vec<FabricBranch>> = runs.iter().map(|r| r.summarize()).collect();
    let traces: Vec<Vec<Vec<f64>>> = runs.iter().map(|r| r.beta_traces()).collect();
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
                tip: crate::varrates::tip_name(&runs[0].taxa, &n.taxa),
                pp: summaries.iter().map(|s| s[j].beta_prop).collect(),
                node_pp: summaries.iter().map(|s| s[j].node_prop).collect(),
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
        quantity: "LS_Beta",
    }
}

pub fn scaled_tree(tree: &Tree) -> Tree {
    let mut t = tree.compact();
    for n in t.preorder() {
        let s = t.nodes[n].attrs.get("fabric_node_scalar").and_then(Attr::as_f64);
        if let (Some(s), Some(l)) = (s, t.nodes[n].length.as_mut()) {
            *l *= s;
        }
    }
    t.name = Some("mean_node_scalar_scaled".into());
    t
}

/// Column names, values per branch and decimals per column.
pub fn columns(rows: &[FabricBranch]) -> (Vec<String>, Vec<Vec<f64>>, Vec<usize>) {
    let header = [
        "Node ID",
        "No Taxa",
        "Original BL",
        "Mean (β * BL)",
        "Median (β * BL)",
        "β 95% HPD Lower",
        "β 95% HPD Upper",
        "β Scalar PP",
        "Mean Node Scalar",
        "Median Node Scalar",
        "Node Scalar 95% HPD Lower",
        "Node Scalar 95% HPD Upper",
        "Node Scalar PP",
    ];
    let values = rows
        .iter()
        .map(|r| vec![r.id as f64, r.n_taxa as f64, r.original_length, r.beta_mean, r.beta_median, r.beta_hpd.0, r.beta_hpd.1, r.beta_prop, r.node_mean, r.node_median, r.node_hpd.0, r.node_hpd.1, r.node_prop])
        .collect();
    (header.iter().map(|s| s.to_string()).collect(), values, vec![0, 0, 6, 6, 6, 6, 6, 9, 6, 6, 6, 6, 9])
}

/// The summary table as text cells (proportions to 9 decimals, the rest to 6).
pub fn table(rows: &[FabricBranch]) -> (Vec<String>, Vec<Vec<String>>) {
    let (header, values, digits) = columns(rows);
    let cells = values.iter().map(|r| r.iter().zip(&digits).map(|(&x, &d)| crate::tree::format_num(x, d)).collect()).collect();
    (header, cells)
}

#[cfg(test)]
mod tests {
    use super::*;

    const FILE: &str = "9\n0\tHomo_sapiens\n1\tPan_troglodytes\n2\tGorilla_gorilla\n3\tPongo_abelii\n4\tHylobates_lar\n5\tMacaca_mulatta\n6\tPapio_anubis\n7\tColobus_guereza\n8\tNasalis_larvatus\n17\n0\t-1.000000\t9\t0\t1\t2\t3\t4\t5\t6\t7\t8\t\n1\t10.195000\t5\t0\t1\t2\t3\t4\t\n2\t3.851600\t4\t0\t1\t2\t3\t\n3\t7.325500\t3\t0\t1\t2\t\n4\t1.993571\t2\t0\t1\t\n5\t6.591829\t1\t0\t\n6\t6.591829\t1\t1\t\n7\t8.585400\t1\t2\t\n8\t15.910500\t1\t3\t\n9\t19.762500\t1\t4\t\n10\t10.898300\t4\t5\t6\t7\t8\t\n11\t8.786356\t2\t5\t6\t\n12\t10.273344\t1\t5\t\n13\t10.273344\t1\t6\t\n14\t7.599467\t2\t7\t8\t\n15\t11.460033\t1\t7\t\n16\t11.460033\t1\t8\t\nIt\tLh\tLh + Prior\tNo Pram\tAlpha\tSigma^2\tAlpha Scale Prior\tNode ID\tScaler\tCreat It\tNode / Branch\t\n10000\t-45.376750\t-60.097217\t2\t19.207254\t61.958074\t1.200000\t11\t-0.947936\t9893\tLS_Beta\t15\t-0.447276\t9811\tLS_Beta\t\n11000\t-43.416938\t-59.511589\t2\t11.622481\t30.084890\t1.200000\t3\t3.379407\t10942\tNode\t16\t-1.977708\t10884\tLS_Beta\t\n12000\t-45.260738\t-64.123927\t3\t28.959254\t65.340008\t1.200000\t1\t-0.828450\t11863\tLS_Beta\t10\t0.486563\t11998\tLS_Beta\t11\t-2.117662\t11834\tLS_Beta\t\n";

    #[test]
    fn runs_combine_and_compare() {
        let a = parse(FILE).unwrap();
        let b = parse(FILE).unwrap();
        assert_eq!(a.lh.len(), 3);
        assert!((a.lh[0] - FILE.lines().nth(29).unwrap().split('\t').nth(1).unwrap().parse::<f64>().unwrap()).abs() < 1e-12);
        let both = Fabric::combine(&[a, b]).unwrap();
        assert_eq!((both.iterations.len(), both.lh.len()), (6, 6));
        // Identical runs: the same posterior probabilities in each, so no spread.
        let runs = [parse(FILE).unwrap(), parse(FILE).unwrap()];
        let r = convergence(&runs, &["a".to_string(), "b".to_string()]);
        assert_eq!(r.runs.len(), 2);
        assert_eq!(r.asd_pp(), (0.0, 0.0));
        assert_eq!(r.table().1.len(), runs[0].nodes.len());
        // Runs from a different tree can't be pooled.
        let other = FILE.replacen("3\t7.325500\t3\t0\t1\t2", "3\t7.325500\t3\t0\t1\t3", 1);
        assert!(Fabric::combine(&[parse(FILE).unwrap(), parse(&other).unwrap()]).is_err());
    }

    #[test]
    fn matches_target_table() {
        let (header, rows) = table(&parse(FILE).unwrap().summarize());
        assert_eq!(header[..5], ["Node ID", "No Taxa", "Original BL", "Mean (β * BL)", "Median (β * BL)"]);
        assert_eq!(header[5..], ["β 95% HPD Lower", "β 95% HPD Upper", "β Scalar PP", "Mean Node Scalar", "Median Node Scalar", "Node Scalar 95% HPD Lower", "Node Scalar 95% HPD Upper", "Node Scalar PP"]);
        let target = "0\t9\t0\t0\t0\t0\t0
1\t5\t10.195\t-0.82845\t-0.82845\t0.333333333\t0
2\t4\t3.8516\t0\t0\t0\t0
3\t3\t7.3255\t0\t0\t0\t0.333333333
4\t2\t1.993571\t0\t0\t0\t0
5\t1\t6.591829\t0\t0\t0\t0
6\t1\t6.591829\t0\t0\t0\t0
7\t1\t8.5854\t0\t0\t0\t0
8\t1\t15.9105\t0\t0\t0\t0
9\t1\t19.7625\t0\t0\t0\t0
10\t4\t10.8983\t0.486563\t0.486563\t0.333333333\t0
11\t2\t8.786356\t-1.532799\t-1.532799\t0.666666667\t0
12\t1\t10.273344\t0\t0\t0\t0
13\t1\t10.273344\t0\t0\t0\t0
14\t2\t7.599467\t0\t0\t0\t0
15\t1\t11.460033\t-0.447276\t-0.447276\t0.333333333\t0
16\t1\t11.460033\t-1.977708\t-1.977708\t0.333333333\t0";
        // The target's columns: the first five, Beta Scalar PP and Node Scalar PP.
        let got: Vec<String> = rows.iter().map(|r| format!("{}\t{}\t{}", r[..5].join("\t"), r[7], r[12])).collect();
        assert_eq!(got.join("\n"), target);
        // Mean / Median Node Scalar: node 3 had one node scalar (3.379407); the rest none (1).
        for r in &rows {
            let want = if r[0] == "3" { ["3.379407", "3.379407"] } else { ["1", "1"] };
            assert_eq!(r[8..10], want, "node {}", r[0]);
            // Node Scalar 95% HPD: one value or none.
            assert_eq!(r[10..12], want, "node {}", r[0]);
        }
        // Beta 95% HPD: node 11 had LS_Beta -0.947936 and -2.117662; node 1 one value; node 2 none.
        let hpd = |id: &str| rows.iter().find(|r| r[0] == id).map(|r| (r[5].clone(), r[6].clone())).unwrap();
        assert_eq!(hpd("11"), ("-2.117662".to_string(), "-0.947936".to_string()));
        assert_eq!(hpd("1"), ("-0.82845".to_string(), "-0.82845".to_string()));
        assert_eq!(hpd("2"), ("0".to_string(), "0".to_string()));
    }

    #[test]
    fn applies_to_tree_by_species() {
        let mut t = crate::io::newick::parse_newick(
            "(((((Homo_sapiens:6.5918285714,Pan_troglodytes:6.5918285714):1.9935714286,Gorilla_gorilla:8.5854):7.3255,Pongo_abelii:15.9105):3.8516,Hylobates_lar:19.7625):10.195,((Macaca_mulatta:10.2733444444,Papio_anubis:10.2733444444):8.7863555556,(Colobus_guereza:11.4600333333,Nasalis_larvatus:11.4600333333):7.5994666667):10.8983);",
        )
        .unwrap();
        let fab = parse(FILE).unwrap();
        let mut rows = fab.summarize();
        let warnings = fab.apply(&mut rows, &mut t);
        assert!(warnings.is_empty(), "{:?}", warnings);
        let mac = t.find_label("Macaca_mulatta").unwrap();
        let mp = t.nodes[mac].parent.unwrap();
        assert!((t.nodes[mp].attrs["fabric_beta_mean"].as_f64().unwrap() + 1.532799).abs() < 1e-6);
        assert!((t.nodes[mp].attrs["fabric_beta_prop"].as_f64().unwrap() - 2.0 / 3.0).abs() < 1e-9);
        // Node 3 (Homo, Pan, Gorilla) had a 3.379407 node scalar in one of three iterations,
        // which also rescales the branches inside it.
        let homo = t.find_label("Homo_sapiens").unwrap();
        let want = (2.0 + 3.379407) / 3.0;
        assert!((t.nodes[homo].attrs["fabric_node_scalar"].as_f64().unwrap() - want).abs() < 1e-9);
        assert!((t.nodes[mp].attrs["fabric_node_scalar"].as_f64().unwrap() - 1.0).abs() < 1e-12);
        let s = scaled_tree(&t);
        let h = s.find_label("Homo_sapiens").unwrap();
        assert!((s.nodes[h].length.unwrap() - 6.5918285714 * want).abs() < 1e-6);
        assert!(crate::varrates::node_posterior(&t, t.nodes[t.nodes[homo].parent.unwrap()].parent.unwrap()).is_some());
    }

    #[test]
    fn rejects_branch_scalars() {
        // A regular variable-rates run (Branch entries) is not Fabric output.
        let vr = FILE.replace("LS_Beta", "Branch");
        assert!(parse(&vr).is_err());
    }
}
