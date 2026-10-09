//! Summaries of posterior tree samples: clade frequencies, majority-rule
//! (and greedy extended-majority) consensus trees, and the maximum clade
//! credibility (MCC) tree, with posterior support and 95% HPD intervals of
//! node heights in the style of TreeAnnotator.

use crate::tree::{Attr, NodeId, Tree};
use anyhow::{bail, Result};
use std::collections::HashMap;

pub type Bits = Box<[u64]>;

fn bits_new(n: usize) -> Vec<u64> {
    vec![0u64; n.div_ceil(64)]
}

fn bits_set(b: &mut [u64], i: usize) {
    b[i / 64] |= 1 << (i % 64);
}

fn bits_get(b: &[u64], i: usize) -> bool {
    b[i / 64] >> (i % 64) & 1 == 1
}

fn bits_count(b: &[u64]) -> usize {
    b.iter().map(|w| w.count_ones() as usize).sum()
}

fn is_subset(a: &[u64], b: &[u64]) -> bool {
    a.iter().zip(b).all(|(x, y)| x & !y == 0)
}

fn disjoint(a: &[u64], b: &[u64]) -> bool {
    a.iter().zip(b).all(|(x, y)| x & y == 0)
}

fn compatible(a: &[u64], b: &[u64]) -> bool {
    disjoint(a, b) || is_subset(a, b) || is_subset(b, a)
}

#[derive(Clone, Debug, Default)]
pub struct CladeStats {
    pub count: usize,
    pub lengths: Vec<f64>,
    pub heights: Vec<f64>,
}

pub struct PosteriorSummary {
    pub n_trees: usize,
    pub taxa: Vec<String>,
    pub index: HashMap<String, usize>,
    /// Treat trees as rooted (clades) or unrooted (splits).
    pub rooted: bool,
    pub clades: HashMap<Bits, CladeStats>,
    /// Attribute that receives clade frequencies: "posterior" for Bayesian
    /// samples, "bootstrap" for bootstrap replicates.
    pub support_key: &'static str,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HeightMode {
    /// Keep the target tree's own branch lengths.
    Keep,
    /// Mean node height across trees containing the clade.
    Mean,
    /// Median node height across trees containing the clade.
    Median,
}

impl PosteriorSummary {
    fn normalize(&self, mut b: Vec<u64>) -> Vec<u64> {
        if !self.rooted && bits_get(&b, 0) {
            let n = self.taxa.len();
            for (i, w) in b.iter_mut().enumerate() {
                *w = !*w;
                let lo = i * 64;
                if lo + 64 > n {
                    let keep = n - lo;
                    *w &= if keep == 64 { u64::MAX } else { (1u64 << keep) - 1 };
                }
            }
        }
        b
    }

    /// Clade (bitset) for each reachable node, in postorder.
    pub fn clades_of(&self, t: &Tree) -> Result<Vec<(NodeId, Bits)>> {
        let n = self.taxa.len();
        let mut sets: HashMap<NodeId, Vec<u64>> = HashMap::new();
        let mut out = Vec::new();
        for id in t.postorder() {
            let mut b = bits_new(n);
            if t.is_tip(id) {
                let l = t.label(id);
                let Some(&i) = self.index.get(l) else { bail!("taxon '{}' not present in the first tree", l) };
                bits_set(&mut b, i);
            } else {
                for c in &t.nodes[id].children {
                    for (w, x) in b.iter_mut().zip(&sets[c]) {
                        *w |= x;
                    }
                }
            }
            out.push((id, self.normalize(b.clone()).into_boxed_slice()));
            sets.insert(id, b);
        }
        Ok(out)
    }

    fn is_trivial(&self, b: &[u64]) -> bool {
        let c = bits_count(b);
        let n = self.taxa.len();
        c <= 1 || c >= n || (!self.rooted && c >= n - 1)
    }

    pub fn frequency(&self, b: &Bits) -> f64 {
        self.clades.get(b).map(|s| s.count as f64 / self.n_trees as f64).unwrap_or(0.0)
    }
}

/// Count clades across a set of trees (after any burn-in has been removed).
pub fn summarize(trees: &[Tree], rooted: bool) -> Result<PosteriorSummary> {
    let Some(first) = trees.first() else { bail!("no trees") };
    let taxa: Vec<String> = first.tips().into_iter().map(|n| first.label(n).to_string()).collect();
    let index: HashMap<String, usize> = taxa.iter().enumerate().map(|(i, t)| (t.clone(), i)).collect();
    if index.len() != taxa.len() {
        bail!("duplicate tip labels in first tree");
    }
    let mut s = PosteriorSummary { n_trees: trees.len(), taxa, index, rooted, clades: HashMap::new(), support_key: "posterior" };
    for t in trees {
        if t.num_tips() != s.taxa.len() {
            bail!("trees have different numbers of tips ({} vs {})", t.num_tips(), s.taxa.len());
        }
        let heights = t.heights();
        // Unrooted, the two branches at a binary root are one edge and give the
        // same split: count it once per tree, joining their lengths.
        let mut seen: HashMap<Bits, bool> = HashMap::new();
        for (id, b) in s.clades_of(t)? {
            let len = t.nodes[id].length;
            let e = s.clades.entry(b.clone()).or_default();
            match seen.get(&b) {
                None => {
                    e.count += 1;
                    if let Some(l) = len {
                        e.lengths.push(l);
                    }
                    e.heights.push(heights[id]);
                    seen.insert(b, len.is_some());
                }
                Some(&pushed) => {
                    if let (true, Some(l), Some(last)) = (pushed, len, e.lengths.last_mut()) {
                        *last += l;
                    }
                }
            }
        }
    }
    Ok(s)
}

pub fn mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.iter().sum::<f64>() / v.len() as f64
}

pub fn median(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = s.len() / 2;
    if s.len() % 2 == 0 {
        (s[m - 1] + s[m]) / 2.0
    } else {
        s[m]
    }
}

/// Shortest interval containing `prob` of the samples.
pub fn hpd(v: &[f64], prob: f64) -> Option<(f64, f64)> {
    if v.is_empty() {
        return None;
    }
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let k = ((prob * s.len() as f64).ceil() as usize).clamp(1, s.len());
    let mut best = (s[0], s[k - 1]);
    for i in 0..=s.len() - k {
        let (a, b) = (s[i], s[i + k - 1]);
        if b - a < best.1 - best.0 {
            best = (a, b);
        }
    }
    Some(best)
}

fn annotate_node(t: &mut Tree, id: NodeId, st: &CladeStats, key: &str, freq: f64, internal: bool) {
    let a = &mut t.nodes[id].attrs;
    if internal {
        a.insert(key.into(), Attr::Num(freq));
    }
    if !st.heights.is_empty() {
        a.insert("height".into(), Attr::Num(mean(&st.heights)));
        a.insert("height_median".into(), Attr::Num(median(&st.heights)));
        if let Some((lo, hi)) = hpd(&st.heights, 0.95) {
            a.insert("height_95%_HPD".into(), Attr::Range(lo, hi));
        }
    }
    if !st.lengths.is_empty() {
        if let Some((lo, hi)) = hpd(&st.lengths, 0.95) {
            a.insert("length_95%_HPD".into(), Attr::Range(lo, hi));
        }
    }
}

/// Add clade support (`s.support_key`) and HPD intervals to every node of `t`; optionally
/// replace its branch lengths using summarized node heights.
pub fn annotate(t: &mut Tree, s: &PosteriorSummary, heights: HeightMode) -> Result<()> {
    let clades = s.clades_of(t)?;
    let mut h: HashMap<NodeId, f64> = HashMap::new();
    for (id, b) in &clades {
        let internal = !t.is_tip(*id) && *id != t.root;
        if let Some(st) = s.clades.get(b) {
            annotate_node(t, *id, st, s.support_key, s.frequency(b), internal);
            let v = match heights {
                HeightMode::Mean => mean(&st.heights),
                HeightMode::Median => median(&st.heights),
                HeightMode::Keep => f64::NAN,
            };
            h.insert(*id, v);
        } else if internal {
            t.nodes[*id].attrs.insert(s.support_key.into(), Attr::Num(0.0));
        }
    }
    if heights != HeightMode::Keep {
        for id in t.preorder() {
            if let Some(p) = t.nodes[id].parent {
                if let (Some(&hp), Some(&hc)) = (h.get(&p), h.get(&id)) {
                    if hp.is_finite() && hc.is_finite() {
                        t.nodes[id].length = Some((hp - hc).max(0.0));
                    }
                }
            }
        }
    }
    Ok(())
}

/// Majority-rule consensus. With `threshold` < 0.5, clades are added greedily
/// in order of decreasing frequency while compatible (extended majority rule).
pub fn consensus(s: &PosteriorSummary, threshold: f64, heights: HeightMode) -> Tree {
    let keep = |f: f64| {
        if threshold >= 0.5 {
            f > 0.5 && f >= threshold - 1e-12
        } else {
            f >= threshold - 1e-12
        }
    };
    let mut cands: Vec<(&Bits, &CladeStats)> =
        s.clades.iter().filter(|(b, _)| !s.is_trivial(b) && keep(s.frequency(b))).collect();
    cands.sort_by(|a, b| b.1.count.cmp(&a.1.count).then_with(|| a.0.cmp(b.0)));
    let mut accepted: Vec<&Bits> = Vec::new();
    for (b, _) in &cands {
        if accepted.iter().all(|a| compatible(a, b)) {
            accepted.push(b);
        }
    }
    accepted.sort_by_key(|b| std::cmp::Reverse(bits_count(b)));

    let mut t = Tree::new();
    t.rooted = Some(s.rooted);
    let mut placed: Vec<(&Bits, NodeId)> = Vec::new();
    for b in &accepted {
        let parent = placed.iter().rev().find(|(pb, _)| is_subset(b, pb)).map(|(_, id)| *id).unwrap_or(t.root);
        let id = t.add_node(Some(parent));
        placed.push((b, id));
    }
    let n = s.taxa.len();
    for i in 0..n {
        let parent = placed.iter().rev().find(|(pb, _)| bits_get(pb, i)).map(|(_, id)| *id).unwrap_or(t.root);
        let id = t.add_node(Some(parent));
        t.nodes[id].label = Some(s.taxa[i].clone());
    }
    // Mean branch lengths from the sample, then posterior/HPD annotations.
    if let Ok(clades) = s.clades_of(&t) {
        for (id, b) in clades {
            if id == t.root {
                continue;
            }
            if let Some(st) = s.clades.get(&b) {
                if !st.lengths.is_empty() {
                    t.nodes[id].length = Some(mean(&st.lengths));
                }
            }
        }
    }
    let _ = annotate(&mut t, s, heights);
    t.name = Some(format!("consensus_{:.0}", threshold * 100.0));
    t
}

/// Index of the maximum clade credibility tree (max sum of log clade frequencies).
pub fn mcc_index(trees: &[Tree], s: &PosteriorSummary) -> Result<(usize, f64)> {
    let mut best = (0usize, f64::NEG_INFINITY);
    for (i, t) in trees.iter().enumerate() {
        let mut score = 0.0;
        for (_, b) in s.clades_of(t)? {
            if !s.is_trivial(&b) {
                score += s.frequency(&b).ln();
            }
        }
        if score > best.1 {
            best = (i, score);
        }
    }
    Ok(best)
}

pub fn mcc_tree(trees: &[Tree], s: &PosteriorSummary, heights: HeightMode) -> Result<Tree> {
    let (i, score) = mcc_index(trees, s)?;
    let mut t = trees[i].clone().compact();
    annotate(&mut t, s, heights)?;
    t.name = Some(format!("MCC (sample {}, log clade credibility {:.3})", i + 1, score));
    Ok(t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick_multi;

    fn sample() -> Vec<Tree> {
        parse_newick_multi(
            "((A:1,B:1):1,(C:1,D:1):1);\
             ((A:1,B:1):1,(C:1,D:1):1);\
             ((A:1,C:1):1,(B:1,D:1):1);\
             ((A:2,B:2):1,(C:1,D:1):2);",
        )
        .unwrap()
    }

    #[test]
    fn frequencies() {
        let trees = sample();
        let s = summarize(&trees, true).unwrap();
        let t = &trees[0];
        let a = t.find_label("A").unwrap();
        let ab = t.nodes[a].parent.unwrap();
        let clades = s.clades_of(t).unwrap();
        let b = &clades.iter().find(|(id, _)| *id == ab).unwrap().1;
        assert!((s.frequency(b) - 0.75).abs() < 1e-12);
    }

    #[test]
    fn majority_and_mcc() {
        let trees = sample();
        let s = summarize(&trees, true).unwrap();
        let c = consensus(&s, 0.5, HeightMode::Keep);
        assert_eq!(c.num_tips(), 4);
        let a = c.find_label("A").unwrap();
        let p = c.nodes[a].parent.unwrap();
        assert_eq!(c.tips_below(p).len(), 2);
        assert_eq!(c.nodes[p].attrs["posterior"], Attr::Num(0.75));
        let m = mcc_tree(&trees, &s, HeightMode::Mean).unwrap();
        let a = m.find_label("A").unwrap();
        let p = m.nodes[a].parent.unwrap();
        assert!(m.nodes[p].attrs.contains_key("height_95%_HPD"));
        assert_eq!(m.tips_below(p).len(), 2);
    }

    #[test]
    fn unrooted_splits() {
        let trees = parse_newick_multi("(A,B,(C,D));((A,B),C,D);").unwrap();
        let s = summarize(&trees, false).unwrap();
        let c = consensus(&s, 0.5, HeightMode::Keep);
        // AB|CD is shared by both unrooted trees.
        let cc = c.find_label("C").unwrap();
        let p = c.nodes[cc].parent.unwrap();
        assert_eq!(c.nodes[p].attrs["posterior"], Attr::Num(1.0));
    }

    #[test]
    fn hpd_interval() {
        let v: Vec<f64> = (0..100).map(|x| x as f64).collect();
        let (a, b) = hpd(&v, 0.95).unwrap();
        assert_eq!(b - a, 94.0);
    }
}
