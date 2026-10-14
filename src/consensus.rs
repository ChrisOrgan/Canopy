//! Summaries of posterior tree samples: clade frequencies, majority-rule
//! (and greedy extended-majority) consensus trees, and the maximum clade
//! credibility (MCC) tree, with posterior support and mean heights. A summary is
//! a single tree, so it carries no medians or 95% HPD intervals; trees stepped
//! through in a set do (heights and branch lengths across the set).

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

/// 128-bit clade fingerprint: the XOR of random keys of the clade's taxa.
/// Collisions between different clades are vanishingly unlikely.
pub type Key = u128;

/// Hasher for fingerprints, which are already random.
#[derive(Default)]
pub struct KeyHasher(u64);

impl std::hash::Hasher for KeyHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = self.0.rotate_left(8) ^ b as u64;
        }
    }
    fn write_u128(&mut self, x: u128) {
        self.0 = x as u64 ^ (x >> 64) as u64;
    }
}

pub type KeyMap<V> = HashMap<Key, V, std::hash::BuildHasherDefault<KeyHasher>>;

#[derive(Clone, Debug, Default)]
pub struct CladeStats {
    pub count: usize,
    pub lengths: Vec<f64>,
    pub heights: Vec<f64>,
    /// The clade's taxa (unrooted: the side of the split without the first taxon).
    pub bits: Bits,
    /// A tip or the whole tree (or, unrooted, its complement): not a real clade.
    pub trivial: bool,
    /// Computed once by `summarize`: mean, median and 95% HPD of the node
    /// heights, and mean and 95% HPD of the branch lengths (NaN / None if none).
    pub height_mean: f64,
    pub height_median: f64,
    pub height_hpd: Option<(f64, f64)>,
    pub length_mean: f64,
    pub length_median: f64,
    pub length_hpd: Option<(f64, f64)>,
}

impl CladeStats {
    fn finish(&mut self) {
        self.height_mean = mean(&self.heights);
        self.height_median = median(&self.heights);
        self.height_hpd = hpd(&self.heights, 0.95);
        self.length_mean = mean(&self.lengths);
        self.length_median = median(&self.lengths);
        self.length_hpd = hpd(&self.lengths, 0.95);
    }
}

pub struct PosteriorSummary {
    pub n_trees: usize,
    pub taxa: Vec<String>,
    pub index: HashMap<String, usize>,
    /// Treat trees as rooted (clades) or unrooted (splits).
    pub rooted: bool,
    pub clades: KeyMap<CladeStats>,
    /// Random key per taxon, and their XOR (the whole tree).
    keys: Vec<Key>,
    all: Key,
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
    /// Fingerprint (unrooted: of the side without the first taxon) of each
    /// reachable node, in postorder.
    pub fn clades_of(&self, t: &Tree) -> Result<Vec<(NodeId, Key)>> {
        let mut fp: Vec<Key> = vec![0; t.nodes.len()];
        let mut has0 = vec![false; t.nodes.len()];
        let mut out = Vec::new();
        for id in t.postorder() {
            if t.is_tip(id) {
                let l = t.label(id);
                let Some(&i) = self.index.get(l) else { bail!("taxon '{}' not present in the first tree", l) };
                fp[id] = self.keys[i];
                has0[id] = i == 0;
            } else {
                for &c in &t.nodes[id].children {
                    fp[id] ^= fp[c];
                    has0[id] |= has0[c];
                }
            }
            let k = if !self.rooted && has0[id] { self.all ^ fp[id] } else { fp[id] };
            out.push((id, k));
        }
        Ok(out)
    }

    /// The taxon set of the clade at `id` (unrooted: the side without the first taxon).
    fn bits_of(&self, t: &Tree, id: NodeId) -> Bits {
        let n = self.taxa.len();
        let mut b = bits_new(n);
        for tip in t.tips_below(id) {
            if let Some(&i) = self.index.get(t.label(tip)) {
                bits_set(&mut b, i);
            }
        }
        if !self.rooted && bits_get(&b, 0) {
            for (i, w) in b.iter_mut().enumerate() {
                *w = !*w;
                let lo = i * 64;
                if lo + 64 > n {
                    let keep = n - lo;
                    *w &= if keep == 64 { u64::MAX } else { (1u64 << keep) - 1 };
                }
            }
        }
        b.into_boxed_slice()
    }

    fn is_trivial(&self, b: &[u64]) -> bool {
        let c = bits_count(b);
        let n = self.taxa.len();
        c <= 1 || c >= n || (!self.rooted && c >= n - 1)
    }

    pub fn frequency(&self, k: &Key) -> f64 {
        self.clades.get(k).map(|s| s.count as f64 / self.n_trees as f64).unwrap_or(0.0)
    }

    /// Count the clades of some trees (one thread's share).
    fn count(&self, trees: &[Tree]) -> Result<KeyMap<CladeStats>> {
        let mut clades: KeyMap<CladeStats> = KeyMap::default();
        let mut seen: KeyMap<bool> = KeyMap::default();
        for t in trees {
            if t.num_tips() != self.taxa.len() {
                bail!("trees have different numbers of tips ({} vs {})", t.num_tips(), self.taxa.len());
            }
            let heights = t.heights();
            // Unrooted, the two branches at a binary root are one edge and give the
            // same split: count it once per tree, joining their lengths.
            seen.clear();
            for (id, k) in self.clades_of(t)? {
                let len = t.nodes[id].length;
                let e = clades.entry(k).or_insert_with(|| {
                    let bits = self.bits_of(t, id);
                    CladeStats { trivial: self.is_trivial(&bits), bits, ..Default::default() }
                });
                match seen.get(&k) {
                    None => {
                        e.count += 1;
                        if let Some(l) = len {
                            e.lengths.push(l);
                        }
                        e.heights.push(heights[id]);
                        seen.insert(k, len.is_some());
                    }
                    Some(&pushed) => {
                        if let (true, Some(l), Some(last)) = (pushed, len, e.lengths.last_mut()) {
                            *last += l;
                        }
                    }
                }
            }
        }
        Ok(clades)
    }
}

/// Number of worker threads for `n` items.
fn threads_for(n: usize) -> usize {
    std::thread::available_parallelism().map(|p| p.get()).unwrap_or(1).min(n).max(1)
}

/// Count clades across a set of trees (after any burn-in has been removed), on
/// all cores, and compute each clade's height and length summaries once.
pub fn summarize(trees: &[Tree], rooted: bool) -> Result<PosteriorSummary> {
    let Some(first) = trees.first() else { bail!("no trees") };
    let taxa: Vec<String> = first.tips().into_iter().map(|n| first.label(n).to_string()).collect();
    let index: HashMap<String, usize> = taxa.iter().enumerate().map(|(i, t)| (t.clone(), i)).collect();
    if index.len() != taxa.len() {
        bail!("duplicate tip labels in first tree");
    }
    // Deterministic random keys (xorshift).
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let keys: Vec<Key> = (0..taxa.len()).map(|_| (next() as u128) << 64 | next() as u128).collect();
    let all = keys.iter().fold(0, |a, k| a ^ k);
    let mut s = PosteriorSummary { n_trees: trees.len(), taxa, index, rooted, clades: KeyMap::default(), keys, all, support_key: "posterior" };
    let chunk = trees.len().div_ceil(threads_for(trees.len()));
    let parts: Vec<Result<KeyMap<CladeStats>>> = std::thread::scope(|sc| {
        let s = &s;
        let handles: Vec<_> = trees.chunks(chunk).map(|c| sc.spawn(move || s.count(c))).collect();
        handles.into_iter().map(|h| h.join().expect("clade counting thread panicked")).collect()
    });
    for part in parts {
        for (k, st) in part? {
            match s.clades.entry(k) {
                std::collections::hash_map::Entry::Vacant(v) => {
                    v.insert(st);
                }
                std::collections::hash_map::Entry::Occupied(mut o) => {
                    let e = o.get_mut();
                    e.count += st.count;
                    e.lengths.extend(st.lengths);
                    e.heights.extend(st.heights);
                }
            }
        }
    }
    let mut stats: Vec<&mut CladeStats> = s.clades.values_mut().collect();
    let chunk = stats.len().div_ceil(threads_for(stats.len())).max(1);
    std::thread::scope(|sc| {
        for c in stats.chunks_mut(chunk) {
            sc.spawn(move || c.iter_mut().for_each(|st| st.finish()));
        }
    });
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
        // The set's mean, so `height` can stay the tree's own (summary trees rename it back).
        a.insert("height_mean".into(), Attr::Num(st.height_mean));
        a.insert("height_median".into(), Attr::Num(st.height_median));
        if let Some((lo, hi)) = st.height_hpd {
            a.insert("height_95%_HPD".into(), Attr::Range(lo, hi));
        }
    }
    if !st.lengths.is_empty() {
        a.insert("length_mean".into(), Attr::Num(st.length_mean));
        a.insert("length_median".into(), Attr::Num(st.length_median));
    }
    if let Some((lo, hi)) = st.length_hpd {
        a.insert("length_95%_HPD".into(), Attr::Range(lo, hi));
    }
}

/// Label the nodes of `t` with their clade frequency (support) and the height
/// and length summaries from the tree set; optionally set node heights to the
/// mean or median.
pub fn annotate(t: &mut Tree, s: &PosteriorSummary, heights: HeightMode) -> Result<()> {
    let clades = s.clades_of(t)?;
    let mut h: HashMap<NodeId, f64> = HashMap::new();
    for (id, k) in &clades {
        let internal = !t.is_tip(*id) && *id != t.root;
        if let Some(st) = s.clades.get(k) {
            annotate_node(t, *id, st, s.support_key, s.frequency(k), internal);
            let v = match heights {
                HeightMode::Mean => st.height_mean,
                HeightMode::Median => st.height_median,
                HeightMode::Keep => f64::NAN,
            };
            h.insert(*id, v);
        } else if internal {
            t.nodes[*id].attrs.insert(s.support_key.into(), Attr::Num(0.0));
        }
    }
    if heights != HeightMode::Keep {
        // A clade's mean (or median) height is taken over the trees that contain
        // it, so a weakly supported clade can come out older than its parent.
        // Cap each node at its parent's height (a zero-length branch) so that
        // branch lengths stay non-negative and, with tips at height 0, the tree
        // stays ultrametric. Clamping the branch alone would push tips apart.
        for id in t.preorder() {
            if let Some(p) = t.nodes[id].parent {
                if let (Some(&hp), Some(&hc)) = (h.get(&p), h.get(&id)) {
                    if hp.is_finite() && hc.is_finite() {
                        let hc = hc.min(hp);
                        h.insert(id, hc);
                        t.nodes[id].length = Some(hp - hc);
                    }
                }
            }
        }
    }
    Ok(())
}

/// Drop the spread of the sample (medians and 95% HPD intervals of heights and
/// lengths) from a summary tree: it is a single tree. Support stays, and the mean
/// height becomes its `height` (TreeAnnotator's convention).
pub fn strip_spread(t: &mut Tree) {
    for n in 0..t.nodes.len() {
        let a = &mut t.nodes[n].attrs;
        if let Some(h) = a.remove("height_mean") {
            a.insert("height".into(), h);
        }
        a.retain(|k, _| k != "length_mean" && !k.ends_with("_median") && !k.contains("_95%_HPD"));
    }
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
    let n_trees = s.n_trees as f64;
    let mut cands: Vec<&CladeStats> = s.clades.values().filter(|st| !st.trivial && keep(st.count as f64 / n_trees)).collect();
    cands.sort_by(|a, b| b.count.cmp(&a.count).then_with(|| a.bits.cmp(&b.bits)));
    let mut accepted: Vec<&Bits> = Vec::new();
    for st in &cands {
        if accepted.iter().all(|a| compatible(a, &st.bits)) {
            accepted.push(&st.bits);
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
        for (id, k) in clades {
            if id == t.root {
                continue;
            }
            if let Some(st) = s.clades.get(&k) {
                if !st.lengths.is_empty() {
                    t.nodes[id].length = Some(st.length_mean);
                }
            }
        }
    }
    let _ = annotate(&mut t, s, heights);
    strip_spread(&mut t);
    t.name = Some(format!("consensus_{:.0}", threshold * 100.0));
    t
}

/// Index of the maximum clade credibility tree (max sum of log clade
/// frequencies), scoring the trees on all cores.
pub fn mcc_index(trees: &[Tree], s: &PosteriorSummary) -> Result<(usize, f64)> {
    let score = |t: &Tree| -> Result<f64> {
        let mut score = 0.0;
        for (_, k) in s.clades_of(t)? {
            match s.clades.get(&k) {
                Some(st) if st.trivial => {}
                _ => score += s.frequency(&k).ln(),
            }
        }
        Ok(score)
    };
    let chunk = trees.len().div_ceil(threads_for(trees.len())).max(1);
    let parts: Vec<Result<Vec<f64>>> = std::thread::scope(|sc| {
        let handles: Vec<_> = trees.chunks(chunk).map(|c| sc.spawn(move || c.iter().map(score).collect::<Result<Vec<f64>>>())).collect();
        handles.into_iter().map(|h| h.join().expect("MCC scoring thread panicked")).collect()
    });
    let mut best = (0usize, f64::NEG_INFINITY);
    let mut i = 0;
    for part in parts {
        for sc in part? {
            if sc > best.1 {
                best = (i, sc);
            }
            i += 1;
        }
    }
    Ok(best)
}

pub fn mcc_tree(trees: &[Tree], s: &PosteriorSummary, heights: HeightMode) -> Result<Tree> {
    let (i, score) = mcc_index(trees, s)?;
    let mut t = trees[i].clone().compact();
    annotate(&mut t, s, heights)?;
    strip_spread(&mut t);
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
        // A single summary tree: mean heights, no medians or HPD intervals.
        assert!(m.nodes[p].attrs.contains_key("height") && !m.nodes[p].attrs.contains_key("height_95%_HPD") && !m.nodes[p].attrs.contains_key("height_median"));
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
    fn mean_heights_keep_trees_ultrametric() {
        // (A,B) is young in the first tree and old in the others, while (A,B,C) is
        // only in the first, so on that tree the mean height of (A,B) (4.33) is
        // older than its parent's (2).
        let trees = parse_newick_multi("(((A:1,B:1):1,C:2):8,D:10);(((A:6,B:6):3,D:9):1,C:10);(((A:6,B:6):3,D:9):1,C:10);").unwrap();
        let s = summarize(&trees, true).unwrap();
        let mut t = trees[0].clone();
        annotate(&mut t, &s, HeightMode::Mean).unwrap();
        let d = t.depths(true);
        let tips: Vec<f64> = t.tips().into_iter().map(|n| d[n]).collect();
        assert!(tips.iter().all(|x| (x - tips[0]).abs() < 1e-9), "{:?}", tips);
        assert!(t.preorder().iter().all(|&n| t.nodes[n].length.unwrap_or(0.0) >= 0.0));
        // The capped node sits at its parent's height: a zero-length branch.
        let ab = t.nodes[t.find_label("A").unwrap()].parent.unwrap();
        assert!(t.nodes[ab].length.unwrap().abs() < 1e-12);
    }

    #[test]
    fn hpd_interval() {
        let v: Vec<f64> = (0..100).map(|x| x as f64).collect();
        let (a, b) = hpd(&v, 0.95).unwrap();
        assert_eq!(b - a, 94.0);
    }
}
