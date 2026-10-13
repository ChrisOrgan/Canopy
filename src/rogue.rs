//! Entropy-based rogue detection on a set of rooted trees (Klawitter, Bouckaert &
//! Drummond 2024, "Skeletons in the forest").
//!
//! The tree set is summarized by a conditional clade distribution (CCD0: every
//! observed clade, split into any pair of observed clades, with the probability
//! of a tree proportional to the product of its clade frequencies). Its entropy
//! measures the uncertainty in the set. A clade's rogue score is how much the
//! entropy would drop if its taxa were removed from every tree: high scores flag
//! taxa whose unstable placement causes the uncertainty. The skeleton analysis
//! repeats this, setting aside the top-scoring clade each round and re-scoring
//! the rest, to rank rogues in order. Both only flag taxa; the trees are never
//! changed.
//!
//! Klawitter, J., Bouckaert, R. R. & Drummond, A. J. (2024) Skeletons in the forest:
//! using entropy-based rogue detection on Bayesian phylogenetic tree distributions.
//! bioRxiv, doi:10.1101/2024.09.25.615070.

use crate::tree::Tree;
use anyhow::{bail, Result};
use std::collections::HashMap;

/// Size limits: building the CCD0 pairs up clades, so the work grows with the
/// square of their number (about 20 s for 1,000 taxa and 100 trees).
pub const MAX_TAXA: usize = 1000;
pub const MAX_CLADES: usize = 100_000;

/// A set of taxa (indices into `Ccd::taxa`).
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Bits(Vec<u64>);

impl Bits {
    fn empty(n: usize) -> Bits {
        Bits(vec![0; n.div_ceil(64).max(1)])
    }
    fn set(&mut self, i: usize) {
        self.0[i / 64] |= 1 << (i % 64);
    }
    fn count(&self) -> usize {
        self.0.iter().map(|w| w.count_ones() as usize).sum()
    }
    fn is_empty(&self) -> bool {
        self.0.iter().all(|&w| w == 0)
    }
    fn or(&self, o: &Bits) -> Bits {
        Bits(self.0.iter().zip(&o.0).map(|(a, b)| a | b).collect())
    }
    fn minus(&self, o: &Bits) -> Bits {
        Bits(self.0.iter().zip(&o.0).map(|(a, b)| a & !b).collect())
    }
    fn is_subset(&self, o: &Bits) -> bool {
        self.0.iter().zip(&o.0).all(|(a, b)| a & !b == 0)
    }
    fn ones(&self) -> impl Iterator<Item = usize> + '_ {
        self.0.iter().enumerate().flat_map(|(k, &w)| (0..64).filter(move |b| w >> b & 1 == 1).map(move |b| k * 64 + b))
    }
}

/// A conditional clade distribution: clades (smallest first) and, for each, its
/// splits into two child clades with their conditional probabilities.
pub struct Ccd {
    pub taxa: Vec<String>,
    clades: Vec<Bits>,
    splits: Vec<Vec<(usize, usize, f64)>>,
    root: usize,
}

/// One clade's rogue score: the drop in entropy when its taxa are removed.
#[derive(Clone, Debug)]
pub struct RogueScore {
    pub taxa: Vec<String>,
    /// Probability of the clade in the tree distribution.
    pub probability: f64,
    pub score: f64,
}

/// One round of the skeleton analysis (round 0 is the full tree set): the clade
/// set aside this round and the distribution of the taxa still in.
#[derive(Clone, Debug)]
pub struct Round {
    pub flagged: Vec<String>,
    pub taxa_left: usize,
    pub entropy: f64,
    /// Probability of the most probable tree on the taxa still in.
    pub map_probability: f64,
}

/// Canonical order for the two children of a split.
fn ordered(a: Bits, b: Bits) -> (Bits, Bits) {
    if a.0 <= b.0 {
        (a, b)
    } else {
        (b, a)
    }
}

impl Ccd {
    /// CCD0 of a set of rooted, binary trees on the same taxa.
    pub fn from_trees(trees: &[Tree]) -> Result<Ccd> {
        let Some(first) = trees.first() else { bail!("no trees") };
        let mut taxa: Vec<String> = first.tips().into_iter().map(|t| first.label(t).to_string()).collect();
        taxa.sort();
        let index: HashMap<&str, usize> = taxa.iter().enumerate().map(|(i, t)| (t.as_str(), i)).collect();
        let n = taxa.len();
        if n < 3 {
            bail!("rogue detection needs at least three taxa");
        }
        if n > MAX_TAXA {
            bail!("rogue detection is limited to {} taxa (this tree set has {}): its clade distribution (CCD0) pairs up clades, so the work grows with the square of the number of clades", MAX_TAXA, n);
        }
        // Clade frequencies.
        let mut freq: HashMap<Bits, f64> = HashMap::new();
        for t in trees {
            let mut below: Vec<Option<Bits>> = vec![None; t.nodes.len()];
            for v in t.postorder() {
                let bits = if t.is_tip(v) {
                    let Some(&i) = index.get(t.label(v)) else { bail!("taxon '{}' is not in every tree", t.label(v)) };
                    let mut b = Bits::empty(n);
                    b.set(i);
                    b
                } else {
                    t.nodes[v].children.iter().fold(Bits::empty(n), |acc, &c| acc.or(below[c].as_ref().unwrap()))
                };
                *freq.entry(bits.clone()).or_default() += 1.0;
                below[v] = Some(bits);
            }
            if t.num_tips() != n {
                bail!("every tree must have the same {} taxa", n);
            }
        }
        if freq.len() > MAX_CLADES {
            bail!("rogue detection is limited to {} distinct clades (this tree set has {}); try fewer trees or a larger burn-in", MAX_CLADES, freq.len());
        }
        let k = trees.len() as f64;
        let freq: HashMap<Bits, f64> = freq.into_iter().map(|(b, c)| (b, c / k)).collect();
        // Every split of an observed clade into two observed clades.
        let mut clades: Vec<Bits> = freq.keys().cloned().collect();
        clades.sort_by_key(|b| (b.count(), b.0.clone()));
        let id: HashMap<Bits, usize> = clades.iter().cloned().enumerate().map(|(i, b)| (b, i)).collect();
        let mut by_size: Vec<Vec<usize>> = vec![Vec::new(); n + 1];
        for (i, c) in clades.iter().enumerate() {
            by_size[c.count()].push(i);
        }
        let mut child_pairs: Vec<Vec<(usize, usize)>> = vec![Vec::new(); clades.len()];
        for (ci, c) in clades.iter().enumerate() {
            let size = c.count();
            for s in 1..=size / 2 {
                for &ai in &by_size[s] {
                    let a = &clades[ai];
                    if !a.is_subset(c) {
                        continue;
                    }
                    let b = c.minus(a);
                    if let Some(&bi) = id.get(&b) {
                        // Count each pair once when both halves are the same size.
                        if s * 2 < size || a.0 < b.0 {
                            child_pairs[ci].push((ai, bi));
                        }
                    }
                }
            }
        }
        // log g(C) = log f(C) + log Σ g(A) g(B); a tree's probability is ∏ f(C) / g(root).
        let mut lg = vec![f64::NEG_INFINITY; clades.len()];
        let mut splits = vec![Vec::new(); clades.len()];
        for ci in 0..clades.len() {
            if clades[ci].count() == 1 {
                lg[ci] = 0.0;
                continue;
            }
            let terms: Vec<f64> = child_pairs[ci].iter().map(|&(a, b)| lg[a] + lg[b]).collect();
            let m = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            if !m.is_finite() {
                continue;
            }
            let sum = m + terms.iter().map(|t| (t - m).exp()).sum::<f64>().ln();
            lg[ci] = freq[&clades[ci]].ln() + sum;
            splits[ci] = child_pairs[ci].iter().zip(&terms).filter(|(_, t)| t.is_finite()).map(|(&(a, b), t)| (a, b, (t - sum).exp())).collect();
        }
        let root = clades.len() - 1;
        if clades[root].count() != n || splits[root].is_empty() {
            bail!("the trees must be rooted and binary");
        }
        Ok(Ccd { taxa, clades, splits, root })
    }

    /// Probability of each clade under the distribution.
    fn clade_probs(&self) -> Vec<f64> {
        let mut p = vec![0.0; self.clades.len()];
        p[self.root] = 1.0;
        for c in (0..self.clades.len()).rev() {
            if p[c] > 0.0 {
                for &(a, b, q) in &self.splits[c] {
                    p[a] += p[c] * q;
                    p[b] += p[c] * q;
                }
            }
        }
        p
    }

    /// Entropy (in nats) of the tree distribution.
    pub fn entropy(&self) -> f64 {
        let p = self.clade_probs();
        let h: f64 = (0..self.clades.len()).map(|c| p[c] * self.splits[c].iter().map(|s| if s.2 > 0.0 { -s.2 * s.2.ln() } else { 0.0 }).sum::<f64>()).sum();
        h.max(0.0)
    }

    /// Probability of the most probable tree.
    pub fn map_probability(&self) -> f64 {
        let mut best = vec![0.0f64; self.clades.len()];
        for c in 0..self.clades.len() {
            best[c] = if self.splits[c].is_empty() { 0.0 } else { self.splits[c].iter().map(|&(a, b, q)| q.ln() + best[a] + best[b]).fold(f64::NEG_INFINITY, f64::max) };
        }
        best[self.root].exp()
    }

    pub fn num_taxa(&self) -> usize {
        self.clades[self.root].count()
    }

    /// Split weights of the distribution with the taxa `x` removed: for each
    /// remaining clade, the probability of each way it splits.
    fn filtered_weights(&self, probs: &[f64], x: &Bits) -> HashMap<Bits, HashMap<(Bits, Bits), f64>> {
        let mut w: HashMap<Bits, HashMap<(Bits, Bits), f64>> = HashMap::new();
        for c in 0..self.clades.len() {
            if probs[c] <= 0.0 {
                continue;
            }
            for &(a, b, q) in &self.splits[c] {
                let (a2, b2) = (self.clades[a].minus(x), self.clades[b].minus(x));
                // A split that loses one side entirely collapses into the other.
                if !a2.is_empty() && !b2.is_empty() {
                    *w.entry(self.clades[c].minus(x)).or_default().entry(ordered(a2, b2)).or_default() += probs[c] * q;
                }
            }
        }
        w
    }

    /// Entropy after removing the taxa `x`, without building the new distribution.
    fn filtered_entropy(&self, probs: &[f64], x: &Bits) -> f64 {
        self.filtered_weights(probs, x)
            .values()
            .map(|s| {
                let total: f64 = s.values().sum();
                s.values().filter(|&&v| v > 0.0).map(|&v| v * (total / v).ln()).sum::<f64>()
            })
            .sum()
    }

    /// The distribution with the taxa `x` left out of every tree (used only inside
    /// the skeleton analysis).
    fn filter(&self, x: &Bits) -> Ccd {
        let w = self.filtered_weights(&self.clade_probs(), x);
        let mut clades: Vec<Bits> = w.keys().cloned().collect();
        for s in w.values() {
            for (a, b) in s.keys() {
                clades.push(a.clone());
                clades.push(b.clone());
            }
        }
        clades.sort_by_key(|b| (b.count(), b.0.clone()));
        clades.dedup();
        let id: HashMap<&Bits, usize> = clades.iter().enumerate().map(|(i, b)| (b, i)).collect();
        let mut splits = vec![Vec::new(); clades.len()];
        for (c, s) in &w {
            let total: f64 = s.values().sum();
            splits[id[c]] = s.iter().map(|((a, b), v)| (id[a], id[b], v / total)).collect();
        }
        let root = clades.len() - 1;
        Ccd { taxa: self.taxa.clone(), clades, splits, root }
    }

    fn names(&self, b: &Bits) -> Vec<String> {
        b.ones().map(|i| self.taxa[i].clone()).collect()
    }

    /// Clades to score: at most `max_size` taxa, probability at least
    /// `min_prob`, leaving at least three taxa.
    fn candidates(&self, probs: &[f64], max_size: usize, min_prob: f64) -> Vec<usize> {
        let left = self.num_taxa();
        (0..self.clades.len()).filter(|&c| c != self.root && probs[c] >= min_prob - 1e-12 && self.clades[c].count() <= max_size && left - self.clades[c].count() >= 3).collect()
    }

    /// Rogue scores of the candidate clades, highest first.
    pub fn rogue_scores(&self, max_size: usize, min_prob: f64) -> Vec<RogueScore> {
        let probs = self.clade_probs();
        let h = self.entropy();
        let mut out: Vec<RogueScore> = self
            .candidates(&probs, max_size, min_prob)
            .into_iter()
            .map(|c| RogueScore { taxa: self.names(&self.clades[c]), probability: probs[c], score: h - self.filtered_entropy(&probs, &self.clades[c]) })
            .collect();
        out.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));
        out
    }

    /// Skeleton analysis: each round, set aside the clade with the highest rogue
    /// score and re-score the rest, until the entropy is (near) zero, no clade
    /// lowers it, no candidate is left, or three taxa remain. Flags rogues in
    /// order; the trees themselves are not changed.
    pub fn skeleton(&self, max_size: usize, min_prob: f64) -> Vec<Round> {
        let mut ccd = self.filter(&Bits::empty(self.taxa.len()));
        let mut rounds = vec![Round { flagged: Vec::new(), taxa_left: ccd.num_taxa(), entropy: ccd.entropy(), map_probability: ccd.map_probability() }];
        while rounds.last().unwrap().entropy > 1e-9 {
            let probs = ccd.clade_probs();
            let best = ccd
                .candidates(&probs, max_size, min_prob)
                .into_iter()
                .map(|c| (c, ccd.filtered_entropy(&probs, &ccd.clades[c])))
                // Near-ties (sums in a different order) go to the first clade, so runs repeat.
                .fold(None, |best: Option<(usize, f64)>, (c, h)| match best {
                    Some((_, bh)) if h >= bh - 1e-9 => best,
                    _ => Some((c, h)),
                });
            // Stop when no clade lowers the entropy any further.
            let Some((c, h)) = best else { break };
            if h > rounds.last().unwrap().entropy - 1e-6 {
                break;
            }
            let x = ccd.clades[c].clone();
            ccd = ccd.filter(&x);
            rounds.push(Round { flagged: self.names(&x), taxa_left: ccd.num_taxa(), entropy: ccd.entropy(), map_probability: ccd.map_probability() });
        }
        rounds
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick;

    fn trees(nwk: &[&str]) -> Vec<Tree> {
        nwk.iter().map(|s| parse_newick(s).unwrap()).collect()
    }

    #[test]
    fn identical_trees_have_zero_entropy() {
        let ccd = Ccd::from_trees(&trees(&["(((A,B),C),D);", "(((A,B),C),D);"])).unwrap();
        assert!(ccd.entropy().abs() < 1e-12);
        assert!((ccd.map_probability() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn a_rogue_taxon_carries_the_entropy() {
        // X jumps between three places on an otherwise fixed tree ((A,B),(C,D)).
        let set = trees(&["(((A,X),B),(C,D));", "((A,B),((C,X),D));", "(((A,B),X),(C,D));"]);
        let ccd = Ccd::from_trees(&set).unwrap();
        let h = ccd.entropy();
        assert!(h > 0.5, "entropy {}", h);
        let scores = ccd.rogue_scores(1, 0.0);
        assert_eq!(scores[0].taxa, ["X"]);
        // Removing X leaves no uncertainty at all.
        assert!((scores[0].score - h).abs() < 1e-9);
        let rounds = ccd.skeleton(1, 0.0);
        assert_eq!(rounds.len(), 2);
        assert_eq!(rounds[1].flagged, ["X"]);
        assert_eq!(rounds[1].taxa_left, 4);
        assert!(rounds[1].entropy.abs() < 1e-12 && (rounds[1].map_probability - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ccd0_tree_probabilities_sum_to_one() {
        // Two different trees on four taxa: the CCD0 mixes their clades.
        let ccd = Ccd::from_trees(&trees(&["(((A,B),C),D);", "((A,B),(C,D));", "(((A,B),C),D);"])).unwrap();
        let p = ccd.clade_probs();
        assert!((p[ccd.root] - 1.0).abs() < 1e-12);
        let total: f64 = ccd.splits[ccd.root].iter().map(|s| s.2).sum();
        assert!((total - 1.0).abs() < 1e-12);
        // Entropy of a two-tree distribution with probabilities 2/3 and 1/3 (both
        // trees share (A,B); CCD0 weights each tree by the product of clade frequencies).
        let h = ccd.entropy();
        assert!(h > 0.0 && h < 2f64.ln() + 1e-9);
    }
}
