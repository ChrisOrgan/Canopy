//! Ancestral state estimation for mapping comparative data onto branches.
//!
//! * Continuous traits: maximum-likelihood ancestral states under Brownian
//!   motion (equivalent to phytools `fastAnc`), computed with an exact
//!   two-pass algorithm in O(n).
//! * Discrete traits: Fitch parsimony with a deterministic tie-break.

use crate::tree::{Attr, NodeId, Tree};

const EPS: f64 = 1e-9;

/// Combine independent Gaussian estimates (mean, variance) by inverse variance.
fn combine(items: impl Iterator<Item = (f64, f64)>) -> Option<(f64, f64)> {
    let (mut w, mut wx) = (0.0, 0.0);
    for (x, v) in items {
        let wi = 1.0 / v.max(EPS);
        w += wi;
        wx += wi * x;
    }
    if w > 0.0 {
        Some((wx / w, 1.0 / w))
    } else {
        None
    }
}

/// Estimate the continuous attribute `key` at internal nodes (and tips with
/// missing data). Writes the estimates to `key` on those nodes and returns
/// the number of nodes filled in.
pub fn reconstruct_continuous(t: &mut Tree, key: &str) -> usize {
    let n = t.nodes.len();
    let len = |t: &Tree, i: NodeId| t.nodes[i].length.unwrap_or(1.0).max(EPS);
    let observed: Vec<Option<f64>> = (0..n)
        .map(|i| if t.is_tip(i) { t.nodes[i].attrs.get(key).and_then(Attr::as_f64).filter(|v| v.is_finite()) } else { None })
        .collect();

    // Down pass: estimate at each node from its own subtree.
    let mut down: Vec<Option<(f64, f64)>> = vec![None; n];
    for u in t.postorder() {
        down[u] = if t.is_tip(u) {
            observed[u].map(|x| (x, 0.0))
        } else {
            combine(t.nodes[u].children.iter().filter_map(|&c| down[c].map(|(x, v)| (x, v + len(t, c)))))
        };
    }
    // Up pass: estimate at each node from everything outside its subtree.
    let mut up: Vec<Option<(f64, f64)>> = vec![None; n];
    for u in t.preorder() {
        let kids = t.nodes[u].children.clone();
        for &c in &kids {
            let mut parts: Vec<(f64, f64)> =
                kids.iter().filter(|&&d| d != c).filter_map(|&d| down[d].map(|(x, v)| (x, v + len(t, d)))).collect();
            if let Some((x, v)) = up[u] {
                parts.push((x, v));
            }
            up[c] = combine(parts.into_iter()).map(|(x, v)| (x, v + len(t, c)));
        }
    }
    let mut filled = 0;
    for u in t.preorder() {
        if observed[u].is_some() {
            continue;
        }
        // Up contributions already include the branch above u.
        let below = if t.is_tip(u) { None } else { down[u] };
        if let Some((x, _)) = combine(below.into_iter().chain(up[u])) {
            t.nodes[u].attrs.insert(key.to_string(), Attr::Num(x));
            filled += 1;
        }
    }
    filled
}

/// Fitch parsimony reconstruction for a categorical attribute.
pub fn reconstruct_discrete(t: &mut Tree, key: &str) -> usize {
    let n = t.nodes.len();
    let mut states: Vec<String> = Vec::new();
    let mut obs: Vec<Option<usize>> = vec![None; n];
    for tip in t.tips() {
        if let Some(a) = t.nodes[tip].attrs.get(key) {
            let s = a.display(6);
            let i = match states.iter().position(|x| *x == s) {
                Some(i) => i,
                None => {
                    states.push(s);
                    states.len() - 1
                }
            };
            obs[tip] = Some(i);
        }
    }
    let k = states.len();
    if k == 0 {
        return 0;
    }
    let full = vec![true; k];
    let mut sets: Vec<Vec<bool>> = vec![full.clone(); n];
    for u in t.postorder() {
        if t.is_tip(u) {
            if let Some(i) = obs[u] {
                sets[u] = (0..k).map(|j| j == i).collect();
            }
        } else {
            let ch = &t.nodes[u].children;
            let inter: Vec<bool> = (0..k).map(|j| ch.iter().all(|&c| sets[c][j])).collect();
            sets[u] = if inter.iter().any(|&b| b) {
                inter
            } else {
                // Union, keeping only states present in the most children.
                let counts: Vec<usize> = (0..k).map(|j| ch.iter().filter(|&&c| sets[c][j]).count()).collect();
                let m = *counts.iter().max().unwrap();
                counts.iter().map(|&c| c == m).collect()
            };
        }
    }
    let mut assigned: Vec<Option<usize>> = vec![None; n];
    let mut filled = 0;
    for u in t.preorder() {
        let choice = match t.nodes[u].parent.and_then(|p| assigned[p]) {
            Some(ps) if sets[u][ps] => ps,
            _ => sets[u].iter().position(|&b| b).unwrap_or(0),
        };
        assigned[u] = Some(choice);
        if obs[u].is_none() {
            t.nodes[u].attrs.insert(key.to_string(), Attr::Text(states[choice].clone()));
            filled += 1;
        }
    }
    filled
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick;

    #[test]
    fn continuous_two_tips() {
        let mut t = parse_newick("(A:1,B:3);").unwrap();
        let a = t.find_label("A").unwrap();
        let b = t.find_label("B").unwrap();
        t.nodes[a].attrs.insert("x".into(), Attr::Num(0.0));
        t.nodes[b].attrs.insert("x".into(), Attr::Num(4.0));
        reconstruct_continuous(&mut t, "x");
        let r = t.nodes[t.root].attrs["x"].as_f64().unwrap();
        assert!((r - 1.0).abs() < 1e-9);
    }

    #[test]
    fn continuous_matches_reroot_estimate() {
        // For node AB in ((A:1,B:1):1,C:2), the ML estimate is the inverse-variance
        // mean of A (v=1), B (v=1) and C through the root (v=3).
        let mut t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        for (l, v) in [("A", 1.0), ("B", 3.0), ("C", 8.0)] {
            let n = t.find_label(l).unwrap();
            t.nodes[n].attrs.insert("x".into(), Attr::Num(v));
        }
        reconstruct_continuous(&mut t, "x");
        let ab = t.nodes[t.find_label("A").unwrap()].parent.unwrap();
        let expected = (1.0 + 3.0 + 8.0 / 3.0) / (1.0 + 1.0 + 1.0 / 3.0);
        let got = t.nodes[ab].attrs["x"].as_f64().unwrap();
        assert!((got - expected).abs() < 1e-9, "{} vs {}", got, expected);
    }

    #[test]
    fn fitch() {
        let mut t = parse_newick("((A,B),(C,D));").unwrap();
        for (l, v) in [("A", "red"), ("B", "red"), ("C", "blue"), ("D", "red")] {
            let n = t.find_label(l).unwrap();
            t.nodes[n].attrs.insert("c".into(), Attr::Text(v.into()));
        }
        reconstruct_discrete(&mut t, "c");
        assert_eq!(t.nodes[t.root].attrs["c"], Attr::Text("red".into()));
    }
}
