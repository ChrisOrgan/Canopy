//! Phylogenetic imbalance ratio (Gardner & Organ 2021): how suitable categorical
//! comparative data are for phylogenetic comparative models. PIR = NIR × CI,
//! where the normalized imbalance ratio NIR = (Tmax − Tmin) / n over the classes
//! (every combination of the characters' states) and CI is the ensemble
//! consistency index (Kluge & Farris 1969) on the tree. Values below 0.1 are
//! recommended.

use crate::io::data::{normalize_name, DataTable};
use crate::tree::Tree;
use std::collections::{HashMap, HashSet};

/// The PIR of a set of categorical characters, with its parts.
#[derive(Clone, Debug)]
pub struct Pir {
    /// Taxa in the tree with a state for every character.
    pub taxa: usize,
    /// Each class (state combination, "a, b") with its number of taxa, observed or not.
    pub classes: Vec<(String, usize)>,
    /// Minimum and parsimony (Fitch) number of state changes, summed over characters.
    pub min_steps: usize,
    pub steps: usize,
    pub ci: f64,
    pub nir: f64,
    pub pir: f64,
}

fn missing(v: &str) -> bool {
    v.is_empty() || v.eq_ignore_ascii_case("na") || v == "?" || v == "-"
}

/// Fitch (Hartigan for polytomies) parsimony length of a character with `k`
/// states (tip → state index in `obs`); tips not in `obs` count as any state.
pub fn fitch_steps(tree: &Tree, obs: &HashMap<usize, usize>, k: usize) -> usize {
    let full: u64 = if k >= 64 { u64::MAX } else { (1u64 << k) - 1 };
    let mut sets = vec![full; tree.nodes.len()];
    let mut steps = 0;
    for u in tree.postorder() {
        let ch = &tree.nodes[u].children;
        if ch.is_empty() {
            if let Some(&s) = obs.get(&u) {
                sets[u] = 1 << s;
            }
            continue;
        }
        let counts: Vec<usize> = (0..k).map(|j| ch.iter().filter(|&&c| sets[c] >> j & 1 == 1).count()).collect();
        let m = counts.iter().copied().max().unwrap_or(0);
        sets[u] = (0..k).filter(|&j| counts[j] == m).fold(0, |a, j| a | 1 << j);
        steps += ch.len() - m;
    }
    steps
}

/// PIR of the categorical columns `cols` of `data` on `tree`, over the tips in
/// `only` (normalized names) if given. None without taxa or with over 64 states.
pub fn pir(tree: &Tree, data: &DataTable, cols: &[usize], only: Option<&HashSet<String>>) -> Option<Pir> {
    let rows: HashMap<String, &Vec<String>> = data.rows.iter().map(|r| (normalize_name(&r[data.key]), r)).collect();
    // Tips with a state for every character.
    let tips: Vec<(usize, Vec<&str>)> = tree
        .tips()
        .into_iter()
        .filter_map(|t| {
            let key = normalize_name(tree.label(t));
            if only.map(|s| !s.contains(&key)).unwrap_or(false) {
                return None;
            }
            let r = rows.get(&key)?;
            let vals: Vec<&str> = cols.iter().map(|&c| r.get(c).map(|v| v.trim()).unwrap_or("")).collect();
            if vals.iter().any(|v| missing(v)) {
                None
            } else {
                Some((t, vals))
            }
        })
        .collect();
    if tips.is_empty() || cols.is_empty() {
        return None;
    }
    let (mut min_steps, mut steps) = (0, 0);
    let mut states: Vec<Vec<&str>> = Vec::new();
    for i in 0..cols.len() {
        let mut s: Vec<&str> = tips.iter().map(|t| t.1[i]).collect();
        s.sort_unstable();
        s.dedup();
        if s.len() > 64 {
            return None;
        }
        let obs: HashMap<usize, usize> = tips.iter().map(|t| (t.0, s.iter().position(|x| *x == t.1[i]).unwrap())).collect();
        min_steps += s.len() - 1;
        steps += fitch_steps(tree, &obs, s.len());
        states.push(s);
    }
    // Every combination of states is a class, including unobserved ones.
    let mut combos: Vec<Vec<&str>> = vec![Vec::new()];
    for s in &states {
        combos = combos.iter().flat_map(|c| s.iter().map(move |x| [c.clone(), vec![*x]].concat())).collect();
    }
    let classes: Vec<(String, usize)> = combos.iter().map(|c| (c.join(", "), tips.iter().filter(|t| t.1 == *c).count())).collect();
    let n = tips.len();
    let (tmax, tmin) = (classes.iter().map(|c| c.1).max().unwrap_or(0), classes.iter().map(|c| c.1).min().unwrap_or(0));
    let nir = (tmax - tmin) as f64 / n as f64;
    let ci = if steps > 0 { min_steps as f64 / steps as f64 } else { f64::NAN };
    Some(Pir { taxa: n, classes, min_steps, steps, ci, nir, pir: nir * ci })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick;

    #[test]
    fn darwins_scenario_and_homoplasy() {
        // Darwin's scenario: both characters change once, on the same branch.
        let t = parse_newick("(((A,B),(C,D)),((E,F),(G,H)));").unwrap();
        let d = DataTable::from_str("taxon,x,y\nA,1,1\nB,1,1\nC,1,1\nD,1,1\nE,0,0\nF,0,0\nG,0,0\nH,0,0\n", "d", None).unwrap();
        let p = pir(&t, &d, &[1, 2], None).unwrap();
        assert_eq!((p.taxa, p.min_steps, p.steps), (8, 2, 2));
        assert_eq!(p.classes.len(), 4);
        assert!((p.ci - 1.0).abs() < 1e-12 && (p.nir - 0.5).abs() < 1e-12 && (p.pir - 0.5).abs() < 1e-12);
        // Alternating states: x changes four times (CI 1/4); y has missing data for H.
        let d = DataTable::from_str("taxon,x,y\nA,1,a\nB,0,a\nC,1,b\nD,0,b\nE,1,b\nF,0,b\nG,1,a\nH,0,?\n", "d", None).unwrap();
        let p = pir(&t, &d, &[1], None).unwrap();
        assert_eq!((p.steps, p.taxa), (4, 8));
        assert!((p.ci - 0.25).abs() < 1e-12 && p.nir.abs() < 1e-12);
        let p = pir(&t, &d, &[1, 2], None).unwrap();
        assert_eq!(p.taxa, 7);
        // Restricted to a clade.
        let only: HashSet<String> = ["a", "b"].iter().map(|s| s.to_string()).collect();
        assert_eq!(pir(&t, &d, &[1], Some(&only)).unwrap().taxa, 2);
    }
}
