//! Lineages through time: the number of lineages in a tree at each time, from
//! the root to the youngest tip (Nee et al. 1992). Fossil tips end lineages
//! where they sit, so non-ultrametric trees are handled.

use crate::tree::Tree;

/// Steps of the lineage count: (depth from the root, lineages from there on),
/// starting at the root with the number of its children. Depths of the tips'
/// ends are included, so the curve stops at the deepest tip.
pub fn ltt(tree: &Tree) -> Vec<(f64, usize)> {
    let d = tree.depths(true);
    // Events: +k-1 at a node with k children, -1 where a lineage ends at a tip
    // younger than the deepest one. Tips within 0.1% of the deepest count as
    // extant: MCC and consensus trees summed from rounded branch lengths can
    // leave extant tips a few thousandths of a percent short.
    let deepest = tree.tips().into_iter().map(|t| d[t]).fold(0.0, f64::max);
    let mut events: Vec<(f64, i64)> = Vec::new();
    for n in tree.preorder() {
        let k = tree.nodes[n].children.len() as i64;
        if k > 0 && n != tree.root {
            events.push((d[n], k - 1));
        } else if k == 0 && d[n] < deepest - 1e-3 * deepest.max(1e-12) {
            events.push((d[n], -1));
        }
    }
    events.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(b.1.cmp(&a.1)));
    let mut count = tree.nodes[tree.root].children.len() as i64;
    let mut out = vec![(0.0, count as usize)];
    for (x, k) in events {
        count += k;
        match out.last_mut() {
            Some(last) if (last.0 - x).abs() < 1e-12 => last.1 = count.max(0) as usize,
            _ => out.push((x, count.max(0) as usize)),
        }
    }
    out.push((deepest, count.max(0) as usize));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick;

    #[test]
    fn counts_lineages() {
        // Splits at 0 (root), 1 and 2; four tips at depth 3.
        let t = parse_newick("((A:2,B:2):1,(C:1,D:1):2);").unwrap();
        assert_eq!(ltt(&t), vec![(0.0, 2), (1.0, 3), (2.0, 4), (3.0, 4)]);
        // A fossil tip (E, ending at 1.5) removes a lineage.
        let f = parse_newick("((A:2,E:0.5):1,(C:1,D:1):2);").unwrap();
        assert_eq!(ltt(&f), vec![(0.0, 2), (1.0, 3), (1.5, 2), (2.0, 3), (3.0, 3)]);
        // Tips a rounding error short of the others are still extant.
        let r = parse_newick("((A:2.0000001,B:2):1,(C:1,D:0.9999999):2);").unwrap();
        assert_eq!(ltt(&r).last(), Some(&(3.0000001, 4)));
        // As in an MCC tree: tips 0.003% short (29.9571 vs 29.958) stay extant.
        let m = parse_newick("((A:29.9571,B:29.958):0.0,(C:29.9575,D:29.958):0.0);").unwrap();
        assert_eq!(ltt(&m).last().map(|x| x.1), Some(4));
        // A fossil 1% short of the present still ends its lineage.
        let g = parse_newick("((A:29.66,B:29.958):0.0,(C:29.958,D:29.958):0.0);").unwrap();
        assert_eq!(ltt(&g).last().map(|x| x.1), Some(3));
    }
}
