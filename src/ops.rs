//! Tree editing operations (ggtree/treeio/ape equivalents): rerooting,
//! midpoint rooting, rotation, ladderizing, tip dropping, subtree
//! pruning-and-regrafting, and clade extraction.
//!
//! Operations keep node ids stable; removed nodes are left detached in the arena.

use crate::tree::{Attr, NodeId, Tree};
use anyhow::{bail, Result};

/// Attributes describing the branch *above* a node rather than the node itself.
const EDGE_KEYS: &[&str] = &["support", "bootstrap", "posterior", "prob", "length_95%_HPD", "rate"];

struct EdgeData {
    length: Option<f64>,
    attrs: Vec<(String, Attr)>,
    label: Option<String>,
}

fn take_edge(t: &mut Tree, n: NodeId) -> EdgeData {
    let node = &mut t.nodes[n];
    let mut attrs = Vec::new();
    for k in EDGE_KEYS {
        if let Some(v) = node.attrs.remove(*k) {
            attrs.push((k.to_string(), v));
        }
    }
    // Numeric internal labels are support values: they travel with the edge.
    let label = if !node.children.is_empty() && node.label.as_deref().map(|l| l.trim().parse::<f64>().is_ok()).unwrap_or(false) {
        node.label.take()
    } else {
        None
    };
    EdgeData { length: node.length.take(), attrs, label }
}

fn put_edge(t: &mut Tree, n: NodeId, e: EdgeData) {
    let node = &mut t.nodes[n];
    node.length = e.length;
    for (k, v) in e.attrs {
        node.attrs.insert(k, v);
    }
    if e.label.is_some() && !node.children.is_empty() {
        node.label = e.label;
    }
}

fn replace_child(t: &mut Tree, parent: NodeId, old: NodeId, new: NodeId) {
    if let Some(slot) = t.nodes[parent].children.iter_mut().find(|c| **c == old) {
        *slot = new;
    }
}

/// Remove a node with a single child, joining its two branches.
fn splice_unary(t: &mut Tree, n: NodeId) {
    if t.nodes[n].children.len() != 1 {
        return;
    }
    let c = t.nodes[n].children[0];
    let extra = t.nodes[n].length;
    match t.nodes[n].parent {
        Some(p) => {
            replace_child(t, p, n, c);
            t.nodes[c].parent = Some(p);
            if let Some(e) = extra {
                t.nodes[c].length = Some(t.nodes[c].length.unwrap_or(0.0) + e);
            }
        }
        None => {
            t.root = c;
            t.nodes[c].parent = None;
            // Keep the joined length as a root edge only if the old root had one.
            t.nodes[c].length = extra.map(|e| e + t.nodes[c].length.unwrap_or(0.0));
            if extra.is_none() {
                t.nodes[c].length = None;
            }
        }
    }
    t.nodes[n].children.clear();
    t.nodes[n].parent = None;
}

/// Reroot on the branch above `n`, placing the new root at fraction `frac`
/// (0 = at `n`, 1 = at its parent) of the branch length.
pub fn reroot(t: &mut Tree, n: NodeId, frac: f64) -> Result<NodeId> {
    let Some(p) = t.nodes[n].parent else { bail!("node is already the root") };
    let frac = frac.clamp(0.0, 1.0);
    let old_root = t.root;
    let r = t.add_node(None);
    t.nodes[p].children.retain(|&c| c != n);
    let n_edge = take_edge(t, n);
    let total = n_edge.length;
    // Child-side half keeps the original edge data.
    put_edge(
        t,
        n,
        EdgeData { length: total.map(|l| l * frac), attrs: n_edge.attrs.clone(), label: n_edge.label.clone() },
    );
    let mut carry = EdgeData { length: total.map(|l| l * (1.0 - frac)), attrs: n_edge.attrs, label: n_edge.label };
    t.nodes[r].children = vec![n, p];
    t.nodes[n].parent = Some(r);

    let mut prev = r;
    let mut cur = p;
    loop {
        let old_parent = t.nodes[cur].parent;
        let old_edge = take_edge(t, cur);
        t.nodes[cur].parent = Some(prev);
        t.nodes[cur].children.retain(|&c| c != prev);
        if let Some(q) = old_parent {
            t.nodes[cur].children.push(q);
        }
        put_edge(t, cur, carry);
        carry = old_edge;
        match old_parent {
            Some(q) => {
                prev = cur;
                cur = q;
            }
            None => break,
        }
    }
    t.root = r;
    t.rooted = Some(true);
    if t.nodes[old_root].children.len() == 1 {
        splice_unary(t, old_root);
    }
    Ok(r)
}

/// Reroot so that `outgroup` tips form a clade beside the root.
pub fn reroot_outgroup(t: &mut Tree, outgroup: &[NodeId]) -> Result<NodeId> {
    let mut m = t.mrca(outgroup).ok_or_else(|| anyhow::anyhow!("empty outgroup"))?;
    if m == t.root {
        // Outgroup straddles the root: reroot on an ingroup tip first.
        let set: std::collections::HashSet<_> = outgroup.iter().copied().collect();
        let Some(other) = t.tips().into_iter().find(|x| !set.contains(x)) else { bail!("outgroup contains all tips") };
        reroot(t, other, 0.5)?;
        m = t.mrca(outgroup).unwrap();
        if m == t.root {
            bail!("outgroup is not monophyletic on any rooting");
        }
    }
    reroot(t, m, 0.5)
}

/// Unrooted path distance between all reachable nodes from `src`.
fn distances_from(t: &Tree, src: NodeId) -> Vec<f64> {
    let mut d = vec![f64::NAN; t.nodes.len()];
    d[src] = 0.0;
    let mut stack = vec![src];
    while let Some(u) = stack.pop() {
        let mut nb: Vec<(NodeId, f64)> = t.nodes[u].children.iter().map(|&c| (c, t.nodes[c].length.unwrap_or(0.0))).collect();
        if let Some(p) = t.nodes[u].parent {
            nb.push((p, t.nodes[u].length.unwrap_or(0.0)));
        }
        for (v, l) in nb {
            if d[v].is_nan() {
                d[v] = d[u] + l;
                stack.push(v);
            }
        }
    }
    d
}

/// Root at the midpoint of the longest tip-to-tip path.
pub fn midpoint_root(t: &mut Tree) -> Result<NodeId> {
    let tips = t.tips();
    if tips.len() < 2 {
        bail!("need at least two tips");
    }
    let far = |d: &[f64]| *tips.iter().max_by(|&&a, &&b| d[a].partial_cmp(&d[b]).unwrap()).unwrap();
    let a = far(&distances_from(t, tips[0]));
    let da = distances_from(t, a);
    let b = far(&da);
    let half = da[b] / 2.0;
    // Walk from b towards a: path b -> mrca -> a.
    let m = t.mrca(&[a, b]).unwrap();
    let depths = t.depths(true);
    let mut path_b: Vec<NodeId> = t.path_to_root(b).into_iter().take_while(|&x| x != m).collect();
    let mut path_a: Vec<NodeId> = t.path_to_root(a).into_iter().take_while(|&x| x != m).collect();
    // Distance from b along its path to m.
    let db_m = depths[b] - depths[m];
    if db_m >= half {
        for &x in &path_b {
            let dist_b_to_parent = depths[b] - depths[x] + t.nodes[x].length.unwrap_or(0.0);
            if dist_b_to_parent >= half {
                let l = t.nodes[x].length.unwrap_or(0.0);
                let into = half - (depths[b] - depths[x]);
                let frac = if l > 0.0 { into / l } else { 0.0 };
                return reroot(t, x, frac);
            }
        }
    } else {
        path_a.reverse();
        path_b.clear();
        let mut acc = db_m;
        for &x in &path_a {
            let l = t.nodes[x].length.unwrap_or(0.0);
            if acc + l >= half {
                let frac = if l > 0.0 { 1.0 - (half - acc) / l } else { 0.0 };
                return reroot(t, x, frac);
            }
            acc += l;
        }
    }
    bail!("could not locate midpoint")
}

/// Reverse the order of children of `n` (ggtree `rotate`).
pub fn rotate(t: &mut Tree, n: NodeId) {
    t.nodes[n].children.reverse();
}

/// Swap the positions of two sibling clades (ggtree `flip`).
pub fn flip(t: &mut Tree, a: NodeId, b: NodeId) -> Result<()> {
    let (Some(pa), Some(pb)) = (t.nodes[a].parent, t.nodes[b].parent) else { bail!("cannot flip the root") };
    if pa != pb {
        bail!("flip requires two sister clades");
    }
    let ch = &mut t.nodes[pa].children;
    let ia = ch.iter().position(|&c| c == a).unwrap();
    let ib = ch.iter().position(|&c| c == b).unwrap();
    ch.swap(ia, ib);
    Ok(())
}

/// Sort children by clade size. `right = true` puts the largest clade last (towards the bottom).
pub fn ladderize(t: &mut Tree, n: NodeId, right: bool) {
    let counts = t.tip_counts();
    for x in t.preorder_from(n) {
        let mut ch = t.nodes[x].children.clone();
        ch.sort_by_key(|&c| counts[c]);
        if !right {
            ch.reverse();
        }
        t.nodes[x].children = ch;
    }
}

/// Detach the subtree at `n` from its parent, cleaning up a resulting unary node.
fn prune(t: &mut Tree, n: NodeId) -> Result<()> {
    let Some(p) = t.nodes[n].parent else { bail!("cannot prune the root") };
    t.nodes[p].children.retain(|&c| c != n);
    t.nodes[n].parent = None;
    if t.nodes[p].children.len() == 1 {
        splice_unary(t, p);
    } else if t.nodes[p].children.is_empty() {
        prune(t, p)?;
    }
    Ok(())
}

/// Remove tips (or whole clades) from the tree.
pub fn drop_nodes(t: &mut Tree, nodes: &[NodeId]) -> Result<()> {
    let reach = t.reachable();
    let remaining = t.tips().into_iter().filter(|x| !nodes.iter().any(|&d| t.is_ancestor(d, *x))).count();
    if remaining < 2 {
        bail!("cannot drop: fewer than two tips would remain");
    }
    for &n in nodes {
        if reach[n] && t.reachable()[n] && n != t.root {
            prune(t, n)?;
        }
    }
    Ok(())
}

/// Keep only the given tips.
pub fn keep_tips(t: &mut Tree, keep: &[NodeId]) -> Result<()> {
    let drop: Vec<NodeId> = t.tips().into_iter().filter(|x| !keep.contains(x)).collect();
    drop_nodes(t, &drop)
}

/// Move the subtree at `s` so it becomes sister to `target` (drag-and-drop SPR).
pub fn regraft(t: &mut Tree, s: NodeId, target: NodeId) -> Result<()> {
    if s == t.root {
        bail!("cannot move the root");
    }
    if t.is_ancestor(s, target) {
        bail!("cannot attach a clade inside itself");
    }
    if t.nodes[s].parent == Some(target) {
        bail!("clade is already attached there");
    }
    let snapshot = t.clone();
    prune(t, s)?;
    if !t.reachable()[target] {
        *t = snapshot;
        bail!("target disappeared while pruning");
    }
    let m = t.add_node(None);
    match t.nodes[target].parent {
        Some(p) => {
            replace_child(t, p, target, m);
            t.nodes[m].parent = Some(p);
            let l = t.nodes[target].length;
            t.nodes[m].length = l.map(|x| x / 2.0);
            t.nodes[target].length = l.map(|x| x / 2.0);
        }
        None => {
            t.root = m;
            if t.nodes[target].length.is_none() && t.has_lengths() {
                t.nodes[target].length = Some(0.0);
            }
        }
    }
    t.nodes[m].children = vec![target, s];
    t.nodes[target].parent = Some(m);
    t.nodes[s].parent = Some(m);
    Ok(())
}

/// Collapse internal branches shorter than `min_len` or with support below `min_support` into polytomies.
pub fn collapse_weak(t: &mut Tree, min_len: f64, min_support: Option<f64>) -> usize {
    let mut count = 0;
    for n in t.postorder() {
        if n == t.root || t.is_tip(n) {
            continue;
        }
        let short = t.nodes[n].length.map(|l| l < min_len).unwrap_or(false);
        let weak = match (min_support, t.nodes[n].attrs.get("support").and_then(Attr::as_f64)) {
            (Some(th), Some(s)) => s < th,
            _ => false,
        };
        if short || weak {
            let p = t.nodes[n].parent.unwrap();
            let add = t.nodes[n].length.unwrap_or(0.0);
            let kids = std::mem::take(&mut t.nodes[n].children);
            let pos = t.nodes[p].children.iter().position(|&c| c == n).unwrap();
            t.nodes[p].children.remove(pos);
            for (i, &k) in kids.iter().enumerate() {
                t.nodes[k].parent = Some(p);
                if let Some(l) = t.nodes[k].length.as_mut() {
                    *l += add;
                }
                t.nodes[p].children.insert(pos + i, k);
            }
            t.nodes[n].parent = None;
            count += 1;
        }
    }
    count
}

/// Branch-length transformations used in comparative methods (as in
/// geiger::rescale).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BranchTransform {
    /// λ multiplies the depth of internal nodes; tips keep their depth.
    Lambda,
    /// κ raises every branch length to the power κ.
    Kappa,
    /// δ raises node depths (scaled to the clade height) to the power δ.
    Delta,
    /// Ornstein–Uhlenbeck with strength α (pull toward an optimum).
    OrnsteinUhlenbeck,
}

impl BranchTransform {
    pub const ALL: [BranchTransform; 4] =
        [BranchTransform::Lambda, BranchTransform::Kappa, BranchTransform::Delta, BranchTransform::OrnsteinUhlenbeck];

    pub fn name(self) -> &'static str {
        match self {
            BranchTransform::Lambda => "λ (lambda)",
            BranchTransform::Kappa => "κ (kappa)",
            BranchTransform::Delta => "δ (delta)",
            BranchTransform::OrnsteinUhlenbeck => "OU (α)",
        }
    }

    /// Parameter value that leaves the tree unchanged.
    pub fn identity(self) -> f64 {
        match self {
            BranchTransform::OrnsteinUhlenbeck => 0.0,
            _ => 1.0,
        }
    }

    pub fn description(self) -> &'static str {
        match self {
            BranchTransform::Lambda => "λ scales internal branches while tips keep their distance from the root. λ = 0 gives a star tree, λ = 1 the original.",
            BranchTransform::Kappa => "κ raises each branch length to the power κ. κ = 0 makes every branch equal (speciational change), κ = 1 the original.",
            BranchTransform::Delta => "δ raises node depths to the power δ. δ < 1 lengthens early branches, δ > 1 lengthens late ones; δ = 1 the original.",
            BranchTransform::OrnsteinUhlenbeck => "α is the strength of pull toward an optimum under an Ornstein–Uhlenbeck process. Larger α erodes deep shared history, lengthening tip branches relative to internal ones; α = 0 is Brownian motion (the original). The result is rescaled to the original height.",
        }
    }

    pub fn reference(self) -> &'static str {
        match self {
            BranchTransform::Lambda | BranchTransform::Delta => {
                "Pagel, M. (1999) Inferring the historical patterns of biological evolution. Nature 401: 877–884.\nPagel, M. (1997) Inferring evolutionary processes from phylogenies. Zoologica Scripta 26: 331–348."
            }
            BranchTransform::Kappa => {
                "Pagel, M. (1994) Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters. Proc. R. Soc. Lond. B 255: 37–45.\nPagel, M. (1999) Nature 401: 877–884."
            }
            BranchTransform::OrnsteinUhlenbeck => {
                "Hansen, T. F. (1997) Stabilizing selection and the comparative analysis of adaptation. Evolution 51: 1341–1351.\nButler, M. A. & King, A. A. (2004) Phylogenetic comparative analysis: a modeling approach for adaptive evolution. Am. Nat. 164: 683–695."
            }
        }
    }
}

/// Reference for the implementation these transforms follow.
pub const TRANSFORM_IMPLEMENTATION_REF: &str =
    "Implemented as in geiger: Pennell, M. W. et al. (2014) geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. Bioinformatics 30: 2216–2218.";

/// Transform the branch lengths inside the clade rooted at `root` (pass the
/// tree root for the whole tree). The clade's stem branch is left unchanged.
pub fn transform_branches(t: &mut Tree, root: NodeId, kind: BranchTransform, p: f64) -> Result<()> {
    if !p.is_finite() || p < 0.0 {
        bail!("parameter must be a non-negative number");
    }
    let nodes = t.preorder_from(root);
    if nodes.iter().skip(1).all(|&n| t.nodes[n].length.is_none()) {
        bail!("the tree has no branch lengths to transform");
    }
    // Depths measured from the clade root.
    let mut d: std::collections::HashMap<NodeId, f64> = std::collections::HashMap::new();
    for &n in &nodes {
        let v = if n == root { 0.0 } else { d[&t.nodes[n].parent.unwrap()] + t.nodes[n].length.unwrap_or(0.0).max(0.0) };
        d.insert(n, v);
    }
    if kind == BranchTransform::Kappa {
        for &n in nodes.iter().skip(1) {
            if let Some(l) = t.nodes[n].length.as_mut() {
                if *l > 0.0 {
                    *l = l.powf(p);
                }
            }
        }
        return Ok(());
    }
    let h = nodes.iter().map(|n| d[n]).fold(0.0, f64::max);
    if h <= 0.0 {
        bail!("the clade has zero height");
    }
    // OU: shared path t maps to (1/2α)·e^(−2α(H−t))·(1 − e^(−2αt)) (Hansen 1997),
    // rescaled so the clade keeps its height.
    let ou = |x: f64| (-2.0 * p * (h - x)).exp() * (1.0 - (-2.0 * p * x).exp()) / (2.0 * p);
    let ou_h = if p * h > 1e-9 { ou(h) } else { h };
    let new_depth = |n: NodeId, t: &Tree| -> f64 {
        let x = d[&n];
        match kind {
            BranchTransform::Lambda if t.is_tip(n) => x,
            BranchTransform::Lambda => p * x,
            BranchTransform::Delta => h * (x / h).powf(p),
            BranchTransform::OrnsteinUhlenbeck if p * h > 1e-9 => ou(x) * h / ou_h,
            _ => x,
        }
    };
    let nd: std::collections::HashMap<NodeId, f64> = nodes.iter().map(|&n| (n, new_depth(n, t))).collect();
    for &n in nodes.iter().skip(1) {
        let par = t.nodes[n].parent.unwrap();
        t.nodes[n].length = Some((nd[&n] - nd[&par]).max(0.0));
    }
    Ok(())
}
/// Labels of the tips below `n`.
pub fn tip_set(t: &Tree, n: NodeId) -> std::collections::BTreeSet<String> {
    t.tips_below(n).into_iter().map(|x| t.label(x).to_string()).collect()
}

/// The node of `to` holding exactly the same tips as node `n` of `from`
/// (None if the clade is not monophyletic in `to`).
pub fn map_node(from: &Tree, n: NodeId, to: &Tree) -> Option<NodeId> {
    if n >= from.nodes.len() || !from.reachable()[n] {
        return None;
    }
    let want = tip_set(from, n);
    let ids: Vec<NodeId> = want.iter().filter_map(|l| to.tips().into_iter().find(|&x| to.label(x) == l)).collect();
    if ids.len() != want.len() {
        return None;
    }
    let m = to.mrca(&ids)?;
    (tip_set(to, m) == want).then_some(m)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PasteMode {
    /// Attach on the branch above the target, as its sister.
    Sister,
    /// Add as an extra child of the target (creates a polytomy).
    Child,
    /// Replace the target clade.
    Replace,
}

/// Copy `clip` into `t` at `target`. Returns the id of the pasted clade's root.
pub fn graft(t: &mut Tree, target: NodeId, clip: &Tree, mode: PasteMode) -> Result<NodeId> {
    if mode == PasteMode::Child && t.is_tip(target) {
        bail!("cannot add a child to a tip; paste as sister instead");
    }
    // Copy nodes into the arena.
    let mut map = std::collections::HashMap::new();
    for n in clip.preorder() {
        let src = &clip.nodes[n];
        let id = t.nodes.len();
        let parent = src.parent.map(|p| map[&p]);
        t.nodes.push(crate::tree::Node { parent, children: Vec::new(), label: src.label.clone(), length: src.length, attrs: src.attrs.clone() });
        if let Some(p) = parent {
            t.nodes[p].children.push(id);
        }
        map.insert(n, id);
    }
    let s = map[&clip.root];
    t.nodes[s].parent = None;
    match mode {
        PasteMode::Child => {
            t.nodes[target].children.push(s);
            t.nodes[s].parent = Some(target);
        }
        PasteMode::Replace => {
            t.nodes[s].length = t.nodes[target].length;
            match t.nodes[target].parent {
                Some(p) => {
                    replace_child(t, p, target, s);
                    t.nodes[s].parent = Some(p);
                }
                None => t.root = s,
            }
            t.nodes[target].parent = None;
        }
        PasteMode::Sister => {
            let m = t.add_node(None);
            match t.nodes[target].parent {
                Some(p) => {
                    replace_child(t, p, target, m);
                    t.nodes[m].parent = Some(p);
                    let l = t.nodes[target].length;
                    t.nodes[m].length = l.map(|x| x / 2.0);
                    t.nodes[target].length = l.map(|x| x / 2.0);
                }
                None => {
                    t.root = m;
                    if t.nodes[target].length.is_none() && t.has_lengths() {
                        t.nodes[target].length = Some(0.0);
                    }
                }
            }
            t.nodes[m].children = vec![target, s];
            t.nodes[target].parent = Some(m);
            t.nodes[s].parent = Some(m);
        }
    }
    if t.has_lengths() && t.nodes[s].length.is_none() {
        t.nodes[s].length = Some(0.0);
    }
    Ok(s)
}

/// Multiply all branch lengths by `factor`.
pub fn scale_lengths(t: &mut Tree, factor: f64) {
    for n in t.preorder() {
        if let Some(l) = t.nodes[n].length.as_mut() {
            *l *= factor;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::{parse_newick, write_newick, WriteOptions};

    fn nwk(t: &Tree) -> String {
        write_newick(&t.compact(), &WriteOptions::default())
    }

    #[test]
    fn reroot_preserves_tips_and_length() {
        let mut t = parse_newick("((A:1,B:2):3,(C:4,D:5):6);").unwrap();
        let total: f64 = t.preorder().iter().filter_map(|&n| t.nodes[n].length).sum();
        let c = t.find_label("C").unwrap();
        reroot(&mut t, c, 0.5).unwrap();
        assert_eq!(t.num_tips(), 4);
        let total2: f64 = t.preorder().iter().filter_map(|&n| t.nodes[n].length).sum();
        assert!((total - total2).abs() < 1e-9, "{} vs {}: {}", total, total2, nwk(&t));
        let root_children = &t.nodes[t.root].children;
        assert!(root_children.contains(&c));
        assert_eq!(t.nodes[c].length, Some(2.0));
    }

    #[test]
    fn reroot_moves_support() {
        let mut t = parse_newick("(((A:1,B:1)80:1,C:1)70:1,D:1);").unwrap();
        let d = t.find_label("D").unwrap();
        let a = t.find_label("A").unwrap();
        reroot(&mut t, a, 0.5).unwrap();
        // Clade (C,D) is the same bipartition as the old (A,B) edge, so it carries 80.
        let cd = t.mrca(&[t.find_label("C").unwrap(), d]).unwrap();
        assert_eq!(t.tips_below(cd).len(), 2);
        assert_eq!(t.nodes[cd].attrs.get("support"), Some(&Attr::Num(80.0)));
    }

    #[test]
    fn midpoint() {
        let mut t = parse_newick("((A:1,B:1):1,C:10);").unwrap();
        midpoint_root(&mut t).unwrap();
        let c = t.find_label("C").unwrap();
        assert_eq!(t.nodes[c].parent, Some(t.root));
        // Longest path A..C is 12, so C sits 6 from the midpoint root.
        assert!((t.nodes[c].length.unwrap() - 6.0).abs() < 1e-9, "{}", nwk(&t));
    }

    #[test]
    fn drop_and_regraft() {
        let mut t = parse_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
        let a = t.find_label("A").unwrap();
        drop_nodes(&mut t, &[a]).unwrap();
        assert_eq!(nwk(&t), "(B:2,(C:1,D:1):1);");
        let b = t.find_label("B").unwrap();
        let c = t.find_label("C").unwrap();
        regraft(&mut t, b, c).unwrap();
        assert_eq!(t.num_tips(), 3);
        let bc = t.mrca(&[b, c]).unwrap();
        assert_eq!(t.tips_below(bc).len(), 2);
    }

    fn len(t: &Tree, l: &str) -> f64 {
        t.nodes[t.find_label(l).unwrap()].length.unwrap()
    }

    #[test]
    fn branch_transforms() {
        let base = parse_newick("((A:1,B:1):2,C:3);").unwrap();
        let mut t = base.clone();
        let r = t.root;
        transform_branches(&mut t, r, BranchTransform::Lambda, 0.0).unwrap();
        assert_eq!(len(&t, "A"), 3.0); // star tree, tips keep their depth
        let ab = t.nodes[t.find_label("A").unwrap()].parent.unwrap();
        assert_eq!(t.nodes[ab].length, Some(0.0));

        let mut t = base.clone();
        transform_branches(&mut t, r, BranchTransform::Lambda, 1.0).unwrap();
        assert_eq!(nwk(&t), nwk(&base));

        let mut t = base.clone();
        transform_branches(&mut t, r, BranchTransform::Kappa, 0.0).unwrap();
        assert_eq!(nwk(&t), "((A:1,B:1):1,C:1);");

        let mut t = base.clone();
        transform_branches(&mut t, r, BranchTransform::Delta, 2.0).unwrap();
        // AB node at depth 2 of 3 -> 3 * (2/3)^2 = 4/3; tips stay at 3.
        assert!((t.nodes[ab].length.unwrap() - 4.0 / 3.0).abs() < 1e-9);
        assert!((len(&t, "A") - 5.0 / 3.0).abs() < 1e-9);

        // OU: α = 0 is the original; strong α shortens internal branches, height kept.
        let mut t = base.clone();
        transform_branches(&mut t, r, BranchTransform::OrnsteinUhlenbeck, 0.0).unwrap();
        assert_eq!(nwk(&t), nwk(&base));
        let mut t = base.clone();
        transform_branches(&mut t, r, BranchTransform::OrnsteinUhlenbeck, 2.0).unwrap();
        assert!(t.nodes[ab].length.unwrap() < 0.1);
        assert!((len(&t, "C") - 3.0).abs() < 1e-9);
        assert!((t.nodes[ab].length.unwrap() + len(&t, "A") - 3.0).abs() < 1e-9);
    }

    #[test]
    fn graft_and_map() {
        let mut t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        let other = parse_newick("((X:1,Y:1):3,Z:4);").unwrap();
        let x = other.find_label("X").unwrap();
        let mut clip = other.extract(other.nodes[x].parent.unwrap());
        clip.nodes[clip.root].length = Some(0.5);
        let c = t.find_label("C").unwrap();
        graft(&mut t, c, &clip, PasteMode::Sister).unwrap();
        assert_eq!(nwk(&t), "((A:1,B:1):1,(C:1,(X:1,Y:1):0.5):1);");
        let a = t.find_label("A").unwrap();
        graft(&mut t, a, &clip, PasteMode::Replace).unwrap();
        assert!(t.find_label("A").is_none());
        // Clade (X,Y) of `other` maps to a node of `t`.
        let xy = other.nodes[x].parent.unwrap();
        assert!(map_node(&other, xy, &t).is_some());
        let z = other.find_label("Z").unwrap();
        assert!(map_node(&other, z, &t).is_none());
    }

    #[test]
    fn ladderize_sorts() {
        let mut t = parse_newick("(((A,B),C),D);").unwrap();
        let root = t.root;
        ladderize(&mut t, root, true);
        assert_eq!(nwk(&t), "(D,(C,(A,B)));");
    }
}
