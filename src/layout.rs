//! Logical tree layouts (ggtree `layout=` options).
//!
//! A layout assigns each visible node a depth `x` (branch-length or
//! cladogram units) and a slot position `y` (tips occupy integer slots,
//! internal nodes the mean of their children). Screen projection
//! (rectangular, circular, ...) happens in `scene`.

use crate::style::ViewState;
use crate::tree::{NodeId, Tree};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum LayoutKind {
    Rectangular,
    Slanted,
    Roundrect,
    Ellipse,
    Dendrogram,
    Circular,
    Fan,
    InwardCircular,
    Radial,
}

impl LayoutKind {
    pub const ALL: [LayoutKind; 9] = [
        LayoutKind::Rectangular,
        LayoutKind::Slanted,
        LayoutKind::Roundrect,
        LayoutKind::Ellipse,
        LayoutKind::Dendrogram,
        LayoutKind::Circular,
        LayoutKind::Fan,
        LayoutKind::InwardCircular,
        LayoutKind::Radial,
    ];

    /// Layouts offered in the Layout menu (the others remain available from
    /// the command line).
    pub const MENU: [LayoutKind; 6] =
        [LayoutKind::Rectangular, LayoutKind::Slanted, LayoutKind::Roundrect, LayoutKind::Circular, LayoutKind::Fan, LayoutKind::Radial];

    pub fn name(self) -> &'static str {
        match self {
            LayoutKind::Rectangular => "Rectangular",
            LayoutKind::Slanted => "Slanted",
            LayoutKind::Roundrect => "Round rect",
            LayoutKind::Ellipse => "Ellipse",
            LayoutKind::Dendrogram => "Dendrogram",
            LayoutKind::Circular => "Circular",
            LayoutKind::Fan => "Fan",
            LayoutKind::InwardCircular => "Inward circular",
            LayoutKind::Radial => "Unrooted (equal angle)",
        }
    }

    pub fn is_polar(self) -> bool {
        matches!(self, LayoutKind::Circular | LayoutKind::Fan | LayoutKind::InwardCircular)
    }
}

/// Side of the plot the tips face in rectangular-type layouts.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum TipSide {
    #[default]
    Right,
    Left,
    Top,
    Bottom,
}

impl TipSide {
    pub const ALL: [TipSide; 4] = [TipSide::Right, TipSide::Left, TipSide::Top, TipSide::Bottom];

    pub fn name(self) -> &'static str {
        match self {
            TipSide::Right => "Right",
            TipSide::Left => "Left",
            TipSide::Top => "Top",
            TipSide::Bottom => "Bottom",
        }
    }

    pub fn vertical(self) -> bool {
        matches!(self, TipSide::Top | TipSide::Bottom)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LayoutOptions {
    pub kind: LayoutKind,
    /// Use branch lengths (false = cladogram, ggtree `branch.length="none"`).
    pub use_lengths: bool,
    /// Gap in degrees for fan layouts (centred at the bottom; 180 = upper half circle).
    pub open_angle: f32,
    /// Rotation in degrees for polar layouts.
    pub rotate: f32,
    /// Reverse the tip order.
    pub flip_y: bool,
    /// Where the tips face in rectangular-type layouts: Left is ggtree's
    /// `scale_x_reverse`; Top/Bottom rotate the tree.
    pub tips: TipSide,
    /// Draw the root edge if the root has a length.
    pub root_edge: bool,
    /// Slanted layout as a V-shaped cladogram: each node sits back from its
    /// tips by half their spread, so branches never cross (ignores lengths).
    pub slanted_cladogram: bool,
}

impl Default for LayoutOptions {
    fn default() -> Self {
        LayoutOptions {
            kind: LayoutKind::Rectangular,
            use_lengths: true,
            open_angle: 180.0,
            rotate: 0.0,
            flip_y: false,
            tips: TipSide::Right,
            root_edge: true,
            slanted_cladogram: true,
        }
    }
}

pub struct Layout {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Reachable and not hidden inside a collapsed clade.
    pub visible: Vec<bool>,
    /// Collapsed node -> max depth of its hidden clade.
    pub collapsed_depth: Vec<Option<f64>>,
    /// Visible leaves (tips and collapsed nodes) in slot order.
    pub leaves: Vec<NodeId>,
    pub min_x: f64,
    pub max_x: f64,
    /// Length of the drawn root edge (0 if none).
    pub root_len: f64,
    /// Unrooted coordinates for the radial layout.
    pub radial: Option<Vec<(f64, f64)>>,
}

impl Layout {
    pub fn n_slots(&self) -> usize {
        self.leaves.len().max(1)
    }

    /// Visible leaf slots below `n`: (min slot, max slot).
    pub fn slot_span(&self, tree: &Tree, n: NodeId) -> (f64, f64) {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        let mut stack = vec![n];
        while let Some(u) = stack.pop() {
            if !self.visible[u] {
                continue;
            }
            if self.collapsed_depth[u].is_some() || tree.is_tip(u) {
                lo = lo.min(self.y[u]);
                hi = hi.max(self.y[u]);
            } else {
                stack.extend(tree.nodes[u].children.iter().copied());
            }
        }
        (lo, hi)
    }

    /// Max depth reached in the (visible or collapsed) clade below `n`.
    pub fn clade_max_x(&self, tree: &Tree, n: NodeId) -> f64 {
        let mut m = self.x[n];
        let mut stack = vec![n];
        while let Some(u) = stack.pop() {
            if !self.visible[u] {
                continue;
            }
            m = m.max(self.collapsed_depth[u].unwrap_or(self.x[u]));
            if self.collapsed_depth[u].is_none() {
                stack.extend(tree.nodes[u].children.iter().copied());
            }
        }
        m
    }
}

pub fn compute(tree: &Tree, view: &ViewState) -> Layout {
    let opts = &view.layout;
    let n = tree.nodes.len();
    let use_len = opts.use_lengths && tree.has_lengths();
    let mut visible = vec![false; n];
    let mut collapsed_depth = vec![None; n];
    let mut leaves = Vec::new();
    let full_depth = tree.depths(use_len);

    // Visit in preorder, skipping the interior of collapsed clades.
    let mut stack = vec![tree.root];
    while let Some(u) = stack.pop() {
        visible[u] = true;
        if tree.is_tip(u) {
            leaves.push(u);
        } else if view.collapsed.contains(&u) && u != tree.root {
            let m = tree.preorder_from(u).into_iter().map(|d| full_depth[d]).fold(full_depth[u], f64::max);
            collapsed_depth[u] = Some(m);
            leaves.push(u);
        } else {
            for &c in tree.nodes[u].children.iter().rev() {
                stack.push(c);
            }
        }
    }

    let mut x = full_depth;
    if !use_len {
        // Cladogram: align all tips at the same depth (ggtree branch.length="none").
        let mut h = vec![0.0f64; n];
        for u in tree.postorder() {
            if !visible[u] || collapsed_depth[u].is_some() || tree.is_tip(u) {
                h[u] = 0.0;
            } else {
                h[u] = tree.nodes[u].children.iter().filter(|&&c| visible[c]).map(|&c| h[c] + 1.0).fold(0.0, f64::max);
            }
        }
        let total = h[tree.root];
        for u in 0..n {
            x[u] = total - h[u];
            if collapsed_depth[u].is_some() {
                collapsed_depth[u] = Some(total);
                x[u] = total - 1.0;
            }
        }
    }

    let mut y = vec![0.0f64; n];
    for (i, &l) in leaves.iter().enumerate() {
        y[l] = i as f64;
    }
    for u in tree.postorder() {
        if visible[u] && collapsed_depth[u].is_none() && !tree.is_tip(u) {
            let ch: Vec<f64> = tree.nodes[u].children.iter().filter(|&&c| visible[c]).map(|&c| y[c]).collect();
            if !ch.is_empty() {
                // Midpoint of the outermost children, as in ggtree.
                let lo = ch.iter().cloned().fold(f64::INFINITY, f64::min);
                let hi = ch.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                y[u] = (lo + hi) / 2.0;
            }
        }
    }

    // V-shaped slanted cladogram: x = (root spread − node spread) / 2 in tip
    // slots, so every child lies on its parent's arms and nothing crosses.
    let v_cladogram = opts.kind == LayoutKind::Slanted && opts.slanted_cladogram;
    if v_cladogram {
        let mut span = vec![(f64::INFINITY, f64::NEG_INFINITY); n];
        for u in tree.postorder() {
            if !visible[u] {
                continue;
            }
            span[u] = if tree.is_tip(u) || collapsed_depth[u].is_some() {
                (y[u], y[u])
            } else {
                tree.nodes[u].children.iter().filter(|&&c| visible[c]).fold((f64::INFINITY, f64::NEG_INFINITY), |a, &c| (a.0.min(span[c].0), a.1.max(span[c].1)))
            };
        }
        let total = (span[tree.root].1 - span[tree.root].0).max(1.0);
        for u in 0..n {
            if visible[u] {
                if collapsed_depth[u].is_some() {
                    x[u] = total / 2.0 - 0.5;
                    collapsed_depth[u] = Some(total / 2.0);
                } else {
                    // Midpoint of the tip span, matching the V geometry.
                    y[u] = (span[u].0 + span[u].1) / 2.0;
                    x[u] = (total - (span[u].1 - span[u].0)) / 2.0;
                }
            }
        }
    }
    let root_len = if opts.root_edge && use_len && !v_cladogram { tree.nodes[tree.root].length.unwrap_or(0.0) } else { 0.0 };
    let mut max_x = 0.0f64;
    for u in 0..n {
        if visible[u] {
            max_x = max_x.max(collapsed_depth[u].unwrap_or(x[u]));
        }
    }
    if max_x <= 0.0 {
        max_x = 1.0;
    }

    let radial = if opts.kind == LayoutKind::Radial { Some(equal_angle(tree, &visible, &collapsed_depth, &x, &leaves)) } else { None };

    Layout { x, y, visible, collapsed_depth, leaves, min_x: -root_len, max_x, root_len, radial }
}

/// Equal-angle unrooted layout (Felsenstein's algorithm as used in ggtree `layout="equal_angle"`).
fn equal_angle(tree: &Tree, visible: &[bool], collapsed: &[Option<f64>], x: &[f64], leaves: &[NodeId]) -> Vec<(f64, f64)> {
    let n = tree.nodes.len();
    let mut count = vec![0usize; n];
    for u in tree.postorder() {
        if !visible[u] {
            continue;
        }
        count[u] = if tree.is_tip(u) || collapsed[u].is_some() {
            1
        } else {
            tree.nodes[u].children.iter().filter(|&&c| visible[c]).map(|&c| count[c]).sum()
        };
    }
    let total = leaves.len().max(1) as f64;
    let mut pos = vec![(0.0, 0.0); n];
    let mut start = vec![0.0f64; n];
    for u in tree.preorder() {
        if !visible[u] || collapsed[u].is_some() {
            continue;
        }
        let mut a = start[u];
        for &c in &tree.nodes[u].children {
            if !visible[c] {
                continue;
            }
            let w = std::f64::consts::TAU * count[c] as f64 / total;
            let mid = a + w / 2.0;
            let len = (x[c] - x[u]).max(0.0);
            pos[c] = (pos[u].0 + len * mid.cos(), pos[u].1 + len * mid.sin());
            start[c] = a;
            a += w;
        }
    }
    pos
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick;

    #[test]
    fn slots_and_collapse() {
        let t = parse_newick("((A:1,B:2):1,(C:1,D:1):1);").unwrap();
        let mut v = ViewState::default();
        let l = compute(&t, &v);
        assert_eq!(l.leaves.len(), 4);
        assert_eq!(l.max_x, 3.0);
        let ab = t.nodes[t.find_label("A").unwrap()].parent.unwrap();
        assert_eq!(l.y[ab], 0.5);
        v.collapsed.insert(ab);
        let l = compute(&t, &v);
        assert_eq!(l.leaves.len(), 3);
        assert_eq!(l.collapsed_depth[ab], Some(3.0));
    }

    #[test]
    fn slanted_v_cladogram() {
        let t = parse_newick("(((A:5,B:1):1,C:9):1,D:1);").unwrap();
        let mut v = ViewState::default();
        v.layout.kind = LayoutKind::Slanted;
        let l = compute(&t, &v);
        let id = |s: &str| t.find_label(s).unwrap();
        let ab = t.nodes[id("A")].parent.unwrap();
        let abc = t.nodes[ab].parent.unwrap();
        // Tips aligned; each child lies on the line from its parent to the parent's outermost tip.
        assert_eq!(l.x[id("A")], l.x[id("D")]);
        for (p, c, edge_tip) in [(t.root, abc, id("A")), (abc, ab, id("A"))] {
            let slope = (l.y[edge_tip] - l.y[p]) / (l.x[edge_tip] - l.x[p]);
            let expect = l.y[p] + slope * (l.x[c] - l.x[p]);
            assert!((l.y[c] - expect).abs() < 1e-9);
        }
    }

    #[test]
    fn cladogram() {
        let t = parse_newick("((A:1,B:2):1,C:5);").unwrap();
        let mut v = ViewState::default();
        v.layout.use_lengths = false;
        let l = compute(&t, &v);
        let a = t.find_label("A").unwrap();
        let c = t.find_label("C").unwrap();
        assert_eq!(l.x[a], l.x[c]);
    }
}
