//! Core tree data structure.
//!
//! Nodes live in an arena (`Vec<Node>`) and refer to each other by index.
//! Structural edits never remove nodes from the arena; detached nodes are
//! simply unreachable from the root. This keeps `NodeId`s stable across
//! edits, so selections, highlights and other per-node annotations survive.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type NodeId = usize;

/// A node or branch annotation (e.g. from BEAST/MrBayes `[&...]` comments).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Attr {
    Num(f64),
    Text(String),
    Range(f64, f64),
    List(Vec<Attr>),
}

impl Attr {
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Attr::Num(v) => Some(*v),
            Attr::Text(s) => s.trim().parse().ok(),
            _ => None,
        }
    }

    pub fn as_range(&self) -> Option<(f64, f64)> {
        match self {
            Attr::Range(a, b) => Some((*a, *b)),
            Attr::List(v) if v.len() == 2 => Some((v[0].as_f64()?, v[1].as_f64()?)),
            _ => None,
        }
    }

    pub fn display(&self, digits: usize) -> String {
        match self {
            Attr::Num(v) => format_num(*v, digits),
            Attr::Text(s) => s.clone(),
            Attr::Range(a, b) => format!("[{}, {}]", format_num(*a, digits), format_num(*b, digits)),
            Attr::List(v) => {
                let parts: Vec<String> = v.iter().map(|a| a.display(digits)).collect();
                format!("{{{}}}", parts.join(","))
            }
        }
    }
}

/// Format a number with up to `digits` decimals, trimming trailing zeros.
pub fn format_num(v: f64, digits: usize) -> String {
    if v.is_nan() {
        return "NA".into();
    }
    if v != 0.0 && (v.abs() >= 1e6 || v.abs() < 1e-4) {
        return format!("{:.*e}", digits.min(4), v);
    }
    let s = format!("{:.*}", digits, v);
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Node {
    pub parent: Option<NodeId>,
    pub children: Vec<NodeId>,
    pub label: Option<String>,
    /// Length of the branch leading to this node.
    pub length: Option<f64>,
    pub attrs: BTreeMap<String, Attr>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Tree {
    pub nodes: Vec<Node>,
    pub root: NodeId,
    pub name: Option<String>,
    /// `Some(true)` for `[&R]`, `Some(false)` for `[&U]`, `None` if unspecified.
    pub rooted: Option<bool>,
}

impl Default for Tree {
    fn default() -> Self {
        Self::new()
    }
}

impl Tree {
    pub fn new() -> Self {
        Tree { nodes: vec![Node::default()], root: 0, name: None, rooted: None }
    }

    pub fn add_node(&mut self, parent: Option<NodeId>) -> NodeId {
        let id = self.nodes.len();
        self.nodes.push(Node { parent, ..Default::default() });
        if let Some(p) = parent {
            self.nodes[p].children.push(id);
        }
        id
    }

    pub fn node(&self, id: NodeId) -> &Node {
        &self.nodes[id]
    }

    pub fn is_tip(&self, id: NodeId) -> bool {
        self.nodes[id].children.is_empty()
    }

    pub fn parent(&self, id: NodeId) -> Option<NodeId> {
        self.nodes[id].parent
    }

    pub fn label(&self, id: NodeId) -> &str {
        self.nodes[id].label.as_deref().unwrap_or("")
    }

    /// Nodes reachable from the root, parents before children, children in order.
    pub fn preorder(&self) -> Vec<NodeId> {
        self.preorder_from(self.root)
    }

    pub fn preorder_from(&self, start: NodeId) -> Vec<NodeId> {
        let mut out = Vec::new();
        let mut stack = vec![start];
        while let Some(n) = stack.pop() {
            out.push(n);
            for &c in self.nodes[n].children.iter().rev() {
                stack.push(c);
            }
        }
        out
    }

    /// Children before parents.
    pub fn postorder(&self) -> Vec<NodeId> {
        let mut v = self.preorder();
        v.reverse();
        v
    }

    /// Tips in display (preorder) order.
    pub fn tips(&self) -> Vec<NodeId> {
        self.preorder().into_iter().filter(|&n| self.is_tip(n)).collect()
    }

    pub fn tips_below(&self, id: NodeId) -> Vec<NodeId> {
        self.preorder_from(id).into_iter().filter(|&n| self.is_tip(n)).collect()
    }

    pub fn num_tips(&self) -> usize {
        self.tips().len()
    }

    /// Boolean mask of nodes reachable from the root.
    pub fn reachable(&self) -> Vec<bool> {
        let mut m = vec![false; self.nodes.len()];
        for n in self.preorder() {
            m[n] = true;
        }
        m
    }

    /// True if any non-root reachable node has a branch length.
    pub fn has_lengths(&self) -> bool {
        self.preorder().into_iter().any(|n| n != self.root && self.nodes[n].length.is_some())
    }

    /// Root-to-node distances (indexed by NodeId). If `use_lengths` is false,
    /// or the tree has no lengths, every edge counts as 1.
    pub fn depths(&self, use_lengths: bool) -> Vec<f64> {
        let use_len = use_lengths && self.has_lengths();
        let mut d = vec![0.0; self.nodes.len()];
        for n in self.preorder() {
            if let Some(p) = self.nodes[n].parent {
                let l = if use_len { self.nodes[n].length.unwrap_or(0.0).max(0.0) } else { 1.0 };
                d[n] = d[p] + l;
            }
        }
        d
    }

    /// Node heights: distance from the deepest tip (as in time-calibrated trees).
    pub fn heights(&self) -> Vec<f64> {
        let d = self.depths(true);
        let max = self.preorder().iter().map(|&n| d[n]).fold(0.0, f64::max);
        d.iter().map(|x| max - x).collect()
    }

    pub fn is_ancestor(&self, anc: NodeId, mut n: NodeId) -> bool {
        loop {
            if n == anc {
                return true;
            }
            match self.nodes[n].parent {
                Some(p) => n = p,
                None => return false,
            }
        }
    }

    pub fn path_to_root(&self, mut n: NodeId) -> Vec<NodeId> {
        let mut v = vec![n];
        while let Some(p) = self.nodes[n].parent {
            v.push(p);
            n = p;
        }
        v
    }

    /// Most recent common ancestor of a set of nodes.
    pub fn mrca(&self, ids: &[NodeId]) -> Option<NodeId> {
        let first = *ids.first()?;
        let mut path = self.path_to_root(first);
        for &id in &ids[1..] {
            let other: std::collections::HashSet<NodeId> = self.path_to_root(id).into_iter().collect();
            path.retain(|n| other.contains(n));
        }
        path.first().copied()
    }

    pub fn find_label(&self, label: &str) -> Option<NodeId> {
        self.preorder().into_iter().find(|&n| self.nodes[n].label.as_deref() == Some(label))
    }

    /// Number of tips below each node (0 for unreachable nodes).
    pub fn tip_counts(&self) -> Vec<usize> {
        let mut c = vec![0usize; self.nodes.len()];
        for n in self.postorder() {
            if self.is_tip(n) {
                c[n] = 1;
            } else {
                c[n] = self.nodes[n].children.iter().map(|&k| c[k]).sum();
            }
        }
        c
    }

    /// Copy of the subtree rooted at `start`, with a fresh compact arena.
    pub fn extract(&self, start: NodeId) -> Tree {
        let mut t = Tree { nodes: Vec::new(), root: 0, name: self.name.clone(), rooted: self.rooted };
        let mut map = std::collections::HashMap::new();
        for n in self.preorder_from(start) {
            let src = &self.nodes[n];
            let parent = if n == start { None } else { src.parent.map(|p| map[&p]) };
            let id = t.nodes.len();
            t.nodes.push(Node {
                parent,
                children: Vec::new(),
                label: src.label.clone(),
                length: if n == start { None } else { src.length },
                attrs: src.attrs.clone(),
            });
            if let Some(p) = parent {
                t.nodes[p].children.push(id);
            }
            map.insert(n, id);
        }
        t
    }

    /// Copy of the reachable tree with a fresh compact arena.
    pub fn compact(&self) -> Tree {
        let mut t = self.extract(self.root);
        t.nodes[0].length = self.nodes[self.root].length;
        t
    }

    /// Colless's imbalance index for the clade at `n`: the sum over binary
    /// nodes of |tips(left) − tips(right)|. Returns (index, normalized index,
    /// whether polytomies were skipped). Normalized = I / ((k−1)(k−2)/2) for
    /// k tips, so 0 is perfectly balanced and 1 a fully pectinate (caterpillar) clade.
    pub fn colless(&self, n: NodeId) -> (usize, Option<f64>, bool) {
        let counts = self.tip_counts();
        let mut sum = 0usize;
        let mut polytomy = false;
        for u in self.preorder_from(n) {
            match self.nodes[u].children.as_slice() {
                [a, b] => sum += counts[*a].abs_diff(counts[*b]),
                [] => {}
                _ => polytomy = true,
            }
        }
        let k = counts[n] as f64;
        let norm = (k > 2.0).then(|| sum as f64 / ((k - 1.0) * (k - 2.0) / 2.0));
        (sum, norm, polytomy)
    }

    /// Number of cherries (pairs of sister tips) in the clade at `n`.
    pub fn cherries(&self, n: NodeId) -> usize {
        self.preorder_from(n)
            .into_iter()
            .filter(|&u| {
                let c = &self.nodes[u].children;
                c.len() == 2 && c.iter().all(|&x| self.is_tip(x))
            })
            .count()
    }

    /// All attribute keys present on reachable nodes, sorted.
    pub fn attr_keys(&self) -> Vec<String> {
        let mut keys = std::collections::BTreeSet::new();
        for n in self.preorder() {
            for k in self.nodes[n].attrs.keys() {
                keys.insert(k.clone());
            }
        }
        keys.into_iter().collect()
    }

    /// Value used for an attribute lookup that also understands "label" and "length".
    pub fn value(&self, id: NodeId, key: &str) -> Option<Attr> {
        match key {
            "label" => self.nodes[id].label.clone().map(Attr::Text),
            "branch.length" | "length" => self.nodes[id].length.map(Attr::Num),
            _ => self.nodes[id].attrs.get(key).cloned(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::io::newick::parse_newick;

    #[test]
    fn mrca_and_counts() {
        let t = parse_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
        let a = t.find_label("A").unwrap();
        let b = t.find_label("B").unwrap();
        let c = t.find_label("C").unwrap();
        let ab = t.mrca(&[a, b]).unwrap();
        assert_eq!(t.tips_below(ab).len(), 2);
        assert_eq!(t.mrca(&[a, c]), Some(t.root));
        assert_eq!(t.tip_counts()[t.root], 4);
        let h = t.heights();
        assert!((h[t.root] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn colless_and_cherries() {
        let bal = parse_newick("((A,B),(C,D));").unwrap();
        assert_eq!(bal.colless(bal.root), (0, Some(0.0), false));
        assert_eq!(bal.cherries(bal.root), 2);
        let cat = parse_newick("(((A,B),C),D);").unwrap();
        // |3-1| + |2-1| + |1-1| = 3, the maximum for 4 tips.
        assert_eq!(cat.colless(cat.root), (3, Some(1.0), false));
        assert_eq!(cat.cherries(cat.root), 1);
        let poly = parse_newick("((A,B,C),D);").unwrap();
        assert!(poly.colless(poly.root).2);
    }
}
