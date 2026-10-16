//! Pre- and post-order tree traversal, drawn step by step over a tree figure:
//! visited nodes are marked and numbered in visiting order, the current one
//! highlighted (GUI: Tree › Tree traversal…; animated GIF export).

use crate::scene::{Prim, Scene};
use crate::style::Color;
use crate::tree::{NodeId, Tree};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Order {
    /// A node, then its children: the root first.
    Pre,
    /// The children, then the node: the root last.
    Post,
}

impl Order {
    pub fn name(self) -> &'static str {
        match self {
            Order::Pre => "Pre-order",
            Order::Post => "Post-order",
        }
    }

    /// One-line description of the visiting rule.
    pub fn rule(self) -> &'static str {
        match self {
            Order::Pre => "each node before its children (root first)",
            Order::Post => "each node after its children (root last)",
        }
    }
}

/// Color of nodes already visited, and of the node being visited.
pub const DONE: Color = Color { r: 29, g: 158, b: 117, a: 255 };
pub const NOW: Color = Color { r: 239, g: 159, b: 39, a: 255 };

/// The drawn nodes in visiting order (nodes hidden in collapsed clades are left out).
pub fn visits(tree: &Tree, scene: &Scene, order: Order) -> Vec<NodeId> {
    let all = match order {
        Order::Pre => tree.preorder(),
        Order::Post => postorder_left_to_right(tree),
    };
    all.into_iter().filter(|&n| scene.node_pos.get(n).copied().flatten().is_some()).collect()
}

/// Post-order with children taken in drawing order (top to bottom), as taught;
/// `Tree::postorder` takes them the other way round.
fn postorder_left_to_right(tree: &Tree) -> Vec<NodeId> {
    let mut out = Vec::new();
    let mut stack = vec![(tree.root, false)];
    while let Some((n, expanded)) = stack.pop() {
        if expanded {
            out.push(n);
        } else {
            stack.push((n, true));
            stack.extend(tree.nodes[n].children.iter().rev().map(|&c| (c, false)));
        }
    }
    out
}

/// Mark the first `step` visits on `scene`: earlier ones as done, the last as
/// current, each numbered if `numbers`, with a caption in the top-right corner.
pub fn overlay(scene: &mut Scene, visits: &[NodeId], step: usize, order: Order, numbers: bool, pt: f32) {
    let step = step.min(visits.len());
    let ink = Color::rgb(40, 40, 40);
    for (i, &n) in visits.iter().take(step).enumerate() {
        let Some(p) = scene.node_pos.get(n).copied().flatten() else { continue };
        let now = i + 1 == step;
        let r = if now { 6.5 } else { 4.5 } * pt;
        scene.prims.push(Prim::Circle { c: p, r, fill: if now { NOW } else { DONE }, stroke: Some((Color::WHITE, 1.0 * pt)) });
        if numbers {
            let size = 7.5 * pt;
            scene.prims.push(Prim::Text { pos: [p[0] - r * 0.6, p[1] - r - 1.5 * pt], text: (i + 1).to_string(), size, color: ink, angle: 0.0, halign: 1.0, valign: 1.0, italic: false, bold: now });
        }
    }
    let caption = format!("{}: {} — step {} of {}", order.name(), order.rule(), step, visits.len());
    scene.prims.push(Prim::Text { pos: [scene.width - 8.0 * pt, 8.0 * pt], text: caption, size: 8.0 * pt, color: ink, angle: 0.0, halign: 1.0, valign: 0.0, italic: false, bold: false });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::app::document::smart_layers;
    use crate::io::newick::parse_newick;
    use crate::scene::{build, ApproxEnv, BuildInput};
    use crate::style::ViewState;

    #[test]
    fn orders_and_overlay() {
        let t = parse_newick("((A:1,B:1)X:1,(C:1,(D:1,E:1)Z:1)Y:1)R;").unwrap();
        let layers = smart_layers(&t, false);
        let view = ViewState::default();
        let mut sc = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        let names = |v: &[NodeId]| v.iter().map(|&n| t.label(n).to_string()).collect::<Vec<_>>().join("");
        let pre = visits(&t, &sc, Order::Pre);
        let post = visits(&t, &sc, Order::Post);
        assert_eq!(names(&pre), "RXABYCZDE");
        assert_eq!(names(&post), "ABXCDEZYR");
        let before = sc.prims.len();
        overlay(&mut sc, &post, 4, Order::Post, true, 1.0);
        let added = &sc.prims[before..];
        // Four circles (three done, one current), four numbers and a caption.
        let circles: Vec<Color> = added.iter().filter_map(|p| if let Prim::Circle { fill, .. } = p { Some(*fill) } else { None }).collect();
        assert_eq!(circles, [DONE, DONE, DONE, NOW]);
        assert_eq!(added.iter().filter(|p| matches!(p, Prim::Text { .. })).count(), 5);
    }
}
