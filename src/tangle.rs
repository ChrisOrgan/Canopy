//! Tanglegrams and mirror trees: two trees drawn face to face, the left with
//! its tips to the right and the right with its tips to the left, with lines
//! joining tips that share a name (none for a mirror tree). `untangle` rotates
//! nodes of both trees to reduce the crossings of those lines.

use crate::io::data::normalize_name;
use crate::layout::TipSide;
use crate::scene::{self, BuildInput, Prim, Scene, SceneEnv, P};
use crate::style::{Color, ViewState};
use crate::tree::{NodeId, Tree};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// How the two trees and their connecting lines are drawn.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TangleStyle {
    /// Lines between matching tips (off: a mirror tree).
    pub lines: bool,
    pub line_color: Color,
    pub line_width: f32,
    /// Color each line with the left tree's clade color at that tip, where set.
    pub clade_colors: bool,
    /// Space between the two trees' tip labels, in points.
    pub gap: f32,
}

impl Default for TangleStyle {
    fn default() -> Self {
        TangleStyle { lines: true, line_color: Color::rgb(120, 120, 120), line_width: 0.8, clade_colors: true, gap: 70.0 }
    }
}

/// Tips of `left` and `right` with the same name (matched as data rows are).
pub fn matches(left: &Tree, right: &Tree) -> Vec<(NodeId, NodeId)> {
    let by_name: HashMap<String, NodeId> = right.tips().into_iter().map(|t| (normalize_name(right.label(t)), t)).collect();
    left.tips().into_iter().filter_map(|t| by_name.get(&normalize_name(left.label(t))).map(|&r| (t, r))).collect()
}

/// Number of crossing connecting lines with the trees' current tip orders.
pub fn crossings(left: &Tree, right: &Tree) -> usize {
    let pos = |t: &Tree| -> HashMap<NodeId, usize> { t.tips().into_iter().enumerate().map(|(i, n)| (n, i)).collect() };
    let (pl, pr) = (pos(left), pos(right));
    let mut pairs: Vec<(usize, usize)> = matches(left, right).into_iter().map(|(a, b)| (pl[&a], pr[&b])).collect();
    pairs.sort_unstable();
    let mut n = 0;
    for i in 0..pairs.len() {
        for j in i + 1..pairs.len() {
            if pairs[j].1 < pairs[i].1 {
                n += 1;
            }
        }
    }
    n
}

/// Reorder the children of every node of `t` by the mean target position of
/// their matched tips (barycenter heuristic); clades without matches keep their place.
fn reorder(t: &mut Tree, target: &HashMap<NodeId, f64>) {
    let mut key: Vec<Option<(f64, usize)>> = vec![None; t.nodes.len()];
    for u in t.postorder() {
        if t.is_tip(u) {
            key[u] = target.get(&u).map(|&p| (p, 1));
            continue;
        }
        let ch = t.nodes[u].children.clone();
        let (sum, cnt) = ch.iter().filter_map(|&c| key[c]).fold((0.0, 0), |a, (m, k)| (a.0 + m * k as f64, a.1 + k));
        key[u] = if cnt > 0 { Some((sum / cnt as f64, cnt)) } else { None };
        // Unmatched children sit where they were, between their neighbours' keys.
        let mut keyed: Vec<(f64, NodeId)> = Vec::new();
        let mut last = f64::NEG_INFINITY;
        for &c in &ch {
            let k = key[c].map(|x| x.0).unwrap_or(last);
            keyed.push((k, c));
            last = k;
        }
        keyed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        t.nodes[u].children = keyed.into_iter().map(|k| k.1).collect();
    }
}

/// Rotate nodes of both trees to reduce crossing lines: alternately order one
/// tree by the other's tip positions, keeping the best orders found. Returns
/// the crossings before and after.
pub fn untangle(left: &mut Tree, right: &mut Tree) -> (usize, usize) {
    let before = crossings(left, right);
    let mut best = (before, left.clone(), right.clone());
    let pairs = matches(left, right);
    for round in 0..20 {
        let pos = |t: &Tree| -> HashMap<NodeId, f64> { t.tips().into_iter().enumerate().map(|(i, n)| (n, i as f64)).collect() };
        if round % 2 == 0 {
            let pl = pos(left);
            let target: HashMap<NodeId, f64> = pairs.iter().map(|&(a, b)| (b, pl[&a])).collect();
            reorder(right, &target);
        } else {
            let pr = pos(right);
            let target: HashMap<NodeId, f64> = pairs.iter().map(|&(a, b)| (a, pr[&b])).collect();
            reorder(left, &target);
        }
        let c = crossings(left, right);
        if c < best.0 {
            best = (c, left.clone(), right.clone());
        }
        if c == 0 {
            break;
        }
    }
    *left = best.1;
    *right = best.2;
    (before, best.0)
}

/// Move a primitive right by `dx`.
fn shift(p: &mut Prim, dx: f32) {
    let s = |q: &mut P| q[0] += dx;
    match p {
        Prim::Path { pts, .. } | Prim::Poly { pts, .. } => pts.iter_mut().for_each(s),
        Prim::Circle { c, .. } => s(c),
        Prim::Text { pos, .. } => s(pos),
        Prim::Strip { a, b, .. } => a.iter_mut().chain(b.iter_mut()).for_each(s),
        Prim::Shade { a, b, fade, .. } => {
            a.iter_mut().chain(b.iter_mut()).for_each(s);
            match fade {
                scene::Fade::Linear(p0, p1) => {
                    s(p0);
                    s(p1)
                }
                scene::Fade::Radial(c, ..) => s(c),
            }
        }
        Prim::Image { center, .. } => s(center),
    }
}

/// A tree's view for one side: tips facing `side`, rectangular-type layout,
/// no background (the figure has one) and no title (the figure's goes on top).
fn side_view(v: &ViewState, side: TipSide) -> ViewState {
    let mut v = v.clone();
    v.layout.tips = side;
    v.layout.flip_y = false;
    if v.layout.kind.is_polar() || v.layout.kind == crate::layout::LayoutKind::Radial {
        v.layout.kind = crate::layout::LayoutKind::Rectangular;
    }
    v.background = None;
    v.title.clear();
    v
}

/// Where a tip's connecting line starts: just past its label on the inner side,
/// else at the tip itself.
fn anchor(s: &Scene, n: NodeId, inner_right: bool, pad: f32) -> Option<P> {
    let p = s.node_pos.get(n).copied().flatten()?;
    match s.label_boxes.iter().find(|b| b.0 == n) {
        Some((_, b)) => Some([if inner_right { b[2] + pad } else { b[0] - pad }, p[1]]),
        None => Some([p[0] + if inner_right { pad } else { -pad }, p[1]]),
    }
}

/// The two trees face to face in `width` × `height` pixels (points × `pt`),
/// with connecting lines unless `style.lines` is off. The left input's view
/// gives the background and title.
pub fn build(left: &BuildInput, right: &BuildInput, style: &TangleStyle, env: &dyn SceneEnv, width: f32, height: f32, pt: f32) -> Scene {
    let gap = style.gap * pt;
    let title = left.view.title.trim().to_string();
    let title_h = if title.is_empty() { 0.0 } else { left.view.title_size * pt * 2.0 };
    let w = ((width - gap) / 2.0).max(10.0);
    let h = (height - title_h).max(10.0);
    let (lv, rv) = (side_view(left.view, TipSide::Right), side_view(right.view, TipSide::Left));
    let mut ls = scene::build(&BuildInput { view: &lv, ..*left }, env, w, h, pt);
    let mut rs = scene::build(&BuildInput { view: &rv, ..*right }, env, w, h, pt);
    let mut out = Scene { width, height, background: left.view.background, ..Default::default() };
    if let Some(bg) = left.view.background {
        out.prims.push(Prim::Poly { pts: vec![[0.0, 0.0], [width, 0.0], [width, height], [0.0, height]], fill: bg, stroke: None });
    }
    // Shift both halves below the title, the right half past the gap.
    for p in ls.prims.iter_mut() {
        shift_y(p, title_h);
    }
    for p in rs.prims.iter_mut() {
        shift_y(p, title_h);
        shift(p, w + gap);
    }
    let lift = |q: &mut Option<P>, dx: f32| {
        if let Some(p) = q {
            p[0] += dx;
            p[1] += title_h;
        }
    };
    ls.node_pos.iter_mut().for_each(|q| lift(q, 0.0));
    rs.node_pos.iter_mut().for_each(|q| lift(q, w + gap));
    for b in ls.label_boxes.iter_mut() {
        b.1[1] += title_h;
        b.1[3] += title_h;
    }
    for b in rs.label_boxes.iter_mut() {
        b.1 = [b.1[0] + w + gap, b.1[1] + title_h, b.1[2] + w + gap, b.1[3] + title_h];
    }
    out.prims.append(&mut ls.prims);
    out.prims.append(&mut rs.prims);
    if style.lines {
        let colors = left.view.resolved_clade_colors(left.tree);
        let pad = 3.0 * pt;
        for (a, b) in matches(left.tree, right.tree) {
            if let (Some(p), Some(q)) = (anchor(&ls, a, true, pad), anchor(&rs, b, false, pad)) {
                let color = if style.clade_colors { colors[a].unwrap_or(style.line_color) } else { style.line_color };
                out.prims.push(Prim::Path { pts: vec![p, q], color, width: style.line_width * pt, dashed: false });
            }
        }
    }
    if !title.is_empty() {
        let size = left.view.title_size * pt;
        out.prims.push(Prim::Text { pos: [width / 2.0, title_h * 0.5], text: title, size, color: Color::BLACK, angle: 0.0, halign: 0.5, valign: 0.5, italic: false, bold: true });
    }
    out
}

/// Move a primitive down by `dy`.
pub(crate) fn shift_y(p: &mut Prim, dy: f32) {
    if dy == 0.0 {
        return;
    }
    let s = |q: &mut P| q[1] += dy;
    match p {
        Prim::Path { pts, .. } | Prim::Poly { pts, .. } => pts.iter_mut().for_each(s),
        Prim::Circle { c, .. } => s(c),
        Prim::Text { pos, .. } => s(pos),
        Prim::Strip { a, b, .. } => a.iter_mut().chain(b.iter_mut()).for_each(s),
        Prim::Shade { a, b, fade, .. } => {
            a.iter_mut().chain(b.iter_mut()).for_each(s);
            match fade {
                scene::Fade::Linear(p0, p1) => {
                    s(p0);
                    s(p1)
                }
                scene::Fade::Radial(c, ..) => s(c),
            }
        }
        Prim::Image { center, .. } => s(center),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::app::document::smart_layers;
    use crate::io::newick::parse_newick;
    use crate::scene::ApproxEnv;

    #[test]
    fn untangle_removes_crossings() {
        let mut l = parse_newick("((A,B),(C,D));").unwrap();
        let mut r = parse_newick("((D,C),(B,A));").unwrap();
        assert_eq!(crossings(&l, &r), 6);
        let (before, after) = untangle(&mut l, &mut r);
        assert_eq!((before, after), (6, 0));
        // Unmatched tips are fine.
        let mut l = parse_newick("((A,B),(C,X));").unwrap();
        let mut r = parse_newick("((C,Y),(B,A));").unwrap();
        assert_eq!(untangle(&mut l, &mut r).1, 0);
        assert_eq!(matches(&l, &r).len(), 3);
    }

    #[test]
    fn tanglegram_scene() {
        let l = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        let r = parse_newick("(C:1,(B:2,A:2):1);").unwrap();
        let (ll, rl) = (smart_layers(&l, false), smart_layers(&r, false));
        let view = ViewState::default();
        let inp = |t, layers| BuildInput { tree: t, view: &view, layers, overlay: &[] };
        let mut style = TangleStyle { clade_colors: false, ..Default::default() };
        let lines = |s: &Scene| s.prims.iter().filter(|p| matches!(p, Prim::Path { color, .. } if *color == style_color())).count();
        fn style_color() -> Color {
            TangleStyle::default().line_color
        }
        let s = build(&inp(&l, &ll), &inp(&r, &rl), &style, &ApproxEnv, 800.0, 400.0, 1.0);
        assert_eq!(lines(&s), 3);
        // Left labels end before the gap; right ones start after it.
        let texts: Vec<f32> = s.prims.iter().filter_map(|p| if let Prim::Text { pos, text, .. } = p { (text == "A").then_some(pos[0]) } else { None }).collect();
        assert_eq!(texts.len(), 2);
        assert!(texts[0] < 400.0 && texts[1] > 400.0);
        style.lines = false;
        assert_eq!(lines(&build(&inp(&l, &ll), &inp(&r, &rl), &style, &ApproxEnv, 800.0, 400.0, 1.0)), 0);
    }
}
