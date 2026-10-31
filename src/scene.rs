//! Scene building: turns a tree, its view state and a stack of ggtree-style
//! layers into renderer-independent drawing primitives in pixel space.
//! The same scene feeds the interactive canvas and the PNG/TIFF/SVG exporters,
//! so what you see is what you export.

use crate::layout::{self, Layout, LayoutKind, TipSide};
use crate::style::*;
use crate::tree::{format_num, Attr, NodeId, Tree};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

pub type P = [f32; 2];

#[derive(Clone, Debug)]
pub enum Prim {
    Path { pts: Vec<P>, color: Color, width: f32, dashed: bool },
    /// Convex polygon.
    Poly { pts: Vec<P>, fill: Color, stroke: Option<(Color, f32)> },
    Circle { c: P, r: f32, fill: Color, stroke: Option<(Color, f32)> },
    /// `angle` in radians, clockwise on screen. `halign`/`valign` in 0..1
    /// locate the anchor within the text box (0,0 = top-left).
    Text { pos: P, text: String, size: f32, color: Color, angle: f32, halign: f32, valign: f32, italic: bool, bold: bool },
    /// Band between two polylines of equal length (annular sectors etc.).
    Strip { a: Vec<P>, b: Vec<P>, fill: Color, stroke: Option<(Color, f32)> },
    /// Band between two polylines (as `Strip`) filled with a gradient from `c0` to `c1`.
    Shade { a: Vec<P>, b: Vec<P>, c0: Color, c1: Color, fade: Fade },
    /// A raster image (PhyloPic silhouette) looked up by key at render time.
    Image { key: String, center: P, height: f32, tint: Option<Color> },
}

/// Direction of a `Prim::Shade` gradient: `c0` at the first point (radius), `c1` at the second.
#[derive(Clone, Copy, Debug)]
pub enum Fade {
    Linear(P, P),
    /// Center, then the radii of `c0` and `c1`.
    Radial(P, f32, f32),
}

impl Fade {
    /// Position of a point along the gradient: 0 at `c0`, 1 at `c1` (clamped).
    pub fn t(&self, p: P) -> f32 {
        let t = match *self {
            Fade::Linear(a, b) => {
                let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
                ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy).max(1e-9)
            }
            Fade::Radial(c, r0, r1) if (r1 - r0).abs() > 1e-9 => (((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt() - r0) / (r1 - r0),
            Fade::Radial(..) => 1.0,
        };
        t.clamp(0.0, 1.0)
    }
}

#[derive(Default)]
pub struct Scene {
    pub width: f32,
    pub height: f32,
    pub background: Option<Color>,
    pub prims: Vec<Prim>,
    /// Screen position of every visible node (for hit-testing and selection).
    pub node_pos: Vec<Option<P>>,
    /// Hit-test segments for branches: (child node, a, b).
    pub branches: Vec<(NodeId, P, P)>,
    pub label_boxes: Vec<(NodeId, [f32; 4])>,
    /// Image keys requested by PhyloPic layers.
    pub wanted_images: Vec<String>,
    /// Pixels per branch-length unit.
    pub px_per_unit: f32,
    /// Unit screen direction in which a longer branch moves each visible node
    /// (for dragging a node to change its branch length).
    pub depth_dirs: Vec<Option<P>>,
    /// Time-axis ticks (where grid lines go): depth and displayed value, and
    /// whether the axis counts back from the present (for snapping drags).
    pub axis_ticks: Vec<(f64, f64)>,
    pub axis_before_present: bool,
    /// Value the axis shows at the youngest tip when counting back.
    pub axis_youngest: f64,
}

/// Services the scene builder needs from the renderer.
pub trait SceneEnv {
    fn text_width(&self, text: &str, size: f32, italic: bool) -> f32;
    /// Width / height of a loaded image, if available.
    fn image_aspect(&self, key: &str) -> Option<f32>;
}

/// Character-count text measurement for tests and headless use.
pub struct ApproxEnv;

impl SceneEnv for ApproxEnv {
    fn text_width(&self, text: &str, size: f32, _italic: bool) -> f32 {
        text.chars().count() as f32 * size * 0.55
    }
    fn image_aspect(&self, _key: &str) -> Option<f32> {
        None
    }
}

/// Lookup key for a taxon's silhouette (label with underscores as spaces).
pub fn image_key(label: &str) -> String {
    label.replace('_', " ").split_whitespace().collect::<Vec<_>>().join(" ")
}

pub struct BuildInput<'a> {
    pub tree: &'a Tree,
    pub view: &'a ViewState,
    pub layers: &'a [LayerEntry],
    /// Posterior sample trees for DensiTree layers (empty if none).
    pub overlay: &'a [Tree],
}

/// Up to `max` indices spread evenly over `0..n`.
fn subsample(n: usize, max: usize) -> Vec<usize> {
    let k = max.min(n).max(1);
    if n == 0 {
        return Vec::new();
    }
    (0..k).map(|i| i * n / k).collect()
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Orient {
    Right,
    Left,
    Down,
    Up,
}

enum GeoKind {
    Lin { orient: Orient, w: f32, u0: f32, v0: f32, sp: f32, flip_y: bool },
    Polar { c: P, r0: f32, inward: bool, theta0: f64, step: f64 },
    Free { pos: Vec<P>, dir: Vec<f64>, sp: f32 },
}

struct Geo {
    kind: GeoKind,
    s: f32,
    min_x: f64,
    max_x: f64,
    n: usize,
}

impl Geo {
    fn uv(&self, u: f32, v: f32) -> P {
        match &self.kind {
            GeoKind::Lin { orient, w, .. } => match orient {
                Orient::Right => [u, v],
                Orient::Left => [*w - u, v],
                Orient::Down => [v, u],
                Orient::Up => [v, *w - u],
            },
            _ => [u, v],
        }
    }

    fn u(&self, x: f64) -> f32 {
        match &self.kind {
            GeoKind::Lin { u0, .. } => *u0 + ((x - self.min_x) as f32) * self.s,
            _ => 0.0,
        }
    }

    fn v(&self, slot: f64) -> f32 {
        match &self.kind {
            GeoKind::Lin { v0, sp, flip_y, .. } => {
                let k = if *flip_y { self.n as f64 - 1.0 - slot } else { slot };
                *v0 + (k as f32 + 0.5) * *sp
            }
            _ => 0.0,
        }
    }

    fn r(&self, x: f64) -> f32 {
        match &self.kind {
            GeoKind::Polar { r0, inward, .. } => {
                if *inward {
                    *r0 + ((self.max_x - x) as f32) * self.s
                } else {
                    *r0 + ((x - self.min_x) as f32) * self.s
                }
            }
            _ => 0.0,
        }
    }

    fn theta(&self, slot: f64) -> f64 {
        match &self.kind {
            GeoKind::Polar { theta0, step, .. } => *theta0 + slot * *step,
            _ => 0.0,
        }
    }

    fn polar(&self, r: f32, th: f64) -> P {
        match &self.kind {
            GeoKind::Polar { c, .. } => [c[0] + r * th.cos() as f32, c[1] - r * th.sin() as f32],
            _ => [0.0, 0.0],
        }
    }

    /// Slot spacing in pixels (approximate for polar layouts).
    fn spacing(&self, at_x: f64) -> f32 {
        match &self.kind {
            GeoKind::Lin { sp, .. } => *sp,
            GeoKind::Polar { step, .. } => (self.r(at_x) * *step as f32).abs().max(1.0),
            GeoKind::Free { sp, .. } => *sp,
        }
    }

    /// Unit screen direction of increasing depth at node `n` drawn at `p`.
    fn depth_dir(&self, n: NodeId, p: P) -> P {
        match &self.kind {
            GeoKind::Lin { orient, .. } => match orient {
                Orient::Right => [1.0, 0.0],
                Orient::Left => [-1.0, 0.0],
                Orient::Down => [0.0, 1.0],
                Orient::Up => [0.0, -1.0],
            },
            GeoKind::Polar { c, inward, .. } => {
                let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
                let len = (dx * dx + dy * dy).sqrt();
                let s = if *inward { -1.0 } else { 1.0 };
                if len < 1e-3 { [s, 0.0] } else { [s * dx / len, s * dy / len] }
            }
            GeoKind::Free { dir, .. } => [dir[n].cos() as f32, -(dir[n].sin() as f32)],
        }
    }

    fn at(&self, x: f64, slot: f64, node: NodeId) -> P {
        match &self.kind {
            GeoKind::Lin { .. } => self.uv(self.u(x), self.v(slot)),
            GeoKind::Polar { .. } => self.polar(self.r(x), self.theta(slot)),
            GeoKind::Free { pos, .. } => pos[node],
        }
    }

    /// Outward direction (math angle) at a slot / node.
    fn out_angle(&self, slot: f64, node: NodeId) -> f64 {
        match &self.kind {
            GeoKind::Lin { orient, .. } => match orient {
                Orient::Right => 0.0,
                Orient::Left => PI,
                Orient::Down => -FRAC_PI_2,
                Orient::Up => FRAC_PI_2,
            },
            GeoKind::Polar { inward, .. } => self.theta(slot) + if *inward { PI } else { 0.0 },
            GeoKind::Free { dir, .. } => dir[node],
        }
    }

    /// Point `off` pixels outward from depth `base_x` at a slot/node.
    fn ext(&self, slot: f64, base_x: f64, off: f32, node: NodeId) -> P {
        match &self.kind {
            GeoKind::Lin { .. } => self.uv(self.u(base_x) + off, self.v(slot)),
            GeoKind::Polar { inward, .. } => {
                let d = if *inward { -off } else { off };
                self.polar(self.r(base_x) + d, self.theta(slot))
            }
            GeoKind::Free { pos, dir, .. } => {
                let a = dir[node];
                [pos[node][0] + off * a.cos() as f32, pos[node][1] - off * a.sin() as f32]
            }
        }
    }

    /// Region covering slots [lo, hi] between `off0` and `off1` beyond `base_x`.
    fn cell(&self, lo: f64, hi: f64, base_x: f64, off0: f32, off1: f32, node: NodeId) -> Region {
        match &self.kind {
            GeoKind::Lin { sp, .. } => {
                let (a, b) = (self.u(base_x) + off0, self.u(base_x) + off1);
                let (v1, v2) = (self.v(lo) - sp / 2.0, self.v(hi) + sp / 2.0);
                Region::Poly(vec![self.uv(a, v1), self.uv(b, v1), self.uv(b, v2), self.uv(a, v2)])
            }
            GeoKind::Polar { inward, .. } => {
                let d = if *inward { -1.0 } else { 1.0 };
                let r1 = self.r(base_x) + d * off0;
                let r2 = self.r(base_x) + d * off1;
                self.sector(r1, r2, self.theta(lo - 0.5), self.theta(hi + 0.5))
            }
            GeoKind::Free { pos, dir, sp } => {
                let a = dir[node];
                let (ux, uy) = (a.cos() as f32, -a.sin() as f32);
                let (px, py) = (-uy, ux);
                let hh = sp / 2.0 * ((hi - lo) as f32 + 1.0);
                let p = pos[node];
                let q = |o: f32, h: f32| [p[0] + ux * o + px * h, p[1] + uy * o + py * h];
                Region::Poly(vec![q(off0, -hh), q(off1, -hh), q(off1, hh), q(off0, hh)])
            }
        }
    }

    fn sector(&self, r1: f32, r2: f32, t1: f64, t2: f64) -> Region {
        Region::Strip(self.arc(r1, t1, t2), self.arc(r2, t1, t2))
    }

    fn arc(&self, r: f32, t1: f64, t2: f64) -> Vec<P> {
        let steps = (((t2 - t1).abs() / 0.02).ceil() as usize).max(1);
        (0..=steps).map(|i| self.polar(r, t1 + (t2 - t1) * i as f64 / steps as f64)).collect()
    }
}

enum Region {
    Poly(Vec<P>),
    Strip(Vec<P>, Vec<P>),
}

impl Region {
    fn prim(self, fill: Color, stroke: Option<(Color, f32)>) -> Prim {
        match self {
            Region::Poly(pts) => Prim::Poly { pts, fill, stroke },
            Region::Strip(a, b) => Prim::Strip { a, b, fill, stroke },
        }
    }
}

/// Text rotation that keeps text upright while pointing along math angle `a`.
/// Returns (screen angle, flipped) where flipped text should be right-anchored.
fn text_rot(a: f64) -> (f32, bool) {
    let a = a.rem_euclid(TAU);
    if a.cos() < -1e-6 {
        ((-a + PI) as f32, true)
    } else {
        (-a as f32, false)
    }
}

fn nice_step(range: f64, target: usize) -> f64 {
    if range <= 0.0 || !range.is_finite() {
        return 1.0;
    }
    let raw = range / target as f64;
    let mag = 10f64.powf(raw.log10().floor());
    let f = raw / mag;
    let nice = if f < 1.5 {
        1.0
    } else if f < 3.0 {
        2.0
    } else if f < 7.0 {
        5.0
    } else {
        10.0
    };
    nice * mag
}

fn convex_hull(mut pts: Vec<P>) -> Vec<P> {
    pts.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap().then(a[1].partial_cmp(&b[1]).unwrap()));
    pts.dedup();
    if pts.len() < 3 {
        return pts;
    }
    let cross = |o: P, a: P, b: P| (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    let mut lower: Vec<P> = Vec::new();
    for &p in &pts {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 0.0 {
            lower.pop();
        }
        lower.push(p);
    }
    let mut upper: Vec<P> = Vec::new();
    for &p in pts.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 0.0 {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

fn shape_prim(shape: Shape, c: P, size: f32, fill: Color) -> Prim {
    let r = size / 2.0;
    let stroke = None;
    match shape {
        Shape::Circle => Prim::Circle { c, r, fill, stroke },
        Shape::Square => Prim::Poly {
            pts: vec![[c[0] - r, c[1] - r], [c[0] + r, c[1] - r], [c[0] + r, c[1] + r], [c[0] - r, c[1] + r]],
            fill,
            stroke,
        },
        Shape::Triangle => Prim::Poly {
            pts: vec![[c[0], c[1] - r * 1.15], [c[0] + r * 1.1, c[1] + r * 0.8], [c[0] - r * 1.1, c[1] + r * 0.8]],
            fill,
            stroke,
        },
        Shape::Diamond => Prim::Poly {
            pts: vec![[c[0], c[1] - r * 1.2], [c[0] + r, c[1]], [c[0], c[1] + r * 1.2], [c[0] - r, c[1]]],
            fill,
            stroke,
        },
    }
}

struct Ctx<'a> {
    tree: &'a Tree,
    lay: Layout,
    geo: Geo,
    pt: f32,
    env: &'a dyn SceneEnv,
    scene: Scene,
    clade_colors: Vec<Option<Color>>,
    /// BEAST `height` of the youngest tip (above 0 with fossil tips): heights are
    /// measured back from it.
    height0: f64,
}

impl<'a> Ctx<'a> {
    #[allow(clippy::too_many_arguments)]
    fn text(&mut self, pos: P, text: String, size: f32, color: Color, angle: f32, halign: f32, valign: f32, italic: bool, bold: bool) {
        self.scene.prims.push(Prim::Text { pos, text, size, color, angle, halign, valign, italic, bold });
    }

    fn leaf_x(&self, n: NodeId) -> f64 {
        self.lay.collapsed_depth[n].unwrap_or(self.lay.x[n])
    }

    fn pos(&self, n: NodeId) -> P {
        self.geo.at(self.lay.x[n], self.lay.y[n], n)
    }
}

/// The text before the first underscore or space ("Homo_sapiens" → "Homo").
pub fn first_word(label: &str) -> &str {
    let l = label.trim();
    l.split(|c: char| c == '_' || c.is_whitespace()).find(|w| !w.is_empty()).unwrap_or(l)
}

/// The text of tip `n`'s label as the tip-label style shows it.
fn tip_text(tree: &Tree, lay: &Layout, n: NodeId, s: &TipLabelStyle) -> String {
    let base = if lay.collapsed_depth[n].is_some() {
        match tree.nodes[n].label.as_deref().filter(|l| l.parse::<f64>().is_err()) {
            Some(l) => l.to_string(),
            None => format!("{} tips", tree.tips_below(n).len()),
        }
    } else {
        let l = tree.label(n);
        if s.first_word { first_word(l).to_string() } else { l.to_string() }
    };
    if s.underscores_as_spaces {
        base.replace('_', " ")
    } else {
        base
    }
}

pub fn build(inp: &BuildInput, env: &dyn SceneEnv, width: f32, height: f32, pt: f32) -> Scene {
    let tree = inp.tree;
    let view = inp.view;
    let opts = &view.layout;
    let mut lay = layout::compute(tree, view);
    if opts.flip_y && opts.kind.is_polar() {
        // Reverse tip order runs the tips the other way around the circle.
        let last = (lay.n_slots() - 1) as f64;
        for y in lay.y.iter_mut() {
            *y = last - *y;
        }
        lay.leaves.reverse();
    }
    let layers: Vec<&Layer> = inp.layers.iter().filter(|l| l.enabled).map(|l| &l.layer).collect();
    // BEAST heights (absolute, so above 0 at the youngest of fossil tips) are
    // measured back from the youngest tip.
    let height0 = tree.youngest_tip_age().unwrap_or(0.0);
    // Make room for node-height intervals and fossil ranges that reach back past the root.
    let fossil_frame = time_frame(layers.iter().copied(), tree);
    for l in &layers {
        match l {
            Layer::NodeBars(s) if s.attr.contains("height") => {
                for n in tree.preorder() {
                    if let Some((_, hi)) = tree.value(n, &s.attr).and_then(|v| v.as_range()) {
                        lay.min_x = lay.min_x.min(lay.max_x - (hi - height0));
                    }
                }
            }
            Layer::FossilRanges(_) => {
                let mpu = if fossil_frame.1 > 0.0 { fossil_frame.1 } else { 1.0 };
                for n in tree.preorder() {
                    if let Some((_, oldest)) = crate::fossil::range(tree, n) {
                        lay.min_x = lay.min_x.min(lay.max_x - (oldest - fossil_frame.0) / mpu);
                    }
                }
            }
            _ => {}
        }
    }
    // DensiTree: make room for sample trees older than the displayed one.
    if let Some(d) = layers.iter().find_map(|l| if let Layer::DensiTree(d) = l { Some(d) } else { None }) {
        for i in subsample(inp.overlay.len(), d.max_trees) {
            let t = &inp.overlay[i];
            if t.has_lengths() {
                lay.min_x = lay.min_x.min(lay.max_x - t.heights()[t.root]);
            }
        }
    }
    let pad = 10.0 * pt;

    let tiplab = layers.iter().find_map(|l| if let Layer::TipLabels(s) = l { Some(s.clone()) } else { None });
    let label_widths: Vec<f32> = match &tiplab {
        Some(s) => lay
            .leaves
            .iter()
            .map(|&n| env.text_width(&tip_text(tree, &lay, n, s), s.size * pt, s.italic))
            .collect(),
        None => vec![0.0; lay.leaves.len()],
    };
    let max_label_w = label_widths.iter().cloned().fold(0.0, f32::max);

    // Tip points moved out from the tips push the tip labels out past them.
    let points_room = layers
        .iter()
        .filter_map(|l| if let Layer::TipPoints(p) = l { Some(p) } else { None })
        .filter(|p| p.offset > 0.0)
        .map(|p| (p.offset + p.size / 2.0) * pt)
        .fold(0.0f32, f32::max);

    // ---- Columns beyond the tips (labels, silhouettes, heatmaps, bars, clade labels).
    let mut starts: Vec<f32> = vec![0.0; layers.len()];
    let mut cursor = points_room;
    let mut clade_ext = 0.0f32;
    let mut has_extras = false;
    for (i, l) in layers.iter().enumerate() {
        match l {
            Layer::TipLabels(s) => {
                starts[i] = s.offset * pt + points_room;
                cursor = starts[i] + max_label_w + 4.0 * pt;
            }
            Layer::Phylopic(s) => {
                has_extras = true;
                starts[i] = cursor + s.offset * pt;
                let aspect = lay
                    .leaves
                    .iter()
                    .filter_map(|&n| env.image_aspect(&image_key(tree.label(n))))
                    .fold(1.0f32, f32::max)
                    .min(2.5);
                cursor = starts[i] + s.size * pt * aspect + 4.0 * pt;
            }
            Layer::Heatmap(s) => {
                has_extras = true;
                starts[i] = cursor + s.offset * pt;
                cursor = starts[i] + s.columns.len() as f32 * s.cell_width * pt + 4.0 * pt;
            }
            Layer::Bars(s) => {
                has_extras = true;
                starts[i] = cursor + s.offset * pt;
                // Extra room so the scale's end labels aren't cut off.
                cursor = starts[i] + s.max_width * pt + if s.show_scale { 10.0 } else { 4.0 } * pt;
            }
            Layer::CladeLabel(s) => {
                has_extras = true;
                starts[i] = cursor + s.offset * pt;
                let w = env.text_width(&s.text, s.size * pt, false);
                // Room for a silhouette after the text (square until its image loads).
                let pic = clade_pic_key(s).map(|k| 5.0 * pt + s.phylopic_size * pt * env.image_aspect(&k).unwrap_or(1.0).min(2.5)).unwrap_or(0.0);
                clade_ext = clade_ext.max(starts[i] + s.bar_width * pt + 4.0 * pt + w + pic);
            }
            _ => {}
        }
    }
    let ext_total = cursor.max(clade_ext);

    // ---- Legends.
    let visible_nodes: Vec<NodeId> = tree.preorder().into_iter().filter(|&n| lay.visible[n]).collect();
    let mut legends: Vec<ColorMap> = Vec::new();
    let add_legend = |m: &Option<ColorMap>, legends: &mut Vec<ColorMap>| {
        if let Some(m) = m {
            if !legends.iter().any(|l| l.title() == m.title()) {
                legends.push(m.clone());
            }
        }
    };
    let mut maps: Vec<Option<ColorMap>> = Vec::with_capacity(layers.len());
    for l in &layers {
        let m = match l {
            Layer::Tree(s) => s.color_by.as_ref().and_then(|k| ColorMap::build(tree, k, &visible_nodes, s.palette)).map(|m| m.reversed(s.reverse_palette)),
            Layer::TipLabels(s) => s.color_by.as_ref().and_then(|k| ColorMap::build(tree, k, &lay.leaves, s.palette)).map(|m| m.reversed(s.reverse_palette)),
            Layer::TipPoints(s) => point_map(tree, s, &lay.leaves),
            Layer::NodePoints(s) => point_map(tree, s, &visible_nodes),
            Layer::NodeScalars(s) => node_scalar_map(tree, s, &visible_nodes),
            Layer::Phylopic(s) => phylopic_map(tree, s, &lay.leaves, &visible_nodes),
            Layer::Heatmap(s) => heatmap_map(tree, &lay.leaves, s),
            _ => None,
        };
        add_legend(&m, &mut legends);
        maps.push(m);
    }
    // Size legends, one per trait and mark, below the color legends.
    let internal_points: Vec<NodeId> = visible_nodes.iter().copied().filter(|&n| !tree.is_tip(n) && lay.collapsed_depth[n].is_none()).collect();
    let mut sizes: Vec<SizeKey> = Vec::new();
    for l in &layers {
        if let Some(k) = size_key(tree, env, l, &visible_nodes, &lay.leaves, &internal_points) {
            if !sizes.iter().any(|s| s.title == k.title && s.mark == k.mark) {
                sizes.push(k);
            }
        }
    }
    let legend_w = if legends.is_empty() && sizes.is_empty() {
        0.0
    } else {
        let lp = pt * view.legend_scale.max(0.1);
        let mut w = 0.0f32;
        for k in &sizes {
            w = w.max(env.text_width(&k.title, 9.0 * lp, false));
            for (v, _) in k.rows() {
                w = w.max(size_mark_w(k, lp) + 4.0 * lp + env.text_width(&format_num(v, 3), 8.0 * lp, false));
            }
        }
        for m in &legends {
            w = w.max(env.text_width(m.title(), 9.0 * lp, false));
            match m {
                ColorMap::Discrete { levels, .. } => {
                    for (l, _) in levels.iter().take(30) {
                        w = w.max(14.0 * lp + env.text_width(l, 8.0 * lp, false));
                    }
                }
                ColorMap::Continuous { min, max, .. } => {
                    w = w.max(18.0 * lp + env.text_width(&format_num(*min, 3), 8.0 * lp, false));
                    w = w.max(18.0 * lp + env.text_width(&format_num(*max, 3), 8.0 * lp, false));
                }
            }
        }
        w + 16.0 * lp
    };

    let title_h = if view.title.is_empty() { 0.0 } else { (view.title_size + 6.0) * pt };
    let heat_names_h = layers
        .iter()
        .filter_map(|l| if let Layer::Heatmap(h) = l { Some(h) } else { None })
        .filter(|h| h.show_names)
        .map(|h| h.columns.iter().map(|c| env.text_width(c, h.name_size * pt, false)).fold(0.0, f32::max) + 4.0 * pt)
        .fold(0.0, f32::max);
    let has_scale = layers.iter().any(|l| matches!(l, Layer::ScaleBar(_)));
    // Room below the tree for the time axis: ticks and labels.
    let axis_h = layers
        .iter()
        .find_map(|l| if let Layer::TimeAxis(a) = l { Some(a) } else { None })
        .map(|a| (6.0 + a.size + 3.0).max(30.0) * pt)
        .unwrap_or(0.0);
    let top = pad + title_h;
    let rect_like = matches!(opts.kind, LayoutKind::Rectangular | LayoutKind::Slanted | LayoutKind::Roundrect | LayoutKind::Ellipse);
    // Depth runs down the page: the scale bar goes in a left gutter and the axis on the right.
    let vertical = opts.kind == LayoutKind::Dendrogram || (rect_like && opts.tips.vertical());
    // Geologic timescale rows run alongside the depth axis: below the tree, or to its right when vertical.
    let geo_beside = rect_like || vertical;
    let geo_h = layers
        .iter()
        .find_map(|l| if let Layer::Geoscale(g) = l { Some(g) } else { None })
        .filter(|_| geo_beside)
        .map(|g| [g.eras, g.periods, g.epochs, g.ages].iter().filter(|b| **b).count() as f32 * g.row_height * pt + 6.0 * pt)
        .unwrap_or(0.0);
    // Room past the last tip for a bar chart's scale (it sits beside the tree's
    // own axis and timescale, so it only needs that much room overall).
    let bars_scale = layers.iter().filter_map(|l| if let Layer::Bars(b) = l { b.show_scale.then_some((10.0 + 2.0 * b.label_size) * pt) } else { None }).fold(0.0, f32::max);
    let below = (if has_scale { 22.0 * pt } else { 0.0 } + axis_h + geo_h).max(bars_scale);
    let bottom = if vertical { pad } else { pad + below };

    let n_slots = lay.n_slots();
    let range = (lay.max_x - lay.min_x).max(1e-12);
    let unaligned_labels = tiplab.as_ref().map(|s| !s.align).unwrap_or(false) && !has_extras;
    let label_off = tiplab.as_ref().map(|s| s.offset * pt + points_room).unwrap_or(0.0);
    let fit_scale = |avail: f32| -> f32 {
        let mut s = (avail - ext_total) / range as f32;
        if unaligned_labels {
            s = f32::INFINITY;
            for (i, &n) in lay.leaves.iter().enumerate() {
                let x = (lay.collapsed_depth[n].unwrap_or(lay.x[n]) - lay.min_x) as f32;
                let need = label_off + label_widths[i] + 4.0 * pt;
                if x > 0.0 {
                    s = s.min((avail - need) / x);
                }
            }
            if !s.is_finite() {
                s = avail / range as f32;
            }
        }
        s.max(1e-4)
    };

    let plot_w = (width - legend_w).max(50.0);
    let kind = match opts.kind {
        LayoutKind::Circular | LayoutKind::Fan | LayoutKind::InwardCircular => {
            let cx = plot_w / 2.0;
            let cy = (top + height - bottom) / 2.0;
            let r_avail = (plot_w.min(height - top - bottom) / 2.0 - pad).max(20.0);
            let inward = opts.kind == LayoutKind::InwardCircular;
            // A circular tree with a timescale leaves a 12° opening for the age axis.
            let has_geoscale = layers.iter().any(|l| matches!(l, Layer::Geoscale(_)));
            let open = match opts.kind {
                LayoutKind::Fan => Some(opts.open_angle as f64),
                LayoutKind::Circular if has_geoscale => Some(12.0),
                _ => None,
            };
            let fan = opts.kind == LayoutKind::Fan;
            let (theta0, step) = match open {
                Some(gap) => {
                    let sweep = (360.0 - gap).to_radians();
                    // A fan's opening faces down unless rotated (180° = upper half circle).
                    let face = if fan { 1.5 * PI } else { 0.0 };
                    (opts.rotate as f64 * PI / 180.0 + face + gap.to_radians() / 2.0, if n_slots > 1 { sweep / (n_slots - 1) as f64 } else { 0.0 })
                }
                None => (opts.rotate as f64 * PI / 180.0, TAU / n_slots as f64),
            };
            if fan && !inward {
                // Largest scale whose bounding box (tips, labels, rings and the
                // root) fits the canvas; then centre that box.
                let avail = [plot_w - 2.0 * pad, height - top - bottom - 2.0 * pad];
                let (t_lo, t_hi) = (theta0 - step / 2.0, theta0 + (n_slots as f64 - 0.5) * step);
                let edge_room = if has_geoscale { 18.0 * pt } else { 0.0 };
                let bbox = |s: f32| -> [f32; 4] {
                    let mut b = [0.0f32; 4]; // min x, max x, min y, max y around the centre
                    let mut add = |r: f32, th: f64, off: f32| {
                        let (c, sn) = (th.cos() as f32, -(th.sin() as f32));
                        // `off` pushes clockwise of the ray, outside the fan's first edge.
                        let p = [c * r - sn * off, sn * r + c * off];
                        b = [b[0].min(p[0]), b[1].max(p[0]), b[2].min(p[1]), b[3].max(p[1])];
                    };
                    for (i, &n) in lay.leaves.iter().enumerate() {
                        let th = theta0 + lay.y[n] * step;
                        if unaligned_labels {
                            let x = (lay.collapsed_depth[n].unwrap_or(lay.x[n]) - lay.min_x) as f32;
                            add(x * s + label_off + label_widths[i] + 4.0 * pt, th, 0.0);
                        } else {
                            add(range as f32 * s + ext_total, th, 0.0);
                        }
                    }
                    let outer = range as f32 * s + if unaligned_labels { 0.0 } else { ext_total };
                    for k in 0..=48 {
                        add(outer, t_lo + (t_hi - t_lo) * k as f64 / 48.0, 0.0);
                    }
                    add(range as f32 * s, t_lo, edge_room);
                    add(0.0, t_lo, edge_room);
                    b
                };
                let fits = |s: f32| {
                    let b = bbox(s);
                    b[1] - b[0] <= avail[0] && b[3] - b[2] <= avail[1]
                };
                let (mut lo, mut hi) = (1e-4f32, 1e-4f32);
                while fits(hi * 2.0) && hi < 1e9 {
                    hi *= 2.0;
                }
                hi *= 2.0;
                for _ in 0..40 {
                    let mid = (lo + hi) / 2.0;
                    if fits(mid) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let b = bbox(lo);
                let c = [pad + (avail[0] - (b[1] - b[0])) / 2.0 - b[0], top + pad + (avail[1] - (b[3] - b[2])) / 2.0 - b[2]];
                (GeoKind::Polar { c, r0: 0.0, inward, theta0, step }, lo)
            } else {
                let (r0, s) = if inward {
                    let r0 = ext_total + pad;
                    (r0, ((r_avail - r0) / range as f32).max(1e-4))
                } else {
                    // Fit each tip and its label against the canvas width and height
                    // separately, in the direction the label points.
                    let half_w = plot_w / 2.0 - pad;
                    let half_h = (height - top - bottom) / 2.0 - pad;
                    let mut s = f32::INFINITY;
                    for (i, &n) in lay.leaves.iter().enumerate() {
                        let th = theta0 + lay.y[n] * step;
                        let (c, sn) = (th.cos().abs() as f32, th.sin().abs() as f32);
                        let (x, e) = if unaligned_labels {
                            ((lay.collapsed_depth[n].unwrap_or(lay.x[n]) - lay.min_x) as f32, label_off + label_widths[i] + 4.0 * pt)
                        } else {
                            (range as f32, ext_total)
                        };
                        if x <= 0.0 {
                            continue;
                        }
                        if c > 1e-3 {
                            s = s.min((half_w - e * c) / (x * c));
                        }
                        if sn > 1e-3 {
                            s = s.min((half_h - e * sn) / (x * sn));
                        }
                    }
                    if !s.is_finite() || s <= 0.0 {
                        s = fit_scale(r_avail);
                    }
                    (0.0, s.max(1e-4))
                };
                (GeoKind::Polar { c: [cx, cy], r0, inward, theta0, step }, s)
            }
        }
        LayoutKind::Radial => {
            let raw = lay.radial.clone().unwrap_or_default();
            let mut lo = [f64::INFINITY; 2];
            let mut hi = [f64::NEG_INFINITY; 2];
            for &n in &visible_nodes {
                let p = raw[n];
                let ends = match lay.collapsed_depth[n] {
                    Some(_) => vec![p],
                    None => vec![p],
                };
                for q in ends {
                    lo[0] = lo[0].min(q.0);
                    lo[1] = lo[1].min(q.1);
                    hi[0] = hi[0].max(q.0);
                    hi[1] = hi[1].max(q.1);
                }
            }
            let bw = (hi[0] - lo[0]).max(1e-9) as f32;
            let bh = (hi[1] - lo[1]).max(1e-9) as f32;
            let aw = plot_w - 2.0 * pad - 2.0 * ext_total;
            let ah = height - top - bottom - 2.0 * ext_total;
            let s = (aw / bw).min(ah / bh).max(1e-4);
            let cx = plot_w / 2.0;
            let cy = (top + height - bottom) / 2.0;
            let mx = ((lo[0] + hi[0]) / 2.0) as f32;
            let my = ((lo[1] + hi[1]) / 2.0) as f32;
            let pos: Vec<P> = raw.iter().map(|p| [cx + (p.0 as f32 - mx) * s, cy - (p.1 as f32 - my) * s]).collect();
            let mut dir = vec![0.0f64; tree.nodes.len()];
            for &n in &visible_nodes {
                if let Some(p) = tree.nodes[n].parent {
                    dir[n] = (raw[n].1 - raw[p].1).atan2(raw[n].0 - raw[p].0);
                }
            }
            let sp = (TAU as f32 * (bw.max(bh) * s / 2.0) / n_slots as f32).clamp(2.0, 30.0 * pt);
            (GeoKind::Free { pos, dir, sp }, s)
        }
        _ => {
            let orient = if opts.kind == LayoutKind::Dendrogram {
                Orient::Down
            } else {
                match opts.tips {
                    TipSide::Right => Orient::Right,
                    TipSide::Left => Orient::Left,
                    TipSide::Bottom => Orient::Down,
                    TipSide::Top => Orient::Up,
                }
            };
            let (u0, u_end, v0, v_end) = match orient {
                Orient::Down | Orient::Up => (
                    top + heat_names_h.min(0.0),
                    height - bottom,
                    pad + if has_scale { 22.0 * pt } else { 0.0 },
                    plot_w - pad - (axis_h + geo_h).max(bars_scale * 1.4),
                ),
                _ => (pad, plot_w - pad, top + heat_names_h, height - bottom),
            };
            let s = fit_scale(u_end - u0);
            let sp = ((v_end - v0) / n_slots as f32).max(0.5);
            // `w` mirrors depth: across the plot width (Left) or the vertical span (Up).
            let w = if orient == Orient::Up { u0 + u_end } else { plot_w };
            (GeoKind::Lin { orient, w, u0, v0, sp, flip_y: opts.flip_y }, s)
        }
    };
    let geo = Geo { kind: kind.0, s: kind.1, min_x: lay.min_x, max_x: lay.max_x, n: n_slots };

    let mut node_pos = vec![None; tree.nodes.len()];
    for &n in &visible_nodes {
        node_pos[n] = Some(geo.at(lay.x[n], lay.y[n], n));
    }
    let mut depth_dirs = vec![None; tree.nodes.len()];
    for &n in &visible_nodes {
        depth_dirs[n] = node_pos[n].map(|p| geo.depth_dir(n, p));
    }
    let scene = Scene {
        width,
        height,
        background: view.background,
        node_pos,
        px_per_unit: geo.s,
        depth_dirs,
        ..Default::default()
    };
    let clade_colors = view.resolved_clade_colors(tree);
    let mut cx = Ctx { tree, lay, geo, pt, env, scene, clade_colors, height0 };
    if let Some(bg) = view.background {
        cx.scene.prims.push(Prim::Poly { pts: vec![[0.0, 0.0], [width, 0.0], [width, height], [0.0, height]], fill: bg, stroke: None });
    }

    // Branch-length-based decorations make no sense on a cladogram.
    let lengths_shown = opts.use_lengths && tree.has_lengths() && !(opts.kind == LayoutKind::Slanted && opts.slanted_cladogram);
    // The timescale always sits behind the tree.
    for l in &layers {
        if let Layer::Geoscale(s) = l {
            if lengths_shown {
                draw_geoscale(&mut cx, s, axis_h);
            }
        }
    }
    // Fossil ranges are in Ma: placed with the timescale layer's settings, else the time axis's.
    let fossil_time = time_frame(layers.iter().copied(), tree);
    // Clade-label gradients go behind everything, branches included.
    for (i, l) in layers.iter().enumerate() {
        if let Layer::CladeLabel(s) = l {
            draw_clade_gradient(&mut cx, s, starts[i]);
        }
    }
    for (i, l) in layers.iter().enumerate() {
        let start = starts[i];
        let map = maps[i].as_ref();
        match l {
            Layer::Tree(s) => draw_tree(&mut cx, s, map, opts.kind, &visible_nodes),
            Layer::TipLabels(s) => draw_tip_labels(&mut cx, s, map, start, &label_widths),
            Layer::NodeLabels(s) => draw_node_labels(&mut cx, s, &visible_nodes),
            Layer::BranchLengths(s) => draw_branch_lengths(&mut cx, s, &visible_nodes),
            Layer::NodeScalars(s) => draw_node_scalars(&mut cx, s, map, &visible_nodes),
            Layer::DensiTree(s) => draw_densitree(&mut cx, s, opts.kind, inp.overlay),
            Layer::FossilRanges(s) => draw_fossil_ranges(&mut cx, s, &visible_nodes, fossil_time),
            Layer::Geoscale(_) => {}
            Layer::TipPoints(s) => {
                let leaves = cx.lay.leaves.clone();
                draw_points(&mut cx, s, map, &leaves, true)
            }
            Layer::NodePoints(s) => {
                let internal: Vec<NodeId> = visible_nodes
                    .iter()
                    .copied()
                    .filter(|&n| !tree.is_tip(n) && cx.lay.collapsed_depth[n].is_none())
                    .collect();
                draw_points(&mut cx, s, map, &internal, false)
            }
            Layer::NodeBars(s) if lengths_shown => draw_ranges(&mut cx, s, &visible_nodes),
            Layer::NodeBars(_) => {}
            Layer::Highlight(s) => draw_highlight(&mut cx, s),
            Layer::CladeLabel(s) => draw_clade_label(&mut cx, s, start),
            Layer::ScaleBar(s) if lengths_shown => draw_scale_bar(&mut cx, s, bottom, axis_h, geo_h),
            Layer::ScaleBar(_) => {}
            Layer::TimeAxis(s) if lengths_shown => draw_axis(&mut cx, s),
            Layer::TimeAxis(_) => {}
            Layer::Heatmap(s) => draw_heatmap(&mut cx, s, map, start),
            Layer::Bars(s) => draw_bars(&mut cx, s, start),
            Layer::Phylopic(s) => draw_phylopic(&mut cx, s, map, &visible_nodes, start, label_off, &label_widths, tiplab.as_ref().map(|t| t.align).unwrap_or(false)),
        }
    }

    if !view.title.is_empty() {
        cx.text([plot_w / 2.0, pad], view.title.clone(), view.title_size * pt, Color::BLACK, 0.0, 0.5, 0.0, false, true);
    }
    draw_legends(&mut cx, &legends, &sizes, width - legend_w + 6.0 * pt, top, view.legend_scale);
    cx.scene
}

/// Fossil ranges: a wide bar along each branch from the oldest to the youngest
/// age, with ages placed by `time` = (age at the deepest tip, Ma per unit).
fn draw_fossil_ranges(cx: &mut Ctx, s: &FossilRangeStyle, visible: &[NodeId], time: (f64, f64)) {
    if matches!(cx.geo.kind, GeoKind::Free { .. }) {
        return;
    }
    let (young, mpu) = (time.0, if time.1 > 0.0 { time.1 } else { 1.0 });
    let x_of = |age: f64| cx.lay.max_x - (age - young) / mpu;
    for &n in visible {
        if let Some((lo, hi)) = crate::fossil::range(cx.tree, n) {
            let y = cx.lay.y[n];
            let pts = vec![cx.geo.at(x_of(hi), y, n), cx.geo.at(x_of(lo), y, n)];
            cx.scene.prims.push(Prim::Path { pts, color: s.color, width: s.width * cx.pt, dashed: false });
        }
    }
}

/// Nodes shown by a node-scalar layer: internal nodes whose node-scalar
/// posterior is above zero and at least the threshold.
fn node_scalar_nodes(tree: &Tree, s: &NodeScalarStyle, nodes: &[NodeId]) -> Vec<(NodeId, f64)> {
    nodes
        .iter()
        .filter_map(|&n| crate::varrates::node_posterior(tree, n).map(|p| (n, p)))
        .filter(|&(_, p)| p > 0.0 && p >= s.threshold - 1e-12)
        .collect()
}

/// Posterior (0–1) color scale over the node-scalar circles that are drawn.
fn node_scalar_map(tree: &Tree, s: &NodeScalarStyle, nodes: &[NodeId]) -> Option<ColorMap> {
    if s.single_color.is_some() {
        return None;
    }
    let shown = node_scalar_nodes(tree, s, nodes);
    let min = shown.iter().map(|x| x.1).fold(f64::INFINITY, f64::min);
    let max = shown.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max);
    min.is_finite().then(|| ColorMap::Continuous { title: "Node scalar posterior".into(), min, max, palette: s.palette, reversed: s.reverse_palette })
}

fn draw_node_scalars(cx: &mut Ctx, s: &NodeScalarStyle, map: Option<&ColorMap>, visible: &[NodeId]) {
    for (n, p) in node_scalar_nodes(cx.tree, s, visible) {
        let color = s.single_color.or_else(|| map.and_then(|m| m.color(&Attr::Num(p)))).unwrap_or(Color::BLACK);
        // Sized by posterior: from `size_min` at 0 to `size` at 1.
        let size = if s.size_by_pp { s.size_min + (s.size - s.size_min) * p.clamp(0.0, 1.0) as f32 } else { s.size };
        let pos = cx.pos(n);
        cx.scene.prims.push(shape_prim(Shape::Circle, pos, size * cx.pt, color));
    }
}

/// Mark a size legend shows: a branch line, a point shape or a silhouette
/// (image `key`, width / height `aspect`, `tint` as the layer draws it).
#[derive(Clone, Debug, PartialEq)]
enum SizeMark {
    Line,
    Shape(Shape),
    Image { key: String, aspect: f32, tint: Option<Color> },
}

/// Legend for sizes scaled by a numeric trait: values from `lo` to `hi`
/// drawn at sizes (pt) from `small` to `large`.
#[derive(Clone, Debug)]
struct SizeKey {
    title: String,
    lo: f64,
    hi: f64,
    small: f32,
    large: f32,
    mark: SizeMark,
    color: Color,
}

impl SizeKey {
    /// Range of trait `key` over `nodes`; None when no node has a value.
    fn new(tree: &Tree, key: &str, nodes: &[NodeId], (small, large): (f32, f32), mark: SizeMark, color: Color) -> Option<SizeKey> {
        let (lo, hi) = nodes
            .iter()
            .filter_map(|&n| tree.value(n, key).and_then(|v| v.as_f64()).filter(|v| v.is_finite()))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
        lo.is_finite().then(|| SizeKey { title: key.to_string(), lo, hi, small, large, mark, color })
    }

    /// Rows of (value, size in pt): smallest, middle and largest, or one when
    /// every value is the same (drawn at the full size, as the layers do).
    fn rows(&self) -> Vec<(f64, f32)> {
        if self.hi > self.lo {
            (0..3).map(|i| i as f32 / 2.0).map(|t| (self.lo + (self.hi - self.lo) * t as f64, self.small + (self.large - self.small) * t)).collect()
        } else {
            vec![(self.lo, self.large)]
        }
    }
}

/// Size legend of a layer that scales branch widths or points by a trait.
fn size_key(tree: &Tree, env: &dyn SceneEnv, l: &Layer, visible: &[NodeId], leaves: &[NodeId], internal: &[NodeId]) -> Option<SizeKey> {
    // A layer also colored by a trait shows its sizes in grey.
    let ink = |c: Color, by: &Option<String>| if by.is_some() { Color::rgb(110, 110, 110) } else { c };
    match l {
        Layer::Tree(s) => SizeKey::new(tree, s.width_by.as_ref()?, visible, (s.width_min, s.width), SizeMark::Line, ink(s.color, &s.color_by)),
        Layer::TipPoints(s) => SizeKey::new(tree, s.size_by.as_ref()?, leaves, (s.size_min, s.size), SizeMark::Shape(s.shape), ink(s.color, &s.color_by)),
        Layer::NodePoints(s) => SizeKey::new(tree, s.size_by.as_ref()?, internal, (s.size_min, s.size), SizeMark::Shape(s.shape), ink(s.color, &s.color_by)),
        Layer::Phylopic(s) => {
            // The silhouette of the largest value whose image is loaded.
            let key = s.size_by.as_ref()?;
            let nodes = silhouette_nodes(tree, s, leaves, visible);
            let mut by_value: Vec<(f64, NodeId)> = nodes.iter().filter_map(|&n| Some((tree.value(n, key)?.as_f64().filter(|v| v.is_finite())?, n))).collect();
            by_value.sort_by(|a, b| b.0.total_cmp(&a.0));
            let (image, aspect) = by_value.iter().find_map(|&(_, n)| {
                let k = silhouette_key(tree, n)?;
                let a = env.image_aspect(&k)?;
                Some((k, a))
            })?;
            let tint = if s.color_by.is_some() { Some(Color::rgb(110, 110, 110)) } else { s.tint };
            SizeKey::new(tree, key, &nodes, (s.size_min, s.size), SizeMark::Image { key: image, aspect, tint }, Color::BLACK)
        }
        _ => None,
    }
}

/// Color mapping for a point layer, spanning only the points it draws (those
/// passing its filter), under its legend title if it has one.
fn point_map(tree: &Tree, s: &PointStyle, nodes: &[NodeId]) -> Option<ColorMap> {
    let key = s.color_by.as_ref()?;
    let shown: Vec<NodeId> = match &s.filter {
        Some((k, th)) => nodes.iter().copied().filter(|&n| tree.value(n, k).and_then(|v| v.as_f64()).map(|v| v >= *th).unwrap_or(false)).collect(),
        None => nodes.to_vec(),
    };
    let m = ColorMap::build(tree, key, &shown, s.palette)?.reversed(s.reverse_palette);
    Some(match &s.legend_title {
        Some(t) if !t.is_empty() => m.with_title(t),
        _ => m,
    })
}

fn heatmap_map(tree: &Tree, leaves: &[NodeId], s: &HeatmapStyle) -> Option<ColorMap> {
    let mut vals = Vec::new();
    for c in &s.columns {
        for &n in leaves {
            if let Some(v) = tree.value(n, c) {
                vals.push(v);
            }
        }
    }
    if vals.is_empty() {
        return None;
    }
    let title = if s.columns.len() == 1 { s.columns[0].clone() } else { "heatmap".to_string() };
    if vals.iter().all(|v| matches!(v, Attr::Num(_))) {
        let nums: Vec<f64> = vals.iter().filter_map(Attr::as_f64).filter(|x| x.is_finite()).collect();
        let min = nums.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = nums.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        return min.is_finite().then_some(ColorMap::Continuous { title, min, max, palette: s.palette, reversed: s.reverse_palette });
    }
    let mut levels: Vec<String> = vals.iter().map(|v| v.display(4)).collect();
    levels.sort();
    levels.dedup();
    let n = levels.len();
    Some(ColorMap::Discrete { title, levels: levels.into_iter().enumerate().map(|(i, l)| (l, s.palette.discrete(i, n))).collect() }.reversed(s.reverse_palette))
}

fn branch_color(cx: &Ctx, n: NodeId, s: &TreeStyle, map: Option<&ColorMap>) -> Color {
    if let Some(c) = cx.clade_colors[n] {
        return c;
    }
    if let (Some(m), Some(k)) = (map, &s.color_by) {
        if !crate::varrates::branch_passes(cx.tree, n, k, s.min_pp) {
            return s.color;
        }
        // The value cut-off is a BayesTraits branch filter (its slider shows only then).
        if let Some(min) = s.min_value.filter(|_| crate::varrates::branch_pp_key(k).is_some()) {
            if cx.tree.value(n, k).and_then(|v| v.as_f64()).map(|v| if s.min_abs { v.abs() } else { v } < min - 1e-12).unwrap_or(true) {
                return s.color;
            }
        }
        if let Some(c) = cx.tree.value(n, k).and_then(|v| m.color(&v)) {
            return c;
        }
    }
    s.color
}

/// Points of the branch from parent (xp, yp) to child (xc, yc) in the
/// layout's style (rectangular elbow, slanted, curved, polar arc...).
#[allow(clippy::too_many_arguments)]
fn edge_path(g: &Geo, kind: LayoutKind, pt: f32, xp: f64, yp: f64, xc: f64, yc: f64) -> Vec<P> {
    match &g.kind {
        GeoKind::Lin { .. } => {
            let (up, vp, uc, vc) = (g.u(xp), g.v(yp), g.u(xc), g.v(yc));
            match kind {
                LayoutKind::Slanted => vec![g.uv(up, vp), g.uv(uc, vc)],
                LayoutKind::Ellipse => (0..=16)
                    .map(|i| {
                        let t = i as f32 / 16.0 * std::f32::consts::FRAC_PI_2;
                        g.uv(up + (uc - up) * t.sin(), vc + (vp - vc) * t.cos())
                    })
                    .collect(),
                LayoutKind::Roundrect => {
                    let sgn = if vc >= vp { 1.0 } else { -1.0 };
                    let r = (vc - vp).abs().min(uc - up).min(8.0 * pt).max(0.0);
                    let mut v = vec![g.uv(up, vp)];
                    for i in 0..=8 {
                        let t = i as f32 / 8.0 * std::f32::consts::FRAC_PI_2;
                        v.push(g.uv(up + r - r * t.cos(), vc - sgn * r + sgn * r * t.sin()));
                    }
                    v.push(g.uv(uc, vc));
                    v
                }
                _ => vec![g.uv(up, vp), g.uv(up, vc), g.uv(uc, vc)],
            }
        }
        GeoKind::Polar { .. } => {
            let mut v = g.arc(g.r(xp), g.theta(yp), g.theta(yc));
            v.push(g.polar(g.r(xc), g.theta(yc)));
            v
        }
        GeoKind::Free { .. } => Vec::new(),
    }
}

/// Topology key: the sorted set of clades (as sorted tip-label lists).
fn topology_key(t: &Tree) -> String {
    let mut sets: Vec<Vec<String>> = vec![Vec::new(); t.nodes.len()];
    let mut clades: Vec<String> = Vec::new();
    for n in t.postorder() {
        if t.is_tip(n) {
            sets[n] = vec![t.label(n).to_string()];
        } else {
            let mut v: Vec<String> = t.nodes[n].children.iter().flat_map(|&c| sets[c].clone()).collect();
            v.sort();
            clades.push(v.join(","));
            sets[n] = v;
        }
    }
    clades.sort();
    clades.join("|")
}

/// DensiTree: posterior trees drawn translucently over one another, aligned
/// at the tips and using the displayed tree's tip order.
fn draw_densitree(cx: &mut Ctx, s: &DensiTreeStyle, kind: LayoutKind, trees: &[Tree]) {
    if trees.is_empty() || matches!(cx.geo.kind, GeoKind::Free { .. }) {
        return;
    }
    let slot: std::collections::HashMap<&str, f64> =
        cx.lay.leaves.iter().filter(|&&n| cx.tree.is_tip(n)).map(|&n| (cx.tree.label(n), cx.lay.y[n])).collect();
    let picks: Vec<&Tree> = subsample(trees.len(), s.max_trees).into_iter().map(|i| &trees[i]).filter(|t| t.has_lengths()).collect();
    // Rank topologies by frequency among the drawn trees.
    let palette = [Color::hex("#1F5FBF"), Color::hex("#D62728"), Color::hex("#2CA02C")];
    let keys: Vec<String> = if s.by_topology { picks.iter().map(|t| topology_key(t)).collect() } else { Vec::new() };
    let mut ranked: Vec<(&String, usize)> = Vec::new();
    for k in &keys {
        match ranked.iter_mut().find(|(x, _)| *x == k) {
            Some(e) => e.1 += 1,
            None => ranked.push((k, 1)),
        }
    }
    ranked.sort_by(|a, b| b.1.cmp(&a.1));
    let width = s.width * cx.pt;
    let max_x = cx.lay.max_x;
    for (ti, t) in picks.iter().enumerate() {
        let base = if s.by_topology {
            ranked.iter().position(|(k, _)| **k == keys[ti]).and_then(|r| palette.get(r).copied()).unwrap_or(s.color)
        } else {
            s.color
        };
        let color = base.with_alpha(s.alpha.max(1));
        let h = t.heights();
        let mut y = vec![f64::NAN; t.nodes.len()];
        for n in t.postorder() {
            y[n] = if t.is_tip(n) {
                slot.get(t.label(n)).copied().unwrap_or(f64::NAN)
            } else {
                let ys: Vec<f64> = t.nodes[n].children.iter().map(|&c| y[c]).filter(|v| v.is_finite()).collect();
                if ys.is_empty() {
                    f64::NAN
                } else {
                    (ys.iter().cloned().fold(f64::INFINITY, f64::min) + ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max)) / 2.0
                }
            };
        }
        for c in t.preorder() {
            let Some(p) = t.nodes[c].parent else { continue };
            if !(y[p].is_finite() && y[c].is_finite()) {
                continue;
            }
            let pts = edge_path(&cx.geo, kind, cx.pt, max_x - h[p], y[p], max_x - h[c], y[c]);
            cx.scene.prims.push(Prim::Path { pts, color, width, dashed: false });
        }
    }
}

/// Geologic timescale: colored interval bars under the tree (beside it when
/// the tips face up or down) or translucent rings behind it (circular layouts).
fn draw_geoscale(cx: &mut Ctx, s: &GeoscaleStyle, axis_h: f32) {
    use crate::geotime::{overlapping, Level};
    let pt = cx.pt;
    let mpu = if s.ma_per_unit > 0.0 { s.ma_per_unit } else { 1.0 };
    let max_x = cx.lay.max_x;
    let x_of = |age: f64| max_x - (age - s.youngest_age) / mpu;
    let young = s.youngest_age;
    let old = s.youngest_age + (max_x - cx.lay.min_x) * mpu;
    // Finest level first: nearest the tree.
    let rows: Vec<Level> = [(s.ages, Level::Age), (s.epochs, Level::Epoch), (s.periods, Level::Period), (s.eras, Level::Era)].iter().filter(|(on, _)| *on).map(|(_, l)| *l).collect();
    let Some(&finest) = rows.first() else { return };
    match cx.geo.kind {
        GeoKind::Lin { orient, v0, sp, .. } => {
            // Labels run along each interval: rotated to read upward in vertical trees.
            let angle = if matches!(orient, Orient::Down | Orient::Up) { -FRAC_PI_2 as f32 } else { 0.0 };
            let n = cx.geo.n as f32;
            let v_tree_end = v0 + n * sp;
            let mut v = v_tree_end + 6.0 * pt + axis_h;
            let rh = s.row_height * pt;
            if s.boundaries {
                // Boundaries of the finest level shown (ages, then epochs, periods, eras),
                // plus the present when the youngest tip is at 0 Ma.
                let mut lines: Vec<f64> = overlapping(finest, young, old).map(|i| i.start).filter(|a| *a < old && *a > young).collect();
                if young == 0.0 {
                    lines.push(0.0);
                }
                // Boundaries closer than a line width (the Holocene stages near the present)
                // are drawn once, so overlapping translucent lines don't look darker.
                lines.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let mut drawn: Vec<f32> = Vec::new();
                for a in lines {
                    let u = cx.geo.u(x_of(a));
                    if drawn.iter().any(|d| (d - u).abs() < (s.boundary_width * pt).max(1.0) * 1.5) {
                        continue;
                    }
                    drawn.push(u);
                    let pts = vec![cx.geo.uv(u, v0), cx.geo.uv(u, v_tree_end)];
                    cx.scene.prims.push(Prim::Path { pts, color: s.boundary_color, width: s.boundary_width * pt, dashed: false });
                }
            }
            for level in rows {
                for iv in overlapping(level, young, old) {
                    let (a0, a1) = (iv.start.min(old), iv.end.max(young));
                    let (u0, u1) = (cx.geo.u(x_of(a0)), cx.geo.u(x_of(a1)));
                    let pts = vec![cx.geo.uv(u0, v), cx.geo.uv(u1, v), cx.geo.uv(u1, v + rh), cx.geo.uv(u0, v + rh)];
                    cx.scene.prims.push(Prim::Poly { pts, fill: Color::hex(iv.color), stroke: Some((Color::WHITE, 0.5 * pt)) });
                    if s.labels {
                        let w = (u1 - u0).abs() - 2.0 * pt;
                        let size = s.label_size * pt;
                        let text = [iv.name, iv.abbr].into_iter().find(|t| cx.env.text_width(t, size, false) <= w);
                        if let Some(t) = text {
                            let p = cx.geo.uv((u0 + u1) / 2.0, v + rh / 2.0);
                            cx.text(p, t.to_string(), size, Color::BLACK, angle, 0.5, 0.5, false, false);
                        }
                    }
                }
                v += rh;
            }
        }
        GeoKind::Polar { .. } => {
            // Solid concentric rings, one per interval of the finest level
            // selected, behind the whole tree.
            let n = cx.geo.n as f64;
            let (t0, t1) = (cx.geo.theta(-0.5), cx.geo.theta(n - 0.5));
            // A fan (wide opening) gets sectors over its span and the age axis
            // along its first edge; a circle gets full rings and the axis in its
            // 12° opening.
            let fan = t1 - t0 < TAU - 0.3;
            let (ring0, ring1, label_at) = if fan { (t0, t1, (t0 + t1) / 2.0) } else { (0.0, TAU, FRAC_PI_2) };
            let axis = if fan { t0 } else { (t1 + t0 + TAU) / 2.0 };
            let size = s.label_size * pt;
            for iv in overlapping(finest, young, old) {
                let (a0, a1) = (iv.start.min(old), iv.end.max(young));
                let (r_in, r_out) = (cx.geo.r(x_of(a0)), cx.geo.r(x_of(a1)));
                let region = cx.geo.sector(r_in, r_out, ring0, ring1);
                cx.scene.prims.push(region.prim(Color::hex(iv.color), Some((Color::WHITE, 0.6 * pt))));
                // Name along the ring at 12 o'clock, if the ring is thick enough.
                if s.labels && (r_out - r_in).abs() >= size * 1.1 {
                    let mid_r = (r_in + r_out) / 2.0;
                    let room = mid_r * 1.2;
                    if let Some(t) = [iv.name, iv.abbr].into_iter().find(|t| cx.env.text_width(t, size, false) <= room) {
                        let p = cx.geo.polar(mid_r, label_at);
                        cx.text(p, t.to_string(), size, Color::rgb(60, 60, 60), 0.0, 0.5, 0.5, false, false);
                    }
                }
            }
            // Age axis (Ma) along the opening between the first and last tips.
            let c = cx.geo.polar(0.0, axis);
            let (dx, dy) = (axis.cos() as f32, -(axis.sin() as f32));
            let (px, py) = (-dy, dx);
            let at = |r: f32, off: f32| [c[0] + dx * r + px * off, c[1] + dy * r + py * off];
            let (angle, _) = text_rot(axis);
            let r_old = cx.geo.r(x_of(old));
            let r_young = cx.geo.r(x_of(young));
            cx.scene.prims.push(Prim::Path { pts: vec![at(r_old, 0.0), at(r_young, 0.0)], color: Color::BLACK, width: 0.7 * pt, dashed: false });
            let step_age = nice_step(old - young, 5);
            let mut age = (young / step_age).ceil() * step_age;
            while age <= old + 1e-9 {
                let r = cx.geo.r(x_of(age));
                cx.scene.prims.push(Prim::Path { pts: vec![at(r, -2.5 * pt), at(r, 2.5 * pt)], color: Color::BLACK, width: 0.7 * pt, dashed: false });
                cx.text(at(r, 8.0 * pt), format_num(age, 3), size, Color::BLACK, angle, 0.5, 0.5, false, false);
                age += step_age;
            }
            cx.text(at(r_young + 4.0 * pt, 8.0 * pt), "Ma".to_string(), size, Color::BLACK, angle, 0.0, 0.5, false, false);
        }        _ => {}
    }
}

fn draw_tree(cx: &mut Ctx, s: &TreeStyle, map: Option<&ColorMap>, kind: LayoutKind, visible: &[NodeId]) {
    let tree = cx.tree;
    // Width by a trait: linear over the values on the drawn branches.
    let value = |n: NodeId| s.width_by.as_ref().and_then(|k| tree.value(n, k)).and_then(|v| v.as_f64()).filter(|v| v.is_finite());
    let (lo, hi) = visible.iter().filter_map(|&n| value(n)).fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
    let width_of = |n: NodeId| -> f32 {
        let w = match (&s.width_by, value(n)) {
            (None, _) => s.width,
            (Some(_), Some(v)) if hi > lo => s.width_min + (s.width - s.width_min) * ((v - lo) / (hi - lo)) as f32,
            (Some(_), Some(_)) => s.width,
            (Some(_), None) => s.width_min,
        };
        w * cx.pt
    };
    for &c in visible {
        let width = width_of(c);
        let color = branch_color(cx, c, s, map);
        // No value for the color trait (e.g. a data column not reconstructed at
        // internal nodes): grey and dashed, so it can't be read as a low value.
        let missing = s.color_by.as_ref().map(|k| cx.clade_colors[c].is_none() && tree.value(c, k).is_none()).unwrap_or(false);
        let (color, dashed) = if missing { (Color::rgb(160, 160, 160), true) } else { (color, false) };
        let Some(p) = tree.nodes[c].parent else {
            // Root edge.
            if cx.lay.root_len > 0.0 {
                let a = cx.geo.at(-cx.lay.root_len, cx.lay.y[c], c);
                let b = cx.pos(c);
                if !matches!(cx.geo.kind, GeoKind::Free { .. }) {
                    cx.scene.prims.push(Prim::Path { pts: vec![a, b], color, width, dashed: false });
                }
            }
            continue;
        };
        let (xp, yp, xc, yc) = (cx.lay.x[p], cx.lay.y[p], cx.lay.x[c], cx.lay.y[c]);
        let pts: Vec<P> = match &cx.geo.kind {
            GeoKind::Free { pos, .. } => vec![pos[p], pos[c]],
            _ => edge_path(&cx.geo, kind, cx.pt, xp, yp, xc, yc),
        };
        let n = pts.len();
        cx.scene.branches.push((c, pts[n.saturating_sub(2)], pts[n - 1]));
        cx.scene.prims.push(Prim::Path { pts, color, width, dashed });

        if let Some(d) = cx.lay.collapsed_depth[c] {
            let g = &cx.geo;
            let half = 0.4f64;
            let tri: Vec<P> = match &g.kind {
                GeoKind::Free { .. } => {
                    let len = ((d - xc) as f32 * g.s).max(6.0 * cx.pt);
                    let sp = g.spacing(d);
                    let mid = g.ext(yc, xc, len, c);
                    let a = g.out_angle(yc, c);
                    let (px, py) = (a.sin() as f32 * sp, a.cos() as f32 * sp);
                    vec![cx.pos(c), [mid[0] + px, mid[1] + py], [mid[0] - px, mid[1] - py]]
                }
                _ => vec![g.at(xc, yc, c), g.at(d, yc - half, c), g.at(d, yc + half, c)],
            };
            cx.scene.prims.push(Prim::Poly { pts: tri, fill: color.with_alpha(70), stroke: Some((color, width)) });
        }
    }
}

fn draw_tip_labels(cx: &mut Ctx, s: &TipLabelStyle, map: Option<&ColorMap>, start: f32, widths: &[f32]) {
    let size = s.size * cx.pt;
    let leaves = cx.lay.leaves.clone();
    // Only the user's choice aligns labels; extra columns (silhouettes,
    // heatmaps) are laid out after the widest label either way.
    let aligned = s.align;
    for (i, &n) in leaves.iter().enumerate() {
        let text = tip_text(cx.tree, &cx.lay, n, s);
        let x_tip = cx.leaf_x(n);
        let slot = cx.lay.y[n];
        let base = if aligned && !matches!(cx.geo.kind, GeoKind::Free { .. }) { cx.lay.max_x } else { x_tip };
        let pos = cx.geo.ext(slot, base, start, n);
        if aligned && x_tip < cx.lay.max_x - 1e-12 && !matches!(cx.geo.kind, GeoKind::Free { .. }) {
            let a = cx.geo.ext(slot, x_tip, 2.0 * cx.pt, n);
            let b = cx.geo.ext(slot, cx.lay.max_x, start - 2.0 * cx.pt, n);
            cx.scene.prims.push(Prim::Path { pts: vec![a, b], color: Color::GREY, width: 0.5 * cx.pt, dashed: true });
        }
        let mut color = cx.clade_colors[n].unwrap_or(s.color);
        if let (Some(m), Some(k)) = (map, &s.color_by) {
            if let Some(c) = cx.tree.value(n, k).and_then(|v| m.color(&v)) {
                color = c;
            }
        }
        let (angle, flipped) = text_rot(cx.geo.out_angle(slot, n));
        let halign = if flipped { 1.0 } else { 0.0 };
        // Hit box (axis-aligned bounding box of the rotated text box).
        let w = widths[i];
        let (ca, sa) = (angle.cos(), angle.sin());
        let dir = if flipped { -1.0 } else { 1.0 };
        let far = [pos[0] + dir * w * ca, pos[1] + dir * w * sa];
        let h = size * 0.6;
        let bx = [pos[0].min(far[0]) - h, pos[1].min(far[1]) - h, pos[0].max(far[0]) + h, pos[1].max(far[1]) + h];
        cx.scene.label_boxes.push((n, bx));
        cx.text(pos, text, size, color, angle, halign, 0.5, s.italic, false);
    }
}

fn draw_node_labels(cx: &mut Ctx, s: &NodeLabelStyle, visible: &[NodeId]) {
    let size = s.size * cx.pt;
    for &n in visible {
        if cx.tree.is_tip(n) || cx.lay.collapsed_depth[n].is_some() {
            continue;
        }
        let Some(v) = cx.tree.value(n, &s.attr) else { continue };
        if let (Some(min), Some(x)) = (s.min_value, v.as_f64()) {
            if x < min {
                continue;
            }
        }
        let text = match &v {
            Attr::Num(x) => format!("{}{}", format_num(*x, s.digits), s.suffix),
            other => other.display(s.digits),
        };
        let color = match (s.threshold, v.as_f64()) {
            (Some((th, below, above)), Some(x)) => {
                if x < th {
                    below
                } else {
                    above
                }
            }
            _ => s.color,
        };
        let (x, y) = (cx.lay.x[n], cx.lay.y[n]);
        let gap = 2.0 * cx.pt;
        match (&cx.geo.kind, cx.tree.nodes[n].parent) {
            (GeoKind::Lin { orient, .. }, Some(p)) if s.on_branch => {
                let um = (cx.geo.u(cx.lay.x[p]) + cx.geo.u(x)) / 2.0;
                let v = cx.geo.v(y);
                match orient {
                    Orient::Down | Orient::Up => {
                        let pos = cx.geo.uv(um, v - gap);
                        cx.text(pos, text, size, color, 0.0, 1.0, 0.5, false, false)
                    }
                    _ => {
                        let pos = cx.geo.uv(um, v - gap);
                        cx.text(pos, text, size, color, 0.0, 0.5, 1.0, false, false)
                    }
                }
            }
            (GeoKind::Polar { .. }, Some(p)) if s.on_branch => {
                let r = (cx.geo.r(cx.lay.x[p]) + cx.geo.r(x)) / 2.0;
                let pos = cx.geo.polar(r, cx.geo.theta(y));
                cx.text(pos, text, size, color, 0.0, 0.5, 1.0, false, false)
            }
            _ => {
                let p = cx.pos(n);
                cx.text([p[0] - gap, p[1] - gap], text, size, color, 0.0, 1.0, 1.0, false, false)
            }
        }
    }
}

/// Branch-length values (or a per-branch attribute) at the middle of each branch.
fn draw_branch_lengths(cx: &mut Ctx, s: &BranchLengthStyle, visible: &[NodeId]) {
    let size = s.size * cx.pt;
    let gap = 1.5 * cx.pt;
    for &n in visible {
        let Some(p) = cx.tree.nodes[n].parent else { continue };
        let value = match &s.attr {
            Some(a) => cx.tree.value(n, a).and_then(|v| v.as_f64()),
            None => cx.tree.nodes[n].length,
        };
        let Some(value) = value.filter(|v| v.is_finite()) else { continue };
        if s.hide_value.map(|h| (value - h).abs() < 1e-9).unwrap_or(false) {
            continue;
        }
        if !s.attr.as_deref().map(|a| crate::varrates::branch_passes(cx.tree, n, a, s.min_pp)).unwrap_or(true) {
            continue;
        }
        let bt = s.attr.as_deref().and_then(crate::varrates::branch_pp_key).is_some();
        if bt && s.min_value.map(|m| if s.min_abs { value.abs() } else { value } < m - 1e-12).unwrap_or(false) {
            continue;
        }
        let text = format_num(value, s.digits);
        let (xp, xc, y) = (cx.lay.x[p], cx.lay.x[n], cx.lay.y[n]);
        let (pos, angle, halign) = match &cx.geo.kind {
            GeoKind::Lin { orient: Orient::Down | Orient::Up, .. } => {
                let um = (cx.geo.u(xp) + cx.geo.u(xc)) / 2.0;
                (cx.geo.uv(um, cx.geo.v(y) - gap), 0.0, 1.0)
            }
            GeoKind::Lin { .. } => {
                let um = (cx.geo.u(xp) + cx.geo.u(xc)) / 2.0;
                (cx.geo.uv(um, cx.geo.v(y) - gap), 0.0, 0.5)
            }
            GeoKind::Polar { .. } => {
                let th = cx.geo.theta(y);
                let r = (cx.geo.r(xp) + cx.geo.r(xc)) / 2.0;
                (cx.geo.polar(r, th), text_rot(th).0, 0.5)
            }
            GeoKind::Free { pos, .. } => {
                let (a, b) = (pos[p], pos[n]);
                let th = (-(b[1] - a[1]) as f64).atan2((b[0] - a[0]) as f64);
                ([(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0], text_rot(th).0, 0.5)
            }
        };
        cx.text(pos, text, size, s.color, angle, halign, 1.0, false, false);
    }
}

fn draw_points(cx: &mut Ctx, s: &PointStyle, map: Option<&ColorMap>, nodes: &[NodeId], tips: bool) {
    // Size by a numeric trait: linear from `size_min` (smallest) to `size` (largest).
    let size_of = |tree: &Tree, n: NodeId| s.size_by.as_ref().and_then(|k| tree.value(n, k)).and_then(|v| v.as_f64()).filter(|v| v.is_finite());
    let (lo, hi) = nodes.iter().filter_map(|&n| size_of(cx.tree, n)).fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
    for &n in nodes {
        if let Some((k, th)) = &s.filter {
            match cx.tree.value(n, k).and_then(|v| v.as_f64()) {
                Some(v) if v >= *th => {}
                _ => continue,
            }
        }
        let mut color = s.color;
        if let (Some(m), Some(k)) = (map, &s.color_by) {
            match cx.tree.value(n, k).and_then(|v| m.color(&v)) {
                Some(c) => color = c,
                None => continue,
            }
        }
        let p = if tips { cx.geo.ext(cx.lay.y[n], cx.leaf_x(n), s.offset * cx.pt, n) } else { cx.pos(n) };
        let size = match size_of(cx.tree, n) {
            Some(v) if hi > lo => s.size_min + (s.size - s.size_min) * ((v - lo) / (hi - lo)) as f32,
            Some(_) => s.size,
            // Sized by a trait this node lacks: no point, like a missing color.
            None if s.size_by.is_some() => continue,
            None => s.size,
        };
        cx.scene.prims.push(shape_prim(s.shape, p, size * cx.pt, color));
    }
}

/// The range attribute to draw at node `n`. TreeAnnotator `-heights ca` puts
/// each node at its common-ancestor height (`CAheight_mean`), while
/// `height_95%_HPD` covers only the trees containing the clade, so for weakly
/// supported clades it can sit far from the node. When the node's drawn height
/// matches `CAheight_mean` better than `height`, use `CAheight_95%_HPD` instead.
fn range_attr_at<'a>(cx: &Ctx, n: NodeId, attr: &'a str) -> &'a str {
    if attr != "height_95%_HPD" {
        return attr;
    }
    let num = |k: &str| cx.tree.nodes[n].attrs.get(k).and_then(|v| v.as_f64());
    let (Some(ca), Some(mean)) = (num("CAheight_mean"), num("height")) else { return attr };
    if !cx.tree.nodes[n].attrs.contains_key("CAheight_95%_HPD") {
        return attr;
    }
    let drawn = cx.lay.max_x - cx.lay.x[n];
    if (drawn - ca).abs() < (drawn - mean).abs() {
        "CAheight_95%_HPD"
    } else {
        attr
    }
}

fn draw_ranges(cx: &mut Ctx, s: &RangeStyle, visible: &[NodeId]) {
    if matches!(cx.geo.kind, GeoKind::Free { .. }) {
        return;
    }
    let width = s.width * cx.pt;
    for &n in visible {
        let attr = range_attr_at(cx, n, &s.attr);
        let Some((lo, hi)) = cx.tree.value(n, attr).and_then(|v| v.as_range()) else { continue };
        // Heights are measured back from the youngest tip.
        let (x1, x2) = if attr.contains("height") { (cx.lay.max_x - (hi - cx.height0), cx.lay.max_x - (lo - cx.height0)) } else { (lo, hi) };
        let y = cx.lay.y[n];
        let a = cx.geo.at(x1, y, n);
        let b = cx.geo.at(x2, y, n);
        cx.scene.prims.push(Prim::Path { pts: vec![a, b], color: s.color, width, dashed: false });
    }
}

fn draw_highlight(cx: &mut Ctx, s: &HighlightStyle) {
    let n = s.node;
    if n >= cx.tree.nodes.len() || !cx.lay.visible[n] {
        return;
    }
    let (lo, hi) = cx.lay.slot_span(cx.tree, n);
    if !lo.is_finite() {
        return;
    }
    let xmax = cx.lay.clade_max_x(cx.tree, n);
    let x = cx.lay.x[n];
    let back = match cx.tree.nodes[n].parent {
        Some(p) => (((x - cx.lay.x[p]) as f32 * cx.geo.s) / 2.0).min(6.0 * cx.pt),
        None => 4.0 * cx.pt,
    };
    let ext = s.extend * cx.pt;
    let region = match &cx.geo.kind {
        GeoKind::Lin { .. } => {
            let g = &cx.geo;
            let sp = g.spacing(x);
            let (u1, u2) = (g.u(x) - back, g.u(xmax) + ext + 2.0 * cx.pt);
            let (v1, v2) = (g.v(lo) - sp / 2.0, g.v(hi) + sp / 2.0);
            let (v1, v2) = (v1.min(v2), v1.max(v2));
            Region::Poly(vec![g.uv(u1, v1), g.uv(u2, v1), g.uv(u2, v2), g.uv(u1, v2)])
        }
        GeoKind::Polar { inward, .. } => {
            let g = &cx.geo;
            let d = if *inward { -1.0 } else { 1.0 };
            g.sector(g.r(x) - d * back, g.r(xmax) + d * (ext + 2.0 * cx.pt), g.theta(lo - 0.5), g.theta(hi + 0.5))
        }
        GeoKind::Free { pos, .. } => {
            let pad = 6.0 * cx.pt + ext;
            let mut pts = Vec::new();
            for u in cx.tree.preorder_from(n) {
                if cx.lay.visible[u] {
                    let p = pos[u];
                    for (dx, dy) in [(pad, 0.0), (-pad, 0.0), (0.0, pad), (0.0, -pad), (pad * 0.7, pad * 0.7), (-pad * 0.7, pad * 0.7), (pad * 0.7, -pad * 0.7), (-pad * 0.7, -pad * 0.7)] {
                        pts.push([p[0] + dx, p[1] + dy]);
                    }
                }
            }
            Region::Poly(convex_hull(pts))
        }
    };
    cx.scene.prims.push(region.prim(s.fill, None));
}

/// A clade label's gradient: a band behind the clade from just before its node
/// out to the label bar, transparent at the node and `gradient_color` at the bar.
fn draw_clade_gradient(cx: &mut Ctx, s: &CladeLabelStyle, start: f32) {
    let n = s.node;
    if !s.gradient || n >= cx.tree.nodes.len() || !cx.lay.visible[n] {
        return;
    }
    let (lo, hi) = cx.lay.slot_span(cx.tree, n);
    if !lo.is_finite() {
        return;
    }
    let x = cx.lay.x[n];
    let back = match cx.tree.nodes[n].parent {
        Some(p) => (((x - cx.lay.x[p]) as f32 * cx.geo.s) / 2.0).min(6.0 * cx.pt),
        None => 4.0 * cx.pt,
    };
    let max_x = cx.lay.max_x;
    let g = &cx.geo;
    let (a, b, fade) = match &g.kind {
        GeoKind::Lin { .. } => {
            let sp = g.spacing(max_x);
            let (u0, u1) = (g.u(x) - back, g.u(max_x) + start);
            let (v1, v2) = (g.v(lo), g.v(hi));
            let (v1, v2) = (v1.min(v2) - sp * 0.4, v1.max(v2) + sp * 0.4);
            (vec![g.uv(u0, v1), g.uv(u1, v1)], vec![g.uv(u0, v2), g.uv(u1, v2)], Fade::Linear(g.uv(u0, v1), g.uv(u1, v1)))
        }
        GeoKind::Polar { inward, .. } => {
            let d = if *inward { -1.0 } else { 1.0 };
            let (r0, r1) = ((g.r(x) - d * back).max(0.0), g.r(max_x) + d * start);
            let (t1, t2) = (g.theta(lo - 0.4), g.theta(hi + 0.4));
            (g.arc(r0, t1, t2), g.arc(r1, t1, t2), Fade::Radial(g.polar(0.0, 0.0), r0, r1))
        }
        GeoKind::Free { .. } => return,
    };
    cx.scene.prims.push(Prim::Shade { a, b, c0: s.gradient_color.with_alpha(0), c1: s.gradient_color, fade });
}

/// A clade label's silhouette: its key (requested from PhyloPic), width and
/// height in pixels; None without one, or until its image has loaded.
fn clade_pic(cx: &mut Ctx, s: &CladeLabelStyle) -> Option<(String, f32, f32)> {
    let key = clade_pic_key(s)?;
    cx.scene.wanted_images.push(key.clone());
    let aspect = cx.env.image_aspect(&key)?;
    let h = s.phylopic_size * cx.pt;
    Some((key, h * aspect.min(2.5), h))
}

/// Image key of a clade label's silhouette, if it has one.
fn clade_pic_key(s: &CladeLabelStyle) -> Option<String> {
    s.phylopic.as_deref().map(image_key).filter(|k| !k.is_empty())
}

fn draw_clade_label(cx: &mut Ctx, s: &CladeLabelStyle, start: f32) {
    let n = s.node;
    if n >= cx.tree.nodes.len() || !cx.lay.visible[n] {
        return;
    }
    let (lo, hi) = cx.lay.slot_span(cx.tree, n);
    if !lo.is_finite() {
        return;
    }
    let size = s.size * cx.pt;
    let bw = s.bar_width * cx.pt;
    let max_x = cx.lay.max_x;
    match &cx.geo.kind {
        GeoKind::Lin { .. } => {
            let g = &cx.geo;
            let sp = g.spacing(max_x) * 0.35;
            let u = g.u(max_x) + start;
            let (v1, v2) = (g.v(lo), g.v(hi));
            let (v1, v2) = (v1.min(v2) - sp, v1.max(v2) + sp);
            let bar = vec![g.uv(u, v1), g.uv(u, v2)];
            let tp = g.uv(u + bw + 3.0 * cx.pt, (v1 + v2) / 2.0);
            let (angle, flipped) = text_rot(g.out_angle(0.0, n));
            cx.scene.prims.push(Prim::Path { pts: bar, color: s.color, width: bw, dashed: false });
            cx.text(tp, s.text.clone(), size, s.color, angle, if flipped { 1.0 } else { 0.0 }, 0.5, false, false);
            // Silhouette after the text, along the outward direction.
            if let Some((key, w, h)) = clade_pic(cx, s) {
                let g = &cx.geo;
                let foot = if matches!(g.kind, GeoKind::Lin { orient: Orient::Down | Orient::Up, .. }) { h } else { w };
                let tw = cx.env.text_width(&s.text, size, false);
                let center = g.uv(u + bw + 3.0 * cx.pt + tw + 5.0 * cx.pt + foot / 2.0, (v1 + v2) / 2.0);
                cx.scene.prims.push(Prim::Image { key, center, height: h, tint: Some(s.color) });
            }
        }
        GeoKind::Polar { inward, .. } => {
            let g = &cx.geo;
            let d = if *inward { -1.0 } else { 1.0 };
            let r = g.r(max_x) + d * start;
            let bar = g.arc(r, g.theta(lo - 0.35), g.theta(hi + 0.35));
            let mid = (lo + hi) / 2.0;
            let tp = g.polar(r + d * (bw + 3.0 * cx.pt), g.theta(mid));
            let (angle, flipped) = text_rot(g.out_angle(mid, n));
            cx.scene.prims.push(Prim::Path { pts: bar, color: s.color, width: bw, dashed: false });
            cx.text(tp, s.text.clone(), size, s.color, angle, if flipped { 1.0 } else { 0.0 }, 0.5, false, false);
            if let Some((key, w, h)) = clade_pic(cx, s) {
                let g = &cx.geo;
                let a = g.out_angle(mid, n);
                let foot = w * a.cos().abs() as f32 + h * a.sin().abs() as f32;
                let tw = cx.env.text_width(&s.text, size, false);
                let center = g.polar(r + d * (bw + 3.0 * cx.pt + tw + 5.0 * cx.pt + foot / 2.0), g.theta(mid));
                cx.scene.prims.push(Prim::Image { key, center, height: h, tint: Some(s.color) });
            }
        }
        GeoKind::Free { pos, .. } => {
            let tips: Vec<NodeId> = cx.tree.tips_below(n).into_iter().filter(|&t| cx.lay.visible[t]).collect();
            if tips.is_empty() {
                return;
            }
            let c = cx.pos(n);
            let (mut sx, mut sy) = (0.0, 0.0);
            for &t in &tips {
                sx += pos[t][0];
                sy += pos[t][1];
            }
            let m = [sx / tips.len() as f32, sy / tips.len() as f32];
            let a = (-(m[1] - c[1]) as f64).atan2((m[0] - c[0]) as f64);
            let far = tips.iter().map(|&t| ((pos[t][0] - c[0]).powi(2) + (pos[t][1] - c[1]).powi(2)).sqrt()).fold(0.0, f32::max);
            let r = far + start;
            let tp = [c[0] + r * a.cos() as f32, c[1] - r * a.sin() as f32];
            let (angle, flipped) = text_rot(a);
            cx.text(tp, s.text.clone(), size, s.color, angle, if flipped { 1.0 } else { 0.0 }, 0.5, false, true);
            if let Some((key, w, h)) = clade_pic(cx, s) {
                let foot = w * a.cos().abs() as f32 + h * a.sin().abs() as f32;
                let rr = r + cx.env.text_width(&s.text, size, false) + 5.0 * cx.pt + foot / 2.0;
                cx.scene.prims.push(Prim::Image { key, center: [c[0] + rr * a.cos() as f32, c[1] - rr * a.sin() as f32], height: h, tint: Some(s.color) });
            }
        }
    }
}

fn draw_scale_bar(cx: &mut Ctx, s: &ScaleBarStyle, bottom: f32, axis_h: f32, geo_h: f32) {
    let range = cx.lay.max_x - cx.lay.min_x;
    let len = s.length.filter(|l| *l > 0.0).unwrap_or_else(|| nice_step(range, 5));
    let px = len as f32 * cx.geo.s;
    let pt = cx.pt;
    if let GeoKind::Lin { orient: Orient::Down | Orient::Up, .. } = cx.geo.kind {
        // Vertical bar in the left gutter, starting level with the root.
        let u = cx.geo.u(cx.lay.min_x);
        let (ya, yb) = (cx.geo.uv(u, 0.0)[1], cx.geo.uv(u + px, 0.0)[1]);
        let x = 16.0 * pt;
        cx.scene.prims.push(Prim::Path { pts: vec![[x, ya], [x, yb]], color: s.color, width: s.width * pt, dashed: false });
        for yy in [ya, yb] {
            cx.scene.prims.push(Prim::Path { pts: vec![[x - 3.0 * pt, yy], [x + 3.0 * pt, yy]], color: s.color, width: s.width * pt, dashed: false });
        }
        cx.text([x - 4.0 * pt, (ya + yb) / 2.0], format_num(len, 4), s.size * pt, s.color, -FRAC_PI_2 as f32, 0.5, 1.0, false, false);
        return;
    }
    let y = cx.scene.height - bottom + 10.0 * pt + axis_h + geo_h + 8.0 * pt;
    let x = 14.0 * pt;
    cx.scene.prims.push(Prim::Path { pts: vec![[x, y], [x + px, y]], color: s.color, width: s.width * pt, dashed: false });
    for xx in [x, x + px] {
        cx.scene.prims.push(Prim::Path { pts: vec![[xx, y - 3.0 * pt], [xx, y + 3.0 * pt]], color: s.color, width: s.width * pt, dashed: false });
    }
    cx.text([x + px / 2.0, y - 4.0 * pt], format_num(len, 4), s.size * pt, s.color, 0.0, 0.5, 1.0, false, false);
}

fn draw_axis(cx: &mut Ctx, s: &AxisStyle) {
    let GeoKind::Lin { orient, v0, sp, .. } = cx.geo.kind else { return };
    let pt = cx.pt;
    let n = cx.geo.n as f32;
    let v_axis = v0 + n * sp + 6.0 * pt;
    let (x0, x1) = (0.0f64.max(cx.lay.min_x), cx.lay.max_x);
    let step = nice_step(x1 - x0, 6);
    let line = vec![cx.geo.uv(cx.geo.u(x0), v_axis), cx.geo.uv(cx.geo.u(x1), v_axis)];
    cx.scene.prims.push(Prim::Path { pts: line, color: s.line_color, width: s.line_width * pt, dashed: false });
    // Ticks are placed on round values of the displayed quantity.
    cx.scene.axis_before_present = s.time_before_present;
    // Counting back, the youngest tip sits at `youngest_age` rather than 0.
    let young = if s.time_before_present { s.youngest_age } else { 0.0 };
    cx.scene.axis_youngest = young;
    let disp = |x: f64| if s.time_before_present { young + x1 - x } else { x };
    let (d0, d1) = (disp(x0).min(disp(x1)), disp(x0).max(disp(x1)));
    let mut t = (d0 / step).ceil() * step;
    while t <= d1 + step * 1e-6 {
        let x = if s.time_before_present { x1 - (t - young) } else { t };
        let u = cx.geo.u(x);
        cx.scene.axis_ticks.push((x, t));
        if s.grid {
            let g = vec![cx.geo.uv(u, v0), cx.geo.uv(u, v_axis)];
            cx.scene.prims.push(Prim::Path { pts: g, color: s.grid_color, width: s.grid_width * pt, dashed: false });
        }
        let tick = vec![cx.geo.uv(u, v_axis), cx.geo.uv(u, v_axis + 4.0 * pt)];
        cx.scene.prims.push(Prim::Path { pts: tick, color: s.line_color, width: s.line_width * pt, dashed: false });
        let p = cx.geo.uv(u, v_axis + 6.0 * pt);
        let label = format_num(t.abs(), 4);
        match orient {
            Orient::Down | Orient::Up => cx.text(p, label, s.size * pt, s.line_color, 0.0, 0.0, 0.5, false, false),
            _ => cx.text(p, label, s.size * pt, s.line_color, 0.0, 0.5, 0.0, false, false),
        }
        t += step;
    }
}

fn draw_heatmap(cx: &mut Ctx, s: &HeatmapStyle, map: Option<&ColorMap>, start: f32) {
    let cw = s.cell_width * cx.pt;
    let leaves = cx.lay.leaves.clone();
    let max_x = cx.lay.max_x;
    let free = matches!(cx.geo.kind, GeoKind::Free { .. });
    for (ci, col) in s.columns.iter().enumerate() {
        let o0 = start + ci as f32 * cw;
        for &n in &leaves {
            let y = cx.lay.y[n];
            let base = if free { cx.leaf_x(n) } else { max_x };
            let fill = cx.tree.value(n, col).and_then(|v| map.and_then(|m| m.color(&v))).unwrap_or(Color::rgb(235, 235, 235));
            let region = cx.geo.cell(y, y, base, o0, o0 + cw, n);
            cx.scene.prims.push(region.prim(fill, Some((Color::WHITE, 0.3 * cx.pt))));
        }
        if s.show_names {
            if let GeoKind::Lin { orient: Orient::Right | Orient::Left, v0, .. } = cx.geo.kind {
                let u = cx.geo.u(max_x) + o0 + cw / 2.0;
                let p = cx.geo.uv(u, v0 - 3.0 * cx.pt);
                cx.text(p, col.clone(), s.name_size * cx.pt, Color::BLACK, -std::f32::consts::FRAC_PI_2, 0.0, 0.5, false, false);
            }
        }
    }
}

fn draw_bars(cx: &mut Ctx, s: &BarStyle, start: f32) {
    let pt = cx.pt;
    let leaves = cx.lay.leaves.clone();
    let vals: Vec<Option<f64>> = leaves.iter().map(|&n| cx.tree.value(n, &s.column).and_then(|v| v.as_f64()).filter(|v| v.is_finite())).collect();
    // The value range always includes zero, where the bars start.
    let lo = vals.iter().flatten().fold(0.0f64, |a, &v| a.min(v));
    let hi = vals.iter().flatten().fold(0.0f64, |a, &v| a.max(v));
    if hi - lo <= 0.0 {
        return;
    }
    let width = s.max_width * pt;
    let off = |v: f64| start + ((v - lo) / (hi - lo)) as f32 * width;
    let zero = off(0.0);
    let free = matches!(cx.geo.kind, GeoKind::Free { .. });
    for (i, &n) in leaves.iter().enumerate() {
        let Some(v) = vals[i] else { continue };
        let (a, b) = (zero.min(off(v)), zero.max(off(v)));
        let y = cx.lay.y[n];
        let base = if free { cx.leaf_x(n) } else { cx.lay.max_x };
        let sp = cx.geo.spacing(base);
        let shrink = 0.15 * sp;
        let fill = if v < 0.0 { s.negative_color } else { s.color };
        // Narrow the cell to leave gaps between bars.
        let mut region = cx.geo.cell(y, y, base, a, b, n);
        if let (GeoKind::Lin { orient, .. }, Region::Poly(pts)) = (&cx.geo.kind, &mut region) {
            for (j, p) in pts.iter_mut().enumerate() {
                let sign = if j < 2 { 1.0 } else { -1.0 };
                match orient {
                    Orient::Down | Orient::Up => p[0] += sign * shrink,
                    _ => p[1] += sign * shrink,
                }
            }
        }
        cx.scene.prims.push(region.prim(fill, None));
    }
    let GeoKind::Lin { orient, v0, sp, .. } = cx.geo.kind else { return };
    let ub = cx.geo.u(cx.lay.max_x);
    let v_end = v0 + cx.geo.n as f32 * sp;
    let grey = Color::rgb(120, 120, 120);
    if lo < 0.0 {
        // Zero line through the bars.
        cx.scene.prims.push(Prim::Path { pts: vec![cx.geo.uv(ub + zero, v0), cx.geo.uv(ub + zero, v_end)], color: grey, width: 0.6 * pt, dashed: false });
    }
    if !s.show_scale {
        return;
    }
    // Scale: an axis past the last tip with round tick values, and the column name.
    let size = s.label_size * pt;
    let (line, lw, ink) = (s.scale_color, s.scale_width * pt, s.label_color);
    let vertical = matches!(orient, Orient::Down | Orient::Up);
    let va = v_end + 4.0 * pt;
    cx.scene.prims.push(Prim::Path { pts: vec![cx.geo.uv(ub + off(lo), va), cx.geo.uv(ub + off(hi), va)], color: line, width: lw, dashed: false });
    // Round steps, fine enough for at least three ticks.
    // Fewer ticks for larger labels.
    let mut target = ((width / (5.0 * size)) as usize).clamp(2, 5);
    let mut step = nice_step(hi - lo, target);
    while ((hi / step).floor() - (lo / step).ceil()) < 2.0 && target < 10 {
        target += 1;
        step = nice_step(hi - lo, target);
    }
    let mut t = (lo / step).ceil() * step;
    while t <= hi + step * 1e-6 {
        let u = ub + off(t);
        cx.scene.prims.push(Prim::Path { pts: vec![cx.geo.uv(u, va), cx.geo.uv(u, va + 3.0 * pt)], color: line, width: lw, dashed: false });
        let label = format_num(if t.abs() < step * 1e-9 { 0.0 } else { t }, 4);
        let p = cx.geo.uv(u, va + 4.0 * pt);
        if vertical {
            cx.text(p, label, size, ink, 0.0, 0.0, 0.5, false, false);
        } else {
            cx.text(p, label, size, ink, 0.0, 0.5, 0.0, false, false);
        }
        t += step;
    }
    let mid = ub + (off(lo) + off(hi)) / 2.0;
    if vertical {
        // Past the far end of the bars, along the axis.
        let p = cx.geo.uv(ub + off(hi) + 4.0 * pt, va);
        let valign = if orient == Orient::Down { 0.0 } else { 1.0 };
        cx.text(p, s.column.clone(), size, ink, 0.0, 0.5, valign, false, false);
    } else {
        cx.text(cx.geo.uv(mid, va + size + 6.0 * pt), s.column.clone(), size, ink, 0.0, 0.5, 0.0, false, false);
    }
}

/// Attribute naming the taxon (or image) whose silhouette a node shows.
pub const PHYLOPIC_KEY: &str = "phylopic";

/// Image key of the silhouette node `n` shows: its `phylopic` attribute if set,
/// else a tip's label; internal nodes need the attribute.
fn silhouette_key(tree: &Tree, n: NodeId) -> Option<String> {
    match tree.nodes[n].attrs.get(PHYLOPIC_KEY) {
        Some(Attr::Text(name)) if !name.trim().is_empty() => Some(image_key(name)),
        _ if tree.is_tip(n) => Some(image_key(tree.label(n))),
        _ => None,
    }
}

/// Nodes a PhyloPic layer draws: tips (if shown) and internal nodes with a
/// `phylopic` attribute (if shown).
fn silhouette_nodes(tree: &Tree, s: &PhylopicStyle, leaves: &[NodeId], visible: &[NodeId]) -> Vec<NodeId> {
    let mut v: Vec<NodeId> = if s.tips { leaves.iter().copied().filter(|&n| tree.is_tip(n)).collect() } else { Vec::new() };
    if s.nodes {
        v.extend(visible.iter().copied().filter(|&n| !tree.is_tip(n) && silhouette_key(tree, n).is_some()));
    }
    v
}

/// Color scale of a PhyloPic layer colored by a trait, over the silhouettes it draws.
fn phylopic_map(tree: &Tree, s: &PhylopicStyle, leaves: &[NodeId], visible: &[NodeId]) -> Option<ColorMap> {
    let key = s.color_by.as_ref()?;
    Some(ColorMap::build(tree, key, &silhouette_nodes(tree, s, leaves, visible), s.palette)?.reversed(s.reverse_palette))
}

#[allow(clippy::too_many_arguments)]
fn draw_phylopic(cx: &mut Ctx, s: &PhylopicStyle, map: Option<&ColorMap>, visible: &[NodeId], start: f32, label_off: f32, widths: &[f32], labels_aligned: bool) {
    // Never taller than the space between neighbouring tips.
    let room = match cx.geo.kind {
        GeoKind::Polar { step, .. } => ((cx.geo.r(cx.lay.max_x) + start) * step as f32).abs(),
        _ => cx.geo.spacing(cx.lay.max_x),
    };
    let leaves = cx.lay.leaves.clone();
    // Height by a numeric trait: linear from `size_min` (smallest) to `size` (largest).
    let drawn = silhouette_nodes(cx.tree, s, &leaves, visible);
    let size_of = |tree: &Tree, n: NodeId| s.size_by.as_ref().and_then(|k| tree.value(n, k)).and_then(|v| v.as_f64()).filter(|v| v.is_finite());
    let (lo, hi) = drawn.iter().filter_map(|&n| size_of(cx.tree, n)).fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
    let height = |tree: &Tree, n: NodeId| -> f32 {
        match size_of(tree, n) {
            Some(v) if hi > lo => (s.size_min + (s.size - s.size_min) * ((v - lo) / (hi - lo)) as f32) * cx.pt,
            _ => s.size * cx.pt,
        }
    };
    let tint = |tree: &Tree, n: NodeId| s.color_by.as_ref().and_then(|k| tree.value(n, k)).and_then(|v| map.and_then(|m| m.color(&v))).or(s.tint);
    let free = matches!(cx.geo.kind, GeoKind::Free { .. });
    if s.tips {
        for (i, &n) in leaves.iter().enumerate() {
            if !cx.tree.is_tip(n) {
                continue;
            }
            let Some(key) = silhouette_key(cx.tree, n) else { continue };
            cx.scene.wanted_images.push(key.clone());
            let Some(aspect) = cx.env.image_aspect(&key) else { continue };
            let h = height(cx.tree, n);
            let h = if s.fit_spacing { h.min(room * 0.92).max(2.0) } else { h };
            let w = h * aspect.min(2.5);
            let y = cx.lay.y[n];
            // Rotated layouts: offset along the radial direction by the image's
            // footprint in that direction.
            let a = cx.geo.out_angle(y, n);
            let foot = (w * a.cos().abs() as f32 + h * a.sin().abs() as f32) / 2.0;
            // Unaligned silhouettes follow their own label, wherever it ends.
            let after_label = label_off + widths.get(i).copied().unwrap_or(0.0) + 4.0 * cx.pt + s.offset * cx.pt + foot;
            let (base, off) = if s.align && !free {
                (cx.lay.max_x, start + foot)
            } else if labels_aligned && !free {
                (cx.lay.max_x, after_label)
            } else {
                (cx.leaf_x(n), after_label)
            };
            let center = cx.geo.ext(y, base, off, n);
            cx.scene.prims.push(Prim::Image { key, center, height: h, tint: tint(cx.tree, n) });
        }
    }
    if s.nodes {
        // Internal nodes with a silhouette: just before the node, above its stem
        // branch (clear of the lines to its children).
        for &n in visible {
            if cx.tree.is_tip(n) {
                continue;
            }
            let Some(key) = silhouette_key(cx.tree, n) else { continue };
            cx.scene.wanted_images.push(key.clone());
            let Some(aspect) = cx.env.image_aspect(&key) else { continue };
            let h = height(cx.tree, n);
            let w = h * aspect.min(2.5);
            let y = cx.lay.y[n];
            let back = ((w / 2.0 + 3.0 * cx.pt) / cx.geo.s.max(1e-6)) as f64;
            let p = if free { cx.pos(n) } else { cx.geo.at(cx.lay.x[n] - back, y, n) };
            let center = [p[0], p[1] - h / 2.0 - 2.0 * cx.pt];
            cx.scene.prims.push(Prim::Image { key, center, height: h, tint: tint(cx.tree, n) });
        }
    }
}

/// Width of a size legend's mark column at legend scale `lp`: lines are a
/// short stroke, shapes and silhouettes as wide as the largest.
fn size_mark_w(k: &SizeKey, lp: f32) -> f32 {
    match k.mark {
        SizeMark::Line => 18.0 * lp,
        SizeMark::Shape(_) => (k.large.max(k.small) * lp).max(8.0 * lp),
        SizeMark::Image { aspect, .. } => (k.large.max(k.small) * lp * aspect.min(2.5)).max(8.0 * lp),
    }
}

fn draw_legends(cx: &mut Ctx, legends: &[ColorMap], sizes: &[SizeKey], x: f32, top: f32, scale: f32) {
    // Legend sizes scale with the view's legend size.
    let pt = cx.pt * scale.max(0.1);
    let mut y = top;
    for m in legends {
        cx.text([x, y], m.title().to_string(), 9.0 * pt, Color::BLACK, 0.0, 0.0, 0.0, false, true);
        y += 14.0 * pt;
        match m {
            ColorMap::Discrete { levels, .. } => {
                for (l, c) in levels.iter().take(30) {
                    let r = 4.0 * pt;
                    cx.scene.prims.push(Prim::Poly {
                        pts: vec![[x, y], [x + 2.0 * r, y], [x + 2.0 * r, y + 2.0 * r], [x, y + 2.0 * r]],
                        fill: *c,
                        stroke: None,
                    });
                    cx.text([x + 2.0 * r + 4.0 * pt, y + r], l.clone(), 8.0 * pt, Color::BLACK, 0.0, 0.0, 0.5, false, false);
                    y += 11.0 * pt;
                }
                if levels.len() > 30 {
                    cx.text([x, y], format!("… {} more", levels.len() - 30), 8.0 * pt, Color::GREY, 0.0, 0.0, 0.0, false, false);
                    y += 11.0 * pt;
                }
            }
            ColorMap::Continuous { min, max, palette, reversed, .. } => {
                let h = 70.0 * pt;
                let steps = 30;
                for i in 0..steps {
                    let t0 = i as f32 / steps as f32;
                    let t1 = (i + 1) as f32 / steps as f32;
                    // Top of the bar is the maximum.
                    let t = 1.0 - (t0 + t1) / 2.0;
                    let c = palette.continuous(if *reversed { 1.0 - t } else { t });
                    cx.scene.prims.push(Prim::Poly {
                        pts: vec![[x, y + t0 * h], [x + 10.0 * pt, y + t0 * h], [x + 10.0 * pt, y + t1 * h + 0.5], [x, y + t1 * h + 0.5]],
                        fill: c,
                        stroke: None,
                    });
                }
                cx.text([x + 14.0 * pt, y], format_num(*max, 3), 8.0 * pt, Color::BLACK, 0.0, 0.0, 0.0, false, false);
                cx.text([x + 14.0 * pt, y + h], format_num(*min, 3), 8.0 * pt, Color::BLACK, 0.0, 0.0, 1.0, false, false);
                y += h + 6.0 * pt;
            }
        }
        y += 8.0 * pt;
    }
    for k in sizes {
        cx.text([x, y], k.title.clone(), 9.0 * pt, Color::BLACK, 0.0, 0.0, 0.0, false, true);
        y += 14.0 * pt;
        let mark_w = size_mark_w(k, pt);
        // Largest value on top, as in the color gradients.
        for (v, size) in k.rows().into_iter().rev() {
            // Marks at their size in the tree times the legend size; rows at
            // least a text line apart.
            let row = (size * pt).max(11.0 * pt);
            let mid = y + row / 2.0;
            match &k.mark {
                SizeMark::Line => cx.scene.prims.push(Prim::Path { pts: vec![[x, mid], [x + mark_w, mid]], color: k.color, width: size * pt, dashed: false }),
                SizeMark::Shape(sh) => cx.scene.prims.push(shape_prim(*sh, [x + mark_w / 2.0, mid], size * pt, k.color)),
                SizeMark::Image { key, tint, .. } => cx.scene.prims.push(Prim::Image { key: key.clone(), center: [x + mark_w / 2.0, mid], height: size * pt, tint: *tint }),
            }
            cx.text([x + mark_w + 4.0 * pt, mid], format_num(v, 3), 8.0 * pt, Color::BLACK, 0.0, 0.0, 0.5, false, false);
            y += row + 2.0 * pt;
        }
        y += 8.0 * pt;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick;

    #[test]
    fn axis_boundary_and_legend_styles() {
        let mut t = parse_newick("((A:1,B:2):1,(C:1,D:1):1);").unwrap();
        for (i, n) in t.tips().into_iter().enumerate() {
            t.nodes[n].attrs.insert("diet".into(), Attr::Text(if i % 2 == 0 { "a".into() } else { "b".into() }));
        }
        let red = Color::rgb(200, 0, 0);
        let blue = Color::rgb(0, 0, 200);
        let green = Color::rgb(0, 150, 0);
        let mut axis = Layer::default_axis();
        if let Layer::TimeAxis(s) = &mut axis {
            (s.grid, s.line_color, s.line_width, s.grid_color, s.grid_width) = (true, red, 3.0, blue, 2.0);
        }
        let mut geo = Layer::default_geoscale();
        if let Layer::Geoscale(g) = &mut geo {
            (g.boundaries, g.boundary_color, g.boundary_width, g.ma_per_unit) = (true, green, 1.5, 100.0);
        }
        let mut points = Layer::default_tip_points();
        if let Layer::TipPoints(p) = &mut points {
            p.color_by = Some("diet".into());
        }
        let layers: Vec<LayerEntry> = [Layer::default_tree(), axis, geo, points].into_iter().map(LayerEntry::new).collect();
        let mut view = ViewState::default();
        let build_with = |view: &ViewState| build(&BuildInput { tree: &t, view, layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 600.0, 1.0);
        let s = build_with(&view);
        let has = |c: Color, w: f32| s.prims.iter().any(|p| matches!(p, Prim::Path { color, width, .. } if *color == c && (*width - w).abs() < 1e-6));
        assert!(has(red, 3.0) && has(blue, 2.0) && has(green, 1.5));
        // Axis ticks are reported for snapping: counting back, 0 sits at the deepest tip (depth 3).
        assert!(s.axis_before_present && s.axis_ticks.iter().any(|&(x, t)| t == 0.0 && (x - 3.0).abs() < 1e-9), "{:?}", s.axis_ticks);
        // With the youngest tip at 5, the axis reads 5 there and counts up toward the root.
        let mut aged = layers.clone();
        for e in aged.iter_mut() {
            if let Layer::TimeAxis(a) = &mut e.layer {
                a.youngest_age = 5.0;
            }
        }
        let s5 = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &aged, overlay: &[] }, &ApproxEnv, 600.0, 600.0, 1.0);
        assert!(s5.axis_ticks.iter().any(|&(x, v)| v == 5.0 && (x - 3.0).abs() < 1e-9), "{:?}", s5.axis_ticks);
        assert!(s5.axis_ticks.iter().all(|&(_, v)| v >= 5.0 - 1e-9));
        // Legend title text grows with the legend size.
        let title_size = |s: &Scene| s.prims.iter().find_map(|p| match p { Prim::Text { text, size, .. } if text == "diet" => Some(*size), _ => None }).unwrap();
        let normal = title_size(&s);
        view.legend_scale = 2.0;
        assert!((title_size(&build_with(&view)) - 2.0 * normal).abs() < 1e-4);
    }

    #[test]
    fn clade_label_gradient_behind_branches() {
        let t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        let ab = t.nodes[t.find_label("A").unwrap()].parent.unwrap();
        let mut s = CladeLabelStyle::new(ab, "AB".into());
        s.gradient = true;
        let layers: Vec<LayerEntry> = [Layer::default_tree(), Layer::CladeLabel(s)].into_iter().map(LayerEntry::new).collect();
        let sc = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        // One gradient band, before the first branch, transparent at the clade's node.
        let first_path = sc.prims.iter().position(|p| matches!(p, Prim::Path { .. })).unwrap();
        let shades: Vec<usize> = sc.prims.iter().enumerate().filter(|(_, p)| matches!(p, Prim::Shade { .. })).map(|(i, _)| i).collect();
        assert_eq!(shades.len(), 1);
        assert!(shades[0] < first_path);
        let Prim::Shade { c0, c1, fade, a, .. } = &sc.prims[shades[0]] else { unreachable!() };
        assert_eq!((c0.a, c1.a), (0, 200));
        assert!(fade.t(a[0]) < 1e-6 && (fade.t(a[1]) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn branch_labels_by_absolute_value() {
        let mut t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        for (name, v) in [("A", -3.0), ("B", 0.5), ("C", 2.0)] {
            let n = t.find_label(name).unwrap();
            t.nodes[n].attrs.insert("fabric_beta_mean".into(), Attr::Num(v));
        }
        let labels = |min: f64, abs: bool| {
            let l = Layer::BranchLengths(BranchLengthStyle { size: 7.0, color: Color::BLACK, digits: 1, attr: Some("fabric_beta_mean".into()), hide_value: None, min_pp: 0.0, min_value: Some(min), min_abs: abs });
            let layers = vec![LayerEntry::new(l)];
            let s = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
            let mut v: Vec<String> = s.prims.iter().filter_map(|p| if let Prim::Text { text, .. } = p { Some(text.clone()) } else { None }).collect();
            v.sort();
            v
        };
        assert_eq!(labels(1.0, false), ["2"]);
        assert_eq!(labels(1.0, true), ["-3", "2"]);
    }

    #[test]
    fn root_edge_with_or_without_a_root_length() {
        let paths = |t: &Tree, view: &ViewState| {
            let layers = vec![LayerEntry::new(Layer::default_tree())];
            let s = build(&BuildInput { tree: t, view, layers: &layers, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
            s.prims.iter().filter(|p| matches!(p, Prim::Path { .. })).count()
        };
        let t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        let mut view = ViewState::default();
        let base = paths(&t, &view);
        // No root length: nothing until a length is set.
        view.layout.root_edge = true;
        assert_eq!(paths(&t, &view), base);
        view.layout.root_edge_len = 0.2;
        assert_eq!(paths(&t, &view), base + 1);
        // A root with its own length shows it when no length is set by hand.
        let t2 = parse_newick("((A:1,B:1):1,C:2):0.5;").unwrap();
        view.layout.root_edge_len = 0.0;
        assert_eq!(crate::layout::root_edge_length(&t2, &view.layout), 0.5);
        assert_eq!(paths(&t2, &view), base + 1);
    }

    #[test]
    fn clade_label_silhouette() {
        struct Square;
        impl SceneEnv for Square {
            fn text_width(&self, text: &str, size: f32, _italic: bool) -> f32 {
                text.chars().count() as f32 * size * 0.55
            }
            fn image_aspect(&self, _key: &str) -> Option<f32> {
                Some(1.0)
            }
        }
        let t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        let ab = t.nodes[t.find_label("A").unwrap()].parent.unwrap();
        let mut s = CladeLabelStyle::new(ab, "Hominidae".into());
        s.phylopic = Some("Hominidae".into());
        s.color = Color::rgb(200, 0, 0);
        let layers = vec![LayerEntry::new(Layer::default_tree()), LayerEntry::new(Layer::CladeLabel(s))];
        let sc = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &Square, 600.0, 400.0, 1.0);
        assert!(sc.wanted_images.contains(&"Hominidae".to_string()));
        // One silhouette in the label's color, past the end of the label text, inside the figure.
        let text_x = sc.prims.iter().find_map(|p| if let Prim::Text { text, pos, .. } = p { (text == "Hominidae").then_some(pos[0]) } else { None }).unwrap();
        let pics: Vec<(P, Option<Color>)> = sc.prims.iter().filter_map(|p| if let Prim::Image { center, tint, .. } = p { Some((*center, *tint)) } else { None }).collect();
        assert_eq!(pics.len(), 1);
        assert_eq!(pics[0].1, Some(Color::rgb(200, 0, 0)));
        assert!(pics[0].0[0] > text_x + 9.0 * 11.0 * 0.55 && pics[0].0[0] + 12.0 <= 600.0);
    }

    #[test]
    fn branch_widths_by_trait() {
        let mut t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        for (name, v) in [("A", 10.0), ("B", 20.0)] {
            let n = t.find_label(name).unwrap();
            t.nodes[n].attrs.insert("mass".into(), Attr::Num(v));
        }
        let mut l = Layer::default_tree();
        if let Layer::Tree(s) = &mut l {
            s.width_by = Some("mass".into());
            s.width = 4.0;
            s.width_min = 1.0;
        }
        let layers = vec![LayerEntry::new(l)];
        let s = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
        let mut widths: Vec<f32> = s.prims.iter().filter_map(|p| if let Prim::Path { width, color, .. } = p { (*color == Color::BLACK).then_some(*width) } else { None }).collect();
        widths.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // A (smallest) thinnest, B (largest) thickest; C and the AB branch have no
        // value: thinnest. The size legend adds lines of 1, 2.5 and 4.
        assert_eq!(widths.first(), Some(&1.0));
        assert_eq!(widths.last(), Some(&4.0));
        assert_eq!(widths.iter().filter(|w| **w == 1.0).count(), 4);
        assert_eq!(widths.iter().filter(|w| **w == 2.5).count(), 1);
        let texts: Vec<&str> = s.prims.iter().filter_map(|p| if let Prim::Text { text, .. } = p { Some(text.as_str()) } else { None }).collect();
        for want in ["mass", "10", "15", "20"] {
            assert!(texts.contains(&want), "{want} missing from {texts:?}");
        }
    }

    #[test]
    fn point_sizes_get_a_legend() {
        let mut t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        for (name, v) in [("A", 1.0), ("B", 3.0), ("C", 5.0)] {
            let n = t.find_label(name).unwrap();
            t.nodes[n].attrs.insert("size".into(), Attr::Num(v));
        }
        let mut l = Layer::default_tip_points();
        if let Layer::TipPoints(s) = &mut l {
            (s.size_by, s.size_min, s.size, s.shape) = (Some("size".into()), 2.0, 10.0, Shape::Square);
        }
        let layers = vec![LayerEntry::new(Layer::default_tree()), LayerEntry::new(l)];
        let s = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
        let texts: Vec<&str> = s.prims.iter().filter_map(|p| if let Prim::Text { text, .. } = p { Some(text.as_str()) } else { None }).collect();
        for want in ["size", "1", "3", "5"] {
            assert!(texts.contains(&want), "{want} missing from {texts:?}");
        }
        // Three tip points and three legend marks.
        let marks = s.prims.iter().filter(|p| matches!(p, Prim::Poly { stroke: None, fill, .. } if *fill == Color::rgb(0, 114, 178))).count();
        assert_eq!(marks, 6);
        // The legend size scales the size legend's marks too: the largest
        // (10 pt) becomes 20 pt wide, while the tree's points stay at 10.
        let widest = |s: &Scene| {
            s.prims
                .iter()
                .filter_map(|p| match p {
                    Prim::Poly { pts, stroke: None, .. } => Some(pts.iter().map(|q| q[0]).fold(f32::MIN, f32::max) - pts.iter().map(|q| q[0]).fold(f32::MAX, f32::min)),
                    _ => None,
                })
                .fold(0.0, f32::max)
        };
        let view = ViewState { legend_scale: 2.0, ..ViewState::default() };
        let big = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
        assert!((widest(&s) - 10.0).abs() < 1e-3 && (widest(&big) - 20.0).abs() < 1e-3);
    }

    #[test]
    fn leftover_bayestraits_filters_do_not_hide_trait_colors() {
        // A value cut-off set for rate scalars must not apply once branches are
        // colored by an ordinary trait (its slider is hidden then).
        let mut t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        for n in t.preorder() {
            t.nodes[n].attrs.insert("mass".into(), Attr::Num(n as f64));
        }
        let mut l = Layer::default_tree();
        if let Layer::Tree(s) = &mut l {
            s.color_by = Some("mass".into());
            s.min_value = Some(1e9);
            s.min_abs = true;
            s.min_pp = 0.9;
        }
        let layers = vec![LayerEntry::new(l)];
        let s = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
        let black = s.prims.iter().filter(|p| matches!(p, Prim::Path { color, .. } if *color == Color::BLACK)).count();
        assert_eq!(black, 0, "every branch takes its trait color");
    }

    #[test]
    fn branches_without_a_value_are_grey_and_dashed() {
        // Tip values only (a data column not reconstructed): the internal branch is grey and dashed.
        let mut t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        for (name, v) in [("A", 1.0), ("B", 2.0), ("C", 3.0)] {
            let n = t.find_label(name).unwrap();
            t.nodes[n].attrs.insert("mass".into(), Attr::Num(v));
        }
        let mut l = Layer::default_tree();
        if let Layer::Tree(s) = &mut l {
            s.color_by = Some("mass".into());
        }
        let layers = vec![LayerEntry::new(l)];
        let s = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
        let paths: Vec<(Color, bool)> = s.prims.iter().filter_map(|p| if let Prim::Path { color, dashed, pts, .. } = p { (pts.len() > 1).then_some((*color, *dashed)) } else { None }).collect();
        let dashed: Vec<&(Color, bool)> = paths.iter().filter(|p| p.1).collect();
        assert_eq!(dashed.len(), 1);
        assert_eq!(dashed[0].0, Color::rgb(160, 160, 160));
        assert_eq!(paths.iter().filter(|p| !p.1 && p.0 != Color::BLACK).count(), 3);
    }

    #[test]
    fn tip_points_sized_by_trait() {
        let mut t = parse_newick("((A:1,B:1):1,(C:2,D:2):0);").unwrap();
        for (name, v) in [("A", 10.0), ("B", 20.0), ("C", 30.0)] {
            let n = t.find_label(name).unwrap();
            t.nodes[n].attrs.insert("mass".into(), Attr::Num(v));
        }
        let mut l = Layer::default_tip_points();
        if let Layer::TipPoints(s) = &mut l {
            (s.size_by, s.size_min, s.size) = (Some("mass".into()), 2.0, 10.0);
        }
        let layers: Vec<LayerEntry> = [Layer::default_tree(), l].into_iter().map(LayerEntry::new).collect();
        let s = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        let mut radii: Vec<f32> = s.prims.iter().filter_map(|p| match p { Prim::Circle { r, .. } => Some(*r), _ => None }).collect();
        radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // Three points (D has no mass): smallest, middle and largest size, and
        // the size legend's three circles at the same sizes.
        radii.dedup();
        assert_eq!(radii.len(), 3);
        assert!(radii[0] < radii[1] && radii[1] < radii[2]);
        assert!((radii[2] / radii[0] - 5.0).abs() < 1e-3, "{:?}", radii);
    }

    #[test]
    fn silhouettes_sized_and_colored_by_trait_and_at_nodes() {
        // Every image is "loaded" and square, so silhouettes are drawn.
        struct Square;
        impl SceneEnv for Square {
            fn text_width(&self, text: &str, size: f32, _italic: bool) -> f32 {
                text.chars().count() as f32 * size * 0.55
            }
            fn image_aspect(&self, _key: &str) -> Option<f32> {
                Some(1.0)
            }
        }
        let mut t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        for (name, v) in [("A", 1.0), ("B", 2.0), ("C", 3.0)] {
            let n = t.find_label(name).unwrap();
            t.nodes[n].attrs.insert("mass".into(), Attr::Num(v));
        }
        let ab = t.nodes[t.find_label("A").unwrap()].parent.unwrap();
        t.nodes[ab].attrs.insert(PHYLOPIC_KEY.into(), Attr::Text("Hominidae".into()));
        let mut l = Layer::default_phylopic();
        if let Layer::Phylopic(s) = &mut l {
            (s.size_by, s.size_min, s.size, s.fit_spacing, s.color_by) = (Some("mass".into()), 6.0, 18.0, false, Some("mass".into()));
        }
        let layers: Vec<LayerEntry> = [Layer::default_tree(), l].into_iter().map(LayerEntry::new).collect();
        let s = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &Square, 600.0, 400.0, 1.0);
        let images: Vec<(&String, f32, Option<Color>)> = s.prims.iter().filter_map(|p| match p { Prim::Image { key, height, tint, .. } => Some((key, *height, *tint)), _ => None }).collect();
        let tip = |k: &str| images.iter().find(|i| i.0 == k).copied().unwrap();
        // Heights run from 6 pt (smallest mass) to 18 pt (largest); colors differ.
        assert!((tip("A").1 - 6.0).abs() < 1e-4 && (tip("B").1 - 12.0).abs() < 1e-4 && (tip("C").1 - 18.0).abs() < 1e-4);
        assert!(tip("A").2.is_some() && tip("A").2 != tip("C").2);
        // The node given a taxon gets its silhouette (no mass, so the full height).
        assert_eq!(tip("Hominidae").1, 18.0);
        assert!(s.wanted_images.contains(&"Hominidae".to_string()));
        // Size legend: C's silhouette (the largest mass) at 18, 12 and 6 pt, in grey.
        let grey = Some(Color::rgb(110, 110, 110));
        let legend: Vec<f32> = images.iter().filter(|i| i.0 == "C" && i.2 == grey).map(|i| i.1).collect();
        assert_eq!(legend, [18.0, 12.0, 6.0]);
        let texts: Vec<&str> = s.prims.iter().filter_map(|p| if let Prim::Text { text, .. } = p { Some(text.as_str()) } else { None }).collect();
        assert!(["1", "2", "3"].iter().all(|v| texts.contains(v)));
    }

    #[test]
    fn heights_count_back_from_a_fossil_youngest_tip() {
        // BEAST heights are absolute: the youngest tip is 10 Ma old, its parent 13 Ma
        // with a 95% HPD of 12–14 Ma, i.e. 2–4 units before the tips.
        let t = parse_newick("((A[&height=10]:3,B[&height=10]:3)[&height=13,height_95%_HPD={12,14}]:2,C[&height=10]:5)[&height=15];").unwrap();
        assert_eq!(t.youngest_tip_age(), Some(10.0));
        let layers: Vec<LayerEntry> = [Layer::default_tree(), Layer::default_node_bars()].into_iter().map(LayerEntry::new).collect();
        let s = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        let bar = s.prims.iter().find_map(|p| match p { Prim::Path { pts, width, .. } if (*width - 5.0).abs() < 1e-6 => Some(pts.clone()), _ => None }).unwrap();
        let tip = s.node_pos[t.find_label("A").unwrap()].unwrap();
        assert!(((tip[0] - bar[0][0]) / s.px_per_unit - 4.0).abs() < 1e-3, "{:?} {:?}", bar, tip);
        assert!(((tip[0] - bar[1][0]) / s.px_per_unit - 2.0).abs() < 1e-3);
    }

    #[test]
    fn fossil_range_bars() {
        // Deepest tip at the present, root at 10 Ma (1 Ma per unit).
        let mut t = parse_newick("((A:4,B:10):0,C:10);").unwrap();
        let a = t.find_label("A").unwrap();
        crate::fossil::set_range(&mut t, a, 8.0, 6.0);
        let layers: Vec<LayerEntry> = [Layer::default_tree(), Layer::default_fossil_ranges()].into_iter().map(LayerEntry::new).collect();
        let s = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        let bars: Vec<&Vec<P>> = s.prims.iter().filter_map(|p| match p { Prim::Path { pts, width, .. } if (*width - 6.0).abs() < 1e-6 => Some(pts), _ => None }).collect();
        assert_eq!(bars.len(), 1);
        // From 8 Ma to 6 Ma: two units, at A's height.
        let (p0, p1) = (bars[0][0], bars[0][1]);
        assert!(((p1[0] - p0[0]) / s.px_per_unit - 2.0).abs() < 1e-3);
        assert!((p0[1] - s.node_pos[a].unwrap()[1]).abs() < 1e-3);
        // With a time axis whose youngest tip is at 20 Ma, a 28–26 Ma range sits
        // where 8–6 did: ages follow the axis's time frame.
        let mut t2 = t.clone();
        crate::fossil::set_range(&mut t2, a, 28.0, 26.0);
        let mut axis = Layer::default_axis();
        if let Layer::TimeAxis(x) = &mut axis {
            x.youngest_age = 20.0;
        }
        let layers2: Vec<LayerEntry> = [Layer::default_tree(), axis, Layer::default_fossil_ranges()].into_iter().map(LayerEntry::new).collect();
        let s2 = build(&BuildInput { tree: &t2, view: &ViewState::default(), layers: &layers2, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        let bar = s2.prims.iter().find_map(|p| match p { Prim::Path { pts, width, .. } if (*width - 6.0).abs() < 1e-6 => Some(pts.clone()), _ => None }).unwrap();
        let tip = s2.node_pos[t2.find_label("B").unwrap()].unwrap();
        // 26 Ma is 6 units before the youngest tip (B, at 20 Ma).
        assert!(((tip[0] - bar[1][0]) / s2.px_per_unit - 6.0).abs() < 1e-3);
    }

    #[test]
    fn geoscale_boundaries_follow_the_finest_level() {
        use crate::geotime::{overlapping, Level};
        // Root at 60 Ma, youngest tip at the present.
        let t = parse_newick("((A:20,B:20):40,C:60);").unwrap();
        let green = Color::rgb(0, 150, 0);
        let count = |epochs: bool, youngest: f64| {
            let mut geo = Layer::default_geoscale();
            if let Layer::Geoscale(g) = &mut geo {
                (g.boundaries, g.boundary_color, g.epochs, g.periods, g.youngest_age) = (true, green, epochs, true, youngest);
            }
            let layers: Vec<LayerEntry> = [Layer::default_tree(), geo].into_iter().map(LayerEntry::new).collect();
            let s = build(&BuildInput { tree: &t, view: &ViewState::default(), layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 600.0, 1.0);
            s.prims.iter().filter(|p| matches!(p, Prim::Path { color, dashed: false, .. } if *color == green)).count()
        };
        // Boundaries more than 0.1 Ma from the youngest tip (about 1 px here): the
        // Holocene starts 0.0117 Ma ago, so its line merges with the present's.
        let starts = |level: Level, young: f64| overlapping(level, young, young + 60.0).filter(|i| i.start < young + 60.0 && i.start > young + 0.1).count();
        // Epoch boundaries plus the present; periods when epochs are off.
        assert_eq!(count(true, 0.0), starts(Level::Epoch, 0.0) + 1);
        assert_eq!(count(false, 0.0), starts(Level::Period, 0.0) + 1);
        // No present line when the youngest tip is older than 0 Ma.
        assert_eq!(count(true, 10.0), starts(Level::Epoch, 10.0));
        assert!(starts(Level::Epoch, 0.0) > starts(Level::Period, 0.0));
    }

    #[test]
    fn depth_directions_for_dragging_branch_lengths() {
        let t = parse_newick("((A:1,B:2):1,(C:1,D:1):1);").unwrap();
        let layers = default_layers();
        let mut view = ViewState::default();
        let a = t.find_label("A").unwrap();
        let dirs = |view: &ViewState| build(&BuildInput { tree: &t, view, layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 600.0, 1.0);
        assert_eq!(dirs(&view).depth_dirs[a], Some([1.0, 0.0]));
        view.layout.tips = TipSide::Top;
        assert_eq!(dirs(&view).depth_dirs[a], Some([0.0, -1.0]));
        // Circular: away from the root.
        view.layout.kind = LayoutKind::Circular;
        let s = dirs(&view);
        let (p, r, d) = (s.node_pos[a].unwrap(), s.node_pos[t.root].unwrap(), s.depth_dirs[a].unwrap());
        let (dx, dy) = (p[0] - r[0], p[1] - r[1]);
        assert!(((dx * d[0] + dy * d[1]) / (dx * dx + dy * dy).sqrt() - 1.0).abs() < 1e-3);
    }

    #[test]
    fn reverse_tip_order_in_circular_layouts() {
        let t = parse_newick("((A:1,B:2):1,(C:1,D:1):1);").unwrap();
        let layers = default_layers();
        let mut view = ViewState::default();
        view.layout.kind = LayoutKind::Circular;
        let pos = |view: &ViewState| build(&BuildInput { tree: &t, view, layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 600.0, 1.0).node_pos;
        let plain = pos(&view);
        view.layout.flip_y = true;
        let flipped = pos(&view);
        let (a, d) = (t.find_label("A").unwrap(), t.find_label("D").unwrap());
        // A and D trade places around the circle (angles about the root).
        let angle = |p: &[Option<P>], n: NodeId| {
            let (q, r) = (p[n].unwrap(), p[t.root].unwrap());
            (q[1] - r[1]).atan2(q[0] - r[0])
        };
        let same = |x: f32, y: f32| (x - y).abs() < 1e-3;
        assert!(!same(angle(&plain, a), angle(&flipped, a)));
        assert!(same(angle(&plain, d), angle(&flipped, a)));
        assert!(same(angle(&plain, a), angle(&flipped, d)));
    }

    #[test]
    fn builds_all_layouts() {
        let mut t = parse_newick("((A:1,B:2)0.9:1,(C:1,D:1):1);").unwrap();
        let a = t.find_label("A").unwrap();
        t.nodes[a].attrs.insert("trait".into(), Attr::Num(3.0));
        let mut layers = default_layers();
        layers.push(LayerEntry::new(Layer::default_node_labels("support")));
        layers.push(LayerEntry::new(Layer::default_axis()));
        layers.push(LayerEntry::new(Layer::Heatmap(HeatmapStyle {
            columns: vec!["trait".into()],
            offset: 2.0,
            cell_width: 10.0,
            palette: Palette::Viridis,
            reverse_palette: false,
            show_names: true,
            name_size: 8.0,
        })));
        let ab = t.nodes[a].parent.unwrap();
        layers.insert(0, LayerEntry::new(Layer::Highlight(HighlightStyle { node: ab, fill: Color::rgb(255, 0, 0).with_alpha(60), extend: 0.0 })));
        for kind in LayoutKind::ALL {
            let mut view = ViewState::default();
            view.layout.kind = kind;
            let s = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 800.0, 600.0, 1.0);
            assert!(!s.prims.is_empty(), "{:?}", kind);
            for n in t.preorder() {
                let p = s.node_pos[n].unwrap();
                assert!(p[0].is_finite() && p[1].is_finite(), "{:?}", kind);
                assert!(p[0] >= -1.0 && p[0] <= 801.0 && p[1] >= -1.0 && p[1] <= 601.0, "{:?} {:?}", kind, p);
            }
        }
    }

    #[test]
    fn branch_lengths_and_threshold_colors() {
        let mut t = parse_newick("((A:1,B:2):1,(C:1,D:1):1);").unwrap();
        let a = t.find_label("A").unwrap();
        let ab = t.nodes[a].parent.unwrap();
        let cd = t.nodes[t.find_label("C").unwrap()].parent.unwrap();
        t.nodes[ab].attrs.insert("posterior".into(), Attr::Num(0.3));
        t.nodes[cd].attrs.insert("posterior".into(), Attr::Num(0.9));
        let layers = vec![LayerEntry::new(Layer::default_branch_lengths()), LayerEntry::new(Layer::default_node_labels("posterior"))];
        let view = ViewState::default();
        let s = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        let texts: Vec<(&str, Color)> = s.prims.iter().filter_map(|p| if let Prim::Text { text, color, .. } = p { Some((text.as_str(), *color)) } else { None }).collect();
        // Six branches, each labelled with its length.
        assert_eq!(texts.iter().filter(|(x, _)| ["1", "2"].contains(x)).count(), 6);
        assert!(texts.contains(&("0.3", Color::hex("#C62828"))));
        assert!(texts.contains(&("0.9", Color::hex("#2E7D32"))));
    }

    #[test]
    fn densitree_and_geoscale() {
        let main = parse_newick("((A:10,B:10):20,C:30);").unwrap();
        let samples = crate::io::newick::parse_newick_multi("((A:12,B:12):20,C:32);((A:8,C:8):25,B:33);").unwrap();
        let mut layers = default_layers();
        layers.insert(0, LayerEntry::new(Layer::default_densitree()));
        layers.push(LayerEntry::new(Layer::default_geoscale()));
        let view = ViewState::default();
        let s = build(&BuildInput { tree: &main, view: &view, layers: &layers, overlay: &samples }, &ApproxEnv, 800.0, 500.0, 1.0);
        // 2 sample trees x 4 branches, drawn translucent.
        let faint = s.prims.iter().filter(|p| matches!(p, Prim::Path { color, .. } if color.a == 25)).count();
        assert_eq!(faint, 8);
        // Root at 33 Ma reaches back into the Paleogene (Oligocene / Eocene).
        let fills: Vec<Color> = s.prims.iter().filter_map(|p| if let Prim::Poly { fill, .. } = p { Some(*fill) } else { None }).collect();
        assert!(fills.contains(&Color::hex("#FD9A52")), "Paleogene");
        assert!(fills.contains(&Color::hex("#FFFF00")), "Miocene");
    }

    #[test]
    fn labels_fit_inside_canvas() {
        let t = parse_newick("((A_very_long_taxon_name:1,B:2):1,C:3);").unwrap();
        let view = ViewState::default();
        let layers = default_layers();
        let s = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
        for (_, b) in &s.label_boxes {
            assert!(b[2] <= 400.0 + 12.0, "{:?}", b);
        }
    }

    #[test]
    fn node_bars_follow_common_ancestor_heights() {
        // TreeAnnotator -heights ca: the AB node sits at CAheight_mean = 3, while
        // height / height_95%_HPD (from trees containing the clade) say ~1.
        let t = parse_newick(
            "((A:3,B:3)[&height=1,height_95%_HPD={0.8,1.2},CAheight_mean=3,CAheight_95%_HPD={2.5,3.5}]:1,C:4)[&height=4,height_95%_HPD={3.9,4.1}];",
        )
        .unwrap();
        let mut layers = default_layers();
        layers.push(LayerEntry::new(Layer::default_node_bars()));
        let view = ViewState::default();
        let s = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
        let ab = t.nodes[t.find_label("A").unwrap()].parent.unwrap();
        let p = s.node_pos[ab].unwrap();
        let bar = s
            .prims
            .iter()
            .find_map(|q| match q {
                Prim::Path { pts, width, .. } if *width == 5.0 && (pts[0][1] - p[1]).abs() < 0.5 => Some((pts[0][0], pts[1][0])),
                _ => None,
            })
            .expect("no bar at the AB node");
        assert!(bar.0 < p[0] && p[0] < bar.1, "bar {:?} misses node at x = {}", bar, p[0]);
    }

    #[test]
    fn tip_points_move_out_and_labels_follow() {
        let t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        let a = t.find_label("A").unwrap();
        let mut layers = default_layers();
        layers.push(LayerEntry::new(Layer::default_tip_points()));
        let view = ViewState::default();
        let build_with = |offset: f32| {
            let mut ls = layers.clone();
            if let Some(LayerEntry { layer: Layer::TipPoints(p), .. }) = ls.last_mut() {
                p.offset = offset;
            }
            build(&BuildInput { tree: &t, view: &view, layers: &ls, overlay: &[] }, &ApproxEnv, 500.0, 300.0, 1.0)
        };
        // Label boxes carry a hit margin of 0.6 × the text size (11 pt by default).
        let label_x = |s: &Scene| s.label_boxes.iter().map(|(_, b)| b[0] + 0.6 * 11.0).fold(f32::MAX, f32::min);
        let s1 = build_with(20.0);
        let tip = s1.node_pos[a].unwrap();
        // The point sits 20 pt beyond its tip, and the labels start beyond the point.
        let p = s1.prims.iter().find_map(|q| if let Prim::Circle { c, .. } = q { ((c[1] - tip[1]).abs() < 0.5).then_some(c[0]) } else { None }).unwrap();
        assert!((p - tip[0] - 20.0).abs() < 0.5, "point at {} for tip at {}", p, tip[0]);
        assert!(label_x(&s1) > p, "labels {} should start past the points {}", label_x(&s1), p);
    }

    #[test]
    fn tip_labels_first_word() {
        assert_eq!(first_word("Homo_sapiens"), "Homo");
        assert_eq!(first_word(" Pan troglodytes"), "Pan");
        assert_eq!(first_word("_Gorilla_gorilla"), "Gorilla");
        assert_eq!(first_word("Pongo"), "Pongo");
        let t = parse_newick("(Homo_sapiens:1,Pan_troglodytes:1);").unwrap();
        let mut layers = default_layers();
        for e in layers.iter_mut() {
            if let Layer::TipLabels(s) = &mut e.layer {
                s.first_word = true;
            }
        }
        let view = ViewState::default();
        let s = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 500.0, 300.0, 1.0);
        let texts: Vec<&str> = s.prims.iter().filter_map(|p| if let Prim::Text { text, .. } = p { Some(text.as_str()) } else { None }).collect();
        assert!(texts.contains(&"Homo") && texts.contains(&"Pan") && !texts.iter().any(|t| t.contains("sapiens")), "{:?}", texts);
    }

    #[test]
    fn bars_take_negative_values_and_a_scale() {
        let mut t = parse_newick("((A:1,B:1):1,C:2);").unwrap();
        for (l, v) in [("A", 4.0), ("B", -2.0), ("C", 1.0)] {
            let n = t.find_label(l).unwrap();
            t.nodes[n].attrs.insert("x".into(), Attr::Num(v));
        }
        let mut layers = default_layers();
        layers.push(LayerEntry::new(Layer::default_bars("x")));
        let view = ViewState::default();
        let s = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 500.0, 300.0, 1.0);
        let Layer::Bars(b) = Layer::default_bars("x") else { unreachable!() };
        let x_range = |c: Color| -> Vec<(f32, f32)> {
            s.prims
                .iter()
                .filter_map(|p| match p {
                    Prim::Poly { pts, fill, stroke: None } if *fill == c => Some(pts.iter().fold((f32::MAX, f32::MIN), |(a, z), q| (a.min(q[0]), z.max(q[0])))),
                    _ => None,
                })
                .collect()
        };
        let pos = x_range(b.color);
        let neg = x_range(b.negative_color);
        assert_eq!((pos.len(), neg.len()), (2, 1));
        // The negative bar ends where the positive ones start (the zero line).
        assert!((neg[0].1 - pos[0].0).abs() < 0.01 && (neg[0].1 - pos[1].0).abs() < 0.01, "{:?} {:?}", neg, pos);
        // A 4-unit bar is twice as long as the 2-unit negative one.
        let len = |r: (f32, f32)| r.1 - r.0;
        assert!((len(*pos.iter().max_by(|a, b| len(**a).partial_cmp(&len(**b)).unwrap()).unwrap()) - 2.0 * len(neg[0])).abs() < 0.5);
        // The scale: tick labels and the column name.
        let texts: Vec<&str> = s.prims.iter().filter_map(|p| if let Prim::Text { text, .. } = p { Some(text.as_str()) } else { None }).collect();
        assert!(texts.contains(&"x") && texts.contains(&"0") && texts.iter().any(|t| t.starts_with('-')), "{:?}", texts);
    }

    #[test]
    fn fan_fills_the_canvas() {
        let t = parse_newick("(((A:1,B:1):1,C:2):1,(D:2,E:2):1);").unwrap();
        let layers = default_layers();
        let mut view = ViewState::default();
        view.layout.kind = LayoutKind::Fan;
        let s = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 600.0, 400.0, 1.0);
        let root = s.node_pos[t.root].unwrap();
        let tips: Vec<P> = t.tips().iter().map(|&n| s.node_pos[n].unwrap()).collect();
        // A 180° fan sits above its root and spreads across most of the width.
        assert!(tips.iter().all(|p| p[1] <= root[1] + 1.0), "{:?} {:?}", root, tips);
        let (lo, hi) = tips.iter().fold((f32::MAX, f32::MIN), |(lo, hi), p| (lo.min(p[0]), hi.max(p[0])));
        assert!(hi - lo > 200.0, "fan too narrow: {}", hi - lo);
        for (_, b) in &s.label_boxes {
            assert!(b[0] >= -1.0 && b[1] >= -1.0 && b[2] <= 601.0 && b[3] <= 401.0, "label outside: {:?}", b);
        }
    }

    #[test]
    fn tips_face_the_chosen_side() {
        let t = parse_newick("((A:1,B:2):1,C:3);").unwrap();
        let layers = default_layers();
        let a = t.find_label("A").unwrap();
        for side in TipSide::ALL {
            let mut view = ViewState::default();
            view.layout.tips = side;
            let s = build(&BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
            let (r, tip) = (s.node_pos[t.root].unwrap(), s.node_pos[a].unwrap());
            let ok = match side {
                TipSide::Right => tip[0] > r[0],
                TipSide::Left => tip[0] < r[0],
                TipSide::Top => tip[1] < r[1],
                TipSide::Bottom => tip[1] > r[1],
            };
            assert!(ok, "{:?}: root {:?}, tip {:?}", side, r, tip);
            // The timescale draws its colored interval boxes on every side.
            let mut with_geo = layers.clone();
            with_geo.push(LayerEntry::new(Layer::default_geoscale()));
            let g = build(&BuildInput { tree: &t, view: &view, layers: &with_geo, overlay: &[] }, &ApproxEnv, 400.0, 300.0, 1.0);
            let boxes = g.prims.iter().filter(|p| matches!(p, Prim::Poly { stroke: Some(_), .. })).count();
            assert!(boxes > 0, "{:?}: no timescale drawn", side);
            for (_, b) in &s.label_boxes {
                assert!(b[0] >= -12.0 && b[1] >= -12.0 && b[2] <= 412.0 && b[3] <= 312.0, "{:?} label outside: {:?}", side, b);
            }
        }
    }
}
