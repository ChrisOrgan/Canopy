//! Scene building: turns a tree, its view state and a stack of ggtree-style
//! layers into renderer-independent drawing primitives in pixel space.
//! The same scene feeds the interactive canvas and the PNG/TIFF/SVG exporters,
//! so what you see is what you export.

use crate::layout::{self, Layout, LayoutKind};
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
    /// A raster image (PhyloPic silhouette) looked up by key at render time.
    Image { key: String, center: P, height: f32, tint: Option<Color> },
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

fn tip_text(tree: &Tree, lay: &Layout, n: NodeId, underscores: bool) -> String {
    let base = if lay.collapsed_depth[n].is_some() {
        match tree.nodes[n].label.as_deref().filter(|l| l.parse::<f64>().is_err()) {
            Some(l) => l.to_string(),
            None => format!("{} tips", tree.tips_below(n).len()),
        }
    } else {
        tree.label(n).to_string()
    };
    if underscores {
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
    let layers: Vec<&Layer> = inp.layers.iter().filter(|l| l.enabled).map(|l| &l.layer).collect();
    // Make room for node-height intervals that reach back past the root.
    for l in &layers {
        if let Layer::NodeBars(s) = l {
            if s.attr.contains("height") {
                for n in tree.preorder() {
                    if let Some((_, hi)) = tree.value(n, &s.attr).and_then(|v| v.as_range()) {
                        lay.min_x = lay.min_x.min(lay.max_x - hi);
                    }
                }
            }
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
            .map(|&n| env.text_width(&tip_text(tree, &lay, n, s.underscores_as_spaces), s.size * pt, s.italic))
            .collect(),
        None => vec![0.0; lay.leaves.len()],
    };
    let max_label_w = label_widths.iter().cloned().fold(0.0, f32::max);

    // ---- Columns beyond the tips (labels, silhouettes, heatmaps, bars, clade labels).
    let mut starts: Vec<f32> = vec![0.0; layers.len()];
    let mut cursor = 0.0f32;
    let mut clade_ext = 0.0f32;
    let mut has_extras = false;
    for (i, l) in layers.iter().enumerate() {
        match l {
            Layer::TipLabels(s) => {
                starts[i] = s.offset * pt;
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
                cursor = starts[i] + s.max_width * pt + 4.0 * pt;
            }
            Layer::CladeLabel(s) => {
                has_extras = true;
                starts[i] = cursor + s.offset * pt;
                let w = env.text_width(&s.text, s.size * pt, false);
                clade_ext = clade_ext.max(starts[i] + s.bar_width * pt + 4.0 * pt + w);
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
            Layer::Tree(s) => s.color_by.as_ref().and_then(|k| ColorMap::build(tree, k, &visible_nodes, s.palette)),
            Layer::TipLabels(s) => s.color_by.as_ref().and_then(|k| ColorMap::build(tree, k, &lay.leaves, s.palette)),
            Layer::TipPoints(s) => s.color_by.as_ref().and_then(|k| ColorMap::build(tree, k, &lay.leaves, s.palette)),
            Layer::NodePoints(s) => s.color_by.as_ref().and_then(|k| ColorMap::build(tree, k, &visible_nodes, s.palette)),
            Layer::Heatmap(s) => heatmap_map(tree, &lay.leaves, s),
            _ => None,
        };
        add_legend(&m, &mut legends);
        maps.push(m);
    }
    let legend_w = if legends.is_empty() {
        0.0
    } else {
        let mut w = 0.0f32;
        for m in &legends {
            w = w.max(env.text_width(m.title(), 9.0 * pt, false));
            match m {
                ColorMap::Discrete { levels, .. } => {
                    for (l, _) in levels.iter().take(30) {
                        w = w.max(14.0 * pt + env.text_width(l, 8.0 * pt, false));
                    }
                }
                ColorMap::Continuous { min, max, .. } => {
                    w = w.max(18.0 * pt + env.text_width(&format_num(*min, 3), 8.0 * pt, false));
                    w = w.max(18.0 * pt + env.text_width(&format_num(*max, 3), 8.0 * pt, false));
                }
            }
        }
        w + 16.0 * pt
    };

    let title_h = if view.title.is_empty() { 0.0 } else { 20.0 * pt };
    let heat_names_h = layers
        .iter()
        .filter_map(|l| if let Layer::Heatmap(h) = l { Some(h) } else { None })
        .filter(|h| h.show_names)
        .map(|h| h.columns.iter().map(|c| env.text_width(c, h.name_size * pt, false)).fold(0.0, f32::max) + 4.0 * pt)
        .fold(0.0, f32::max);
    let has_scale = layers.iter().any(|l| matches!(l, Layer::ScaleBar(_)));
    let has_axis = layers.iter().any(|l| matches!(l, Layer::TimeAxis(_)));
    let top = pad + title_h;
    // Geologic timescale rows sit under the tree (horizontal layouts only).
    let horizontal = matches!(opts.kind, LayoutKind::Rectangular | LayoutKind::Slanted | LayoutKind::Roundrect | LayoutKind::Ellipse);
    let geo_h = layers
        .iter()
        .find_map(|l| if let Layer::Geoscale(g) = l { Some(g) } else { None })
        .filter(|_| horizontal)
        .map(|g| [g.eras, g.periods, g.epochs].iter().filter(|b| **b).count() as f32 * g.row_height * pt + 6.0 * pt)
        .unwrap_or(0.0);
    let bottom = pad + if has_scale { 22.0 * pt } else { 0.0 } + if has_axis { 30.0 * pt } else { 0.0 } + geo_h;

    let n_slots = lay.n_slots();
    let range = (lay.max_x - lay.min_x).max(1e-12);
    let unaligned_labels = tiplab.as_ref().map(|s| !s.align).unwrap_or(false) && !has_extras;
    let label_off = tiplab.as_ref().map(|s| s.offset * pt).unwrap_or(0.0);
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
            let (theta0, step) = match open {
                Some(gap) => {
                    let sweep = (360.0 - gap).to_radians();
                    (opts.rotate as f64 * PI / 180.0 + gap.to_radians() / 2.0, if n_slots > 1 { sweep / (n_slots - 1) as f64 } else { 0.0 })
                }
                None => (opts.rotate as f64 * PI / 180.0, TAU / n_slots as f64),
            };
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
            } else if opts.flip_x {
                Orient::Left
            } else {
                Orient::Right
            };
            let (u0, u_end, v0, v_end) = match orient {
                Orient::Down => (top + heat_names_h.min(0.0), height - bottom, pad, plot_w - pad),
                _ => (pad, plot_w - pad, top + heat_names_h, height - bottom),
            };
            let s = fit_scale(u_end - u0);
            let sp = ((v_end - v0) / n_slots as f32).max(0.5);
            (GeoKind::Lin { orient, w: plot_w, u0, v0, sp, flip_y: opts.flip_y }, s)
        }
    };
    let geo = Geo { kind: kind.0, s: kind.1, min_x: lay.min_x, max_x: lay.max_x, n: n_slots };

    let mut node_pos = vec![None; tree.nodes.len()];
    for &n in &visible_nodes {
        node_pos[n] = Some(geo.at(lay.x[n], lay.y[n], n));
    }
    let scene = Scene {
        width,
        height,
        background: view.background,
        node_pos,
        px_per_unit: geo.s,
        ..Default::default()
    };
    let clade_colors = view.resolved_clade_colors(tree);
    let mut cx = Ctx { tree, lay, geo, pt, env, scene, clade_colors };
    if let Some(bg) = view.background {
        cx.scene.prims.push(Prim::Poly { pts: vec![[0.0, 0.0], [width, 0.0], [width, height], [0.0, height]], fill: bg, stroke: None });
    }

    // Branch-length-based decorations make no sense on a cladogram.
    let lengths_shown = opts.use_lengths && tree.has_lengths() && !(opts.kind == LayoutKind::Slanted && opts.slanted_cladogram);
    // The timescale always sits behind the tree.
    for l in &layers {
        if let Layer::Geoscale(s) = l {
            if lengths_shown {
                draw_geoscale(&mut cx, s, has_axis);
            }
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
            Layer::DensiTree(s) => draw_densitree(&mut cx, s, opts.kind, inp.overlay),
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
            Layer::ScaleBar(s) if lengths_shown => draw_scale_bar(&mut cx, s, bottom, has_axis, geo_h),
            Layer::ScaleBar(_) => {}
            Layer::TimeAxis(s) if lengths_shown => draw_axis(&mut cx, s),
            Layer::TimeAxis(_) => {}
            Layer::Heatmap(s) => draw_heatmap(&mut cx, s, map, start),
            Layer::Bars(s) => draw_bars(&mut cx, s, start),
            Layer::Phylopic(s) => draw_phylopic(&mut cx, s, start, label_off, &label_widths, tiplab.as_ref().map(|t| t.align).unwrap_or(false)),
        }
    }

    if !view.title.is_empty() {
        cx.text([plot_w / 2.0, pad], view.title.clone(), 14.0 * pt, Color::BLACK, 0.0, 0.5, 0.0, false, true);
    }
    draw_legends(&mut cx, &legends, width - legend_w + 6.0 * pt, top);
    cx.scene
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
        let palette = match s.palette {
            Palette::OkabeIto | Palette::Set1 | Palette::Dark2 => Palette::Viridis,
            p => p,
        };
        return min.is_finite().then_some(ColorMap::Continuous { title, min, max, palette });
    }
    let mut levels: Vec<String> = vals.iter().map(|v| v.display(4)).collect();
    levels.sort();
    levels.dedup();
    let n = levels.len();
    Some(ColorMap::Discrete { title, levels: levels.into_iter().enumerate().map(|(i, l)| (l, s.palette.discrete(i, n))).collect() })
}

fn branch_color(cx: &Ctx, n: NodeId, s: &TreeStyle, map: Option<&ColorMap>) -> Color {
    if let Some(c) = cx.clade_colors[n] {
        return c;
    }
    if let (Some(m), Some(k)) = (map, &s.color_by) {
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

/// Geologic timescale: colored interval bars under the tree (horizontal
/// layouts) or translucent rings behind it (circular layouts).
fn draw_geoscale(cx: &mut Ctx, s: &GeoscaleStyle, has_axis: bool) {
    use crate::geotime::{overlapping, Level};
    let pt = cx.pt;
    let mpu = if s.ma_per_unit > 0.0 { s.ma_per_unit } else { 1.0 };
    let max_x = cx.lay.max_x;
    let x_of = |age: f64| max_x - (age - s.youngest_age) / mpu;
    let young = s.youngest_age;
    let old = s.youngest_age + (max_x - cx.lay.min_x) * mpu;
    // Finest level first: nearest the tree.
    let rows: Vec<Level> = [(s.epochs, Level::Epoch), (s.periods, Level::Period), (s.eras, Level::Era)].iter().filter(|(on, _)| *on).map(|(_, l)| *l).collect();
    let Some(&finest) = rows.first() else { return };
    match cx.geo.kind {
        GeoKind::Lin { orient: Orient::Right | Orient::Left, v0, sp, .. } => {
            let n = cx.geo.n as f32;
            let v_tree_end = v0 + n * sp;
            let mut v = v_tree_end + 6.0 * pt + if has_axis { 30.0 * pt } else { 0.0 };
            let rh = s.row_height * pt;
            if s.boundaries {
                let lines: Vec<f64> = overlapping(*rows.last().unwrap(), young, old).map(|i| i.start).filter(|a| *a < old).collect();
                for a in lines {
                    let u = cx.geo.u(x_of(a));
                    let pts = vec![cx.geo.uv(u, v0), cx.geo.uv(u, v_tree_end)];
                    cx.scene.prims.push(Prim::Path { pts, color: Color::rgb(150, 150, 150).with_alpha(140), width: 0.6 * pt, dashed: true });
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
                            cx.text(p, t.to_string(), size, Color::BLACK, 0.0, 0.5, 0.5, false, false);
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
            let axis = (t1 + t0 + TAU) / 2.0;
            let size = s.label_size * pt;
            for iv in overlapping(finest, young, old) {
                let (a0, a1) = (iv.start.min(old), iv.end.max(young));
                let (r_in, r_out) = (cx.geo.r(x_of(a0)), cx.geo.r(x_of(a1)));
                let region = cx.geo.sector(r_in, r_out, 0.0, TAU);
                cx.scene.prims.push(region.prim(Color::hex(iv.color), Some((Color::WHITE, 0.6 * pt))));
                // Name along the ring at 12 o'clock, if the ring is thick enough.
                if s.labels && (r_out - r_in).abs() >= size * 1.1 {
                    let mid_r = (r_in + r_out) / 2.0;
                    let room = mid_r * 1.2;
                    if let Some(t) = [iv.name, iv.abbr].into_iter().find(|t| cx.env.text_width(t, size, false) <= room) {
                        let p = cx.geo.polar(mid_r, FRAC_PI_2);
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
    let width = s.width * cx.pt;
    let tree = cx.tree;
    for &c in visible {
        let color = branch_color(cx, c, s, map);
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
        cx.scene.prims.push(Prim::Path { pts, color, width, dashed: false });

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
        let text = tip_text(cx.tree, &cx.lay, n, s.underscores_as_spaces);
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
            Attr::Num(x) => format_num(*x, s.digits),
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
                    Orient::Down => {
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

/// Branch-length values at the middle of each branch.
fn draw_branch_lengths(cx: &mut Ctx, s: &BranchLengthStyle, visible: &[NodeId]) {
    let size = s.size * cx.pt;
    let gap = 1.5 * cx.pt;
    for &n in visible {
        let (Some(p), Some(len)) = (cx.tree.nodes[n].parent, cx.tree.nodes[n].length) else { continue };
        let text = format_num(len, s.digits);
        let (xp, xc, y) = (cx.lay.x[p], cx.lay.x[n], cx.lay.y[n]);
        let (pos, angle, halign) = match &cx.geo.kind {
            GeoKind::Lin { orient: Orient::Down, .. } => {
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
        let p = if tips { cx.geo.at(cx.leaf_x(n), cx.lay.y[n], n) } else { cx.pos(n) };
        cx.scene.prims.push(shape_prim(s.shape, p, s.size * cx.pt, color));
    }
}

fn draw_ranges(cx: &mut Ctx, s: &RangeStyle, visible: &[NodeId]) {
    if matches!(cx.geo.kind, GeoKind::Free { .. }) {
        return;
    }
    let width = s.width * cx.pt;
    for &n in visible {
        let Some((lo, hi)) = cx.tree.value(n, &s.attr).and_then(|v| v.as_range()) else { continue };
        // Heights are measured back from the youngest tip.
        let (x1, x2) = if s.attr.contains("height") { (cx.lay.max_x - hi, cx.lay.max_x - lo) } else { (lo, hi) };
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
        }
    }
}

fn draw_scale_bar(cx: &mut Ctx, s: &ScaleBarStyle, bottom: f32, has_axis: bool, geo_h: f32) {
    let range = cx.lay.max_x - cx.lay.min_x;
    let len = s.length.filter(|l| *l > 0.0).unwrap_or_else(|| nice_step(range, 5));
    let px = len as f32 * cx.geo.s;
    let pt = cx.pt;
    let y = cx.scene.height - bottom + 10.0 * pt + if has_axis { 30.0 * pt } else { 0.0 } + geo_h + 8.0 * pt;
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
    cx.scene.prims.push(Prim::Path { pts: line, color: Color::BLACK, width: 0.8 * pt, dashed: false });
    // Ticks are placed on round values of the displayed quantity.
    let disp = |x: f64| if s.time_before_present { x1 - x } else { x };
    let (d0, d1) = (disp(x0).min(disp(x1)), disp(x0).max(disp(x1)));
    let mut t = (d0 / step).ceil() * step;
    while t <= d1 + step * 1e-6 {
        let x = if s.time_before_present { x1 - t } else { t };
        let u = cx.geo.u(x);
        if s.grid {
            let g = vec![cx.geo.uv(u, v0), cx.geo.uv(u, v_axis)];
            cx.scene.prims.push(Prim::Path { pts: g, color: Color::rgb(225, 225, 225), width: 0.6 * pt, dashed: false });
        }
        let tick = vec![cx.geo.uv(u, v_axis), cx.geo.uv(u, v_axis + 4.0 * pt)];
        cx.scene.prims.push(Prim::Path { pts: tick, color: Color::BLACK, width: 0.8 * pt, dashed: false });
        let p = cx.geo.uv(u, v_axis + 6.0 * pt);
        let label = format_num(t.abs(), 4);
        match orient {
            Orient::Down => cx.text(p, label, s.size * pt, Color::BLACK, 0.0, 0.0, 0.5, false, false),
            _ => cx.text(p, label, s.size * pt, Color::BLACK, 0.0, 0.5, 0.0, false, false),
        }
        t += step;
    }
    if !s.title.is_empty() {
        let um = (cx.geo.u(x0) + cx.geo.u(x1)) / 2.0;
        let p = cx.geo.uv(um, v_axis + 18.0 * pt);
        cx.text(p, s.title.clone(), s.size * pt, Color::BLACK, 0.0, 0.5, 0.0, false, false);
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
    let leaves = cx.lay.leaves.clone();
    let vals: Vec<Option<f64>> = leaves.iter().map(|&n| cx.tree.value(n, &s.column).and_then(|v| v.as_f64()).filter(|v| v.is_finite())).collect();
    let max = vals.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
    if max <= 0.0 {
        return;
    }
    let free = matches!(cx.geo.kind, GeoKind::Free { .. });
    for (i, &n) in leaves.iter().enumerate() {
        let Some(v) = vals[i] else { continue };
        let w = (v.abs() / max) as f32 * s.max_width * cx.pt;
        let y = cx.lay.y[n];
        let base = if free { cx.leaf_x(n) } else { cx.lay.max_x };
        let sp = cx.geo.spacing(base);
        let shrink = 0.15 * sp;
        // Narrow the cell to leave gaps between bars.
        let mut region = cx.geo.cell(y, y, base, start, start + w, n);
        if let (GeoKind::Lin { orient, .. }, Region::Poly(pts)) = (&cx.geo.kind, &mut region) {
            for (j, p) in pts.iter_mut().enumerate() {
                let sign = if j < 2 { 1.0 } else { -1.0 };
                match orient {
                    Orient::Down => p[0] += sign * shrink,
                    _ => p[1] += sign * shrink,
                }
            }
        }
        cx.scene.prims.push(region.prim(s.color, None));
    }
    let tp = cx.geo.ext(-1.0, cx.lay.max_x, start, 0);
    if let GeoKind::Lin { orient: Orient::Right, .. } = cx.geo.kind {
        cx.text([tp[0], tp[1] + cx.geo.spacing(0.0) * 0.5], format!("{} (max {})", s.column, format_num(max, 3)), 8.0 * cx.pt, Color::BLACK, 0.0, 0.0, 1.0, false, false);
    }
}

fn draw_phylopic(cx: &mut Ctx, s: &PhylopicStyle, start: f32, label_off: f32, widths: &[f32], labels_aligned: bool) {
    // Never taller than the space between neighbouring tips.
    let room = match cx.geo.kind {
        GeoKind::Polar { step, .. } => ((cx.geo.r(cx.lay.max_x) + start) * step as f32).abs(),
        _ => cx.geo.spacing(cx.lay.max_x),
    };
    let h = if s.fit_spacing { (s.size * cx.pt).min(room * 0.92).max(2.0) } else { s.size * cx.pt };
    let leaves = cx.lay.leaves.clone();
    let free = matches!(cx.geo.kind, GeoKind::Free { .. });
    for (i, &n) in leaves.iter().enumerate() {
        if !cx.tree.is_tip(n) {
            continue;
        }
        let key = image_key(cx.tree.label(n));
        cx.scene.wanted_images.push(key.clone());
        let Some(aspect) = cx.env.image_aspect(&key) else { continue };
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
        cx.scene.prims.push(Prim::Image { key, center, height: h, tint: s.tint });
    }
}

fn draw_legends(cx: &mut Ctx, legends: &[ColorMap], x: f32, top: f32) {
    let pt = cx.pt;
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
            ColorMap::Continuous { min, max, palette, .. } => {
                let h = 70.0 * pt;
                let steps = 30;
                for i in 0..steps {
                    let t0 = i as f32 / steps as f32;
                    let t1 = (i + 1) as f32 / steps as f32;
                    // Top of the bar is the maximum.
                    let c = palette.continuous(1.0 - (t0 + t1) / 2.0);
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick;

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
}
