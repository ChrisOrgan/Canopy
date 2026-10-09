//! Visual styling: colors, palettes, and ggtree-style layers.

use crate::layout::LayoutOptions;
use crate::tree::{Attr, NodeId, Tree};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}

impl Color {
    pub const BLACK: Color = Color::rgb(0, 0, 0);
    pub const WHITE: Color = Color::rgb(255, 255, 255);
    pub const GREY: Color = Color::rgb(128, 128, 128);
    pub const LIGHT_GREY: Color = Color::rgb(200, 200, 200);
    pub const TRANSPARENT: Color = Color { r: 0, g: 0, b: 0, a: 0 };

    pub const fn rgb(r: u8, g: u8, b: u8) -> Color {
        Color { r, g, b, a: 255 }
    }

    pub fn with_alpha(self, a: u8) -> Color {
        Color { a, ..self }
    }

    pub fn hex(s: &str) -> Color {
        let s = s.trim_start_matches('#');
        let p = |i: usize| u8::from_str_radix(s.get(i..i + 2).unwrap_or("00"), 16).unwrap_or(0);
        let a = if s.len() >= 8 { p(6) } else { 255 };
        Color { r: p(0), g: p(2), b: p(4), a }
    }

    pub fn to_hex(self) -> String {
        format!("#{:02X}{:02X}{:02X}", self.r, self.g, self.b)
    }

    pub fn lerp(self, o: Color, t: f32) -> Color {
        let t = t.clamp(0.0, 1.0);
        let f = |a: u8, b: u8| (a as f32 + (b as f32 - a as f32) * t).round() as u8;
        Color { r: f(self.r, o.r), g: f(self.g, o.g), b: f(self.b, o.b), a: f(self.a, o.a) }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Palette {
    OkabeIto,
    Set1,
    Dark2,
    Viridis,
    Magma,
    BlueRed,
    Grey,
}

impl Palette {
    pub const ALL: [Palette; 7] =
        [Palette::OkabeIto, Palette::Set1, Palette::Dark2, Palette::Viridis, Palette::Magma, Palette::BlueRed, Palette::Grey];

    pub fn name(self) -> &'static str {
        match self {
            Palette::OkabeIto => "Okabe-Ito",
            Palette::Set1 => "Set1",
            Palette::Dark2 => "Dark2",
            Palette::Viridis => "Viridis",
            Palette::Magma => "Magma",
            Palette::BlueRed => "Blue-Red",
            Palette::Grey => "Greys",
        }
    }

    fn stops(self) -> Vec<Color> {
        let h = |v: &[&str]| v.iter().map(|s| Color::hex(s)).collect();
        match self {
            Palette::OkabeIto => h(&["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]),
            Palette::Set1 => h(&["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"]),
            Palette::Dark2 => h(&["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"]),
            Palette::Viridis => h(&["#440154", "#482878", "#3E4A89", "#31688E", "#26828E", "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725"]),
            Palette::Magma => h(&["#000004", "#1C1044", "#4F127B", "#812581", "#B5367A", "#E55064", "#FB8761", "#FEC287", "#FCFDBF"]),
            Palette::BlueRed => h(&["#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B"]),
            Palette::Grey => h(&["#F0F0F0", "#000000"]),
        }
    }

    /// Color for category `i`.
    pub fn discrete(self, i: usize, n: usize) -> Color {
        let s = self.stops();
        match self {
            Palette::OkabeIto | Palette::Set1 | Palette::Dark2 => s[i % s.len()],
            _ => self.continuous(if n <= 1 { 0.5 } else { i as f32 / (n - 1) as f32 }),
        }
    }

    /// Color at position `t` in [0,1].
    pub fn continuous(self, t: f32) -> Color {
        let s = match self {
            // Qualitative palettes fall back to viridis for gradients.
            Palette::OkabeIto | Palette::Set1 | Palette::Dark2 => Palette::Viridis.stops(),
            _ => self.stops(),
        };
        let t = t.clamp(0.0, 1.0) * (s.len() - 1) as f32;
        let i = (t.floor() as usize).min(s.len() - 2);
        s[i].lerp(s[i + 1], t - i as f32)
    }
}

/// A resolved mapping from an attribute's values to colors.
#[derive(Clone, Debug)]
pub enum ColorMap {
    Discrete { title: String, levels: Vec<(String, Color)> },
    Continuous { title: String, min: f64, max: f64, palette: Palette },
}

impl ColorMap {
    /// Build a mapping for `key` over the given nodes. Numeric attributes give
    /// a gradient; everything else gives discrete levels.
    pub fn build(tree: &Tree, key: &str, nodes: &[NodeId], palette: Palette) -> Option<ColorMap> {
        let vals: Vec<Attr> = nodes.iter().filter_map(|&n| tree.value(n, key)).collect();
        if vals.is_empty() {
            return None;
        }
        let numeric = vals.iter().all(|v| matches!(v, Attr::Num(_)));
        if numeric {
            let nums: Vec<f64> = vals.iter().filter_map(Attr::as_f64).filter(|v| v.is_finite()).collect();
            let min = nums.iter().cloned().fold(f64::INFINITY, f64::min);
            let max = nums.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            if min.is_finite() {
                let pal = match palette {
                    Palette::OkabeIto | Palette::Set1 | Palette::Dark2 => Palette::Viridis,
                    p => p,
                };
                return Some(ColorMap::Continuous { title: key.to_string(), min, max, palette: pal });
            }
        }
        let mut levels: Vec<String> = vals.iter().map(|v| v.display(4)).collect();
        levels.sort();
        levels.dedup();
        let n = levels.len();
        Some(ColorMap::Discrete {
            title: key.to_string(),
            levels: levels.into_iter().enumerate().map(|(i, l)| (l, palette.discrete(i, n))).collect(),
        })
    }

    pub fn color(&self, v: &Attr) -> Option<Color> {
        match self {
            ColorMap::Discrete { levels, .. } => {
                let s = v.display(4);
                levels.iter().find(|(l, _)| *l == s).map(|(_, c)| *c)
            }
            ColorMap::Continuous { min, max, palette, .. } => {
                let x = v.as_f64()?;
                if !x.is_finite() {
                    return None;
                }
                let t = if max > min { (x - min) / (max - min) } else { 0.5 };
                Some(palette.continuous(t as f32))
            }
        }
    }

    pub fn title(&self) -> &str {
        match self {
            ColorMap::Discrete { title, .. } | ColorMap::Continuous { title, .. } => title,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Shape {
    Circle,
    Square,
    Triangle,
    Diamond,
}

impl Shape {
    pub const ALL: [Shape; 4] = [Shape::Circle, Shape::Square, Shape::Triangle, Shape::Diamond];
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TreeStyle {
    pub width: f32,
    pub color: Color,
    pub color_by: Option<String>,
    pub palette: Palette,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TipLabelStyle {
    pub size: f32,
    pub color: Color,
    pub italic: bool,
    pub align: bool,
    pub offset: f32,
    pub underscores_as_spaces: bool,
    pub color_by: Option<String>,
    pub palette: Palette,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NodeLabelStyle {
    /// Attribute to display ("label", "support", "posterior", ...).
    pub attr: String,
    pub size: f32,
    pub color: Color,
    pub digits: usize,
    /// Only show values >= this (numeric attributes).
    pub min_value: Option<f64>,
    /// Draw on the branch (above it) instead of beside the node.
    pub on_branch: bool,
    /// Color numeric values by a cut-off: (threshold, color below, color at or above).
    #[serde(default)]
    pub threshold: Option<(f64, Color, Color)>,
}

/// Overlay of posterior sample trees (DensiTree).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DensiTreeStyle {
    /// Maximum number of trees drawn (evenly spaced through the post-burn-in sample).
    pub max_trees: usize,
    pub color: Color,
    /// Opacity of each tree, 1–255.
    pub alpha: u8,
    pub width: f32,
    /// Color the most frequent topologies blue, red and green, as in DensiTree.
    pub by_topology: bool,
}

/// Geologic timescale under a time-calibrated tree (deeptime `coord_geo`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GeoscaleStyle {
    pub eras: bool,
    pub periods: bool,
    pub epochs: bool,
    /// Height of each row, in points.
    pub row_height: f32,
    pub labels: bool,
    pub label_size: f32,
    /// Faint lines across the plot at boundaries.
    pub boundaries: bool,
    /// Millions of years per branch-length unit.
    pub ma_per_unit: f64,
    /// Age of the youngest tip, in Ma.
    pub youngest_age: f64,
}

/// Branch-length values printed along each branch.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BranchLengthStyle {
    pub size: f32,
    pub color: Color,
    pub digits: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PointStyle {
    pub size: f32,
    pub shape: Shape,
    pub color: Color,
    pub color_by: Option<String>,
    pub palette: Palette,
    /// Only draw where attribute `filter.0` >= `filter.1` (e.g. posterior >= 0.95).
    pub filter: Option<(String, f64)>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RangeStyle {
    pub attr: String,
    pub color: Color,
    pub width: f32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HighlightStyle {
    pub node: NodeId,
    pub fill: Color,
    /// Extend the highlight past the tips (points).
    pub extend: f32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CladeLabelStyle {
    pub node: NodeId,
    pub text: String,
    pub color: Color,
    pub size: f32,
    pub offset: f32,
    pub bar_width: f32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ScaleBarStyle {
    pub length: Option<f64>,
    pub width: f32,
    pub size: f32,
    pub color: Color,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AxisStyle {
    /// Show time before present (root at left = oldest).
    pub time_before_present: bool,
    pub title: String,
    pub size: f32,
    pub grid: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HeatmapStyle {
    pub columns: Vec<String>,
    pub offset: f32,
    pub cell_width: f32,
    pub palette: Palette,
    pub show_names: bool,
    pub name_size: f32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BarStyle {
    pub column: String,
    pub offset: f32,
    pub max_width: f32,
    pub color: Color,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PhylopicStyle {
    pub size: f32,
    pub offset: f32,
    pub tint: Option<Color>,
    pub align: bool,
    /// Shrink silhouettes so they never exceed the gap between tips.
    #[serde(default = "default_true")]
    pub fit_spacing: bool,
}

fn default_true() -> bool {
    true
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Layer {
    Tree(TreeStyle),
    TipLabels(TipLabelStyle),
    NodeLabels(NodeLabelStyle),
    BranchLengths(BranchLengthStyle),
    DensiTree(DensiTreeStyle),
    Geoscale(GeoscaleStyle),
    TipPoints(PointStyle),
    NodePoints(PointStyle),
    NodeBars(RangeStyle),
    Highlight(HighlightStyle),
    CladeLabel(CladeLabelStyle),
    ScaleBar(ScaleBarStyle),
    TimeAxis(AxisStyle),
    Heatmap(HeatmapStyle),
    Bars(BarStyle),
    Phylopic(PhylopicStyle),
}

impl Layer {
    /// ggtree function this layer corresponds to.
    pub fn name(&self) -> &'static str {
        match self {
            Layer::Tree(_) => "geom_tree",
            Layer::TipLabels(_) => "geom_tiplab",
            Layer::NodeLabels(_) => "geom_nodelab",
            Layer::BranchLengths(_) => "geom_text (branch.length)",
            Layer::DensiTree(_) => "densiTree (phangorn)",
            Layer::Geoscale(_) => "coord_geo (deeptime)",
            Layer::TipPoints(_) => "geom_tippoint",
            Layer::NodePoints(_) => "geom_nodepoint",
            Layer::NodeBars(_) => "geom_range",
            Layer::Highlight(_) => "geom_hilight",
            Layer::CladeLabel(_) => "geom_cladelab",
            Layer::ScaleBar(_) => "geom_treescale",
            Layer::TimeAxis(_) => "theme_tree2",
            Layer::Heatmap(_) => "gheatmap",
            Layer::Bars(_) => "geom_facet (bar)",
            Layer::Phylopic(_) => "geom_phylopic",
        }
    }

    pub fn description(&self) -> String {
        match self {
            Layer::Tree(_) => "Branches".into(),
            Layer::TipLabels(_) => "Tip labels".into(),
            Layer::NodeLabels(s) => format!("Node labels ({})", s.attr),
            Layer::BranchLengths(_) => "Branch lengths".into(),
            Layer::DensiTree(d) => format!("DensiTree ({} trees)", d.max_trees),
            Layer::Geoscale(_) => "Geologic timescale".into(),
            Layer::TipPoints(_) => "Tip points".into(),
            Layer::NodePoints(_) => "Node points".into(),
            Layer::NodeBars(s) => format!("Node bars ({})", s.attr),
            Layer::Highlight(h) => format!("Highlight clade #{}", h.node),
            Layer::CladeLabel(c) => format!("Clade label \"{}\"", c.text),
            Layer::ScaleBar(_) => "Scale bar".into(),
            Layer::TimeAxis(_) => "Time axis".into(),
            Layer::Heatmap(h) => format!("Heatmap ({} cols)", h.columns.len()),
            Layer::Bars(b) => format!("Bar chart ({})", b.column),
            Layer::Phylopic(_) => "PhyloPic silhouettes".into(),
        }
    }

    pub fn default_tree() -> Layer {
        Layer::Tree(TreeStyle { width: 1.0, color: Color::BLACK, color_by: None, palette: Palette::Viridis })
    }
    pub fn default_tip_labels() -> Layer {
        Layer::TipLabels(TipLabelStyle {
            size: 11.0,
            color: Color::BLACK,
            italic: true,
            align: false,
            offset: 4.0,
            underscores_as_spaces: true,
            color_by: None,
            palette: Palette::OkabeIto,
        })
    }
    pub fn default_node_labels(attr: &str) -> Layer {
        Layer::NodeLabels(NodeLabelStyle {
            attr: attr.into(),
            size: 8.0,
            color: Color::rgb(80, 80, 80),
            digits: 2,
            min_value: None,
            on_branch: false,
            // Posterior support: red below 0.5, green at or above.
            threshold: (attr == "posterior").then_some((0.5, Color::hex("#C62828"), Color::hex("#2E7D32"))),
        })
    }
    pub fn default_densitree() -> Layer {
        Layer::DensiTree(DensiTreeStyle { max_trees: 100, color: Color::rgb(30, 60, 120), alpha: 25, width: 1.0, by_topology: false })
    }
    pub fn default_geoscale() -> Layer {
        Layer::Geoscale(GeoscaleStyle {
            eras: false,
            periods: true,
            epochs: true,
            row_height: 12.0,
            labels: true,
            label_size: 7.0,
            boundaries: false,
            ma_per_unit: 1.0,
            youngest_age: 0.0,
        })
    }
    pub fn default_branch_lengths() -> Layer {
        Layer::BranchLengths(BranchLengthStyle { size: 7.0, color: Color::rgb(90, 90, 90), digits: 3 })
    }
    pub fn default_tip_points() -> Layer {
        Layer::TipPoints(PointStyle {
            size: 4.0,
            shape: Shape::Circle,
            color: Color::rgb(0, 114, 178),
            color_by: None,
            palette: Palette::OkabeIto,
            filter: None,
        })
    }
    pub fn default_node_points() -> Layer {
        Layer::NodePoints(PointStyle {
            size: 4.0,
            shape: Shape::Circle,
            color: Color::BLACK,
            color_by: None,
            palette: Palette::Viridis,
            filter: None,
        })
    }
    pub fn default_node_bars() -> Layer {
        Layer::NodeBars(RangeStyle { attr: "height_95%_HPD".into(), color: Color::rgb(86, 180, 233).with_alpha(160), width: 5.0 })
    }
    pub fn default_scale_bar() -> Layer {
        Layer::ScaleBar(ScaleBarStyle { length: None, width: 1.0, size: 9.0, color: Color::BLACK })
    }
    pub fn default_axis() -> Layer {
        Layer::TimeAxis(AxisStyle { time_before_present: true, title: String::new(), size: 9.0, grid: false })
    }
    pub fn default_phylopic() -> Layer {
        Layer::Phylopic(PhylopicStyle { size: 18.0, offset: 4.0, tint: None, align: true, fit_spacing: true })
    }
}

/// A layer plus its visibility toggle.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LayerEntry {
    pub enabled: bool,
    pub layer: Layer,
}

impl LayerEntry {
    pub fn new(layer: Layer) -> Self {
        LayerEntry { enabled: true, layer }
    }
}

pub fn default_layers() -> Vec<LayerEntry> {
    vec![LayerEntry::new(Layer::default_tree()), LayerEntry::new(Layer::default_tip_labels()), LayerEntry::new(Layer::default_scale_bar())]
}

/// Per-document view state that is not part of the tree itself.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ViewState {
    pub layout: LayoutOptions,
    /// Collapsed clades, drawn as triangles.
    pub collapsed: HashSet<NodeId>,
    /// Branch color overrides applied to a node's whole subtree (ggtree `groupClade`).
    pub clade_colors: HashMap<NodeId, Color>,
    /// Custom silhouette images assigned to tips (label -> image file path).
    pub tip_images: BTreeMap<String, String>,
    pub background: Option<Color>,
    pub title: String,
}

impl ViewState {
    /// Resolve the clade-color override in effect at every node.
    pub fn resolved_clade_colors(&self, tree: &Tree) -> Vec<Option<Color>> {
        let mut out = vec![None; tree.nodes.len()];
        for n in tree.preorder() {
            out[n] = self.clade_colors.get(&n).copied().or_else(|| tree.nodes[n].parent.and_then(|p| out[p]));
        }
        out
    }
}
