//! Visual styling: colors, palettes, and ggtree-style layers.

use crate::layout::LayoutOptions;
use crate::tree::{Attr, NodeId, Tree};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};

/// An sRGB color with alpha. Saved as "#RRGGBB", or "#RRGGBBAA" when not opaque.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}

impl Serialize for Color {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let alpha = if self.a == 255 { String::new() } else { format!("{:02X}", self.a) };
        s.serialize_str(&format!("{}{}", self.to_hex(), alpha))
    }
}

impl<'de> Deserialize<'de> for Color {
    /// Reads "#RRGGBB[AA]", or the {r, g, b, a} form of earlier saves.
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Repr {
            Hex(String),
            Rgba { r: u8, g: u8, b: u8, a: u8 },
        }
        Ok(match Repr::deserialize(d)? {
            Repr::Hex(s) => Color::hex(&s),
            Repr::Rgba { r, g, b, a } => Color { r, g, b, a },
        })
    }
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

/// Color palettes, grouped by kind. Variant names are stored in saved
/// projects, so existing ones keep their names (`BlueRed` is RdBu).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Palette {
    OkabeIto,
    Set1,
    Dark2,
    Viridis,
    Magma,
    BlueRed,
    Grey,
    Plasma,
    Cividis,
    BlackGreen,
    BlackRed,
    BlackBlue,
    PuOr,
    BrBG,
    PiYG,
    RdYlBu,
    Set2,
    Tableau10,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PaletteKind {
    Sequential,
    Diverging,
    Qualitative,
}

impl PaletteKind {
    pub fn name(self) -> &'static str {
        match self {
            PaletteKind::Sequential => "Sequential",
            PaletteKind::Diverging => "Diverging",
            PaletteKind::Qualitative => "Qualitative",
        }
    }
}

impl Palette {
    /// Every palette, in menu order (sequential, diverging, qualitative).
    pub const ALL: [Palette; 18] = [
        Palette::Viridis,
        Palette::Magma,
        Palette::Plasma,
        Palette::Cividis,
        Palette::Grey,
        Palette::BlackGreen,
        Palette::BlackRed,
        Palette::BlackBlue,
        Palette::BlueRed,
        Palette::PuOr,
        Palette::BrBG,
        Palette::PiYG,
        Palette::RdYlBu,
        Palette::OkabeIto,
        Palette::Set1,
        Palette::Dark2,
        Palette::Set2,
        Palette::Tableau10,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Palette::OkabeIto => "Okabe-Ito",
            Palette::Set1 => "Set1",
            Palette::Dark2 => "Dark2",
            Palette::Set2 => "Set2",
            Palette::Tableau10 => "Tableau 10",
            Palette::Viridis => "Viridis",
            Palette::Magma => "Magma",
            Palette::Plasma => "Plasma",
            Palette::Cividis => "Cividis",
            Palette::Grey => "Greys",
            Palette::BlackGreen => "Black–Green",
            Palette::BlackRed => "Black–Red",
            Palette::BlackBlue => "Black–Blue",
            Palette::BlueRed => "Blue–Red (RdBu)",
            Palette::PuOr => "Orange–Purple (PuOr)",
            Palette::BrBG => "Brown–Teal (BrBG)",
            Palette::PiYG => "Pink–Green (PiYG)",
            Palette::RdYlBu => "Red–Yellow–Blue (RdYlBu)",
        }
    }

    pub fn kind(self) -> PaletteKind {
        match self {
            Palette::OkabeIto | Palette::Set1 | Palette::Dark2 | Palette::Set2 | Palette::Tableau10 => PaletteKind::Qualitative,
            Palette::BlueRed | Palette::PuOr | Palette::BrBG | Palette::PiYG | Palette::RdYlBu => PaletteKind::Diverging,
            _ => PaletteKind::Sequential,
        }
    }

    fn stops(self) -> Vec<Color> {
        let h = |v: &[&str]| v.iter().map(|s| Color::hex(s)).collect();
        match self {
            // Okabe & Ito (2008); ColorBrewer (Brewer); Tableau 10.
            Palette::OkabeIto => h(&["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]),
            Palette::Set1 => h(&["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"]),
            Palette::Dark2 => h(&["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"]),
            Palette::Set2 => h(&["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"]),
            Palette::Tableau10 => h(&["#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"]),
            // Matplotlib's perceptually uniform maps (van der Walt & Smith; Nuñez et al. for cividis).
            Palette::Viridis => h(&["#440154", "#482878", "#3E4A89", "#31688E", "#26828E", "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725"]),
            Palette::Magma => h(&["#000004", "#1C1044", "#4F127B", "#812581", "#B5367A", "#E55064", "#FB8761", "#FEC287", "#FCFDBF"]),
            Palette::Plasma => h(&["#0D0887", "#46039F", "#7201A8", "#9C179E", "#BD3786", "#D8576B", "#ED7953", "#FB9F3A", "#FDCA26", "#F0F921"]),
            Palette::Cividis => h(&["#00204D", "#00336F", "#39486B", "#575C6D", "#707173", "#8A8779", "#A69D75", "#C4B56C", "#E4CF5B", "#FFEA46"]),
            Palette::Grey => h(&["#F0F0F0", "#000000"]),
            Palette::BlackGreen => h(&["#000000", "#00FF00"]),
            Palette::BlackRed => h(&["#000000", "#FF0000"]),
            Palette::BlackBlue => h(&["#000000", "#0066FF"]),
            // Diverging, low values first, through black in the middle (ColorBrewer
            // colors for the arms, as in dark-centred maps such as Crameri's "berlin").
            Palette::BlueRed => h(&["#92C5DE", "#4393C3", "#000000", "#D6604D", "#F4A582"]),
            Palette::PuOr => h(&["#FDB863", "#E08214", "#000000", "#8073AC", "#B2ABD2"]),
            Palette::BrBG => h(&["#DFC27D", "#BF812D", "#000000", "#35978F", "#80CDC1"]),
            Palette::PiYG => h(&["#F1B6DA", "#DE77AE", "#000000", "#7FBC41", "#B8E186"]),
            Palette::RdYlBu => h(&["#FDAE61", "#F46D43", "#000000", "#74ADD1", "#ABD9E9"]),
        }
    }

    /// Color for category `i` of `n`.
    pub fn discrete(self, i: usize, n: usize) -> Color {
        match self.kind() {
            PaletteKind::Qualitative => {
                let s = self.stops();
                s[i % s.len()]
            }
            _ => self.continuous(if n <= 1 { 0.5 } else { i as f32 / (n - 1) as f32 }),
        }
    }

    /// Color at position `t` in [0,1]. Qualitative palettes give their colors
    /// in equal steps (a binned scale) rather than blending them.
    pub fn continuous(self, t: f32) -> Color {
        let s = self.stops();
        let t = t.clamp(0.0, 1.0);
        if self.kind() == PaletteKind::Qualitative {
            return s[((t * s.len() as f32) as usize).min(s.len() - 1)];
        }
        let t = t * (s.len() - 1) as f32;
        let i = (t.floor() as usize).min(s.len() - 2);
        s[i].lerp(s[i + 1], t - i as f32)
    }
}

/// A resolved mapping from an attribute's values to colors.
#[derive(Clone, Debug)]
pub enum ColorMap {
    Discrete { title: String, levels: Vec<(String, Color)> },
    /// `reversed`: the palette runs from its last color (minimum) to its first.
    Continuous { title: String, min: f64, max: f64, palette: Palette, reversed: bool },
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
                return Some(ColorMap::Continuous { title: key.to_string(), min, max, palette, reversed: false });
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
            ColorMap::Continuous { min, max, palette, reversed, .. } => {
                let x = v.as_f64()?;
                if !x.is_finite() {
                    return None;
                }
                let t = if max > min { (x - min) / (max - min) } else { 0.5 };
                Some(palette.continuous(if *reversed { 1.0 - t as f32 } else { t as f32 }))
            }
        }
    }

    /// The same mapping with the palette's direction flipped if `flip`.
    pub fn reversed(mut self, flip: bool) -> Self {
        if flip {
            match &mut self {
                ColorMap::Continuous { reversed, .. } => *reversed = !*reversed,
                ColorMap::Discrete { levels, .. } => {
                    let colors: Vec<Color> = levels.iter().rev().map(|x| x.1).collect();
                    for (l, c) in levels.iter_mut().zip(colors) {
                        l.1 = c;
                    }
                }
            }
        }
        self
    }

    /// The same mapping under a different legend title.
    pub fn with_title(mut self, t: &str) -> Self {
        match &mut self {
            ColorMap::Discrete { title, .. } | ColorMap::Continuous { title, .. } => *title = t.to_string(),
        }
        self
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
    /// Run the palette backwards.
    #[serde(default)]
    pub reverse_palette: bool,
    /// Coloring by a BayesTraits branch value (rate scalar, Fabric beta): only
    /// branches whose posterior is at least this; the rest keep `color`.
    #[serde(default)]
    pub min_pp: f64,
    /// Color only branches whose value is at least this (the rest keep `color`).
    #[serde(default)]
    pub min_value: Option<f64>,
    /// Compare the value's size (absolute value) with `min_value`, for values of either sign.
    #[serde(default)]
    pub min_abs: bool,
    /// Line width by a numeric trait: from `width_min` (smallest value) to
    /// `width` (largest); branches without a value get `width_min`.
    #[serde(default)]
    pub width_by: Option<String>,
    #[serde(default = "branch_width_min")]
    pub width_min: f32,
}

fn branch_width_min() -> f32 {
    0.5
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TipLabelStyle {
    pub size: f32,
    pub color: Color,
    pub italic: bool,
    pub align: bool,
    pub offset: f32,
    pub underscores_as_spaces: bool,
    /// Show only the first word of each label (the text before the first
    /// underscore or space), e.g. the genus of "Homo_sapiens".
    #[serde(default)]
    pub first_word: bool,
    pub color_by: Option<String>,
    pub palette: Palette,
    /// Run the palette backwards.
    #[serde(default)]
    pub reverse_palette: bool,
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
    /// Text appended to numeric values (e.g. "%").
    #[serde(default)]
    pub suffix: String,
}

/// BayesTraits node scalars: circles at nodes whose node-scalar posterior (the
/// proportion of the sample with a node scalar placed there) is at least `threshold`.
/// Fossil (stratigraphic) ranges as wide bars along the branches of taxa with
/// `fossil_min_age` and `fossil_max_age` (Ma; placed with the timescale layer's
/// settings, else branch lengths in Ma with the youngest tip at the present).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FossilRangeStyle {
    /// Bar width, in points.
    pub width: f32,
    pub color: Color,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NodeScalarStyle {
    /// Minimum posterior (0–1) for a circle.
    pub threshold: f64,
    pub size: f32,
    pub palette: Palette,
    /// Run the palette backwards.
    #[serde(default)]
    pub reverse_palette: bool,
    /// One color for every circle instead of the palette (no color key).
    #[serde(default)]
    pub single_color: Option<Color>,
    /// Circle size from `size_min` (posterior 0) to `size` (posterior 1).
    #[serde(default)]
    pub size_by_pp: bool,
    #[serde(default = "scalar_min")]
    pub size_min: f32,
}

fn scalar_min() -> f32 {
    2.0
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
    /// Ages (stages), the finest level.
    #[serde(default)]
    pub ages: bool,
    /// Height of each row, in points.
    pub row_height: f32,
    pub labels: bool,
    pub label_size: f32,
    /// Faint lines across the plot at boundaries.
    pub boundaries: bool,
    /// Color and width (points) of the boundary lines.
    #[serde(default = "boundary_color")]
    pub boundary_color: Color,
    #[serde(default = "boundary_width")]
    pub boundary_width: f32,
    /// Millions of years per branch-length unit.
    pub ma_per_unit: f64,
    /// Age of the youngest tip, in Ma.
    pub youngest_age: f64,
}

/// Values printed along each branch: the branch length, or a per-branch
/// attribute such as a rate scalar.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BranchLengthStyle {
    pub size: f32,
    pub color: Color,
    pub digits: usize,
    /// Attribute to show (None = the branch length).
    #[serde(default)]
    pub attr: Option<String>,
    /// Don't label branches whose value equals this (e.g. 1 = unscaled).
    #[serde(default)]
    pub hide_value: Option<f64>,
    /// For a BayesTraits branch value: only label branches whose posterior is at least this.
    #[serde(default)]
    pub min_pp: f64,
    /// Only label branches whose value is at least this.
    #[serde(default)]
    pub min_value: Option<f64>,
    /// Compare the value's size (absolute value) with `min_value`.
    #[serde(default)]
    pub min_abs: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PointStyle {
    pub size: f32,
    pub shape: Shape,
    pub color: Color,
    pub color_by: Option<String>,
    pub palette: Palette,
    /// Run the palette backwards.
    #[serde(default)]
    pub reverse_palette: bool,
    /// Only draw where attribute `filter.0` >= `filter.1` (e.g. posterior >= 0.95).
    pub filter: Option<(String, f64)>,
    /// Tip points only: distance (pt) out from the tip; tip labels move out to make room.
    pub offset: f32,
    /// Legend title (default: the attribute name).
    #[serde(default)]
    pub legend_title: Option<String>,
    /// Scale the size by a numeric trait: `size_min` for the smallest value,
    /// `size` for the largest.
    #[serde(default)]
    pub size_by: Option<String>,
    #[serde(default = "point_min")]
    pub size_min: f32,
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
    /// A band of color behind the clade's branches, fading to transparent from the
    /// label back to the clade's node.
    #[serde(default)]
    pub gradient: bool,
    /// Color at the label end (its alpha is the strongest opacity).
    #[serde(default = "gradient_color")]
    pub gradient_color: Color,
    /// A PhyloPic silhouette after the text: the taxon looked up on phylopic.org
    /// (or the name an image file was loaded under).
    #[serde(default)]
    pub phylopic: Option<String>,
    /// Height of the silhouette, in points.
    #[serde(default = "clade_pic_size")]
    pub phylopic_size: f32,
}

fn clade_pic_size() -> f32 {
    24.0
}

fn gradient_color() -> Color {
    Color::hex("#56B4E9").with_alpha(200)
}

impl CladeLabelStyle {
    /// A clade label with the default look.
    pub fn new(node: NodeId, text: String) -> Self {
        CladeLabelStyle { node, text, color: Color::BLACK, size: 11.0, offset: 6.0, bar_width: 2.0, gradient: false, gradient_color: gradient_color(), phylopic: None, phylopic_size: clade_pic_size() }
    }
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
    /// Age of the youngest (deepest) tip, counting back (e.g. a tree of fossils).
    #[serde(default)]
    pub youngest_age: f64,
    pub size: f32,
    pub grid: bool,
    /// Color and width (points) of the axis line and ticks.
    #[serde(default = "axis_color")]
    pub line_color: Color,
    #[serde(default = "axis_width")]
    pub line_width: f32,
    /// Color and width (points) of the grid lines.
    #[serde(default = "grid_color")]
    pub grid_color: Color,
    #[serde(default = "grid_width")]
    pub grid_width: f32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HeatmapStyle {
    pub columns: Vec<String>,
    pub offset: f32,
    pub cell_width: f32,
    pub palette: Palette,
    /// Run the palette backwards.
    #[serde(default)]
    pub reverse_palette: bool,
    pub show_names: bool,
    pub name_size: f32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BarStyle {
    pub column: String,
    pub offset: f32,
    /// Width of the whole value range, from the smallest (or 0) to the largest (or 0).
    pub max_width: f32,
    pub color: Color,
    /// Fill for bars below zero.
    pub negative_color: Color,
    /// Axis with tick values under the bars.
    pub show_scale: bool,
    /// The scale's line width (pt) and color.
    #[serde(default = "bar_scale_width")]
    pub scale_width: f32,
    #[serde(default = "black")]
    pub scale_color: Color,
    /// Size (pt) and color of the tick values and the column name.
    #[serde(default = "bar_label_size")]
    pub label_size: f32,
    #[serde(default = "black")]
    pub label_color: Color,
}

fn bar_scale_width() -> f32 {
    0.6
}
fn bar_label_size() -> f32 {
    7.0
}
fn black() -> Color {
    Color::BLACK
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
    /// Draw silhouettes at tips, and at nodes given one (the `phylopic` attribute).
    #[serde(default = "default_true")]
    pub tips: bool,
    #[serde(default = "default_true")]
    pub nodes: bool,
    /// Recolor by a trait (instead of `tint`).
    #[serde(default)]
    pub color_by: Option<String>,
    #[serde(default = "viridis")]
    pub palette: Palette,
    #[serde(default)]
    pub reverse_palette: bool,
    /// Scale the height by a numeric trait: `size_min` for the smallest value,
    /// `size` for the largest.
    #[serde(default)]
    pub size_by: Option<String>,
    #[serde(default = "phylopic_min")]
    pub size_min: f32,
}

fn point_min() -> f32 {
    1.5
}
fn viridis() -> Palette {
    Palette::Viridis
}
fn phylopic_min() -> f32 {
    6.0
}
fn boundary_color() -> Color {
    Color::rgb(150, 150, 150).with_alpha(140)
}
fn boundary_width() -> f32 {
    0.6
}
fn axis_color() -> Color {
    Color::BLACK
}
fn axis_width() -> f32 {
    0.8
}
fn grid_color() -> Color {
    Color::rgb(225, 225, 225)
}
fn grid_width() -> f32 {
    0.6
}
fn title_size() -> f32 {
    14.0
}
fn one() -> f32 {
    1.0
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
    NodeScalars(NodeScalarStyle),
    FossilRanges(FossilRangeStyle),
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
    /// Name shown in the layer list.
    pub fn description(&self) -> String {
        match self {
            Layer::Tree(_) => "Branches".into(),
            Layer::TipLabels(_) => "Tip labels".into(),
            Layer::NodeLabels(s) => format!("Node labels ({})", s.attr),
            Layer::BranchLengths(s) => match &s.attr {
                Some(a) => format!("Branch labels ({})", a),
                None => "Branch lengths".into(),
            },
            Layer::DensiTree(_) => "DensiTree".into(),
            Layer::NodeScalars(s) => format!("Node scalars (posterior ≥ {})", crate::tree::format_num(s.threshold, 2)),
            Layer::Geoscale(_) => "Geologic timescale".into(),
            Layer::FossilRanges(_) => "Fossil ranges".into(),
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
        Layer::Tree(TreeStyle { width: 1.0, color: Color::BLACK, color_by: None, palette: Palette::Viridis, reverse_palette: false, min_pp: 0.0, min_value: None, min_abs: false, width_by: None, width_min: branch_width_min() })
    }
    pub fn default_tip_labels() -> Layer {
        Layer::TipLabels(TipLabelStyle {
            size: 11.0,
            color: Color::BLACK,
            italic: true,
            align: false,
            offset: 4.0,
            underscores_as_spaces: true,
            first_word: false,
            color_by: None,
            palette: Palette::OkabeIto,
            reverse_palette: false,
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
            threshold: matches!(attr, "posterior" | "bootstrap").then_some((0.5, Color::hex("#C62828"), Color::hex("#2E7D32"))),
            suffix: String::new(),
        })
    }
    pub fn default_fossil_ranges() -> Layer {
        Layer::FossilRanges(FossilRangeStyle { width: 6.0, color: Color::rgb(70, 70, 70).with_alpha(200) })
    }
    pub fn default_node_scalars() -> Layer {
        // Blue (low) to red (high): distinct from the viridis scale used for branch rates.
        Layer::NodeScalars(NodeScalarStyle { threshold: 0.0, size: 8.0, palette: Palette::BlueRed, reverse_palette: false, single_color: None, size_by_pp: false, size_min: scalar_min() })
    }
    pub fn default_densitree() -> Layer {
        Layer::DensiTree(DensiTreeStyle { max_trees: 100, color: Color::rgb(30, 60, 120), alpha: 25, width: 1.0, by_topology: false })
    }
    pub fn default_geoscale() -> Layer {
        Layer::Geoscale(GeoscaleStyle {
            eras: false,
            periods: true,
            epochs: true,
            ages: false,
            row_height: 12.0,
            labels: true,
            label_size: 7.0,
            boundaries: false,
            boundary_color: boundary_color(),
            boundary_width: boundary_width(),
            ma_per_unit: 1.0,
            youngest_age: 0.0,
        })
    }
    pub fn default_branch_lengths() -> Layer {
        Layer::BranchLengths(BranchLengthStyle { size: 7.0, color: Color::rgb(90, 90, 90), digits: 3, attr: None, hide_value: None, min_pp: 0.0, min_value: None, min_abs: false })
    }
    pub fn default_tip_points() -> Layer {
        Layer::TipPoints(PointStyle {
            size: 4.0,
            shape: Shape::Circle,
            color: Color::rgb(0, 114, 178),
            color_by: None,
            palette: Palette::OkabeIto,
            reverse_palette: false,
            filter: None,
            offset: 0.0,
            legend_title: None,
            size_by: None,
            size_min: point_min(),
        })
    }
    pub fn default_node_points() -> Layer {
        Layer::NodePoints(PointStyle {
            size: 4.0,
            shape: Shape::Circle,
            color: Color::BLACK,
            color_by: None,
            palette: Palette::Viridis,
            reverse_palette: false,
            filter: None,
            offset: 0.0,
            legend_title: None,
            size_by: None,
            size_min: point_min(),
        })
    }
    pub fn default_node_bars() -> Layer {
        Layer::NodeBars(RangeStyle { attr: "height_95%_HPD".into(), color: Color::rgb(86, 180, 233).with_alpha(160), width: 5.0 })
    }
    pub fn default_bars(column: &str) -> Layer {
        Layer::Bars(BarStyle {
            column: column.into(),
            offset: 4.0,
            max_width: 60.0,
            color: Color::hex("#0072B2"),
            negative_color: Color::hex("#D55E00"),
            show_scale: true,
            scale_width: bar_scale_width(),
            scale_color: Color::BLACK,
            label_size: bar_label_size(),
            label_color: Color::BLACK,
        })
    }
    pub fn default_scale_bar() -> Layer {
        Layer::ScaleBar(ScaleBarStyle { length: None, width: 1.0, size: 9.0, color: Color::BLACK })
    }
    /// A timescale or time axis whose youngest tip age is `age` (other layers unchanged).
    pub fn with_youngest_age(mut self, age: f64) -> Layer {
        match &mut self {
            Layer::Geoscale(g) => g.youngest_age = age,
            Layer::TimeAxis(a) => a.youngest_age = age,
            _ => {}
        }
        self
    }

    pub fn default_axis() -> Layer {
        Layer::TimeAxis(AxisStyle {
            time_before_present: true,
            youngest_age: 0.0,
            size: 9.0,
            grid: false,
            line_color: axis_color(),
            line_width: axis_width(),
            grid_color: grid_color(),
            grid_width: grid_width(),
        })
    }
    pub fn default_phylopic() -> Layer {
        Layer::Phylopic(PhylopicStyle {
            size: 18.0,
            offset: 4.0,
            tint: None,
            align: true,
            fit_spacing: true,
            tips: true,
            nodes: true,
            color_by: None,
            palette: Palette::Viridis,
            reverse_palette: false,
            size_by: None,
            size_min: phylopic_min(),
        })
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

/// Time frame for placing ages (Ma) on the tree: (age of the youngest tip, Ma per
/// branch-length unit), from the geologic timescale layer if there is one, else
/// a time axis counting back from its youngest tip age, else the tree's own
/// youngest tip `height` (fossil tips), else (0, 1).
pub fn time_frame<'a>(layers: impl Iterator<Item = &'a Layer> + Clone, tree: &Tree) -> (f64, f64) {
    layers
        .clone()
        .find_map(|l| if let Layer::Geoscale(g) = l { Some((g.youngest_age, g.ma_per_unit)) } else { None })
        .or_else(|| layers.into_iter().find_map(|l| match l {
            Layer::TimeAxis(a) if a.time_before_present => Some((a.youngest_age, 1.0)),
            _ => None,
        }))
        .unwrap_or((tree.youngest_tip_age().unwrap_or(0.0), 1.0))
}

pub fn default_layers() -> Vec<LayerEntry> {
    vec![LayerEntry::new(Layer::default_tree()), LayerEntry::new(Layer::default_tip_labels()), LayerEntry::new(Layer::default_scale_bar())]
}

/// Per-document view state that is not part of the tree itself.
#[derive(Clone, Debug, Serialize, Deserialize)]
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
    /// Size of the figure title, in points.
    #[serde(default = "title_size")]
    pub title_size: f32,
    /// Size of the color keys (legends), relative to the default.
    #[serde(default = "one")]
    pub legend_scale: f32,
}

impl Default for ViewState {
    fn default() -> Self {
        ViewState {
            layout: Default::default(),
            collapsed: Default::default(),
            clade_colors: Default::default(),
            tip_images: Default::default(),
            background: None,
            title: String::new(),
            title_size: title_size(),
            legend_scale: 1.0,
        }
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reversed_palettes() {
        let mut t = crate::io::newick::parse_newick("((A,B),C);").unwrap();
        for (name, x, s) in [("A", 0.0, "x"), ("B", 5.0, "y"), ("C", 10.0, "z")] {
            let n = t.find_label(name).unwrap();
            t.nodes[n].attrs.insert("v".into(), Attr::Num(x));
            t.nodes[n].attrs.insert("s".into(), Attr::Text(s.into()));
        }
        let tips = t.tips();
        let m = ColorMap::build(&t, "v", &tips, Palette::Viridis).unwrap();
        let r = m.clone().reversed(true);
        // The minimum takes the palette's last color, the maximum its first.
        assert_eq!(r.color(&Attr::Num(0.0)), m.color(&Attr::Num(10.0)));
        assert_eq!(r.color(&Attr::Num(10.0)), Some(Palette::Viridis.continuous(0.0)));
        let d = ColorMap::build(&t, "s", &tips, Palette::Set1).unwrap();
        let dr = d.clone().reversed(true);
        assert_eq!(dr.color(&Attr::Text("x".into())), d.color(&Attr::Text("z".into())));
        assert_eq!(dr.color(&Attr::Text("y".into())), d.color(&Attr::Text("y".into())));
        assert_eq!(m.clone().reversed(false).color(&Attr::Num(0.0)), m.color(&Attr::Num(0.0)));
    }

    #[test]
    fn colors_save_as_hex() {
        let c = Color::rgb(70, 70, 70).with_alpha(200);
        assert_eq!(serde_json::to_string(&c).unwrap(), "\"#464646C8\"");
        assert_eq!(serde_json::to_string(&Color::rgb(255, 0, 0)).unwrap(), "\"#FF0000\"");
        assert_eq!(serde_json::from_str::<Color>("\"#464646C8\"").unwrap(), c);
        // The {r, g, b, a} form of earlier saves still reads.
        assert_eq!(serde_json::from_str::<Color>(r#"{"r":70,"g":70,"b":70,"a":200}"#).unwrap(), c);
    }

    #[test]
    fn palette_groups() {
        let count = |k: PaletteKind| Palette::ALL.iter().filter(|p| p.kind() == k).count();
        assert_eq!((count(PaletteKind::Sequential), count(PaletteKind::Diverging), count(PaletteKind::Qualitative)), (8, 5, 5));
        // Qualitative palettes color numbers in steps of their own colors, not viridis.
        assert_eq!(Palette::Set1.continuous(0.0), Color::hex("#E41A1C"));
        assert_eq!(Palette::Set1.continuous(1.0), Color::hex("#999999"));
        assert_eq!(Palette::BlackRed.continuous(1.0), Color::hex("#FF0000"));
        // Diverging palettes pass through black in the middle.
        for p in Palette::ALL.iter().filter(|p| p.kind() == PaletteKind::Diverging) {
            assert_eq!(p.continuous(0.5), Color::BLACK, "{}", p.name());
        }
    }
}
