//! Side and bottom panels: layer stack (drag to reorder), layout options,
//! layer property editors, node inspector, comparative data, and the
//! posterior tree-set tools.

use super::canvas::{c32, Action};
use super::document::{Document, SampleKind};
use crate::ancestral;
use crate::consensus::HeightMode;
use crate::io::data::ColumnKind;
use crate::layout::{LayoutKind, TipSide};
use crate::style::*;
use crate::tree::{format_num, Attr, NodeId};
use egui::{Color32, ComboBox, DragValue, Slider, Stroke};

pub fn color_edit(ui: &mut egui::Ui, label: &str, c: &mut Color) {
    ui.horizontal(|ui| {
        let mut col = c32(*c);
        if egui::color_picker::color_edit_button_srgba(ui, &mut col, egui::color_picker::Alpha::OnlyBlend).changed() {
            let [r, g, b, a] = col.to_srgba_unmultiplied();
            *c = Color { r, g, b, a };
        }
        ui.label(label);
    });
}

fn attr_select(ui: &mut egui::Ui, id: &str, label: &str, value: &mut Option<String>, keys: &[String]) {
    ui.horizontal(|ui| {
        ui.label(label);
        ComboBox::from_id_salt(id).selected_text(value.clone().unwrap_or_else(|| "(none)".into())).show_ui(ui, |ui| {
            ui.selectable_value(value, None, "(none)");
            for k in keys {
                ui.selectable_value(value, Some(k.clone()), k);
            }
        });
    });
}

fn attr_select_req(ui: &mut egui::Ui, id: &str, label: &str, value: &mut String, keys: &[String]) {
    ui.horizontal(|ui| {
        ui.label(label);
        ComboBox::from_id_salt(id).selected_text(value.clone()).show_ui(ui, |ui| {
            for k in keys {
                ui.selectable_value(value, k.clone(), k);
            }
        });
    });
}

fn palette_select(ui: &mut egui::Ui, id: &str, p: &mut Palette, reverse: &mut bool) {
    ui.horizontal(|ui| {
        ui.label("Palette");
        palette_combo(ui, id, p);
        ui.checkbox(reverse, "reverse").on_hover_text("Run the colors in the opposite direction");
    });
}

/// Palette menu grouped into sequential, diverging and qualitative. Returns
/// true if the choice changed.
pub(crate) fn palette_combo(ui: &mut egui::Ui, id: &str, p: &mut Palette) -> bool {
    let before = *p;
    ComboBox::from_id_salt(id).selected_text(p.name()).width(190.0).show_ui(ui, |ui| {
        let mut last = None;
        for q in Palette::ALL {
            if last != Some(q.kind()) {
                if last.is_some() {
                    ui.separator();
                }
                ui.weak(q.kind().name());
                last = Some(q.kind());
            }
            ui.selectable_value(p, q, q.name());
        }
    });
    *p != before
}

fn num(ui: &mut egui::Ui, label: &str, v: &mut f32, range: std::ops::RangeInclusive<f32>) {
    ui.add(Slider::new(v, range).text(label));
}

fn shape_select(ui: &mut egui::Ui, id: &str, s: &mut Shape) {
    ui.horizontal(|ui| {
        ui.label("Shape");
        ComboBox::from_id_salt(id).selected_text(format!("{:?}", s)).show_ui(ui, |ui| {
            for q in Shape::ALL {
                ui.selectable_value(s, q, format!("{:?}", q));
            }
        });
    });
}

fn keys_with_label(doc: &Document) -> Vec<String> {
    let mut k = vec!["label".to_string()];
    k.extend(doc.tree.attr_keys());
    k
}

/// The layer stack. Rows can be dragged to change drawing order (bottom = drawn last).
pub fn layers_panel(ui: &mut egui::Ui, doc: &mut Document) {
    ui.horizontal(|ui| {
        ui.heading("Layers");
        ui.menu_button("➕ Add", |ui| add_layer_menu(ui, doc));
    });
    ui.weak("Drag ☰ to reorder; layers draw top to bottom.");
    let mut from = None;
    let mut to = None;
    let mut delete = None;
    let n = doc.layers.len();
    let available = doc.overlay().len();
    for i in 0..n {
        let selected = doc.selected_layer == Some(i);
        let row = ui
            .horizontal(|ui| {
                ui.dnd_drag_source(egui::Id::new(("layer-drag", i)), i, |ui| {
                    ui.label("☰");
                });
                ui.checkbox(&mut doc.layers[i].enabled, "");
                let mut text = doc.layers[i].layer.description();
                if let Layer::DensiTree(d) = &doc.layers[i].layer {
                    // Trees actually drawn: never more than the post-burn-in sample.
                    text = format!("{} ({} trees)", text, d.max_trees.min(available));
                }
                text.push_str("  ");
                if ui.selectable_label(selected, text).clicked() {
                    doc.selected_layer = if selected { None } else { Some(i) };
                }
                if ui.small_button("🗑").on_hover_text("Remove layer").clicked() {
                    delete = Some(i);
                }
            })
            .response;
        if let (Some(ptr), Some(_)) = (ui.input(|inp| inp.pointer.interact_pos()), row.dnd_hover_payload::<usize>()) {
            let before = ptr.y < row.rect.center().y;
            let y = if before { row.rect.top() } else { row.rect.bottom() };
            ui.painter().hline(row.rect.x_range(), y, Stroke::new(2.0f32, Color32::from_rgb(0, 114, 178)));
            if let Some(dragged) = row.dnd_release_payload::<usize>() {
                from = Some(*dragged);
                to = Some(if before { i } else { i + 1 });
            }
        }
    }
    if let (Some(f), Some(t)) = (from, to) {
        if f != t && f + 1 != t {
            doc.checkpoint();
            let item = doc.layers.remove(f);
            let t = if t > f { t - 1 } else { t };
            doc.layers.insert(t, item);
            doc.selected_layer = Some(t);
        }
    }
    if let Some(i) = delete {
        doc.checkpoint();
        doc.layers.remove(i);
        doc.selected_layer = None;
    }
}

/// "Only values ≥" checkbox and slider spanning `attr`'s values on the tree's
/// branches; with values of both signs (Fabric betas), an "absolute value"
/// option compares each value's size instead.
fn min_value_slider(ui: &mut egui::Ui, tree: &crate::tree::Tree, attr: &str, min: &mut Option<f64>, abs: &mut bool, hover: &str) {
    let vals: Vec<f64> = tree.preorder().into_iter().filter(|&n| tree.nodes[n].parent.is_some()).filter_map(|n| tree.nodes[n].attrs.get(attr).and_then(Attr::as_f64)).filter(|x| x.is_finite()).collect();
    let (lo, hi) = vals.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &x| (a.0.min(x), a.1.max(x)));
    if !(hi > lo) {
        return;
    }
    let signed = lo < 0.0 && hi > 0.0;
    if !signed {
        *abs = false;
    }
    let (lo, hi) = if *abs { (0.0, lo.abs().max(hi.abs())) } else { (lo, hi) };
    let mut on = min.is_some();
    let mut v = min.unwrap_or(lo).clamp(lo, hi);
    ui.horizontal(|ui| {
        ui.checkbox(&mut on, if *abs { "only |values| ≥" } else { "only values ≥" }).on_hover_text(hover);
        ui.add_enabled(on, Slider::new(&mut v, lo..=hi).max_decimals(3));
    });
    if signed {
        let was = *abs;
        ui.checkbox(abs, "Absolute value").on_hover_text("Compare each value's size, so large negative values (e.g. Fabric betas) count as well as large positive ones");
        if *abs != was {
            v = v.abs();
        }
    }
    *min = on.then_some(v);
}

/// For a BayesTraits branch value (rate scalar, Fabric beta), what its branch
/// posterior counts: "a branch scalar" or "an LS_Beta".
fn branch_pp_what(attr: &str) -> Option<&'static str> {
    match crate::varrates::branch_pp_key(attr)? {
        "fabric_beta_prop" => Some("an LS_Beta"),
        _ => Some("a branch scalar"),
    }
}

fn add_layer_menu(ui: &mut egui::Ui, doc: &mut Document) {
    let sel = doc.single_selection().filter(|&n| !doc.tree.is_tip(n));
    let keys = doc.tree.attr_keys();
    let add = |ui: &mut egui::Ui, name: &str, layer: Layer, doc: &mut Document| {
        if ui.button(name).clicked() {
            doc.add_layer(layer);
            ui.close_menu();
        }
    };
    add(ui, "Branches", Layer::default_tree(), doc);
    add(ui, "Tip labels", Layer::default_tip_labels(), doc);
    let node_attr = ["posterior", "support", "label"].iter().find(|k| keys.iter().any(|x| x == *k) || **k == "label").unwrap().to_string();
    add(ui, "Node labels", Layer::default_node_labels(&node_attr), doc);
    add(ui, "Branch lengths", Layer::default_branch_lengths(), doc);
    if doc.tree.preorder().iter().any(|&n| doc.tree.nodes[n].attrs.contains_key("rate_node_pp") || doc.tree.nodes[n].attrs.contains_key("fabric_node_prop")) {
        add(ui, "Node scalars (BayesTraits)", Layer::default_node_scalars(), doc);
    }
    if !doc.overlay().is_empty() {
        if ui.button("DensiTree").clicked() {
            // Underneath the other layers.
            doc.checkpoint();
            doc.layers.insert(0, LayerEntry::new(Layer::default_densitree()));
            doc.selected_layer = Some(0);
            ui.close_menu();
        }
    } else {
        ui.add_enabled(false, egui::Button::new("DensiTree (open a tree set first)"));
    }
    let youngest = doc.tree.youngest_tip_age().unwrap_or(0.0);
    add(ui, "Geologic timescale", Layer::default_geoscale().with_youngest_age(youngest), doc);
    add(ui, "Fossil ranges", Layer::default_fossil_ranges(), doc);
    add(ui, "Tip points", Layer::default_tip_points(), doc);
    add(ui, "Node points", Layer::default_node_points(), doc);
    add(ui, "Node bars / HPD", Layer::default_node_bars(), doc);
    add(ui, "Scale bar", Layer::default_scale_bar(), doc);
    add(ui, "Time axis", Layer::default_axis().with_youngest_age(youngest), doc);
    add(ui, "PhyloPic silhouettes", Layer::default_phylopic(), doc);
    let numeric: Vec<String> = keys.iter().filter(|k| is_numeric_attr(doc, k)).cloned().collect();
    add(
        ui,
        "Heatmap",
        Layer::Heatmap(HeatmapStyle {
            columns: numeric.first().cloned().into_iter().collect(),
            offset: 4.0,
            cell_width: 12.0,
            palette: Palette::Viridis,
            reverse_palette: false,
            show_names: true,
            name_size: 8.0,
        }),
        doc,
    );
    add(
        ui,
        "Bar chart",
        Layer::default_bars(&numeric.first().cloned().unwrap_or_default()),
        doc,
    );
    ui.separator();
    ui.add_enabled_ui(sel.is_some(), |ui| {
        let n = sel.unwrap_or(0);
        add(ui, "Highlight selected clade", Layer::Highlight(HighlightStyle { node: n, fill: Color::hex("#56B4E9").with_alpha(70), extend: 0.0 }), doc);
        add(
            ui,
            "Label selected clade",
            Layer::CladeLabel(CladeLabelStyle::new(n, "Clade".into())),
            doc,
        );
    });
    if sel.is_none() {
        ui.small("Select an internal node to highlight or label a clade.");
    }
}

fn is_numeric_attr(doc: &Document, k: &str) -> bool {
    let mut any = false;
    for n in doc.tree.preorder() {
        if let Some(v) = doc.tree.nodes[n].attrs.get(k) {
            match v {
                Attr::Num(_) => any = true,
                _ => return false,
            }
        }
    }
    any
}

/// Layout settings plus whole-tree arrangement commands (returned as actions).
pub fn layout_panel(ui: &mut egui::Ui, doc: &mut Document) -> Vec<Action> {
    let mut actions = Vec::new();
    let root = doc.tree.root;
    let l = &mut doc.view.layout;
    ui.heading("Layout");
    ComboBox::from_id_salt("layout-kind").selected_text(l.kind.name()).show_ui(ui, |ui| {
        for k in LayoutKind::MENU {
            ui.selectable_value(&mut l.kind, k, k.name());
        }
    });
    if l.kind == LayoutKind::Slanted {
        ui.checkbox(&mut l.slanted_cladogram, "V-shaped cladogram (no crossing branches)")
            .on_hover_text("Untick to place nodes by branch length instead, as ggtree's slanted layout does (branches may cross).");
    }
    if !(l.kind == LayoutKind::Slanted && l.slanted_cladogram) {
        ui.checkbox(&mut l.use_lengths, "Use branch lengths (off = cladogram)");
    }
    if l.kind == LayoutKind::Fan {
        ui.add(Slider::new(&mut l.open_angle, 0.0..=330.0).text("open angle°"));
    }
    if l.kind.is_polar() {
        ui.add(Slider::new(&mut l.rotate, -180.0..=180.0).text("rotate°"));
    }
    if matches!(l.kind, LayoutKind::Rectangular | LayoutKind::Slanted | LayoutKind::Roundrect | LayoutKind::Ellipse) {
        ui.horizontal(|ui| {
            ui.label("Tips at");
            for t in TipSide::ALL {
                ui.selectable_value(&mut l.tips, t, t.name());
            }
        });
    }
    if ui.button("Ladderize").on_hover_text("Order every clade by size, larger clades down (tick Reverse tip order to put them up)").clicked() {
        actions.push(Action::Ladderize(root, true));
    }
    if l.kind != LayoutKind::Radial {
        ui.checkbox(&mut l.flip_y, "Reverse tip order")
            .on_hover_text(if l.kind.is_polar() { "Run the tips the other way around the circle" } else { "Flip the drawing so the tip order runs the other way" });
    }
    if ui.button("Expand all collapsed").clicked() {
        actions.push(Action::ExpandAll);
    }
    if ui.button("Midpoint root").clicked() {
        actions.push(Action::MidpointRoot);
    }
    ui.horizontal(|ui| {
        // The root's own branch length, else a length set here (5% of the tree's
        // height to start), so the edge shows for any tree.
        let height = doc.tree.depths(true).iter().copied().fold(0.0, f64::max);
        if ui.checkbox(&mut l.root_edge, "Show root edge").changed() && l.root_edge && crate::layout::root_edge_length(&doc.tree, l) <= 0.0 {
            l.root_edge_len = (height * 0.05).max(1e-6);
        }
        if l.root_edge && l.use_lengths {
            let mut len = crate::layout::root_edge_length(&doc.tree, l);
            if ui.add(DragValue::new(&mut len).speed((height * 0.005).max(1e-4)).range(0.0..=f64::MAX).max_decimals(4).prefix("length ")).on_hover_text("Length of the root edge, in branch-length units").changed() {
                l.root_edge_len = len;
            }
        }
    });
    ui.horizontal(|ui| {
        ui.label("Title");
        ui.text_edit_singleline(&mut doc.view.title);
    });
    num(ui, "title size", &mut doc.view.title_size, 6.0..=48.0);
    ui.horizontal(|ui| {
        let mut bg = doc.view.background.unwrap_or(Color::WHITE);
        color_edit(ui, "Background", &mut bg);
        doc.view.background = (bg != Color::WHITE).then_some(bg);
        if doc.view.background.is_some() && ui.small_button("White").clicked() {
            doc.view.background = None;
        }
    });
    ui.add(Slider::new(&mut doc.view.legend_scale, 0.5..=3.0).text("legend size"))
        .on_hover_text("Size of the color keys at the upper right (comparative data, variable rates, Fabric, node scalars)");
    actions
}

/// Zoom slider for the canvas (about its centre), with Fit to reset.
pub fn zoom_panel(ui: &mut egui::Ui, doc: &mut Document) {
    ui.horizontal(|ui| {
        let mut z = doc.zoom[0];
        let r = ui.add(Slider::new(&mut z, 0.2..=60.0).logarithmic(true).text("zoom").custom_formatter(|v, _| format!("{:.0}%", v * 100.0)));
        if r.changed() && doc.zoom[0] > 0.0 {
            // Scales both directions, keeping any extra vertical stretch (Alt + / −).
            let f = z / doc.zoom[0];
            doc.zoom_by(f, f, None);
        }
        if ui.button("Fit").clicked() {
            doc.zoom = [1.0, 1.0];
            doc.pan = [0.0, 0.0];
        }
    });
}

pub fn layer_editor(ui: &mut egui::Ui, doc: &mut Document) {
    let Some(i) = doc.selected_layer else {
        ui.weak("Select a layer to edit its properties.");
        return;
    };
    if i >= doc.layers.len() {
        return;
    }
    let keys = keys_with_label(doc);
    let attr_keys = doc.tree.attr_keys();
    let available = doc.overlay().len();
    let data_cols: Vec<String> = doc.data.as_ref().map(|d| d.trait_columns().into_iter().map(|c| c.0).collect()).unwrap_or_default();
    let mut heat_cols: Vec<String> = attr_keys.clone();
    for c in data_cols {
        if !heat_cols.contains(&c) {
            heat_cols.push(c);
        }
    }
    let layer = &mut doc.layers[i].layer;
    let tip_points = matches!(layer, Layer::TipPoints(_));
    ui.strong(layer.description());
    let mut reconstruct: Option<String> = None;
    let id = |s: &str| format!("{}-{}", s, i);
    match layer {
        Layer::Tree(s) => {
            if s.width_by.is_some() {
                num(ui, "thickest", &mut s.width, 0.1..=20.0);
                num(ui, "thinnest", &mut s.width_min, 0.0..=20.0);
            } else {
                num(ui, "line width", &mut s.width, 0.1..=20.0);
            }
            attr_select(ui, &id("twb"), "Width by", &mut s.width_by, &attr_keys);
            if s.width_by.is_some() {
                ui.small("Each branch's width from its node's value: the thinnest at the smallest value, the thickest at the largest. Data columns are on the tips only, so reconstruct ancestral states (the branches button in the Comparative data panel) to give internal branches values; branches without one are drawn thinnest.");
            }
            color_edit(ui, "color", &mut s.color);
            attr_select(ui, &id("tcb"), "Color by", &mut s.color_by, &attr_keys);
            // A trait with no values at internal nodes (a data column not yet
            // reconstructed): those branches would show no color or width.
            let mut asked: Vec<String> = Vec::new();
            for k in [s.color_by.clone(), s.width_by.clone()].into_iter().flatten() {
                if asked.contains(&k) {
                    continue;
                }
                asked.push(k.clone());
                let internal: Vec<NodeId> = doc.tree.preorder().into_iter().filter(|&n| !doc.tree.is_tip(n)).collect();
                let missing = internal.iter().filter(|&&n| doc.tree.value(n, &k).is_none()).count();
                if missing > 0 && missing == internal.len() && doc.tree.tips().iter().any(|&t| doc.tree.value(t, &k).is_some()) {
                    ui.horizontal_wrapped(|ui| {
                        ui.colored_label(Color32::from_rgb(198, 120, 0), format!("{} has values on the tips only; internal branches are drawn grey and dashed.", k));
                        if ui.small_button("Reconstruct").on_hover_text("Estimate it at internal nodes: ML under Brownian motion for numbers, Fitch parsimony for categories").clicked() {
                            reconstruct = Some(k.clone());
                        }
                    });
                }
            }
            palette_select(ui, &id("tp"), &mut s.palette, &mut s.reverse_palette);
            if let Some(what) = s.color_by.as_deref().and_then(branch_pp_what) {
                ui.add(Slider::new(&mut s.min_pp, 0.0..=1.0).text("posterior ≥")).on_hover_text(format!("Color only branches where the proportion of the posterior sample with {} on that branch is at least this; the rest keep the line color", what));
                if let Some(k) = s.color_by.clone() {
                    min_value_slider(ui, &doc.tree, &k, &mut s.min_value, &mut s.min_abs, "Color only branches whose value is at least this; the rest keep the line color");
                }
            }
            ui.small("Clade colors from the context menu override this.");
        }
        Layer::TipLabels(s) => {
            num(ui, "size", &mut s.size, 3.0..=40.0);
            num(ui, "offset", &mut s.offset, 0.0..=60.0);
            color_edit(ui, "color", &mut s.color);
            ui.checkbox(&mut s.italic, "Italic");
            ui.checkbox(&mut s.align, "Align (dotted leaders)");
            ui.checkbox(&mut s.underscores_as_spaces, "Show underscores as spaces");
            ui.checkbox(&mut s.first_word, "First word only").on_hover_text("Drop the first underscore or space and everything after it (Homo_sapiens → Homo)");
            attr_select(ui, &id("lcb"), "Color by", &mut s.color_by, &attr_keys);
            palette_select(ui, &id("lp"), &mut s.palette, &mut s.reverse_palette);
        }
        Layer::NodeLabels(s) => {
            attr_select_req(ui, &id("nla"), "Show", &mut s.attr, &keys);
            num(ui, "size", &mut s.size, 3.0..=30.0);
            color_edit(ui, "color", &mut s.color);
            ui.add(Slider::new(&mut s.digits, 0..=6).text("digits"));
            let mut use_min = s.min_value.is_some();
            ui.horizontal(|ui| {
                ui.checkbox(&mut use_min, "Only ≥");
                let mut v = s.min_value.unwrap_or(0.5);
                ui.add_enabled(use_min, DragValue::new(&mut v).speed(0.01));
                s.min_value = if use_min { Some(v) } else { None };
            });
            ui.checkbox(&mut s.on_branch, "Place on branch");
            ui.horizontal(|ui| {
                ui.label("Suffix");
                ui.add(egui::TextEdit::singleline(&mut s.suffix).desired_width(40.0)).on_hover_text("Text after each number, e.g. %");
            });
            let mut use_th = s.threshold.is_some();
            ui.checkbox(&mut use_th, "Color by threshold");
            if use_th {
                let (mut th, mut below, mut above) = s.threshold.unwrap_or((0.5, Color::hex("#C62828"), Color::hex("#2E7D32")));
                ui.horizontal(|ui| {
                    ui.label("Threshold");
                    ui.add(DragValue::new(&mut th).speed(0.01));
                });
                color_edit(ui, "below threshold", &mut below);
                color_edit(ui, "at or above threshold", &mut above);
                s.threshold = Some((th, below, above));
            } else {
                s.threshold = None;
            }
        }
        Layer::DensiTree(s) => {
            // The slider tops out at the post-burn-in sample size (when one is loaded).
            if available > 0 {
                s.max_trees = s.max_trees.clamp(1, available);
                ui.add(Slider::new(&mut s.max_trees, 1..=available).logarithmic(available > 50).text(format!("trees drawn (of {})", available)));
            } else {
                ui.weak("No tree set loaded: open the posterior or bootstrap trees to draw this layer.");
            }
            let mut a = s.alpha as f32;
            ui.add(Slider::new(&mut a, 1.0..=255.0).logarithmic(true).text("opacity"));
            s.alpha = a as u8;
            num(ui, "line width", &mut s.width, 0.2..=5.0);
            color_edit(ui, "color", &mut s.color);
            ui.checkbox(&mut s.by_topology, "Color by topology (1st blue, 2nd red, 3rd green)");
            ui.small("Trees are spread evenly through the tree set (after any burn-in) and aligned at the tips. Drawing many trees can slow the canvas; keep this below a few hundred while editing.");
        }
        Layer::Geoscale(s) => {
            ui.checkbox(&mut s.ages, "Ages (stages)");
            ui.checkbox(&mut s.epochs, "Epochs");
            ui.checkbox(&mut s.periods, "Periods");
            ui.checkbox(&mut s.eras, "Eras");
            ui.checkbox(&mut s.labels, "Labels");
            ui.checkbox(&mut s.boundaries, "Boundary lines across the tree");
            if s.boundaries {
                color_edit(ui, "boundary color", &mut s.boundary_color);
                num(ui, "boundary width", &mut s.boundary_width, 0.1..=5.0);
            }
            num(ui, "row height", &mut s.row_height, 4.0..=40.0);
            num(ui, "label size", &mut s.label_size, 3.0..=20.0);
            ui.horizontal(|ui| {
                ui.label("Ma per branch-length unit");
                ui.add(DragValue::new(&mut s.ma_per_unit).speed(0.01).range(0.000001..=1e6));
            });
            ui.horizontal(|ui| {
                ui.label("Youngest tip age (Ma)");
                ui.add(DragValue::new(&mut s.youngest_age).speed(0.1).range(0.0..=5000.0));
            });
            ui.small("International Chronostratigraphic Chart v2024/12 (via Macrostrat). Needs a time-scaled tree: branch lengths in Ma (or set the conversion above). Circular layouts show the finest level as rings.");
        }
        Layer::FossilRanges(s) => {
            num(ui, "bar width", &mut s.width, 1.0..=30.0);
            color_edit(ui, "color", &mut s.color);
            ui.small("Bars span each taxon's fossil range (fossil_max_age to fossil_min_age, in Ma). Ages are placed with the geologic timescale layer's settings if there is one, else the time axis's youngest tip age; otherwise branch lengths are read as Ma with the youngest tip at the present. Load ranges with Tree › Import fossil ranges…, or right-click a tip.");
        }
        Layer::NodeScalars(s) => {
            ui.add(Slider::new(&mut s.threshold, 0.0..=1.0).text("posterior ≥"))
                .on_hover_text("Show a circle where the proportion of the posterior sample with a node scalar at that node is at least this");
            ui.checkbox(&mut s.size_by_pp, "Size by posterior").on_hover_text("Larger circles where the node scalar was placed more often");
            if s.size_by_pp {
                num(ui, "size at posterior 1", &mut s.size, 2.0..=30.0);
                num(ui, "size at posterior 0", &mut s.size_min, 0.5..=30.0);
            } else {
                num(ui, "circle size", &mut s.size, 2.0..=30.0);
            }
            let mut single = s.single_color.is_some();
            ui.checkbox(&mut single, "One color (no gradient)").on_hover_text("Draw every circle in one color instead of coloring by posterior");
            if single {
                let mut c = s.single_color.unwrap_or(Color::hex("#C62828"));
                color_edit(ui, "color", &mut c);
                s.single_color = Some(c);
            } else {
                s.single_color = None;
                palette_select(ui, &id("nsp"), &mut s.palette, &mut s.reverse_palette);
            }
            ui.small("Circles mark nodes where BayesTraits placed a node scalar, colored (or sized) by how often (posterior, 0–1). Import a VarRates file first: Tree › Import BayesTraits variable rates or Fabric.");
        }
        Layer::BranchLengths(s) => {
            ui.horizontal(|ui| {
                ui.label("Show");
                ComboBox::from_id_salt(id("bla")).selected_text(s.attr.clone().unwrap_or_else(|| "branch length".into())).show_ui(ui, |ui| {
                    ui.selectable_value(&mut s.attr, None, "branch length");
                    for k in &attr_keys {
                        ui.selectable_value(&mut s.attr, Some(k.clone()), k);
                    }
                });
            });
            let mut hide = s.hide_value.is_some();
            ui.horizontal(|ui| {
                ui.checkbox(&mut hide, "Hide branches equal to");
                let mut v = s.hide_value.unwrap_or(1.0);
                ui.add_enabled(hide, DragValue::new(&mut v).speed(0.01));
                s.hide_value = if hide { Some(v) } else { None };
            });
            if let Some(what) = s.attr.as_deref().and_then(branch_pp_what) {
                ui.add(Slider::new(&mut s.min_pp, 0.0..=1.0).text("posterior ≥")).on_hover_text(format!("Label only branches where the proportion of the posterior sample with {} on that branch is at least this", what));
                if let Some(k) = s.attr.clone() {
                    min_value_slider(ui, &doc.tree, &k, &mut s.min_value, &mut s.min_abs, "Label only branches whose value is at least this");
                }
            }
            num(ui, "size", &mut s.size, 3.0..=30.0);
            color_edit(ui, "color", &mut s.color);
            ui.add(Slider::new(&mut s.digits, 0..=6).text("decimals"));
        }
        Layer::TipPoints(s) | Layer::NodePoints(s) => {
            num(ui, "size", &mut s.size, 1.0..=30.0);
            attr_select(ui, &id("psb"), "Size by", &mut s.size_by, &attr_keys);
            if s.size_by.is_some() {
                num(ui, "smallest size", &mut s.size_min, 0.5..=30.0);
                ui.small("Size runs from the smallest size (lowest value) to the size above (highest value).");
            }
            if tip_points {
                num(ui, "distance from tip", &mut s.offset, 0.0..=100.0);
            }
            shape_select(ui, &id("ps"), &mut s.shape);
            color_edit(ui, "color", &mut s.color);
            attr_select(ui, &id("pcb"), "Color by", &mut s.color_by, &attr_keys);
            palette_select(ui, &id("pp"), &mut s.palette, &mut s.reverse_palette);
            let mut use_f = s.filter.is_some();
            ui.checkbox(&mut use_f, "Only where attribute ≥ threshold");
            if use_f {
                let (mut k, mut v) = s.filter.clone().unwrap_or_else(|| ("posterior".into(), 0.95));
                attr_select_req(ui, &id("pf"), "Attribute", &mut k, &attr_keys);
                ui.add(DragValue::new(&mut v).speed(0.01).prefix("≥ "));
                s.filter = Some((k, v));
            } else {
                s.filter = None;
            }
        }
        Layer::NodeBars(s) => {
            attr_select_req(ui, &id("nb"), "Range attribute", &mut s.attr, &attr_keys);
            num(ui, "width", &mut s.width, 0.5..=20.0);
            color_edit(ui, "color", &mut s.color);
            ui.small("Height ranges (e.g. height_95%_HPD) are drawn back from the youngest tip. For trees summarized with common-ancestor heights (TreeAnnotator -heights ca), height_95%_HPD uses CAheight_95%_HPD at nodes placed by CAheight_mean, so each bar matches its node.");
        }
        Layer::Highlight(s) => {
            color_edit(ui, "fill", &mut s.fill);
            num(ui, "extend", &mut s.extend, 0.0..=400.0);
        }
        Layer::CladeLabel(s) => {
            let mut pic = s.phylopic.is_some();
            if ui.checkbox(&mut pic, "PhyloPic silhouette").on_hover_text("A silhouette after the label, looked up on phylopic.org (starts with the label's text)").changed() {
                s.phylopic = pic.then(|| s.text.clone());
            }
            if let Some(name) = s.phylopic.as_mut() {
                ui.horizontal(|ui| {
                    ui.label("Taxon");
                    ui.text_edit_singleline(name).on_hover_text("Species, genus or a higher group (e.g. Hominidae), looked up on phylopic.org");
                });
                num(ui, "silhouette size", &mut s.phylopic_size, 4.0..=120.0);
                ui.small("Drawn in the label's color. Fetching must be on: Taxa › Show PhyloPic silhouettes for this tree, or Fetch PhyloPic silhouettes automatically.");
            }
            ui.horizontal(|ui| {
                ui.label("Text");
                ui.text_edit_singleline(&mut s.text);
            });
            num(ui, "size", &mut s.size, 4.0..=40.0);
            num(ui, "offset", &mut s.offset, 0.0..=200.0);
            num(ui, "bar width", &mut s.bar_width, 0.0..=10.0);
            color_edit(ui, "color", &mut s.color);
            ui.checkbox(&mut s.gradient, "Gradient behind branches").on_hover_text("A band of color behind the clade, fading to transparent from the label back to the clade's node");
            if s.gradient {
                color_edit(ui, "gradient color", &mut s.gradient_color);
                ui.small("The color's opacity sets how strong the band is at the label.");
            }
        }
        Layer::ScaleBar(s) => {
            let mut auto = s.length.is_none();
            ui.checkbox(&mut auto, "Automatic length");
            if auto {
                s.length = None;
            } else {
                let mut v = s.length.unwrap_or(0.1);
                ui.add(DragValue::new(&mut v).speed(0.001).range(0.0..=f64::MAX).prefix("length "));
                s.length = Some(v);
            }
            num(ui, "line width", &mut s.width, 0.2..=6.0);
            num(ui, "text size", &mut s.size, 4.0..=30.0);
            color_edit(ui, "color", &mut s.color);
        }
        Layer::TimeAxis(s) => {
            ui.checkbox(&mut s.time_before_present, "Time before present (ages)");
            if s.time_before_present {
                ui.horizontal(|ui| {
                    ui.label("Youngest tip age");
                    ui.add(DragValue::new(&mut s.youngest_age).speed(0.1).range(0.0..=f64::MAX));
                })
                .response
                .on_hover_text("Age of the youngest tip, for trees whose tips are not at the present (e.g. fossils)");
            }
            ui.checkbox(&mut s.grid, "Grid lines");
            if s.grid {
                color_edit(ui, "grid color", &mut s.grid_color);
                num(ui, "grid width", &mut s.grid_width, 0.1..=5.0);
            }
            color_edit(ui, "axis color (line, ticks and numbers)", &mut s.line_color);
            num(ui, "axis width", &mut s.line_width, 0.1..=5.0);
            num(ui, "text size", &mut s.size, 4.0..=30.0);
        }
        Layer::Heatmap(s) => {
            ui.label("Columns:");
            for c in &heat_cols {
                let mut on = s.columns.contains(c);
                if ui.checkbox(&mut on, c).changed() {
                    if on {
                        s.columns.push(c.clone());
                    } else {
                        s.columns.retain(|x| x != c);
                    }
                }
            }
            if heat_cols.is_empty() {
                ui.weak("Load trait data first (Data panel).");
            }
            num(ui, "cell width", &mut s.cell_width, 2.0..=60.0);
            num(ui, "offset", &mut s.offset, 0.0..=100.0);
            palette_select(ui, &id("hp"), &mut s.palette, &mut s.reverse_palette);
            ui.checkbox(&mut s.show_names, "Column names");
        }
        Layer::Bars(s) => {
            attr_select_req(ui, &id("bc"), "Column", &mut s.column, &attr_keys);
            num(ui, "width", &mut s.max_width, 10.0..=400.0);
            num(ui, "offset", &mut s.offset, 0.0..=100.0);
            color_edit(ui, "color", &mut s.color);
            color_edit(ui, "color below zero", &mut s.negative_color);
            ui.checkbox(&mut s.show_scale, "Show scale");
            ui.small("Bars start at zero: negative values extend the other way.");
        }
        Layer::Phylopic(s) => {
            ui.horizontal(|ui| {
                ui.checkbox(&mut s.tips, "At tips");
                ui.checkbox(&mut s.nodes, "At nodes").on_hover_text("Nodes given a silhouette: right-click a node › PhyloPic silhouette…");
            });
            num(ui, "height (pt)", &mut s.size, 2.0..=300.0);
            attr_select(ui, &id("phs"), "Size by", &mut s.size_by, &attr_keys);
            if s.size_by.is_some() {
                num(ui, "smallest height (pt)", &mut s.size_min, 1.0..=300.0);
                ui.small("Height runs from the smallest height (lowest value) to the height above (highest value).");
            }
            ui.checkbox(&mut s.fit_spacing, "Limit to space between tips")
                .on_hover_text("Keeps silhouettes from overlapping. Untick to use the height above as-is, or zoom vertically to give them more room.");
            num(ui, "offset", &mut s.offset, 0.0..=100.0);
            ui.checkbox(&mut s.align, "Align in a column");
            attr_select(ui, &id("phc"), "Color by", &mut s.color_by, &attr_keys);
            if s.color_by.is_some() {
                palette_select(ui, &id("php"), &mut s.palette, &mut s.reverse_palette);
            }
            let mut tint = s.tint.is_some();
            ui.checkbox(&mut tint, "Recolor");
            if tint {
                let mut c = s.tint.unwrap_or(Color::BLACK);
                color_edit(ui, "silhouette color", &mut c);
                s.tint = Some(c);
            } else {
                s.tint = None;
            }
            ui.small("Silhouettes are fetched from phylopic.org by tip name (binomial, then genus), or by the taxon given to a node (right-click › PhyloPic silhouette…). Color by a trait overrides Recolor.");
        }
    }
    if let Some(k) = reconstruct {
        // Numbers by ML (Brownian motion), categories by Fitch parsimony.
        let numeric = match doc.data.as_ref().and_then(|d| d.columns.iter().position(|c| *c == k).map(|c| d.kind(c))) {
            Some(kind) => kind == ColumnKind::Numeric,
            None => doc.tree.tips().iter().filter_map(|&t| doc.tree.nodes[t].attrs.get(&k)).all(|a| a.as_f64().is_some()),
        };
        doc.checkpoint();
        if numeric {
            ancestral::reconstruct_continuous(&mut doc.tree, &k);
        } else {
            ancestral::reconstruct_discrete(&mut doc.tree, &k);
        }
    }
}

pub fn inspector(ui: &mut egui::Ui, doc: &mut Document, search: &mut String) -> Vec<Action> {
    let mut actions = Vec::new();
    ui.heading("Selection");
    ui.horizontal(|ui| {
        ui.label("Find");
        let r = ui.text_edit_singleline(search);
        if (r.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter))) || ui.button("Select").clicked() {
            let q = search.to_lowercase();
            doc.selection = doc.tree.tips().into_iter().filter(|&n| !q.is_empty() && doc.tree.label(n).to_lowercase().contains(&q)).collect();
        }
    });
    match doc.selection.len() {
        0 => {
            ui.weak("Click a node, branch or label. Shift-click to add to the selection. Drag a node onto another branch to move it.");
        }
        1 => {
            let n = *doc.selection.iter().next().unwrap();
            node_details(ui, doc, n, &mut actions);
        }
        k => {
            ui.label(format!("{} nodes selected ({} tips)", k, doc.selected_tips().len()));
            let sel: Vec<NodeId> = doc.selection.iter().copied().collect();
            ui.horizontal_wrapped(|ui| {
                if ui.button("MRCA").clicked() {
                    if let Some(m) = doc.tree.mrca(&sel) {
                        actions.push(Action::SelectClade(m));
                    }
                }
                if ui.button("Drop").clicked() {
                    actions.push(Action::Drop(sel.clone()));
                }
                if ui.button("Keep only").clicked() {
                    actions.push(Action::KeepOnly(sel.clone()));
                }
                if ui.button("Reroot on these (outgroup)").clicked() {
                    if let Some(m) = doc.tree.mrca(&sel) {
                        actions.push(Action::Reroot(m));
                    }
                }
            });
            let tips = doc.selected_tips();
            taxa_table(ui, doc, &tips, usize::MAX);
        }
    }
    actions
}

fn node_details(ui: &mut egui::Ui, doc: &mut Document, n: NodeId, actions: &mut Vec<Action>) {
    let is_tip = doc.tree.is_tip(n);
    ui.strong(if is_tip { format!("Tip #{}", n) } else { format!("Internal node #{} — {} tips", n, doc.tree.tips_below(n).len()) });
    egui::Grid::new("node-grid").num_columns(2).striped(true).show(ui, |ui| {
        ui.label("label");
        let mut label = doc.tree.nodes[n].label.clone().unwrap_or_default();
        let r = ui.text_edit_singleline(&mut label);
        if r.gained_focus() {
            doc.checkpoint();
        }
        if r.changed() {
            doc.set_label(n, if label.is_empty() { None } else { Some(label) });
        }
        if r.lost_focus() {
            doc.show_node_label(n);
        }
        ui.end_row();
        ui.label("branch length");
        let mut l = doc.tree.nodes[n].length.unwrap_or(0.0);
        let r = ui.add(DragValue::new(&mut l).speed(0.001).range(0.0..=f64::MAX));
        if r.drag_started() || (r.changed() && !r.dragged()) {
            doc.checkpoint();
        }
        if r.changed() {
            doc.tree.nodes[n].length = Some(l);
        }
        ui.end_row();
        let depth = doc.tree.depths(true)[n];
        ui.label("root distance");
        ui.label(format_num(depth, 6));
        ui.end_row();
        if !doc.tree.nodes[n].attrs.contains_key("height") && doc.tree.has_lengths() {
            let h = doc.tree.heights()[n];
            ui.label("height").on_hover_text("Time before the youngest tip, from branch lengths");
            ui.label(format_num(if h.abs() < 1e-4 { 0.0 } else { h }, 3));
            ui.end_row();
        }
        if !is_tip {
            let (c, norm, poly) = doc.tree.colless(n);
            ui.label("Colless index").on_hover_text(COLLESS_HOVER);
            let mut txt = match norm {
                Some(x) => format!("{} (normalized {})", c, format_num(x, 3)),
                None => c.to_string(),
            };
            if poly {
                txt.push_str(" — polytomies skipped");
            }
            ui.label(txt);
            ui.end_row();
            ui.label("cherries").on_hover_text("Pairs of sister tips in this clade");
            ui.label(doc.tree.cherries(n).to_string());
            ui.end_row();
            ui.label("γ (gamma)").on_hover_text(GAMMA_HOVER);
            match doc.tree.gamma(n) {
                Ok((g, p)) => ui.label(format!("{} (p = {}, two-tailed)", format_num(g, 3), format_num(p, 3))),
                Err(why) => ui.weak("n/a").on_hover_text(format!("γ not computed: {}", why)),
            };
            ui.end_row();
        }
        for (k, v) in &doc.tree.nodes[n].attrs {
            ui.label(k);
            ui.label(crate::io::data::table_value(k, v));
            ui.end_row();
        }
    });
    ui.horizontal_wrapped(|ui| {
        if n != doc.tree.root && ui.button("Reroot").clicked() {
            actions.push(Action::Reroot(n));
        }
        if !is_tip {
            if ui.button("Rotate").clicked() {
                actions.push(Action::Rotate(n));
            }
            if n != doc.tree.root && ui.button(if doc.view.collapsed.contains(&n) { "Expand" } else { "Collapse" }).clicked() {
                actions.push(Action::ToggleCollapse(n));
            }
            if ui.button("Highlight").clicked() {
                actions.push(Action::Highlight(n));
            }
            if ui.button("Label").clicked() {
                actions.push(Action::CladeLabel(n));
            }
            if ui.button("Open in tab").clicked() {
                actions.push(Action::Extract(n));
            }
        }
        if n != doc.tree.root && ui.button("Drop").clicked() {
            actions.push(Action::Drop(vec![n]));
        }
        if ui.button("Copy clade").on_hover_text("Ctrl+C: copy to paste onto another tree").clicked() {
            actions.push(Action::CopyClade(n));
        }
        if n != doc.tree.root && ui.button("Cut clade").on_hover_text("Ctrl+X: copy, then remove it from this tree").clicked() {
            actions.push(Action::CutClade(n));
        }
        ui.menu_button("Add taxon", |ui| {
            if ui.button("As sister").on_hover_text("Splits the branch above this node").clicked() {
                actions.push(Action::AddTaxon(n, crate::ops::PasteMode::Sister));
                ui.close_menu();
            }
            if !is_tip && ui.button("As child (polytomy)").clicked() {
                actions.push(Action::AddTaxon(n, crate::ops::PasteMode::Child));
                ui.close_menu();
            }
        });
        if !is_tip && ui.button("Transform…").on_hover_text("Scale, or λ, κ, δ or Ornstein–Uhlenbeck transform of this clade").clicked() {
            actions.push(Action::Transform(n));
        }
    });
    let tips = doc.tree.tips_below(n);
    taxa_table(ui, doc, &tips, n);
    if !is_tip {
        vcv_matrix(ui, doc, n, tips.len());
    }
}

/// Phylogenetic variance–covariance matrix of a clade, with copy and save buttons.
/// Collapsed by default; built only while open, since it grows with the square of the tip count.
fn vcv_matrix(ui: &mut egui::Ui, doc: &Document, n: NodeId, k: usize) {
    egui::CollapsingHeader::new(format!("Variance–covariance matrix ({} × {})", k, k)).id_salt(("vcv-closed", n)).default_open(false).show(ui, |ui| {
        ui.small("Shared branch length from this node for each pair of tips; the diagonal is each tip's distance from this node (ape vcv).");
        ui.horizontal_wrapped(|ui| {
            if ui.button("Copy matrix").on_hover_text("Tab-separated, full precision; pastes into Excel").clicked() {
                let (h, rows) = crate::io::data::vcv_table(&doc.tree, n, None);
                ui.ctx().copy_text(crate::io::data::table_to_text(&h, &rows, b'\t'));
            }
            if ui.button("Export…").on_hover_text("Full precision, as CSV or tab-delimited text").clicked() {
                let (h, rows) = crate::io::data::vcv_table(&doc.tree, n, None);
                save_table(&h, &rows, "vcv", "the variance–covariance matrix");
            }
        });
        if k > 150 {
            ui.weak("Too many tips to show here; copy or save the matrix instead.");
            return;
        }
        let (header, rows) = crate::io::data::vcv_table(&doc.tree, n, Some(3));
        egui::ScrollArea::both().max_height(300.0).id_salt(("vcv-scroll", n)).show(ui, |ui| {
            egui::Grid::new(("vcv-grid", n)).striped(true).show(ui, |ui| {
                let tid = egui::Id::new(("vcv-grid", n));
                for (i, h) in header.iter().enumerate() {
                    super::table::header_sized(ui, tid, i, &h.replace('_', " "), None);
                }
                ui.end_row();
                for r in super::table::order(ui, tid, &rows).into_iter().map(|k| &rows[k]) {
                    for (j, v) in r.iter().enumerate() {
                        if j == 0 {
                            ui.label(egui::RichText::new(v.replace('_', " ")).italics().strong());
                        } else {
                            ui.label(v);
                        }
                    }
                    ui.end_row();
                }
            });
        });
    });
}

/// One-sentence explanations shown when hovering over the statistics.
pub const COLLESS_HOVER: &str = "Tree balance: higher values mean a more unbalanced, ladder-like clade, and the normalized value runs up to 1 (fully ladder-like); it is 0 only for a perfectly balanced clade, which needs a power-of-two number of tips.";
pub const GAMMA_HOVER: &str = "Tempo of diversification: γ < 0 means branching is concentrated near the root (a slowdown), γ > 0 near the tips (a speed-up), and |γ| > 1.96 departs from a constant rate at p < 0.05.";

/// Species names and their data for a clade, with copy and save buttons.
fn taxa_table(ui: &mut egui::Ui, doc: &Document, tips: &[NodeId], id: NodeId) {
    let (header, rows) = crate::io::data::tips_table(&doc.tree, tips);
    egui::CollapsingHeader::new(format!("Tree and comparative data ({})", rows.len())).id_salt(("taxa-closed", id)).default_open(false).show(ui, |ui| {
        ui.horizontal_wrapped(|ui| {
            if ui.button("Copy names").clicked() {
                let names: Vec<&str> = rows.iter().map(|r| r[0].as_str()).collect();
                ui.ctx().copy_text(names.join("\n"));
            }
            if ui.button("Copy table").on_hover_text("Tab-separated, pastes into Excel").clicked() {
                ui.ctx().copy_text(crate::io::data::table_to_text(&header, &rows, b'\t'));
            }
            if ui.button("Export…").on_hover_text("As CSV or tab-delimited text").clicked() {
                save_table(&header, &rows, "clade", "the table");
            }
        });
        egui::ScrollArea::both().max_height(260.0).id_salt(("taxa-scroll", id)).show(ui, |ui| {
            egui::Grid::new(("taxa-grid", id)).striped(true).show(ui, |ui| {
                let tid = egui::Id::new(("taxa-grid", id));
                super::table::headers(ui, tid, &header);
                ui.end_row();
                for r in super::table::order(ui, tid, &rows).into_iter().map(|k| &rows[k]) {
                    for (j, v) in r.iter().enumerate() {
                        if j == 0 {
                            ui.label(egui::RichText::new(v.replace('_', " ")).italics());
                        } else {
                            ui.label(v);
                        }
                    }
                    ui.end_row();
                }
            });
        });
    });
}

/// "1 taxon" or "N taxa".
fn taxa(n: usize) -> String {
    if n == 1 { "1 taxon".into() } else { format!("{} taxa", n) }
}

/// e.g. "1 taxon in the tree has no data; 3 taxa in the data are not in the tree."
fn mismatch_summary(no_row: usize, no_tip: usize) -> String {
    let verb = |n: usize, one: &'static str, many: &'static str| if n == 1 { one } else { many };
    format!("{} in the tree {} no data; {} in the data {} not in the tree.", taxa(no_row), verb(no_row, "has", "have"), taxa(no_tip), verb(no_tip, "is", "are"))
}

/// Save a table as CSV (comma) or tab-delimited text (.tsv, .txt), chosen in
/// the file dialog. Returns a status message, or None if cancelled.
fn save_table(header: &[String], rows: &[Vec<String>], stem: &str, what: &str) -> Option<String> {
    let mut path = rfd::FileDialog::new().add_filter("CSV", &["csv"]).add_filter("Tab-delimited", &["tsv", "txt"]).set_file_name(format!("{}.csv", stem)).save_file()?;
    if path.extension().is_none() {
        path.set_extension("csv");
    }
    let tab = path.extension().map(|e| !e.eq_ignore_ascii_case("csv")).unwrap_or(false);
    let text = crate::io::data::table_to_text(header, rows, if tab { b'\t' } else { b',' });
    Some(match std::fs::write(&path, text) {
        Ok(()) => format!("Saved {} ({}) to {}", what, super::count(rows.len(), "row", "rows"), path.display()),
        Err(e) => format!("Could not save {}: {}", what, e),
    })
}


pub fn data_panel(ui: &mut egui::Ui, doc: &mut Document, status: &mut String, open_missing: &mut bool, open_editor: &mut bool, open_summary: &mut Option<super::SummaryView>) -> bool {
    let mut want_load = false;
    ui.heading("Trait Data");
    let Some(data) = doc.data.clone() else {
        ui.weak("Load a CSV/TSV with a taxon column (or drop it on the window). Rows are matched to tip labels.");
        if ui.button("Load data table…").clicked() {
            want_load = true;
        }
        return want_load;
    };
    ui.label(format!("{}: {} rows, key column \"{}\"", data.name, data.rows.len(), data.columns[data.key]));
    // Taxa in the tree without data, and rows matching no taxon.
    let (no_row, no_tip) = data.mismatches(&doc.tree);
    if no_row.is_empty() && no_tip.is_empty() {
        ui.weak(format!("All {} in the tree have data, and every taxon in the data is in the tree.", taxa(doc.tree.num_tips())));
    } else {
        let list = |ui: &mut egui::Ui, id: &str, title: String, names: &[String]| {
            if names.is_empty() {
                return;
            }
            egui::CollapsingHeader::new(title).id_salt(id).default_open(false).show(ui, |ui| {
                if ui.small_button("Copy names").clicked() {
                    ui.ctx().copy_text(names.join("\n"));
                }
                egui::ScrollArea::vertical().id_salt((id, "scroll")).max_height(140.0).show(ui, |ui| {
                    for n in names {
                        ui.label(egui::RichText::new(n.replace('_', " ")).italics());
                    }
                });
            });
        };
        ui.colored_label(Color32::from_rgb(198, 120, 0), mismatch_summary(no_row.len(), no_tip.len()));
        list(ui, "no-row", format!("In the tree, not in the data ({})", no_row.len()), &no_row);
        list(ui, "no-tip", format!("In the data, not in the tree ({})", no_tip.len()), &no_tip);
    }
    ui.horizontal_wrapped(|ui| {
        if ui.small_button("Replace…").clicked() {
            want_load = true;
        }
        if ui.small_button("Drop tips with missing data…").on_hover_text("Remove tips lacking values in chosen columns (undoable)").clicked() {
            *open_missing = true;
        }
        if ui
            .add_enabled(!no_tip.is_empty(), egui::Button::new("Drop data with missing tips").small())
            .on_hover_text(format!("Remove the {} for taxa that are not in the tree from the table", super::count(no_tip.len(), "row", "rows")))
            .on_disabled_hover_text("Every row in the table is a taxon in the tree")
            .clicked()
        {
            doc.data = Some(data.without(&no_tip));
            *status = format!("Dropped {} for taxa not in the tree", super::count(no_tip.len(), "data row", "data rows"));
        }
        if ui.small_button("Export…").on_hover_text("Save the comparative data table (with any edits and new columns) as CSV or tab-delimited text").clicked() {
            let stem = std::path::Path::new(&data.name).file_stem().and_then(|s| s.to_str()).unwrap_or("data").to_string();
            if let Some(msg) = save_table(&data.columns, &data.rows, &stem, "the comparative data") {
                *status = msg;
            }
        }
        if ui.small_button("Edit values…").on_hover_text("Edit the table's values in a spreadsheet-style window").clicked() {
            *open_editor = true;
        }
        if ui.small_button("Summary statistics…").on_hover_text("N, missing, mean, SD, quartiles and range of each column (states for categorical ones)").clicked() {
            if open_summary.is_none() {
                *open_summary = Some(super::SummaryView::default());
            }
        }
    });
    egui::Grid::new("data-cols").num_columns(2).striped(true).show(ui, |ui| {
        for (col, kind) in data.trait_columns() {
            ui.horizontal(|ui| {
                ui.label(&col);
                // Numeric or categorical (detected, and changeable: e.g. 0/1 states).
                let mut k = kind;
                ui.selectable_value(&mut k, ColumnKind::Numeric, "num").on_hover_text("Numeric: continuous colors, bars, Brownian-motion reconstruction");
                ui.selectable_value(&mut k, ColumnKind::Categorical, "cat").on_hover_text("Categorical: one color per state, Fitch parsimony reconstruction");
                if k != kind {
                    doc.set_column_kind(&col, k);
                    *status = format!("\"{}\" is now {}", col, if k == ColumnKind::Numeric { "numeric (values that aren't numbers count as missing)" } else { "categorical" });
                }
                // New columns from this one, and deleting it.
                ui.menu_button("⋯", |ui| {
                    use crate::io::data::Transform;
                    let numeric = kind == ColumnKind::Numeric;
                    for (t, label, hover) in [
                        (Transform::Center, "Center", "New column: each value minus the mean"),
                        (Transform::Standardize, "Standardize (z-scores)", "New column: (value − mean) / standard deviation, so it has mean 0 and standard deviation 1"),
                        (Transform::Ln, "Log (ln)", "New column: natural logarithm (values ≤ 0 are left empty)"),
                        (Transform::Log10, "Log₁₀", "New column: base-10 logarithm (values ≤ 0 are left empty)"),
                    ] {
                        if ui.add_enabled(numeric, egui::Button::new(label)).on_hover_text(hover).on_disabled_hover_text("Only for numeric columns").clicked() {
                            if let Some(d) = doc.data.as_mut() {
                                let c = d.columns.iter().position(|x| *x == col).unwrap();
                                match d.add_transform(c, t) {
                                    Ok((name, lost)) => {
                                        doc.rejoin_data(&[]);
                                        *status = if lost > 0 { format!("Added \"{}\"; {} values \u{2264} 0 have no logarithm and were left empty", name, lost) } else { format!("Added \"{}\"", name) };
                                    }
                                    Err(e) => *status = format!("{:#}", e),
                                }
                            }
                            ui.close_menu();
                        }
                    }
                    ui.separator();
                    ui.menu_button("Delete column", |ui| {
                        if ui.button(format!("Delete \"{}\"", col)).on_hover_text("Removes the column from the table and its values from the tree (cannot be undone)").clicked() {
                            if let Some(d) = doc.data.as_mut() {
                                if let Some(c) = d.columns.iter().position(|x| *x == col) {
                                    d.remove_column(c);
                                }
                            }
                            doc.rejoin_data(std::slice::from_ref(&col));
                            *status = format!("Deleted column \"{}\"", col);
                            ui.close_menu();
                        }
                    });
                });
            });
            ui.horizontal(|ui| {
                if ui.small_button("labels").on_hover_text("Color tip labels").clicked() {
                    doc.checkpoint();
                    let mut found = false;
                    for e in doc.layers.iter_mut() {
                        if let Layer::TipLabels(s) = &mut e.layer {
                            s.color_by = Some(col.clone());
                            found = true;
                        }
                    }
                    if !found {
                        doc.layers.push(LayerEntry::new(Layer::default_tip_labels()));
                    }
                }
                if ui.small_button("points").on_hover_text("Add colored tip points").clicked() {
                    let mut l = Layer::default_tip_points();
                    if let Layer::TipPoints(s) = &mut l {
                        s.color_by = Some(col.clone());
                    }
                    doc.add_layer(l);
                }
                if ui.small_button("heatmap").on_hover_text("Add to heatmap").clicked() {
                    doc.checkpoint();
                    let existing = doc.layers.iter_mut().find_map(|e| if let Layer::Heatmap(h) = &mut e.layer { Some(h) } else { None });
                    match existing {
                        Some(h) => {
                            if !h.columns.contains(&col) {
                                h.columns.push(col.clone());
                            }
                        }
                        None => doc.layers.push(LayerEntry::new(Layer::Heatmap(HeatmapStyle {
                            columns: vec![col.clone()],
                            offset: 4.0,
                            cell_width: 12.0,
                            palette: if kind == ColumnKind::Numeric { Palette::Viridis } else { Palette::OkabeIto },
                            reverse_palette: false,
                            show_names: true,
                            name_size: 8.0,
                        }))),
                    }
                }
                if kind == ColumnKind::Numeric && ui.small_button("bars").on_hover_text("Add bar chart").clicked() {
                    doc.add_layer(Layer::default_bars(&col));
                }
                if ui
                    .small_button("branches")
                    .on_hover_text("Reconstruct ancestral states (ML under Brownian motion, or Fitch parsimony for categories) and color branches")
                    .clicked()
                {
                    doc.checkpoint();
                    let filled = match kind {
                        ColumnKind::Numeric => ancestral::reconstruct_continuous(&mut doc.tree, &col),
                        ColumnKind::Categorical => ancestral::reconstruct_discrete(&mut doc.tree, &col),
                    };
                    for e in doc.layers.iter_mut() {
                        if let Layer::Tree(s) = &mut e.layer {
                            s.color_by = Some(col.clone());
                            s.width = s.width.max(1.5);
                            s.palette = if kind == ColumnKind::Numeric { Palette::Viridis } else { Palette::OkabeIto };
                        }
                    }
                    *status = format!("Reconstructed '{}' at {} nodes", col, filled);
                }
            });
            ui.end_row();
        }
    });
    // The values, one row per taxon (click a header to sort).
    egui::CollapsingHeader::new(format!("Table ({})", super::count(data.rows.len(), "row", "rows"))).id_salt("data-table").default_open(true).show(ui, |ui| {
        egui::ScrollArea::both().id_salt("data-table-scroll").max_height(320.0).show(ui, |ui| {
            egui::Grid::new("data-table-grid").striped(true).show(ui, |ui| {
                let tid = egui::Id::new("data-table");
                super::table::headers(ui, tid, &data.columns);
                ui.end_row();
                for k in super::table::order(ui, tid, &data.rows) {
                    for (c, v) in data.rows[k].iter().enumerate() {
                        if c == data.key {
                            ui.label(egui::RichText::new(v.replace('_', " ")).italics());
                        } else {
                            ui.label(v);
                        }
                    }
                    ui.end_row();
                }
            });
        });
    });
    want_load
}

pub enum PosteriorAction {
    Consensus,
    Mcc,
    OpenSample(usize),
    /// Entropy-based rogue taxa and skeleton analysis.
    Rogue,
}

pub fn posterior_panel(ui: &mut egui::Ui, doc: &mut Document) -> Option<PosteriorAction> {
    let mut out = None;
    let mut goto: Option<usize> = None;
    let mut leave = false;
    let mut refresh = false;
    let Some(p) = doc.posterior.as_mut() else { return None };
    let n = p.trees.len();
    ui.horizontal_wrapped(|ui| {
        ui.strong(format!("Tree set: {} trees", n));
        let mut kind = p.kind;
        ComboBox::from_id_salt("sample-kind").selected_text(kind.name()).show_ui(ui, |ui| {
            for k in [SampleKind::Posterior, SampleKind::Bootstrap] {
                ui.selectable_value(&mut kind, k, k.name());
            }
        });
        if kind != p.kind {
            p.set_kind(kind);
            refresh = true;
        }
        ui.separator();
        if p.kind == SampleKind::Posterior {
            let r = ui.add(Slider::new(&mut p.burnin, 0.0..=0.9).text("burn-in").custom_formatter(|v, _| format!("{:.0}%", v * 100.0)));
            // Recount clade support once the slider is released.
            refresh |= r.drag_stopped() || (r.changed() && !r.dragged());
            ui.label(format!("({} kept)", p.kept().len()));
        } else {
            ui.label("all replicates used (no burn-in)");
        }
        ui.separator();
        refresh |= ui.checkbox(&mut p.rooted, "Rooted clades").changed();
        ui.add(Slider::new(&mut p.threshold, 0.0..=1.0).text("consensus threshold"));
        ComboBox::from_id_salt("heights")
            .selected_text(match p.heights {
                HeightMode::Keep => "Keep target heights",
                HeightMode::Mean => "Mean heights",
                HeightMode::Median => "Median heights",
            })
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut p.heights, HeightMode::Keep, "Keep target heights");
                ui.selectable_value(&mut p.heights, HeightMode::Mean, "Mean heights");
                ui.selectable_value(&mut p.heights, HeightMode::Median, "Median heights");
            });
        if ui.button("Consensus").on_hover_text("Threshold < 0.5 adds compatible clades greedily (extended majority rule)").clicked() {
            out = Some(PosteriorAction::Consensus);
        }
        let mcc_help = match p.kind {
            SampleKind::Posterior => "Maximum clade credibility tree: the sampled tree whose clades have the highest product of posterior probabilities, annotated with posteriors and mean node heights",
            SampleKind::Bootstrap => "The replicate whose clades have the highest product of bootstrap proportions, annotated with bootstrap support",
        };
        if ui.button(if p.kind == SampleKind::Posterior { "MCC tree" } else { "Best-supported replicate" }).on_hover_text(mcc_help).clicked() {
            out = Some(PosteriorAction::Mcc);
        }
        if !p.info.is_empty() {
            ui.separator();
            ui.weak(&p.info);
        }
    });
    ui.horizontal_wrapped(|ui| {
        let browsing = p.browsing();
        ui.label("Step through trees:");
        let cur = p.browse;
        if ui.add_enabled(cur > 0 || !browsing, egui::Button::new("⏮")).on_hover_text("First tree").clicked() {
            goto = Some(0);
        }
        if ui.add_enabled(cur > 0, egui::Button::new("◀")).on_hover_text("Previous (←)").clicked() {
            goto = Some(cur - 1);
        }
        let play_label = if p.playing { "⏸ Pause" } else { "▶ Play" };
        if ui.button(play_label).clicked() {
            p.playing = !p.playing;
            p.last_step = None;
            if p.playing && !browsing {
                goto = Some(cur);
            }
        }
        if ui.add_enabled(cur + 1 < n, egui::Button::new("▶|")).on_hover_text("Next (→)").clicked() {
            goto = Some(cur + 1);
        }
        let mut idx = cur;
        let r = ui.add(Slider::new(&mut idx, 0..=n.saturating_sub(1)).custom_formatter(|v, _| format!("{}", v as usize + 1)).text(format!("of {}", n)));
        if r.changed() {
            goto = Some(idx);
        }
        ui.add(Slider::new(&mut p.speed, 0.5..=30.0).logarithmic(true).text("trees/s"));
        if ui.checkbox(&mut p.ladderize, "Ladderize").changed() && browsing {
            goto = Some(cur);
        }
        if browsing {
            if let Some(name) = doc_tree_name(&p.trees, cur) {
                ui.weak(name);
            }
            if ui.button("Back to original tree").clicked() {
                leave = true;
            }
        }
        if ui.button("Open this tree in new tab").clicked() {
            out = Some(PosteriorAction::OpenSample(cur));
        }
    });
    // Playback.
    if p.playing {
        let due = p.last_step.map(|t| t.elapsed().as_secs_f32() >= 1.0 / p.speed).unwrap_or(true);
        if due && goto.is_none() {
            goto = Some(if p.browsing() { (p.browse + 1) % n.max(1) } else { p.browse });
        }
        ui.ctx().request_repaint_after(std::time::Duration::from_secs_f32(1.0 / p.speed));
    }
    if goto.is_some() {
        p.last_step = Some(std::time::Instant::now());
    }
    let browsing = p.browsing();
    let cur = p.browse;
    if leave {
        doc.leave_samples();
    } else if let Some(i) = goto {
        doc.show_sample(i);
    } else if refresh {
        if browsing {
            doc.show_sample(cur);
        } else {
            doc.annotate_posterior();
        }
    }
    out
}

fn doc_tree_name(trees: &[crate::tree::Tree], i: usize) -> Option<String> {
    trees.get(i).map(|t| t.name.clone().unwrap_or_else(|| format!("sample {}", i + 1)))
}
#[cfg(test)]
mod tests {
    #[test]
    fn mismatch_wording() {
        assert_eq!(super::mismatch_summary(1, 1), "1 taxon in the tree has no data; 1 taxon in the data is not in the tree.");
        assert_eq!(super::mismatch_summary(3, 0), "3 taxa in the tree have no data; 0 taxa in the data are not in the tree.");
        assert_eq!(super::taxa(2), "2 taxa");
    }
}
