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

fn palette_select(ui: &mut egui::Ui, id: &str, p: &mut Palette) {
    ui.horizontal(|ui| {
        ui.label("Palette");
        ComboBox::from_id_salt(id).selected_text(p.name()).show_ui(ui, |ui| {
            for q in Palette::ALL {
                ui.selectable_value(p, q, q.name());
            }
        });
    });
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
    add(ui, "Geologic timescale", Layer::default_geoscale(), doc);
    add(ui, "Tip points", Layer::default_tip_points(), doc);
    add(ui, "Node points", Layer::default_node_points(), doc);
    add(ui, "Node bars / HPD", Layer::default_node_bars(), doc);
    add(ui, "Scale bar", Layer::default_scale_bar(), doc);
    add(ui, "Time axis", Layer::default_axis(), doc);
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
            Layer::CladeLabel(CladeLabelStyle { node: n, text: "Clade".into(), color: Color::BLACK, size: 11.0, offset: 6.0, bar_width: 2.0 }),
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

pub fn layout_panel(ui: &mut egui::Ui, doc: &mut Document) {
    let l = &mut doc.view.layout;
    ui.horizontal(|ui| {
        ui.label("Layout");
        ComboBox::from_id_salt("layout-kind").selected_text(l.kind.name()).show_ui(ui, |ui| {
            for k in LayoutKind::MENU {
                ui.selectable_value(&mut l.kind, k, k.name());
            }
        });
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
    if !l.kind.is_polar() && l.kind != LayoutKind::Radial {
        ui.checkbox(&mut l.flip_y, "Reverse tip order");
    }
    ui.checkbox(&mut l.root_edge, "Show root edge");
    ui.horizontal(|ui| {
        ui.label("Title");
        ui.text_edit_singleline(&mut doc.view.title);
    });
    ui.horizontal(|ui| {
        let mut bg = doc.view.background.unwrap_or(Color::WHITE);
        color_edit(ui, "Background", &mut bg);
        doc.view.background = (bg != Color::WHITE).then_some(bg);
        if doc.view.background.is_some() && ui.small_button("White").clicked() {
            doc.view.background = None;
        }
    });
    ui.horizontal(|ui| {
        ui.label("Zoom");
        if ui.button(" − ").on_hover_text("Zoom out (−)").clicked() {
            doc.zoom_by(0.8, 0.8, None);
        }
        if ui.button(" + ").on_hover_text("Zoom in (+)").clicked() {
            doc.zoom_by(1.25, 1.25, None);
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
    let id = |s: &str| format!("{}-{}", s, i);
    match layer {
        Layer::Tree(s) => {
            num(ui, "line width", &mut s.width, 0.1..=8.0);
            color_edit(ui, "color", &mut s.color);
            attr_select(ui, &id("tcb"), "Color by", &mut s.color_by, &attr_keys);
            palette_select(ui, &id("tp"), &mut s.palette);
            ui.small("Clade colors from the context menu override this.");
        }
        Layer::TipLabels(s) => {
            num(ui, "size", &mut s.size, 3.0..=40.0);
            num(ui, "offset", &mut s.offset, 0.0..=60.0);
            color_edit(ui, "color", &mut s.color);
            ui.checkbox(&mut s.italic, "Italic");
            ui.checkbox(&mut s.align, "Align (dotted leaders)");
            ui.checkbox(&mut s.underscores_as_spaces, "Show underscores as spaces");
            attr_select(ui, &id("lcb"), "Color by", &mut s.color_by, &attr_keys);
            palette_select(ui, &id("lp"), &mut s.palette);
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
            ui.checkbox(&mut s.epochs, "Epochs");
            ui.checkbox(&mut s.periods, "Periods");
            ui.checkbox(&mut s.eras, "Eras");
            ui.checkbox(&mut s.labels, "Labels");
            ui.checkbox(&mut s.boundaries, "Boundary lines across the tree");
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
            ui.small("ICS 2023 chart. Needs a time-scaled tree: branch lengths in Ma (or set the conversion above). Circular layouts show the finest level as rings.");
        }
        Layer::BranchLengths(s) => {
            num(ui, "size", &mut s.size, 3.0..=30.0);
            color_edit(ui, "color", &mut s.color);
            ui.add(Slider::new(&mut s.digits, 0..=6).text("decimals"));
        }
        Layer::TipPoints(s) | Layer::NodePoints(s) => {
            num(ui, "size", &mut s.size, 1.0..=30.0);
            if tip_points {
                num(ui, "distance from tip", &mut s.offset, 0.0..=100.0);
            }
            shape_select(ui, &id("ps"), &mut s.shape);
            color_edit(ui, "color", &mut s.color);
            attr_select(ui, &id("pcb"), "Color by", &mut s.color_by, &attr_keys);
            palette_select(ui, &id("pp"), &mut s.palette);
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
            ui.horizontal(|ui| {
                ui.label("Text");
                ui.text_edit_singleline(&mut s.text);
            });
            num(ui, "size", &mut s.size, 4.0..=40.0);
            num(ui, "offset", &mut s.offset, 0.0..=200.0);
            num(ui, "bar width", &mut s.bar_width, 0.0..=10.0);
            color_edit(ui, "color", &mut s.color);
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
            ui.checkbox(&mut s.grid, "Grid lines");
            ui.horizontal(|ui| {
                ui.label("Title");
                ui.text_edit_singleline(&mut s.title);
            });
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
            palette_select(ui, &id("hp"), &mut s.palette);
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
            num(ui, "height (pt)", &mut s.size, 2.0..=300.0);
            ui.checkbox(&mut s.fit_spacing, "Limit to space between tips")
                .on_hover_text("Keeps silhouettes from overlapping. Untick to use the height above as-is, or zoom vertically to give them more room.");
            num(ui, "offset", &mut s.offset, 0.0..=100.0);
            ui.checkbox(&mut s.align, "Align in a column");
            let mut tint = s.tint.is_some();
            ui.checkbox(&mut tint, "Recolor");
            if tint {
                let mut c = s.tint.unwrap_or(Color::BLACK);
                color_edit(ui, "silhouette color", &mut c);
                s.tint = Some(c);
            } else {
                s.tint = None;
            }
            ui.small("Silhouettes are fetched from phylopic.org by tip name (binomial, then genus).");
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
        if ui.button("Transform…").on_hover_text("λ, κ, δ or Ornstein–Uhlenbeck transform of this clade").clicked() {
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
/// Built only while the section is open, since it grows with the square of the tip count.
fn vcv_matrix(ui: &mut egui::Ui, doc: &Document, n: NodeId, k: usize) {
    egui::CollapsingHeader::new(format!("Variance–covariance matrix ({} × {})", k, k)).id_salt(("vcv", n)).default_open(false).show(ui, |ui| {
        ui.small("Shared branch length from this node for each pair of tips; the diagonal is each tip's distance from this node (ape vcv).");
        ui.horizontal_wrapped(|ui| {
            if ui.button("Copy matrix").on_hover_text("Tab-separated, full precision; pastes into Excel").clicked() {
                let (h, rows) = crate::io::data::vcv_table(&doc.tree, n, None);
                ui.ctx().copy_text(crate::io::data::table_to_text(&h, &rows, b'\t'));
            }
            if ui.button("Save CSV…").clicked() {
                if let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).set_file_name("vcv.csv").save_file() {
                    let (h, rows) = crate::io::data::vcv_table(&doc.tree, n, None);
                    let _ = std::fs::write(path, crate::io::data::table_to_text(&h, &rows, b','));
                }
            }
        });
        if k > 150 {
            ui.weak("Too many tips to show here; copy or save the matrix instead.");
            return;
        }
        let (header, rows) = crate::io::data::vcv_table(&doc.tree, n, Some(3));
        egui::ScrollArea::both().max_height(300.0).id_salt(("vcv-scroll", n)).show(ui, |ui| {
            egui::Grid::new(("vcv-grid", n)).striped(true).show(ui, |ui| {
                for h in &header {
                    ui.label(egui::RichText::new(h.replace('_', " ")).italics().strong());
                }
                ui.end_row();
                for r in &rows {
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
pub const COLLESS_HOVER: &str = "Tree balance: higher values mean a more unbalanced, ladder-like clade, and the normalized value runs from 0 (as balanced as possible) to 1 (fully ladder-like).";
pub const GAMMA_HOVER: &str = "Tempo of diversification: γ < 0 means branching is concentrated near the root (a slowdown), γ > 0 near the tips (a speed-up), and |γ| > 1.96 departs from a constant rate at p < 0.05.";

/// Species names and their data for a clade, with copy and save buttons.
fn taxa_table(ui: &mut egui::Ui, doc: &Document, tips: &[NodeId], id: NodeId) {
    let (header, rows) = crate::io::data::tips_table(&doc.tree, tips);
    egui::CollapsingHeader::new(format!("Species & data ({})", rows.len())).id_salt(("taxa", id)).default_open(true).show(ui, |ui| {
        ui.horizontal_wrapped(|ui| {
            if ui.button("Copy names").clicked() {
                let names: Vec<&str> = rows.iter().map(|r| r[0].as_str()).collect();
                ui.ctx().copy_text(names.join("\n"));
            }
            if ui.button("Copy table").on_hover_text("Tab-separated, pastes into Excel").clicked() {
                ui.ctx().copy_text(crate::io::data::table_to_text(&header, &rows, b'\t'));
            }
            if ui.button("Save CSV…").clicked() {
                if let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).set_file_name("clade.csv").save_file() {
                    let _ = std::fs::write(path, crate::io::data::table_to_text(&header, &rows, b','));
                }
            }
        });
        egui::ScrollArea::both().max_height(260.0).id_salt(("taxa-scroll", id)).show(ui, |ui| {
            egui::Grid::new(("taxa-grid", id)).striped(true).show(ui, |ui| {
                for h in &header {
                    ui.strong(h);
                }
                ui.end_row();
                for r in &rows {
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

pub fn data_panel(ui: &mut egui::Ui, doc: &mut Document, status: &mut String) -> bool {
    let mut want_load = false;
    ui.heading("Comparative data");
    let Some(data) = doc.data.clone() else {
        ui.weak("Load a CSV/TSV with a taxon column (or drop it on the window). Rows are matched to tip labels.");
        if ui.button("Load data table…").clicked() {
            want_load = true;
        }
        return want_load;
    };
    ui.label(format!("{}: {} rows, key column \"{}\"", data.name, data.rows.len(), data.columns[data.key]));
    if ui.small_button("Replace…").clicked() {
        want_load = true;
    }
    egui::Grid::new("data-cols").num_columns(2).striped(true).show(ui, |ui| {
        for (col, kind) in data.trait_columns() {
            ui.label(format!("{} {}", col, if kind == ColumnKind::Numeric { "(num)" } else { "(cat)" }));
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
    want_load
}

pub enum PosteriorAction {
    Consensus,
    Mcc,
    OpenSample(usize),
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
        if ui.button("Majority-rule consensus").on_hover_text("Threshold < 0.5 adds compatible clades greedily (extended majority rule)").clicked() {
            out = Some(PosteriorAction::Consensus);
        }
        let mcc_help = match p.kind {
            SampleKind::Posterior => "Maximum clade credibility tree: the sampled tree whose clades have the highest product of posterior probabilities, annotated with posteriors and 95% HPD heights",
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