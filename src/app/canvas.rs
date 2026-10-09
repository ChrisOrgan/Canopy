//! Interactive tree canvas: paints a `Scene` with egui and handles
//! selection, hover, zoom/pan, the context menu and drag-and-drop
//! subtree moves (prune-and-regraft).

use super::document::Document;
use crate::phylopic::{PhyloPic, PicStatus};
use crate::scene::{self, BuildInput, Prim, Scene, SceneEnv, P};
use crate::style::Color;
use crate::tree::NodeId;
use egui::{epaint, pos2, vec2, Color32, FontFamily, FontId, Pos2, Rect, Sense, Shape, Stroke, TextureHandle};
use std::collections::HashMap;

pub const PLOT_FAMILY: &str = "plot";

pub fn c32(c: Color) -> Color32 {
    Color32::from_rgba_unmultiplied(c.r, c.g, c.b, c.a)
}

pub fn plot_font(size: f32) -> FontId {
    FontId::new(size, FontFamily::Name(PLOT_FAMILY.into()))
}

pub struct EguiEnv<'a> {
    pub ctx: &'a egui::Context,
    pub pics: &'a PhyloPic,
}

impl SceneEnv for EguiEnv<'_> {
    fn text_width(&self, text: &str, size: f32, _italic: bool) -> f32 {
        self.ctx.fonts(|f| f.layout_no_wrap(text.to_string(), plot_font(size), Color32::BLACK).size().x)
    }
    fn image_aspect(&self, key: &str) -> Option<f32> {
        match self.pics.entries.get(key) {
            Some(PicStatus::Ready(img, _)) => Some(img.width() as f32 / img.height().max(1) as f32),
            _ => None,
        }
    }
}

/// GPU textures for silhouettes: original colors and a white mask for tinting.
#[derive(Default)]
pub struct TextureCache {
    tex: HashMap<(String, bool), TextureHandle>,
}

impl TextureCache {
    fn get(&mut self, ctx: &egui::Context, pics: &PhyloPic, key: &str, mask: bool) -> Option<egui::TextureId> {
        if let Some(t) = self.tex.get(&(key.to_string(), mask)) {
            return Some(t.id());
        }
        let PicStatus::Ready(img, _) = pics.entries.get(key)? else { return None };
        let mut px = img.as_raw().clone();
        if mask {
            for p in px.chunks_mut(4) {
                p[0] = 255;
                p[1] = 255;
                p[2] = 255;
            }
        }
        let ci = egui::ColorImage::from_rgba_unmultiplied([img.width() as usize, img.height() as usize], &px);
        let h = ctx.load_texture(format!("pic:{}:{}", key, mask), ci, egui::TextureOptions::LINEAR);
        let id = h.id();
        self.tex.insert((key.to_string(), mask), h);
        Some(id)
    }

    pub fn invalidate(&mut self, key: &str) {
        self.tex.retain(|(k, _), _| k != key);
    }
}

/// Operations requested from the canvas (applied by the app with undo support).
#[derive(Clone, Debug)]
pub enum Action {
    Reroot(NodeId),
    MidpointRoot,
    Rotate(NodeId),
    Ladderize(NodeId, bool),
    ToggleCollapse(NodeId),
    ExpandAll,
    Drop(Vec<NodeId>),
    KeepOnly(Vec<NodeId>),
    Extract(NodeId),
    Highlight(NodeId),
    CladeLabel(NodeId),
    ColorClade(NodeId, Color),
    ClearCladeColor(NodeId),
    Regraft(NodeId, NodeId),
    Rename(NodeId),
    SelectClade(NodeId),
    AssignImage(NodeId),
    RefetchImage(NodeId),
    CopyClade(NodeId),
    Paste(NodeId, crate::ops::PasteMode),
    CopyNames(NodeId),
    CopyData(NodeId),
    Transform(NodeId),
}

#[derive(Default)]
pub struct CanvasState {
    pub hovered: Option<NodeId>,
    pub drag_node: Option<NodeId>,
    pub drag_start: Option<Pos2>,
    /// Origin of the last painted scene (for hit-testing dropped files).
    pub origin: Pos2,
    pub last_node_pos: Vec<Option<P>>,
    pub color_pick: Color32,
}

fn to_screen(o: Pos2, p: P) -> Pos2 {
    pos2(o.x + p[0], o.y + p[1])
}

fn dist_to_segment(p: Pos2, a: Pos2, b: Pos2) -> f32 {
    let ab = b - a;
    let t = if ab.length_sq() > 0.0 { ((p - a).dot(ab) / ab.length_sq()).clamp(0.0, 1.0) } else { 0.0 };
    (a + ab * t - p).length()
}

/// Nearest node to `p`: node markers, then label boxes, then branches.
pub fn hit_test(scene: &Scene, o: Pos2, p: Pos2) -> Option<NodeId> {
    let mut best: Option<(NodeId, f32)> = None;
    for (n, q) in scene.node_pos.iter().enumerate() {
        if let Some(q) = q {
            let d = (to_screen(o, *q) - p).length();
            if d < 7.0 && best.map(|b| d < b.1).unwrap_or(true) {
                best = Some((n, d));
            }
        }
    }
    if best.is_some() {
        return best.map(|b| b.0);
    }
    for (n, b) in &scene.label_boxes {
        let r = Rect::from_min_max(pos2(o.x + b[0], o.y + b[1]), pos2(o.x + b[2], o.y + b[3]));
        if r.contains(p) {
            return Some(*n);
        }
    }
    for (n, a, b) in &scene.branches {
        let d = dist_to_segment(p, to_screen(o, *a), to_screen(o, *b));
        if d < 4.0 && best.map(|x| d < x.1).unwrap_or(true) {
            best = Some((*n, d));
        }
    }
    best.map(|b| b.0)
}

pub fn paint_scene(painter: &egui::Painter, scene: &Scene, o: Pos2, textures: &mut TextureCache, pics: &PhyloPic) {
    let ctx = painter.ctx().clone();
    let mut shapes: Vec<Shape> = Vec::with_capacity(scene.prims.len());
    for prim in &scene.prims {
        match prim {
            Prim::Path { pts, color, width, dashed } => {
                let pts: Vec<Pos2> = pts.iter().map(|p| to_screen(o, *p)).collect();
                let stroke = Stroke::new(*width, c32(*color));
                if *dashed {
                    shapes.extend(Shape::dashed_line(&pts, stroke, width * 2.0 + 1.0, width * 2.0 + 1.0));
                } else if pts.len() == 2 {
                    shapes.push(Shape::line_segment([pts[0], pts[1]], stroke));
                } else {
                    shapes.push(Shape::line(pts, stroke));
                }
            }
            Prim::Poly { pts, fill, stroke } => {
                let pts: Vec<Pos2> = pts.iter().map(|p| to_screen(o, *p)).collect();
                let st = stroke.map(|(c, w)| Stroke::new(w, c32(c))).unwrap_or(Stroke::NONE);
                shapes.push(Shape::convex_polygon(pts, c32(*fill), st));
            }
            Prim::Strip { a, b, fill, stroke } => {
                let mut mesh = epaint::Mesh::default();
                let col = c32(*fill);
                for (pa, pb) in a.iter().zip(b) {
                    mesh.colored_vertex(to_screen(o, *pa), col);
                    mesh.colored_vertex(to_screen(o, *pb), col);
                }
                let n = a.len().min(b.len()) as u32;
                for i in 0..n.saturating_sub(1) {
                    let k = i * 2;
                    mesh.add_triangle(k, k + 1, k + 2);
                    mesh.add_triangle(k + 1, k + 3, k + 2);
                }
                shapes.push(Shape::mesh(mesh));
                if let Some((c, w)) = stroke {
                    let mut ring: Vec<Pos2> = a.iter().map(|p| to_screen(o, *p)).collect();
                    ring.extend(b.iter().rev().map(|p| to_screen(o, *p)));
                    shapes.push(Shape::closed_line(ring, Stroke::new(*w, c32(*c))));
                }
            }
            Prim::Circle { c, r, fill, stroke } => {
                let st = stroke.map(|(c, w)| Stroke::new(w, c32(c))).unwrap_or(Stroke::NONE);
                shapes.push(Shape::circle_filled(to_screen(o, *c), *r, c32(*fill)));
                if st.width > 0.0 {
                    shapes.push(Shape::circle_stroke(to_screen(o, *c), *r, st));
                }
            }
            Prim::Text { pos, text, size, color, angle, halign, valign, italic, bold: _ } => {
                let mut job = egui::text::LayoutJob::default();
                job.append(
                    text,
                    0.0,
                    egui::TextFormat { font_id: plot_font(*size), color: c32(*color), italics: *italic, ..Default::default() },
                );
                let galley = ctx.fonts(|f| f.layout_job(job));
                let sz = galley.size();
                let off = vec2(-halign * sz.x, -valign * sz.y);
                let (s, c) = angle.sin_cos();
                let rot = vec2(off.x * c - off.y * s, off.x * s + off.y * c);
                let p = to_screen(o, *pos) + rot;
                shapes.push(Shape::Text(epaint::TextShape::new(p, galley, c32(*color)).with_angle(*angle)));
            }
            Prim::Image { key, center, height, tint } => {
                if let Some(id) = textures.get(&ctx, pics, key, tint.is_some()) {
                    let aspect = match pics.entries.get(key) {
                        Some(PicStatus::Ready(img, _)) => img.width() as f32 / img.height().max(1) as f32,
                        _ => 1.0,
                    };
                    let w = height * aspect;
                    let c = to_screen(o, *center);
                    let rect = Rect::from_center_size(c, vec2(w, *height));
                    let tint = tint.map(c32).unwrap_or(Color32::WHITE);
                    let mut mesh = epaint::Mesh::with_texture(id);
                    mesh.add_rect_with_uv(rect, Rect::from_min_max(pos2(0.0, 0.0), pos2(1.0, 1.0)), tint);
                    shapes.push(Shape::mesh(mesh));
                }
            }
        }
    }
    painter.extend(shapes);
}

/// Show the canvas for `doc`; returns requested actions.
pub fn show(ui: &mut egui::Ui, doc: &mut Document, st: &mut CanvasState, textures: &mut TextureCache, pics: &mut PhyloPic, has_clip: bool) -> Vec<Action> {
    let mut actions = Vec::new();
    let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
    let rect = resp.rect;
    let bg = doc.view.background.map(c32).unwrap_or(Color32::WHITE);
    painter.rect_filled(rect, 0.0, bg);

    let size = rect.size();
    let (w, h) = (size.x * doc.zoom[0], size.y * doc.zoom[1]);
    let scene = {
        let env = EguiEnv { ctx: ui.ctx(), pics };
        scene::build(&BuildInput { tree: &doc.tree, view: &doc.view, layers: &doc.layers, overlay: doc.overlay() }, &env, w, h, 1.0)
    };
    for k in &scene.wanted_images {
        pics.request(k);
    }
    let o = rect.min + vec2(doc.pan[0], doc.pan[1]);
    st.origin = o;
    st.last_node_pos = scene.node_pos.clone();

    // Selection highlights under the tree.
    let sel_col = Color32::from_rgba_unmultiplied(255, 140, 0, 90);
    for &n in &doc.selection {
        for (bn, a, b) in &scene.branches {
            if doc.tree.is_ancestor(n, *bn) && *bn != n {
                painter.line_segment([to_screen(o, *a), to_screen(o, *b)], Stroke::new(6.0f32, sel_col));
            }
        }
        for (ln, b) in &scene.label_boxes {
            if doc.tree.is_ancestor(n, *ln) {
                let r = Rect::from_min_max(pos2(o.x + b[0], o.y + b[1]), pos2(o.x + b[2], o.y + b[3])).shrink(3.0);
                painter.rect_filled(r, 2.0, sel_col);
            }
        }
    }
    paint_scene(&painter, &scene, o, textures, pics);
    for &n in &doc.selection {
        if let Some(Some(p)) = scene.node_pos.get(n) {
            painter.circle(to_screen(o, *p), 5.0, Color32::from_rgb(255, 140, 0), Stroke::new(1.0f32, Color32::BLACK));
        }
    }

    // ---- Hover.
    let pointer = resp.hover_pos();
    st.hovered = pointer.and_then(|p| hit_test(&scene, o, p));
    if let (Some(n), Some(_)) = (st.hovered, pointer) {
        if let Some(Some(p)) = scene.node_pos.get(n) {
            painter.circle_stroke(to_screen(o, *p), 6.0, Stroke::new(1.5f32, Color32::from_rgb(0, 114, 178)));
        }
        if st.drag_node.is_none() {
            let tree = &doc.tree;
            resp.clone().on_hover_ui_at_pointer(|ui| {
                let label = tree.label(n);
                let tips = tree.tips_below(n).len();
                if tree.is_tip(n) {
                    ui.strong(label);
                } else {
                    ui.strong(format!("Clade of {} tips {}", tips, if label.is_empty() { String::new() } else { format!("({})", label) }));
                }
                if let Some(l) = tree.nodes[n].length {
                    ui.label(format!("branch length: {}", crate::tree::format_num(l, 6)));
                }
                for (k, v) in tree.nodes[n].attrs.iter().take(12) {
                    ui.label(format!("{}: {}", k, crate::io::data::table_value(k, v)));
                }
            });
        }
    }

    // ---- Click selection.
    let modifiers = ui.input(|i| i.modifiers);
    if resp.clicked() {
        match st.hovered {
            Some(n) if modifiers.shift || modifiers.command => {
                if !doc.selection.remove(&n) {
                    doc.selection.insert(n);
                }
            }
            Some(n) => {
                doc.selection.clear();
                doc.selection.insert(n);
            }
            None => doc.selection.clear(),
        }
    }
    if resp.double_clicked() {
        if let Some(n) = st.hovered {
            actions.push(Action::Rename(n));
        }
    }
    if resp.secondary_clicked() {
        if let Some(n) = st.hovered {
            if !doc.selection.contains(&n) {
                doc.selection.clear();
                doc.selection.insert(n);
            }
        }
    }

    // ---- Dragging: from a node = move subtree; elsewhere = pan.
    if resp.drag_started_by(egui::PointerButton::Primary) {
        st.drag_start = ui.input(|i| i.pointer.press_origin()).or(resp.interact_pointer_pos());
        st.drag_node = st.drag_start.and_then(|p| hit_test(&scene, o, p)).filter(|&n| n != doc.tree.root);
    }
    if resp.dragged_by(egui::PointerButton::Primary) || resp.dragged_by(egui::PointerButton::Middle) {
        match st.drag_node {
            Some(n) => {
                if let (Some(Some(p)), Some(ptr)) = (scene.node_pos.get(n), pointer) {
                    let from = to_screen(o, *p);
                    painter.line_segment([from, ptr], Stroke::new(2.0f32, Color32::from_rgb(213, 94, 0)));
                    painter.text(
                        ptr + vec2(10.0, -10.0),
                        egui::Align2::LEFT_BOTTOM,
                        match st.hovered {
                            Some(t) if t != n && !doc.tree.is_ancestor(n, t) => "drop to attach here",
                            _ => "drag onto a branch to move this clade",
                        },
                        FontId::proportional(12.0),
                        Color32::from_rgb(213, 94, 0),
                    );
                }
            }
            None => {
                let d = resp.drag_delta();
                doc.pan[0] += d.x;
                doc.pan[1] += d.y;
            }
        }
    }
    if resp.drag_stopped() {
        if let (Some(s), Some(t)) = (st.drag_node, st.hovered) {
            if s != t && !doc.tree.is_ancestor(s, t) {
                actions.push(Action::Regraft(s, t));
            }
        }
        st.drag_node = None;
        st.drag_start = None;
    }

    // ---- Scroll / zoom.
    if resp.hovered() {
        let (scroll, zoom, ptr) = ui.input(|i| (i.smooth_scroll_delta, i.zoom_delta_2d(), i.pointer.hover_pos()));
        if zoom != vec2(1.0, 1.0) {
            let ptr = ptr.unwrap_or(rect.center());
            let (zx, zy) = if modifiers.alt { (1.0, zoom.y) } else { (zoom.x, zoom.y) };
            let nzx = (doc.zoom[0] * zx).clamp(0.2, 60.0);
            let nzy = (doc.zoom[1] * zy).clamp(0.2, 200.0);
            let fx = nzx / doc.zoom[0];
            let fy = nzy / doc.zoom[1];
            doc.pan[0] = ptr.x - rect.min.x - (ptr.x - rect.min.x - doc.pan[0]) * fx;
            doc.pan[1] = ptr.y - rect.min.y - (ptr.y - rect.min.y - doc.pan[1]) * fy;
            doc.zoom = [nzx, nzy];
        } else if scroll != vec2(0.0, 0.0) {
            doc.pan[0] += scroll.x;
            doc.pan[1] += scroll.y;
        }
    }

    // ---- Context menu.
    let sel: Vec<NodeId> = doc.selection.iter().copied().collect();
    let target = doc.single_selection();
    let tree = &doc.tree;
    resp.context_menu(|ui| {
        ui.set_min_width(200.0);
        let Some(n) = target else {
            if sel.len() > 1 {
                ui.label(format!("{} nodes selected", sel.len()));
                if ui.button("Drop selected").clicked() {
                    actions.push(Action::Drop(sel.clone()));
                    ui.close_menu();
                }
                if ui.button("Keep only selected tips").clicked() {
                    actions.push(Action::KeepOnly(sel.clone()));
                    ui.close_menu();
                }
                if let Some(m) = tree.mrca(&sel) {
                    if ui.button("Select MRCA clade").clicked() {
                        actions.push(Action::SelectClade(m));
                        ui.close_menu();
                    }
                }
            } else {
                if ui.button("Midpoint root").clicked() {
                    actions.push(Action::MidpointRoot);
                    ui.close_menu();
                }
                if ui.button("Expand all collapsed clades").clicked() {
                    actions.push(Action::ExpandAll);
                    ui.close_menu();
                }
            }
            return;
        };
        let is_tip = tree.is_tip(n);
        ui.label(egui::RichText::new(if is_tip { tree.label(n).to_string() } else { format!("Clade ({} tips)", tree.tips_below(n).len()) }).strong());
        ui.separator();
        if n != tree.root && ui.button("Reroot here").clicked() {
            actions.push(Action::Reroot(n));
            ui.close_menu();
        }
        if !is_tip {
            if ui.button("Rotate (swap children)").clicked() {
                actions.push(Action::Rotate(n));
                ui.close_menu();
            }
            ui.menu_button("Ladderize", |ui| {
                if ui.button("Larger clades down").clicked() {
                    actions.push(Action::Ladderize(n, true));
                    ui.close_menu();
                }
                if ui.button("Larger clades up").clicked() {
                    actions.push(Action::Ladderize(n, false));
                    ui.close_menu();
                }
            });
            let collapsed = doc.view.collapsed.contains(&n);
            if n != tree.root && ui.button(if collapsed { "Expand clade" } else { "Collapse clade" }).clicked() {
                actions.push(Action::ToggleCollapse(n));
                ui.close_menu();
            }
            if ui.button("Highlight clade").clicked() {
                actions.push(Action::Highlight(n));
                ui.close_menu();
            }
            if ui.button("Label clade…").clicked() {
                actions.push(Action::CladeLabel(n));
                ui.close_menu();
            }
        }
        ui.menu_button("Color branches", |ui| {
            for (name, c) in [
                ("Vermillion", Color::hex("#D55E00")),
                ("Blue", Color::hex("#0072B2")),
                ("Green", Color::hex("#009E73")),
                ("Orange", Color::hex("#E69F00")),
                ("Purple", Color::hex("#CC79A7")),
                ("Sky blue", Color::hex("#56B4E9")),
                ("Grey", Color::GREY),
            ] {
                if ui.add(egui::Button::new(egui::RichText::new(format!("■ {}", name)).color(c32(c)))).clicked() {
                    actions.push(Action::ColorClade(n, c));
                    ui.close_menu();
                }
            }
            ui.horizontal(|ui| {
                ui.color_edit_button_srgba(&mut st.color_pick);
                if ui.button("Custom").clicked() {
                    let c = st.color_pick;
                    actions.push(Action::ColorClade(n, Color { r: c.r(), g: c.g(), b: c.b(), a: 255 }));
                    ui.close_menu();
                }
            });
            if ui.button("Clear color").clicked() {
                actions.push(Action::ClearCladeColor(n));
                ui.close_menu();
            }
        });
        if ui.button(if is_tip { "Rename…" } else { "Set node label…" }).clicked() {
            actions.push(Action::Rename(n));
            ui.close_menu();
        }
        if !is_tip && ui.button("Select all tips in clade").clicked() {
            actions.push(Action::SelectClade(n));
            ui.close_menu();
        }
        if !is_tip && ui.button("Open clade in new tab").clicked() {
            actions.push(Action::Extract(n));
            ui.close_menu();
        }
        ui.separator();
        if ui.button("Copy clade  (Ctrl+C)").clicked() {
            actions.push(Action::CopyClade(n));
            ui.close_menu();
        }
        ui.add_enabled_ui(has_clip, |ui| {
            ui.menu_button("Paste clade", |ui| {
                use crate::ops::PasteMode;
                if ui.button("As sister  (Ctrl+V)").clicked() {
                    actions.push(Action::Paste(n, PasteMode::Sister));
                    ui.close_menu();
                }
                if !is_tip && ui.button("As child (polytomy)").clicked() {
                    actions.push(Action::Paste(n, PasteMode::Child));
                    ui.close_menu();
                }
                if ui.button("Replace this clade").clicked() {
                    actions.push(Action::Paste(n, PasteMode::Replace));
                    ui.close_menu();
                }
            });
        });
        if ui.button("Copy species names").clicked() {
            actions.push(Action::CopyNames(n));
            ui.close_menu();
        }
        if ui.button("Copy species data (table)").clicked() {
            actions.push(Action::CopyData(n));
            ui.close_menu();
        }
        if !is_tip && ui.button("Transform branch lengths…").clicked() {
            actions.push(Action::Transform(n));
            ui.close_menu();
        }
        if is_tip {
            ui.separator();
            if ui.button("Use image file as silhouette…").clicked() {
                actions.push(Action::AssignImage(n));
                ui.close_menu();
            }
            if ui.button("Retry PhyloPic lookup").clicked() {
                actions.push(Action::RefetchImage(n));
                ui.close_menu();
            }
        }
        ui.separator();
        if n != tree.root && ui.button(if is_tip { "Drop tip" } else { "Drop clade" }).clicked() {
            actions.push(Action::Drop(vec![n]));
            ui.close_menu();
        }
    });
    actions
}
