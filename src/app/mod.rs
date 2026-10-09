//! The Canopy desktop application (egui/eframe).

pub mod canvas;
pub mod document;
pub mod export;
pub mod panels;

use crate::consensus;
use crate::io::{self, data::DataTable, newick, nexus, FileKind};
use crate::ops;
use crate::phylopic::{PhyloPic, PicStatus};
use crate::render::fonts::FontBytes;
use crate::render::raster::RasterFonts;
use crate::scene::image_key;
use crate::style::*;
use crate::tree::{NodeId, Tree};
use canvas::{Action, CanvasState, TextureCache};
use document::{Document, Posterior};
use egui::{Color32, Key, KeyboardShortcut, Modifiers};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Bundled example, available from Help › Open example tree.
const EXAMPLE: &str = include_str!("../../examples/primates.nwk");

#[derive(serde::Serialize, serde::Deserialize)]
#[serde(default)]
struct Settings {
    auto_phylopic: bool,
    recent: Vec<PathBuf>,
    /// Stored under a new key so earlier saved defaults don't override "Follow system".
    #[serde(rename = "ui_theme")]
    theme: Theme,
}

/// Interface color theme (the tree canvas stays white, like the exported figure).
#[derive(Clone, Copy, PartialEq, Eq, Debug, serde::Serialize, serde::Deserialize)]
enum Theme {
    Light,
    Dark,
    System,
}

impl Theme {
    fn apply(self, ctx: &egui::Context) {
        ctx.set_theme(match self {
            Theme::Light => egui::ThemePreference::Light,
            Theme::Dark => egui::ThemePreference::Dark,
            Theme::System => egui::ThemePreference::System,
        });
    }
}

impl Default for Settings {
    fn default() -> Self {
        Settings { auto_phylopic: true, recent: Vec::new(), theme: Theme::System }
    }
}

enum Dialog {
    Rename { node: NodeId, text: String },
    CladeLabel { node: NodeId, text: String },
    /// Branch-length transform being previewed; `base` holds the untransformed tree.
    Transform { node: NodeId, kind: ops::BranchTransform, value: f64, base: Tree },
}

pub struct CanopyApp {
    docs: Vec<Document>,
    active: usize,
    fonts: FontBytes,
    raster_fonts: RasterFonts,
    pics: PhyloPic,
    textures: TextureCache,
    canvas: CanvasState,
    export: export::ExportDialog,
    settings: Settings,
    status: String,
    search: String,
    dialog: Option<Dialog>,
    show_credits: bool,
    show_about: bool,
    ctx: egui::Context,
    /// Copied clade, and the Newick text placed on the system clipboard for it.
    clip: Option<Tree>,
    clip_text: String,
    logo: Option<egui::TextureHandle>,
}

fn install_fonts(ctx: &egui::Context, fonts: &FontBytes) {
    let mut defs = egui::FontDefinitions::default();
    defs.font_data.insert("plot-regular".into(), egui::FontData::from_owned(fonts.regular.clone()));
    let mut fam = vec!["plot-regular".to_string()];
    fam.extend(defs.families.get(&egui::FontFamily::Proportional).cloned().unwrap_or_default());
    defs.families.insert(egui::FontFamily::Name(canvas::PLOT_FAMILY.into()), fam);
    ctx.set_fonts(defs);
}

impl CanopyApp {
    pub fn new(cc: &eframe::CreationContext<'_>, files: Vec<PathBuf>) -> Self {
        let fonts = FontBytes::load();
        install_fonts(&cc.egui_ctx, &fonts);
        let raster_fonts = RasterFonts::new(&fonts).expect("no usable font found");
        let settings: Settings = cc.storage.and_then(|s| eframe::get_value(s, eframe::APP_KEY)).unwrap_or_default();
        settings.theme.apply(&cc.egui_ctx);
        let ctx = cc.egui_ctx.clone();
        let mut pics = PhyloPic::new(Arc::new(move || ctx.request_repaint()));
        pics.enabled = settings.auto_phylopic;
        let mut app = CanopyApp {
            docs: Vec::new(),
            active: 0,
            fonts,
            raster_fonts,
            pics,
            textures: TextureCache::default(),
            canvas: CanvasState { color_pick: Color32::from_rgb(213, 94, 0), ..Default::default() },
            export: export::ExportDialog::default(),
            settings,
            status: "Open a tree (Ctrl+O) or drop Newick/NEXUS files, CSV trait tables or images onto the window.".into(),
            search: String::new(),
            dialog: None,
            show_credits: false,
            show_about: false,
            ctx: cc.egui_ctx.clone(),
            clip: None,
            clip_text: String::new(),
            logo: None,
        };
        for f in files {
            app.open_path(&f);
        }
        app
    }

    fn doc(&mut self) -> Option<&mut Document> {
        self.docs.get_mut(self.active)
    }

    fn add_doc(&mut self, d: Document) {
        self.docs.push(d);
        self.active = self.docs.len() - 1;
    }

    fn remember(&mut self, path: &Path) {
        self.settings.recent.retain(|p| p != path);
        self.settings.recent.insert(0, path.to_path_buf());
        self.settings.recent.truncate(10);
    }

    pub fn open_path(&mut self, path: &Path) {
        let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        match io::classify(path) {
            FileKind::Tree | FileKind::Unknown => match io::read_trees(path) {
                Ok(mut trees) => {
                    self.remember(path);
                    let n = trees.len();
                    if n == 1 {
                        let mut d = Document::new(&name, trees.remove(0), self.settings.auto_phylopic);
                        d.path = Some(path.to_path_buf());
                        self.add_doc(d);
                        self.status = format!("Opened {}", name);
                    } else {
                        let first = trees[0].clone();
                        let mut d = Document::new(&name, first, self.settings.auto_phylopic);
                        d.path = Some(path.to_path_buf());
                        d.posterior = Some(Posterior::new(trees));
                        d.annotate_posterior();
                        self.add_doc(d);
                        self.status = format!("Opened {} with {} trees: use the posterior panel for consensus / MCC", name, n);
                    }
                }
                Err(e) => self.status = format!("Could not read {}: {:#}", name, e),
            },
            FileKind::Data => match DataTable::from_path(path) {
                Ok(table) => self.attach_data(table),
                Err(e) => self.status = format!("Could not read table {}: {:#}", name, e),
            },
            FileKind::Project => match Document::load_project(path) {
                Ok(d) => {
                    self.remember(path);
                    self.add_doc(d);
                    self.status = format!("Opened project {}", name);
                }
                Err(e) => self.status = format!("Could not open project: {:#}", e),
            },
            FileKind::Image => self.status = "Drop an image onto a tip to use it as that taxon's silhouette.".into(),
        }
    }

    fn attach_data(&mut self, table: DataTable) {
        let Some(doc) = self.doc() else { return };
        doc.checkpoint();
        let (m, unmatched) = table.join_to_tree(&mut doc.tree);
        let tips = doc.tree.num_tips();
        let name = table.name.clone();
        doc.data = Some(table);
        self.status = format!(
            "Joined {}: {} of {} tips matched{}",
            name,
            m,
            tips,
            if unmatched.is_empty() { String::new() } else { format!("; unmatched rows: {}", unmatched.iter().take(5).cloned().collect::<Vec<_>>().join(", ")) }
        );
    }

    fn assign_image(&mut self, node: NodeId, path: &Path) {
        let Some(doc) = self.docs.get_mut(self.active) else { return };
        match image::open(path) {
            Ok(img) => {
                let key = image_key(doc.tree.label(node));
                self.pics.set_custom(&key, img.to_rgba8());
                self.textures.invalidate(&key);
                doc.view.tip_images.insert(key.clone(), path.display().to_string());
                if !doc.layers.iter().any(|l| matches!(l.layer, Layer::Phylopic(_))) {
                    doc.add_layer(Layer::default_phylopic());
                }
                self.status = format!("Image assigned to {}", key);
            }
            Err(e) => self.status = format!("Could not open image: {}", e),
        }
    }

    fn handle_dropped(&mut self, ctx: &egui::Context) {
        let (files, pos) = ctx.input(|i| (i.raw.dropped_files.clone(), i.pointer.hover_pos()));
        for f in files {
            let Some(path) = f.path else { continue };
            if io::classify(&path) == FileKind::Image {
                // Assign to the tip under the pointer, or the selected tip.
                let target = pos
                    .and_then(|p| {
                        let doc = self.docs.get(self.active)?;
                        let o = self.canvas.origin;
                        doc.tree.tips().into_iter().find(|&t| {
                            self.canvas.last_node_pos.get(t).copied().flatten().map(|q| (egui::pos2(o.x + q[0], o.y + q[1]) - p).length() < 40.0).unwrap_or(false)
                        })
                    })
                    .or_else(|| self.docs.get(self.active).and_then(|d| d.single_selection()).filter(|&n| self.docs[self.active].tree.is_tip(n)));
                match target {
                    Some(n) => self.assign_image(n, &path),
                    None => self.status = "Drop the image onto a tip (or select a tip first).".into(),
                }
            } else {
                self.open_path(&path);
            }
        }
    }

    fn apply(&mut self, a: Action) {
        let auto = self.settings.auto_phylopic;
        if self.docs.get(self.active).is_none() {
            return;
        }
        // Actions that need app-level state besides the document.
        match a {
            Action::Extract(n) => {
                let doc = &self.docs[self.active];
                let mut d = Document::new(&format!("{} (clade)", doc.name), doc.tree.extract(n), auto);
                d.view.layout = doc.view.layout.clone();
                self.add_doc(d);
                return;
            }
            Action::CladeLabel(n) => {
                let doc = &self.docs[self.active];
                let text = doc.tree.nodes[n].label.clone().filter(|l| l.parse::<f64>().is_err()).unwrap_or_default();
                self.dialog = Some(Dialog::CladeLabel { node: n, text });
                return;
            }
            Action::Rename(n) => {
                let text = self.docs[self.active].tree.nodes[n].label.clone().unwrap_or_default();
                self.dialog = Some(Dialog::Rename { node: n, text });
                return;
            }
            Action::AssignImage(n) => {
                if let Some(p) = rfd::FileDialog::new().add_filter("Image", &["png", "jpg", "jpeg"]).pick_file() {
                    self.assign_image(n, &p);
                }
                return;
            }
            Action::CopyClade(n) => {
                let doc = &self.docs[self.active];
                let mut clip = doc.tree.extract(n);
                clip.nodes[clip.root].length = doc.tree.nodes[n].length;
                self.clip_text = newick::write_newick(&clip, &newick::WriteOptions { annotations: true, ..Default::default() });
                self.ctx.copy_text(self.clip_text.clone());
                self.status = format!("Copied clade of {} tips (also on the clipboard as Newick). Select a node in any tree and paste.", clip.num_tips());
                self.clip = Some(clip);
                return;
            }
            Action::Paste(n, mode) => {
                let Some(clip) = self.clip.clone() else { return };
                let doc = &mut self.docs[self.active];
                let existing: std::collections::HashSet<String> = doc.tree.tips().into_iter().map(|t| doc.tree.label(t).to_string()).collect();
                let dups: Vec<String> = clip.tips().into_iter().map(|t| clip.label(t).to_string()).filter(|l| existing.contains(l)).collect();
                doc.checkpoint();
                match ops::graft(&mut doc.tree, n, &clip, mode) {
                    Ok(s) => {
                        doc.selection.clear();
                        doc.selection.insert(s);
                        self.status = if dups.is_empty() {
                            format!("Pasted clade of {} tips", clip.num_tips())
                        } else {
                            format!("Pasted clade; note these tips now appear twice: {}", dups.iter().take(6).cloned().collect::<Vec<_>>().join(", "))
                        };
                    }
                    Err(e) => {
                        doc.undo();
                        self.status = format!("{:#}", e);
                    }
                }
                return;
            }
            Action::CopyNames(n) | Action::CopyData(n) => {
                let doc = &self.docs[self.active];
                let tips = doc.tree.tips_below(n);
                let (header, rows) = io::data::tips_table(&doc.tree, &tips);
                let text = if matches!(a, Action::CopyNames(_)) {
                    rows.iter().map(|r| r[0].clone()).collect::<Vec<_>>().join("\n")
                } else {
                    io::data::table_to_text(&header, &rows, b'\t')
                };
                self.ctx.copy_text(text);
                self.status = format!("Copied {} species to the clipboard", rows.len());
                return;
            }
            Action::Transform(n) => {
                let doc = &mut self.docs[self.active];
                if !doc.tree.has_lengths() {
                    self.status = "This tree has no branch lengths to transform.".into();
                    return;
                }
                doc.checkpoint();
                self.dialog = Some(Dialog::Transform { node: n, kind: ops::BranchTransform::Lambda, value: 1.0, base: doc.tree.clone() });
                return;
            }
            Action::RefetchImage(n) => {
                let key = image_key(self.docs[self.active].tree.label(n));
                self.pics.entries.remove(&key);
                self.textures.invalidate(&key);
                if let Some(d) = dirs::cache_dir() {
                    let stem: String = key.chars().map(|c| if c.is_ascii_alphanumeric() { c.to_ascii_lowercase() } else { '_' }).collect();
                    let _ = std::fs::remove_file(d.join("Canopy").join("phylopic").join(format!("{}.none", stem)));
                }
                self.pics.request(&key);
                return;
            }
            _ => {}
        }
        let doc = &mut self.docs[self.active];
        let result: anyhow::Result<()> = (|| {
            match a {
                Action::Extract(_)
                | Action::CladeLabel(_)
                | Action::Rename(_)
                | Action::AssignImage(_)
                | Action::RefetchImage(_)
                | Action::CopyClade(_)
                | Action::Paste(..)
                | Action::CopyNames(_)
                | Action::CopyData(_)
                | Action::Transform(_) => {}
                Action::Reroot(n) => {
                    doc.checkpoint();
                    let r = ops::reroot(&mut doc.tree, n, 0.5)?;
                    doc.selection.clear();
                    doc.selection.insert(r);
                }
                Action::MidpointRoot => {
                    doc.checkpoint();
                    ops::midpoint_root(&mut doc.tree)?;
                }
                Action::Rotate(n) => {
                    doc.checkpoint();
                    ops::rotate(&mut doc.tree, n);
                }
                Action::Ladderize(n, right) => {
                    doc.checkpoint();
                    ops::ladderize(&mut doc.tree, n, right);
                }
                Action::ToggleCollapse(n) => {
                    doc.checkpoint();
                    if !doc.view.collapsed.remove(&n) {
                        doc.view.collapsed.insert(n);
                    }
                }
                Action::ExpandAll => {
                    doc.checkpoint();
                    doc.view.collapsed.clear();
                }
                Action::Drop(nodes) => {
                    doc.checkpoint();
                    if let Err(e) = ops::drop_nodes(&mut doc.tree, &nodes) {
                        doc.undo();
                        return Err(e);
                    }
                    doc.selection.clear();
                }
                Action::KeepOnly(nodes) => {
                    let tips: Vec<NodeId> = nodes.iter().flat_map(|&n| doc.tree.tips_below(n)).collect();
                    doc.checkpoint();
                    if let Err(e) = ops::keep_tips(&mut doc.tree, &tips) {
                        doc.undo();
                        return Err(e);
                    }
                    doc.selection.clear();
                }
                Action::Highlight(n) => {
                    let k = doc.layers.iter().filter(|l| matches!(l.layer, Layer::Highlight(_))).count();
                    let fill = Palette::OkabeIto.discrete(k + 1, 8).with_alpha(70);
                    doc.add_layer(Layer::Highlight(HighlightStyle { node: n, fill, extend: 0.0 }));
                }
                Action::ColorClade(n, c) => {
                    doc.checkpoint();
                    doc.view.clade_colors.insert(n, c);
                }
                Action::ClearCladeColor(n) => {
                    doc.checkpoint();
                    let sub: Vec<NodeId> = doc.tree.preorder_from(n);
                    for s in sub {
                        doc.view.clade_colors.remove(&s);
                    }
                }
                Action::Regraft(s, t) => {
                    doc.checkpoint();
                    if let Err(e) = ops::regraft(&mut doc.tree, s, t) {
                        doc.undo();
                        return Err(e);
                    }
                }
                Action::SelectClade(n) => {
                    doc.selection.clear();
                    doc.selection.insert(n);
                }
            }
            Ok(())
        })();
        if let Err(e) = result {
            self.status = format!("{:#}", e);
        }
    }

    fn open_dialog(&mut self) {
        let files = rfd::FileDialog::new()
            .add_filter("Trees, data and projects", &["nwk", "newick", "tre", "tree", "trees", "nex", "nexus", "nxs", "t", "con", "treefile", "contree", "txt", "csv", "tsv", "canopy"])
            .add_filter("All files", &["*"])
            .pick_files();
        for f in files.unwrap_or_default() {
            self.open_path(&f);
        }
    }

    fn save_tree(&mut self, nexus_fmt: bool) {
        let Some(doc) = self.docs.get(self.active) else { return };
        let ext = if nexus_fmt { "nex" } else { "nwk" };
        let Some(path) = rfd::FileDialog::new().add_filter(if nexus_fmt { "NEXUS" } else { "Newick" }, &[ext]).set_file_name(format!("tree.{}", ext)).save_file() else { return };
        let t = doc.tree.compact();
        let opts = newick::WriteOptions { annotations: nexus_fmt, ..Default::default() };
        let text = if nexus_fmt { nexus::write_nexus(&[&t], &opts) } else { newick::write_newick(&t, &opts) + "\n" };
        self.status = match std::fs::write(&path, text) {
            Ok(_) => format!("Saved {}", path.display()),
            Err(e) => format!("Save failed: {}", e),
        };
    }

    fn save_project(&mut self) {
        let Some(doc) = self.docs.get_mut(self.active) else { return };
        let default = doc.path.as_ref().and_then(|p| p.file_stem()).map(|s| format!("{}.canopy", s.to_string_lossy())).unwrap_or("tree.canopy".into());
        let Some(path) = rfd::FileDialog::new().add_filter("Canopy project", &["canopy"]).set_file_name(default).save_file() else { return };
        self.status = match doc.save_project(&path) {
            Ok(_) => {
                doc.dirty = false;
                format!("Saved project {}", path.display())
            }
            Err(e) => format!("Save failed: {:#}", e),
        };
    }

    fn posterior_action(&mut self, a: panels::PosteriorAction) {
        let auto = self.settings.auto_phylopic;
        let Some(doc) = self.docs.get_mut(self.active) else { return };
        let Some(p) = doc.posterior.as_mut() else { return };
        let started = std::time::Instant::now();
        let kept = p.kept();
        let res: anyhow::Result<(String, Tree)> = (|| match a {
            panels::PosteriorAction::Consensus => {
                let s = consensus::summarize(kept, p.rooted)?;
                let t = consensus::consensus(&s, p.threshold as f64, p.heights);
                Ok((format!("Consensus {:.0}%", p.threshold * 100.0), t))
            }
            panels::PosteriorAction::Mcc => {
                let s = consensus::summarize(kept, p.rooted)?;
                let t = consensus::mcc_tree(kept, &s, p.heights)?;
                Ok(("MCC".to_string(), t))
            }
            panels::PosteriorAction::OpenSample(i) => Ok((format!("sample {}", i + 1), p.trees[i].clone())),
        })();
        match res {
            Ok((label, t)) => {
                p.info = format!("{} from {} trees in {:.2}s", label, kept.len(), started.elapsed().as_secs_f32());
                let name = format!("{} – {}", doc.name, label);
                let mut d = Document::new(&name, t, auto);
                d.view.layout = doc.view.layout.clone();
                // Keep the sample (shared, not copied) so this tab can show a DensiTree.
                let mut post = Posterior::shared(p.trees.clone());
                post.burnin = p.burnin;
                post.rooted = p.rooted;
                d.posterior = Some(post);
                self.status = format!("Created {}", name);
                self.add_doc(d);
            }
            Err(e) => self.status = format!("{:#}", e),
        }
    }

    fn shortcuts(&mut self, ctx: &egui::Context) {
        let cmd = Modifiers::COMMAND;
        let sc = |m, k| KeyboardShortcut::new(m, k);
        if ctx.input_mut(|i| i.consume_shortcut(&sc(cmd, Key::O))) {
            self.open_dialog();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&sc(cmd, Key::S))) {
            self.save_project();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&sc(cmd, Key::E))) {
            self.export.open = true;
        }
        if ctx.input_mut(|i| i.consume_shortcut(&sc(cmd | Modifiers::SHIFT, Key::Z)) || i.consume_shortcut(&sc(cmd, Key::Y))) {
            if let Some(d) = self.doc() {
                d.redo();
            }
        }
        if ctx.input_mut(|i| i.consume_shortcut(&sc(cmd, Key::Z))) {
            if let Some(d) = self.doc() {
                d.undo();
            }
        }
        let typing = ctx.wants_keyboard_input();
        if !typing {
            if ctx.input(|i| i.key_pressed(Key::Delete)) {
                if let Some(d) = self.docs.get(self.active) {
                    let sel: Vec<NodeId> = d.selection.iter().copied().collect();
                    if !sel.is_empty() {
                        self.apply(Action::Drop(sel));
                    }
                }
            }
            if ctx.input(|i| i.key_pressed(Key::Escape)) {
                if let Some(d) = self.doc() {
                    d.selection.clear();
                }
            }
            // Step through posterior samples with the arrow keys.
            let (left, right) = ctx.input(|i| (i.key_pressed(Key::ArrowLeft), i.key_pressed(Key::ArrowRight)));
            if left || right {
                if let Some(d) = self.doc() {
                    if let Some(p) = &d.posterior {
                        let n = p.trees.len();
                        let i = if !p.browsing() {
                            p.browse
                        } else if right {
                            (p.browse + 1).min(n.saturating_sub(1))
                        } else {
                            p.browse.saturating_sub(1)
                        };
                        d.show_sample(i);
                    }
                }
            }
            // Clipboard: Ctrl+C copies the selected clade; Ctrl+V pastes it
            // (or Newick text from another program) next to the selection.
            let (copy, paste) = ctx.input(|i| {
                let copy = i.events.iter().any(|e| matches!(e, egui::Event::Copy));
                let paste = i.events.iter().find_map(|e| if let egui::Event::Paste(s) = e { Some(s.clone()) } else { None });
                (copy, paste)
            });
            let selected = self.docs.get(self.active).and_then(|d| d.single_selection());
            if copy {
                if let Some(n) = selected {
                    self.apply(Action::CopyClade(n));
                }
            }
            if let Some(text) = paste {
                let text = text.trim().to_string();
                if !text.is_empty() && text != self.clip_text.trim() {
                    match newick::parse_newick(&text) {
                        Ok(t) => {
                            self.clip_text = text.clone();
                            self.clip = Some(t);
                        }
                        Err(_) => {
                            self.status = "The clipboard does not contain a Newick tree.".into();
                            return;
                        }
                    }
                }
                match (selected, &self.clip) {
                    (Some(n), Some(_)) => self.apply(Action::Paste(n, ops::PasteMode::Sister)),
                    (None, Some(t)) => {
                        let d = Document::new("Pasted tree", t.clone(), self.settings.auto_phylopic);
                        self.add_doc(d);
                    }
                    _ => {}
                }
            }
        }
    }

    fn transform_dialog(&mut self, ctx: &egui::Context) {
        let Some(Dialog::Transform { node, kind, value, base }) = &mut self.dialog else { return };
        let whole = Some(*node) == self.docs.get(self.active).map(|d| d.tree.root);
        let (mut changed, mut apply, mut cancel) = (false, false, false);
        egui::Window::new("Transform branch lengths").collapsible(false).resizable(false).show(ctx, |ui| {
            ui.label(if whole { "Scope: whole tree".to_string() } else { format!("Scope: clade of {} tips", base.tips_below(*node).len()) });
            ui.horizontal(|ui| {
                for k in ops::BranchTransform::ALL {
                    if ui.selectable_value(kind, k, k.name()).changed() {
                        *value = k.identity();
                        changed = true;
                    }
                }
            });
            // α is in 1/(branch-length units); offer up to 20 / clade height.
            let height = {
                let d = base.depths(true);
                base.tips_below(*node).iter().map(|&t| d[t] - d[*node]).fold(0.0, f64::max).max(1e-9)
            };
            let slider = match kind {
                ops::BranchTransform::Lambda => egui::Slider::new(value, 0.0..=1.0),
                ops::BranchTransform::Kappa => egui::Slider::new(value, 0.0..=2.0),
                ops::BranchTransform::Delta => egui::Slider::new(value, 0.05..=5.0),
                ops::BranchTransform::OrnsteinUhlenbeck => egui::Slider::new(value, 0.0..=20.0 / height).logarithmic(true),
            };
            changed |= ui.add(slider.text("value")).changed();
            ui.small(kind.description());
            ui.add_space(4.0);
            ui.label(egui::RichText::new("References").strong().small());
            ui.small(kind.reference());
            ui.small(ops::TRANSFORM_IMPLEMENTATION_REF);
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                if ui.button("Apply").clicked() {
                    apply = true;
                }
                if ui.button("Cancel").clicked() {
                    cancel = true;
                }
                if ui.button("Reset").on_hover_text("Back to the untransformed tree").clicked() {
                    *value = kind.identity();
                    changed = true;
                }
            });
        });
        let job = changed.then(|| (*node, *kind, *value, base.clone()));
        if let Some((n, k, v, base)) = job {
            if let Some(doc) = self.docs.get_mut(self.active) {
                doc.tree = base;
                if let Err(e) = ops::transform_branches(&mut doc.tree, n, k, v) {
                    self.status = format!("{:#}", e);
                }
            }
        }
        if cancel {
            if let Some(doc) = self.docs.get_mut(self.active) {
                doc.undo();
            }
            self.dialog = None;
        } else if apply {
            if let Some(Dialog::Transform { kind, value, .. }) = &self.dialog {
                self.status = format!("Applied {} = {}", kind.name(), crate::tree::format_num(*value, 3));
            }
            self.dialog = None;
        }
    }

    fn menu_bar(&mut self, ui: &mut egui::Ui) {
        egui::menu::bar(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.button("Open…  (Ctrl+O)").clicked() {
                    ui.close_menu();
                    self.open_dialog();
                }
                ui.menu_button("Open recent", |ui| {
                    for p in self.settings.recent.clone() {
                        if ui.button(p.display().to_string()).clicked() {
                            ui.close_menu();
                            self.open_path(&p);
                        }
                    }
                });
                if ui.button("Load trait data…").clicked() {
                    ui.close_menu();
                    if let Some(p) = rfd::FileDialog::new().add_filter("Table", &["csv", "tsv", "txt", "tab"]).pick_file() {
                        self.open_path(&p);
                    }
                }
                ui.separator();
                if ui.button("Save project…  (Ctrl+S)").clicked() {
                    ui.close_menu();
                    self.save_project();
                }
                if ui.button("Save tree as Newick…").clicked() {
                    ui.close_menu();
                    self.save_tree(false);
                }
                if ui.button("Save tree as NEXUS (with annotations)…").clicked() {
                    ui.close_menu();
                    self.save_tree(true);
                }
                ui.separator();
                if ui.button("Export figure…  (Ctrl+E)").clicked() {
                    ui.close_menu();
                    self.export.open = true;
                }
                ui.separator();
                if ui.button("Close tab").clicked() {
                    ui.close_menu();
                    if !self.docs.is_empty() {
                        self.docs.remove(self.active);
                        self.active = self.active.min(self.docs.len().saturating_sub(1));
                    }
                }
                if ui.button("Quit").clicked() {
                    ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
                }
            });
            ui.menu_button("Edit", |ui| {
                let (cu, cr) = self.docs.get(self.active).map(|d| (d.can_undo(), d.can_redo())).unwrap_or((false, false));
                if ui.add_enabled(cu, egui::Button::new("Undo  (Ctrl+Z)")).clicked() {
                    if let Some(d) = self.doc() {
                        d.undo();
                    }
                    ui.close_menu();
                }
                if ui.add_enabled(cr, egui::Button::new("Redo  (Ctrl+Y)")).clicked() {
                    if let Some(d) = self.doc() {
                        d.redo();
                    }
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Select all tips").clicked() {
                    if let Some(d) = self.doc() {
                        d.selection = d.tree.tips().into_iter().collect();
                    }
                    ui.close_menu();
                }
                if ui.button("Clear selection (Esc)").clicked() {
                    if let Some(d) = self.doc() {
                        d.selection.clear();
                    }
                    ui.close_menu();
                }
            });
            ui.menu_button("Tree", |ui| {
                if ui.button("Midpoint root").clicked() {
                    self.apply(Action::MidpointRoot);
                    ui.close_menu();
                }
                if ui.button("Ladderize (larger clades down)").clicked() {
                    if let Some(r) = self.docs.get(self.active).map(|d| d.tree.root) {
                        self.apply(Action::Ladderize(r, true));
                    }
                    ui.close_menu();
                }
                if ui.button("Ladderize (larger clades up)").clicked() {
                    if let Some(r) = self.docs.get(self.active).map(|d| d.tree.root) {
                        self.apply(Action::Ladderize(r, false));
                    }
                    ui.close_menu();
                }
                if ui.button("Expand all collapsed").clicked() {
                    self.apply(Action::ExpandAll);
                    ui.close_menu();
                }
                if ui.button("Transform branch lengths (λ, κ, δ, OU)…").clicked() {
                    if let Some(r) = self.docs.get(self.active).map(|d| d.tree.root) {
                        self.apply(Action::Transform(r));
                    }
                    ui.close_menu();
                }
                if ui.button("Collapse weak nodes (support < 50)…").clicked() {
                    if let Some(d) = self.doc() {
                        d.checkpoint();
                        let max_support = d.tree.preorder().iter().filter_map(|&n| d.tree.nodes[n].attrs.get("support").and_then(|a| a.as_f64())).fold(0.0, f64::max);
                        let th = if max_support <= 1.0 { 0.5 } else { 50.0 };
                        let k = ops::collapse_weak(&mut d.tree, 0.0, Some(th));
                        self.status = format!("Collapsed {} nodes with support < {}", k, th);
                    }
                    ui.close_menu();
                }
                if ui.button("Clear all clade colors").clicked() {
                    if let Some(d) = self.doc() {
                        d.checkpoint();
                        d.view.clade_colors.clear();
                    }
                    ui.close_menu();
                }
            });
            ui.menu_button("View", |ui| {
                ui.label("Theme");
                let before = self.settings.theme;
                ui.radio_value(&mut self.settings.theme, Theme::Light, "Light");
                ui.radio_value(&mut self.settings.theme, Theme::Dark, "Dark");
                ui.radio_value(&mut self.settings.theme, Theme::System, "Follow system");
                if self.settings.theme != before {
                    self.settings.theme.apply(ui.ctx());
                    ui.close_menu();
                }
            });
            ui.menu_button("PhyloPic", |ui| {
                if ui.checkbox(&mut self.settings.auto_phylopic, "Fetch silhouettes automatically").changed() {
                    self.pics.enabled = self.settings.auto_phylopic;
                }
                if ui.button("Show silhouettes for this tree").clicked() {
                    if let Some(d) = self.doc() {
                        if !d.layers.iter().any(|l| matches!(l.layer, Layer::Phylopic(_))) {
                            d.add_layer(Layer::default_phylopic());
                        }
                        for e in d.layers.iter_mut() {
                            if matches!(e.layer, Layer::Phylopic(_)) {
                                e.enabled = true;
                            }
                        }
                    }
                    self.pics.enabled = true;
                    ui.close_menu();
                }
                if ui.button("Retry failed lookups").clicked() {
                    self.pics.retry_failed();
                    ui.close_menu();
                }
                if ui.button("Image credits…").clicked() {
                    self.show_credits = true;
                    ui.close_menu();
                }
            });
            ui.menu_button("Help", |ui| {
                if ui.button("Open example tree").clicked() {
                    if let Ok(t) = newick::parse_newick(EXAMPLE) {
                        self.add_doc(Document::new("Example: primates", t, self.settings.auto_phylopic));
                    }
                    ui.close_menu();
                }

                if ui.button("About Canopy").clicked() {
                    self.show_about = true;
                    ui.close_menu();
                }
            });
        });
    }

    fn tabs(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            let mut close = None;
            for (i, d) in self.docs.iter().enumerate() {
                let label = format!("{}{}", d.name, if d.dirty { " •" } else { "" });
                if ui.selectable_label(i == self.active, label).clicked() {
                    self.active = i;
                }
                if ui.small_button("×").clicked() {
                    close = Some(i);
                }
                ui.separator();
            }
            if let Some(i) = close {
                self.docs.remove(i);
                if self.active >= i && self.active > 0 {
                    self.active -= 1;
                }
            }
        });
    }

    fn dialogs(&mut self, ctx: &egui::Context) {
        let mut done = false;
        if matches!(self.dialog, Some(Dialog::Transform { .. })) {
            self.transform_dialog(ctx);
        } else if let Some(dialog) = &mut self.dialog {
            let (title, node, text) = match dialog {
                Dialog::Rename { node, text } => ("Rename", *node, text),
                Dialog::CladeLabel { node, text } => ("Clade label", *node, text),
                Dialog::Transform { .. } => unreachable!(),
            };
            let mut ok = false;
            egui::Window::new(title).collapsible(false).resizable(false).show(ctx, |ui| {
                let r = ui.text_edit_singleline(text);
                r.request_focus();
                if r.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                    ok = true;
                }
                ui.horizontal(|ui| {
                    if ui.button("OK").clicked() {
                        ok = true;
                    }
                    if ui.button("Cancel").clicked() {
                        done = true;
                    }
                });
            });
            if ok {
                let value = text.clone();
                let is_label = matches!(dialog, Dialog::CladeLabel { .. });
                if let Some(doc) = self.docs.get_mut(self.active) {
                    if is_label {
                        doc.add_layer(Layer::CladeLabel(CladeLabelStyle { node, text: value, color: Color::BLACK, size: 11.0, offset: 6.0, bar_width: 2.0 }));
                    } else {
                        doc.checkpoint();
                        doc.tree.nodes[node].label = if value.is_empty() { None } else { Some(value) };
                    }
                }
                done = true;
            }
        }
        if done {
            self.dialog = None;
        }

        if self.show_credits {
            let keys: Vec<String> = self.docs.get(self.active).map(|d| d.tree.tips().into_iter().map(|t| image_key(d.tree.label(t))).collect()).unwrap_or_default();
            let credits = self.pics.credits(&keys);
            let mut open = true;
            egui::Window::new("PhyloPic image credits").open(&mut open).default_width(520.0).show(ctx, |ui| {
                ui.label("Silhouettes from PhyloPic (phylopic.org). Credit contributors as their licenses require:");
                let text: String = credits
                    .iter()
                    .map(|m| format!("{} — {} ({}), {}", m.query, m.attribution.clone().unwrap_or_else(|| "anonymous".into()), m.license_name(), m.page_url()))
                    .collect::<Vec<_>>()
                    .join("\n");
                egui::ScrollArea::vertical().max_height(320.0).show(ui, |ui| {
                    for m in &credits {
                        ui.horizontal_wrapped(|ui| {
                            ui.strong(&m.query);
                            if m.matched_name.to_lowercase() != m.query.to_lowercase() {
                                ui.weak(format!("(as {})", m.matched_name));
                            }
                            ui.label(format!("— {} · {}", m.attribution.clone().unwrap_or_else(|| "anonymous".into()), m.license_name()));
                            ui.hyperlink_to("source", m.page_url());
                        });
                    }
                });
                if ui.button("Copy as text").clicked() {
                    ui.ctx().copy_text(text);
                }
            });
            self.show_credits = open;
        }
        if self.show_about {
            let mut open = true;
            let logo = self
                .logo
                .get_or_insert_with(|| {
                    let img = crate::logo::render(256, false);
                    let ci = egui::ColorImage::from_rgba_unmultiplied([img.width() as usize, img.height() as usize], img.as_raw());
                    ctx.load_texture("canopy-logo", ci, egui::TextureOptions::LINEAR)
                })
                .clone();
            egui::Window::new("About Canopy").open(&mut open).show(ctx, |ui| {
                ui.add(egui::Image::new((logo.id(), egui::vec2(128.0, 128.0))));
                ui.heading(format!("Canopy {}", env!("CARGO_PKG_VERSION")));
                ui.label("Phylogenetic tree visualization in the spirit of ggtree.");
                ui.add_space(6.0);
                ui.strong("Credits");
                ui.label("Created by Chris Organ and Claude Code, 2026.");
                ui.add_space(6.0);
                ui.label(format!("Figure font: {}", self.fonts.name));
                ui.hyperlink_to("Silhouettes: PhyloPic", "https://www.phylopic.org/");
            });
            self.show_about = open;
        }
    }
}

impl eframe::App for CanopyApp {
    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        eframe::set_value(storage, eframe::APP_KEY, &self.settings);
    }

    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if self.pics.poll() {
            ctx.request_repaint();
        }
        self.handle_dropped(ctx);
        self.shortcuts(ctx);

        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            self.menu_bar(ui);
            self.tabs(ui);
        });
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(&self.status);
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let (ready, pending, missing, failed) = self.pics.counts();
                    if pending > 0 {
                        ui.spinner();
                    }
                    if ready + pending + missing + failed > 0 {
                        ui.weak(format!("PhyloPic: {} found, {} pending, {} not found{}", ready, pending, missing, if failed > 0 { format!(", {} failed", failed) } else { String::new() }));
                    }
                });
            });
        });

        let mut post_action = None;
        if let Some(doc) = self.docs.get_mut(self.active) {
            if doc.posterior.is_some() {
                egui::TopBottomPanel::bottom("posterior").show(ctx, |ui| {
                    post_action = panels::posterior_panel(ui, doc);
                });
            }
        }
        if let Some(a) = post_action {
            self.posterior_action(a);
        }

        let mut actions = Vec::new();
        let mut load_data = false;
        if let Some(doc) = self.docs.get_mut(self.active) {
            egui::SidePanel::left("left").resizable(true).default_width(300.0).show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    egui::CollapsingHeader::new("Layout").default_open(true).show(ui, |ui| panels::layout_panel(ui, doc));
                    ui.separator();
                    panels::layers_panel(ui, doc);
                    ui.separator();
                    panels::layer_editor(ui, doc);
                });
            });
            egui::SidePanel::right("right").resizable(true).default_width(300.0).show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    actions.extend(panels::inspector(ui, doc, &mut self.search));
                    ui.separator();
                    load_data = panels::data_panel(ui, doc, &mut self.status);
                    ui.separator();
                    ui.collapsing("Tree summary", |ui| {
                        let t = &doc.tree;
                        ui.label(format!("{} tips, {} internal nodes", t.num_tips(), t.preorder().len() - t.num_tips()));
                        ui.label(format!("Branch lengths: {}", if t.has_lengths() { "yes" } else { "no" }));
                        ui.label(format!("Rooted: {}", match t.rooted { Some(true) => "yes", Some(false) => "no ([&U])", None => "unspecified" }));
                        let h = t.heights();
                        ui.label(format!("Root height: {}", crate::tree::format_num(h[t.root], 5)));
                        let keys = t.attr_keys();
                        if !keys.is_empty() {
                            ui.label(format!("Annotations: {}", keys.join(", ")));
                        }
                    });
                });
            });
        }

        let hovering_files = ctx.input(|i| !i.raw.hovered_files.is_empty());
        egui::CentralPanel::default().frame(egui::Frame::none()).show(ctx, |ui| {
            if let Some(doc) = self.docs.get_mut(self.active) {
                actions.extend(canvas::show(ui, doc, &mut self.canvas, &mut self.textures, &mut self.pics, self.clip.is_some()));
            } else {
                ui.centered_and_justified(|ui| ui.heading("Drop a Newick or NEXUS tree here, or use File › Open"));
            }
            if hovering_files {
                let r = ui.max_rect();
                ui.painter().rect_filled(r, 0.0, Color32::from_rgba_unmultiplied(0, 114, 178, 40));
                ui.painter().text(r.center(), egui::Align2::CENTER_CENTER, "Drop trees, tables or images", egui::FontId::proportional(24.0), Color32::from_rgb(0, 80, 140));
            }
        });

        for a in actions {
            self.apply(a);
        }
        if load_data {
            if let Some(p) = rfd::FileDialog::new().add_filter("Table", &["csv", "tsv", "txt", "tab"]).pick_file() {
                self.open_path(&p);
            }
        }

        self.dialogs(ctx);
        if let Some(doc) = self.docs.get(self.active) {
            let images = self.pics.images();
            if let Some(job) = self.export.show(ctx, doc, &self.raster_fonts, &images) {
                self.status = match export::export(doc, &self.raster_fonts, &self.fonts.name, &images, &job) {
                    Ok(_) => format!("Exported {} ({}×{} px, {} dpi)", job.path.display(), job.width_px, job.height_px, job.dpi),
                    Err(e) => format!("Export failed: {:#}", e),
                };
            }
        }
        // Pending silhouettes: keep polling.
        if self.pics.entries.values().any(|s| matches!(s, PicStatus::Pending)) {
            ctx.request_repaint_after(std::time::Duration::from_millis(250));
        }
    }
}
