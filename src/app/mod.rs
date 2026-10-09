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
use crate::taxonomy;
use crate::tree::{NodeId, Tree};
use canvas::{Action, CanvasState, TextureCache};
use document::{Document, Posterior, SampleKind};
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
    /// The right-hand panel (tree summary, inspector, data).
    show_right_panel: bool,
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
        Settings { auto_phylopic: true, recent: Vec::new(), theme: Theme::System, show_right_panel: true }
    }
}

enum Dialog {
    Rename { node: NodeId, text: String },
    CladeLabel { node: NodeId, text: String },
    /// Branch-length transform being previewed; `base` holds the untransformed tree.
    Transform { node: NodeId, kind: ops::BranchTransform, value: f64, base: Tree },
    /// Collapse nodes whose support attribute `key` is below `threshold`.
    Collapse { key: String, threshold: f64 },
    /// Regex find-and-replace on labels.
    Regex(RegexEdit),
    SaveTree(SaveTreeOptions),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum TreeFormat {
    Newick,
    Nexus,
    PhyloXml,
}

/// What "Save tree as…" writes; remembered for the next save.
#[derive(Clone)]
struct SaveTreeOptions {
    format: TreeFormat,
    lengths: bool,
    internal_labels: bool,
    annotations: bool,
    /// For a tab holding a tree set: which trees to write.
    which: SaveWhich,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SaveWhich {
    /// Every tree in the set.
    All,
    /// The trees kept after burn-in.
    AfterBurnin,
    /// Only the tree on screen.
    Shown,
}

impl Default for SaveTreeOptions {
    fn default() -> Self {
        SaveTreeOptions { format: TreeFormat::Nexus, lengths: true, internal_labels: true, annotations: true, which: SaveWhich::All }
    }
}

#[derive(Default)]
struct RegexEdit {
    pattern: String,
    replace: String,
    tips: bool,
    internal: bool,
    selected_only: bool,
    ignore_case: bool,
}

/// Largest value of support attribute `key` on the tree's nodes.
fn max_support(t: &Tree, key: &str) -> f64 {
    t.preorder().iter().filter_map(|&n| t.nodes[n].attrs.get(key).and_then(|a| a.as_f64())).fold(0.0, f64::max)
}

/// 0.5 for probabilities, 50 for percentages.
fn default_support_threshold(t: &Tree, key: &str) -> f64 {
    if max_support(t, key) <= 1.0 {
        0.5
    } else {
        50.0
    }
}

/// A check of the active tree's tip names against the Open Tree taxonomy.
struct TaxonCheck {
    /// Tips checked (id and label at the time of the check).
    tips: Vec<(NodeId, String)>,
    pending: Option<std::sync::mpsc::Receiver<anyhow::Result<Vec<taxonomy::NameCheck>>>>,
    results: Vec<taxonomy::NameCheck>,
    /// Rename this tip to the suggested name.
    rename: Vec<bool>,
    only_problems: bool,
    error: Option<String>,
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
    /// Open Tree of Life name check in progress or showing results.
    taxon_check: Option<TaxonCheck>,
    /// Error shown in a window until dismissed: (title, message).
    error: Option<(String, String)>,
    show_node_report: bool,
    save_opts: SaveTreeOptions,
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
            taxon_check: None,
            error: None,
            show_node_report: false,
            save_opts: SaveTreeOptions::default(),
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
                        let mut post = Posterior::new(trees);
                        post.set_kind(SampleKind::guess(&name));
                        d.posterior = Some(post);
                        d.annotate_posterior();
                        let kind = d.posterior.as_ref().map(|p| p.kind.name()).unwrap_or_default();
                        self.add_doc(d);
                        self.status = format!("Opened {} with {} trees ({}; change it in the tree-set panel if wrong)", name, n, kind.to_lowercase());
                    }
                }
                Err(e) => self.fail(format!("Could not read {}", name), format!("{:#}", e)),
            },
            FileKind::Data => match DataTable::from_path(path) {
                Ok(table) => self.attach_data(table),
                Err(e) => self.fail(format!("Could not load data table {}", name), format!("{:#}", e)),
            },
            FileKind::Image => self.status = "Drop an image onto a tip to use it as that taxon's silhouette.".into(),
        }
    }

    /// Report an error in a window (and the status bar) until dismissed.
    fn fail(&mut self, title: String, message: String) {
        self.status = format!("{}: {}", title, message);
        self.error = Some((title, message));
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
                // With a tree set, the new tab gets the post-burn-in samples that
                // contain exactly this clade (the ones behind its posterior
                // support), so it keeps its own support, consensus and DensiTree.
                if let Some(p) = &doc.posterior {
                    let want = ops::tip_set(&doc.tree, n);
                    let kept = p.kept();
                    let trees: Vec<Tree> = kept.iter().filter_map(|t| ops::clade_subtree(t, &want, p.rooted)).collect();
                    let found = trees.len();
                    let pct = 100.0 * found as f64 / kept.len().max(1) as f64;
                    self.status = format!(
                        "Opened the clade ({} tips) from the {} of {} post-burn-in trees that contain it ({}%).",
                        want.len(),
                        found,
                        kept.len(),
                        crate::tree::format_num(pct, 1)
                    );
                    if found > 1 {
                        let mut q = Posterior::new(trees);
                        q.burnin = 0.0;
                        q.rooted = p.rooted;
                        q.kind = p.kind;
                        d.posterior = Some(q);
                        d.annotate_posterior();
                    } else {
                        self.status.push_str(" Too few to keep a tree set; the tab shows this tree only.");
                    }
                }
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
            Action::Polytomies(n, hard) => {
                let doc = &mut self.docs[self.active];
                let mut t = doc.tree.clone();
                let k = if hard { ops::hard_polytomies(&mut t, n, 0.0) } else { ops::soft_polytomies(&mut t, n) };
                let scope = if n == doc.tree.root { "the tree" } else { "this clade" };
                if k > 0 {
                    doc.checkpoint();
                    doc.tree = t;
                }
                self.status = match (hard, k) {
                    (_, 0) => format!("No {} polytomies to convert in {}.", if hard { "soft" } else { "hard" }, scope),
                    (true, k) => format!("Collapsed {} zero-length branches in {} into hard polytomies.", k, scope),
                    (false, k) => format!("Resolved {} hard polytomies in {} with zero-length branches.", k, scope),
                };
                return;
            }
            Action::CopyVcv(n) => {
                let (header, rows) = io::data::vcv_table(&self.docs[self.active].tree, n, None);
                self.ctx.copy_text(io::data::table_to_text(&header, &rows, b'\t'));
                self.status = format!("Copied the {0} × {0} variance–covariance matrix to the clipboard", rows.len());
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
                | Action::CopyVcv(_)
                | Action::Polytomies(..)
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
            .add_filter("Trees and data", &["nwk", "newick", "tre", "tree", "trees", "nex", "nexus", "nxs", "t", "con", "treefile", "contree", "xml", "phyloxml", "txt", "csv", "tsv"])
            .add_filter("All files", &["*"])
            .pick_files();
        for f in files.unwrap_or_default() {
            self.open_path(&f);
        }
    }

    /// Format and contents, then a file dialog.
    fn save_tree_dialog(&mut self, ctx: &egui::Context) {
        // Sizes of the active tab's tree set, if it has one: (all, after burn-in).
        let set = self.docs.get(self.active).and_then(|d| d.posterior.as_ref()).map(|p| (p.trees.len(), p.kept().len()));
        let Some(Dialog::SaveTree(o)) = &mut self.dialog else { return };
        let (mut save, mut cancel) = (false, false);
        egui::Window::new(if set.is_some() { "Save trees" } else { "Save tree" }).collapsible(false).resizable(false).show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Format");
                ui.selectable_value(&mut o.format, TreeFormat::Newick, "Newick");
                ui.selectable_value(&mut o.format, TreeFormat::Nexus, "NEXUS");
                ui.selectable_value(&mut o.format, TreeFormat::PhyloXml, "phyloXML");
            });
            if let Some((all, kept)) = set {
                ui.label("Trees");
                ui.radio_value(&mut o.which, SaveWhich::All, format!("All {} trees in the set", all));
                if kept < all {
                    ui.radio_value(&mut o.which, SaveWhich::AfterBurnin, format!("The {} trees after burn-in", kept));
                } else if o.which == SaveWhich::AfterBurnin {
                    o.which = SaveWhich::All;
                }
                ui.radio_value(&mut o.which, SaveWhich::Shown, "Only the tree shown");
                ui.separator();
            }
            ui.checkbox(&mut o.lengths, "Branch lengths");
            ui.checkbox(&mut o.internal_labels, "Internal node labels").on_hover_text("Clade names and numeric support labels on internal nodes; tip labels are always written");
            ui.checkbox(&mut o.annotations, "Annotations").on_hover_text("Posterior, bootstrap, HPD intervals, rates and other node data: [&key=value] comments in Newick/NEXUS, confidence and properties in phyloXML");
            if !o.internal_labels && !o.annotations {
                ui.weak("Only the topology, tip names and (if ticked) branch lengths will be written.");
            }
            ui.horizontal(|ui| {
                if ui.button("Save…").clicked() {
                    save = true;
                }
                if ui.button("Cancel").clicked() {
                    cancel = true;
                }
            });
        });
        if save || cancel {
            let o = o.clone();
            self.dialog = None;
            self.save_opts = o.clone();
            if save {
                self.save_tree(&o);
            }
        }
    }

    fn save_tree(&mut self, o: &SaveTreeOptions) {
        let Some(doc) = self.docs.get(self.active) else { return };
        let fmt = o.format;
        let (name, ext) = match fmt {
            TreeFormat::Newick => ("Newick", "nwk"),
            TreeFormat::Nexus => ("NEXUS", "nex"),
            TreeFormat::PhyloXml => ("phyloXML", "xml"),
        };
        // The whole tree set (or its post-burn-in part), or just the tree shown.
        let shown = doc.tree.compact();
        let trees: Vec<&Tree> = match (&doc.posterior, o.which) {
            (Some(p), SaveWhich::All) => p.trees.iter().collect(),
            (Some(p), SaveWhich::AfterBurnin) => p.kept().iter().collect(),
            _ => vec![&shown],
        };
        let stem = if trees.len() > 1 { "trees" } else { "tree" };
        let Some(path) = rfd::FileDialog::new().add_filter(name, &[ext]).set_file_name(format!("{}.{}", stem, ext)).save_file() else { return };
        let opts = newick::WriteOptions { lengths: o.lengths, internal_labels: o.internal_labels, annotations: o.annotations, ..Default::default() };
        let text = match fmt {
            TreeFormat::Newick => trees.iter().map(|t| newick::write_newick(t, &opts) + "\n").collect(),
            TreeFormat::Nexus => nexus::write_nexus(&trees, &opts),
            TreeFormat::PhyloXml => io::phyloxml::write_phyloxml(&trees, &opts),
        };
        self.status = match std::fs::write(&path, text) {
            Ok(_) if trees.len() > 1 => format!("Saved {} trees to {}", trees.len(), path.display()),
            Ok(_) => format!("Saved {}", path.display()),
            Err(e) => format!("Save failed: {}", e),
        };
    }

    fn posterior_action(&mut self, a: panels::PosteriorAction) {
        let auto = self.settings.auto_phylopic;
        let Some(doc) = self.docs.get_mut(self.active) else { return };
        let Some(p) = doc.posterior.as_mut() else { return };
        let started = std::time::Instant::now();
        let kept = p.kept();
        let n_kept = kept.len();
        let res: anyhow::Result<(String, Tree)> = (|| match a {
            panels::PosteriorAction::Consensus => {
                let mut s = consensus::summarize(kept, p.rooted)?;
                s.support_key = p.kind.support_key();
                let t = consensus::consensus(&s, p.threshold as f64, p.heights);
                Ok((format!("Consensus {:.0}%", p.threshold * 100.0), t))
            }
            panels::PosteriorAction::Mcc => {
                let mut s = consensus::summarize(kept, p.rooted)?;
                s.support_key = p.kind.support_key();
                let t = consensus::mcc_tree(kept, &s, p.heights)?;
                Ok(("MCC".to_string(), t))
            }
            panels::PosteriorAction::OpenSample(i) => Ok((format!("sample {}", i + 1), p.trees[i].clone())),
        })();
        match res {
            Ok((label, t)) => {
                p.info = format!("{} from {} trees in {:.2}s", label, n_kept, started.elapsed().as_secs_f32());
                let name = format!("{} – {}", doc.name, label);
                let mut d = Document::new(&name, t, auto);
                d.view.layout = doc.view.layout.clone();
                // A summary or single sample is one tree. Consensus/MCC tabs keep a
                // reference to the post-burn-in sample (shared, not copied) for DensiTree only.
                if !matches!(a, panels::PosteriorAction::OpenSample(_)) {
                    d.source_sample = Some((p.trees.clone(), p.trees.len() - n_kept));
                }
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
            self.dialog = Some(Dialog::SaveTree(self.save_opts.clone()));
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

    /// Send the active tree's tip names to Open Tree of Life in the background.
    fn start_taxon_check(&mut self) {
        let Some(doc) = self.docs.get(self.active) else { return };
        let tips: Vec<(NodeId, String)> = doc.tree.tips().into_iter().map(|t| (t, doc.tree.label(t).to_string())).collect();
        let names: Vec<String> = tips.iter().map(|(_, l)| taxonomy::query_name(l)).collect();
        let (tx, rx) = std::sync::mpsc::channel();
        let ctx = self.ctx.clone();
        std::thread::spawn(move || {
            let _ = tx.send(taxonomy::check_names(&names));
            ctx.request_repaint();
        });
        self.taxon_check = Some(TaxonCheck { tips, pending: Some(rx), results: Vec::new(), rename: Vec::new(), only_problems: true, error: None });
    }

    /// Progress, then a table of tip names with Open Tree's verdict and
    /// optional renaming to the accepted names.
    fn taxon_dialog(&mut self, ctx: &egui::Context) {
        let Some(tc) = &mut self.taxon_check else { return };
        if let Some(rx) = &tc.pending {
            match rx.try_recv() {
                Ok(Ok(r)) => {
                    tc.rename = vec![false; r.len()];
                    tc.results = r;
                    tc.pending = None;
                }
                Ok(Err(e)) => {
                    tc.error = Some(format!("{:#}", e));
                    tc.pending = None;
                }
                Err(std::sync::mpsc::TryRecvError::Empty) => {}
                Err(std::sync::mpsc::TryRecvError::Disconnected) => tc.pending = None,
            }
        }
        let (mut open, mut apply, mut copy) = (true, false, false);
        egui::Window::new("Check taxon names").open(&mut open).default_width(620.0).show(ctx, |ui| {
            ui.small("Tip names checked against the Open Tree of Life taxonomy (OTT), with fuzzy matching for misspellings.");
            if tc.pending.is_some() {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(format!("Checking {} names with Open Tree of Life…", tc.tips.len()));
                });
                return;
            }
            if let Some(e) = &tc.error {
                ui.colored_label(Color32::from_rgb(198, 40, 40), e);
                return;
            }
            let count = |s: taxonomy::NameStatus| tc.results.iter().filter(|r| r.status == s).count();
            ui.label(format!(
                "{} ok · {} synonyms · {} possible misspellings · {} not found",
                count(taxonomy::NameStatus::Exact),
                count(taxonomy::NameStatus::Synonym),
                count(taxonomy::NameStatus::Approximate),
                count(taxonomy::NameStatus::Unmatched)
            ));
            ui.horizontal(|ui| {
                ui.checkbox(&mut tc.only_problems, "Show only names that need attention");
                if ui.button("Select all suggestions").clicked() {
                    for (i, r) in tc.results.iter().enumerate() {
                        tc.rename[i] = r.suggestion().is_some();
                    }
                }
            });
            egui::ScrollArea::vertical().max_height(380.0).show(ui, |ui| {
                egui::Grid::new("taxon-grid").striped(true).num_columns(4).show(ui, |ui| {
                    ui.strong("Tip");
                    ui.strong("Status");
                    ui.strong("Open Tree name");
                    ui.strong("Rename");
                    ui.end_row();
                    for (i, r) in tc.results.iter().enumerate() {
                        if tc.only_problems && r.status == taxonomy::NameStatus::Exact {
                            continue;
                        }
                        ui.label(egui::RichText::new(&r.query).italics());
                        let color = match r.status {
                            taxonomy::NameStatus::Exact => Color32::from_rgb(46, 125, 50),
                            taxonomy::NameStatus::Synonym => Color32::from_rgb(2, 119, 189),
                            taxonomy::NameStatus::Approximate => Color32::from_rgb(230, 120, 0),
                            taxonomy::NameStatus::Unmatched => Color32::from_rgb(198, 40, 40),
                        };
                        let mut hover = format!("Score {:.2}", r.score);
                        if let Some(id) = r.ott_id {
                            hover.push_str(&format!(" · OTT {}", id));
                        }
                        if let Some(rank) = &r.rank {
                            hover.push_str(&format!(" · {}", rank));
                        }
                        if !r.flags.is_empty() {
                            hover.push_str(&format!(" · {}", r.flags.join(", ")));
                        }
                        if !r.alternatives.is_empty() {
                            hover.push_str(&format!("\nAlso matches: {}", r.alternatives.join(", ")));
                        }
                        ui.colored_label(color, r.status.label()).on_hover_text(hover);
                        ui.label(egui::RichText::new(r.accepted.as_deref().unwrap_or("—")).italics());
                        if r.suggestion().is_some() {
                            ui.checkbox(&mut tc.rename[i], "");
                        } else {
                            ui.label("");
                        }
                        ui.end_row();
                    }
                });
            });
            ui.horizontal(|ui| {
                let n = tc.rename.iter().filter(|b| **b).count();
                if ui.add_enabled(n > 0, egui::Button::new(format!("Rename {} tips", n))).clicked() {
                    apply = true;
                }
                if ui.button("Copy report").on_hover_text("Tab-separated, pastes into Excel").clicked() {
                    copy = true;
                }
            });
        });
        if copy {
            let mut text = String::from("label\tstatus\taccepted_name\tott_id\tflags\n");
            for ((_, l), r) in tc.tips.iter().zip(&tc.results) {
                text.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    l,
                    r.status.label(),
                    r.accepted.as_deref().unwrap_or(""),
                    r.ott_id.map(|i| i.to_string()).unwrap_or_default(),
                    r.flags.join(",")
                ));
            }
            ctx.copy_text(text);
            self.status = "Copied the name check report".into();
        }
        if apply {
            let tc = self.taxon_check.as_mut().unwrap();
            let Some(doc) = self.docs.get_mut(self.active) else { return };
            doc.checkpoint();
            let (mut done, mut skipped) = (0, 0);
            for (i, ((node, label), r)) in tc.tips.iter().zip(&tc.results).enumerate() {
                let Some(new) = r.suggestion().filter(|_| tc.rename[i]) else { continue };
                // Only if this tab still has that tip under the same name.
                if *node >= doc.tree.nodes.len() || doc.tree.nodes[*node].label.as_deref() != Some(label.as_str()) {
                    skipped += 1;
                    continue;
                }
                // Keep the label's style: underscores stay underscores.
                let new = if label.contains('_') && !label.contains(' ') { new.replace(' ', "_") } else { new.to_string() };
                doc.tree.nodes[*node].label = Some(new);
                done += 1;
            }
            self.status = format!("Renamed {} tips to their Open Tree names{}", done, if skipped > 0 { format!("; {} skipped (tree changed)", skipped) } else { String::new() });
            open = false;
        }
        if !open {
            self.taxon_check = None;
        }
    }

    /// Every clade (bipartition) of the active tree with its support values;
    /// clicking a node number selects it.
    fn node_report(&mut self, ctx: &egui::Context) {
        if !self.show_node_report {
            return;
        }
        let Some(doc) = self.docs.get_mut(self.active) else { return };
        let (header, rows, ids) = io::data::node_report(&doc.tree);
        let mut open = true;
        egui::Window::new("Node report").open(&mut open).default_width(640.0).show(ctx, |ui| {
            ui.small(
                "One row per clade (bipartition): its support values, branch length and taxa. For unrooted trees (e.g. RAxML bipartitions) each row is the split between these taxa and the rest. Click a node number to select it.",
            );
            ui.horizontal(|ui| {
                if ui.button("Copy table").on_hover_text("Tab-separated, pastes into Excel").clicked() {
                    ui.ctx().copy_text(io::data::table_to_text(&header, &rows, b'\t'));
                }
                if ui.button("Save CSV…").clicked() {
                    if let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).set_file_name("node_report.csv").save_file() {
                        let _ = std::fs::write(path, io::data::table_to_text(&header, &rows, b','));
                    }
                }
                ui.label(format!("{} nodes", rows.len()));
            });
            egui::ScrollArea::both().max_height(420.0).show(ui, |ui| {
                egui::Grid::new("node-report").striped(true).show(ui, |ui| {
                    for h in &header {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for (r, &n) in rows.iter().zip(&ids) {
                        for (j, v) in r.iter().enumerate() {
                            if j == 0 {
                                if ui.link(v).clicked() {
                                    doc.selection.clear();
                                    doc.selection.insert(n);
                                }
                            } else if j == r.len() - 1 {
                                // Taxa: shortened, full list on hover.
                                let short: String = if v.chars().count() > 60 { format!("{}…", v.chars().take(60).collect::<String>()) } else { v.clone() };
                                ui.label(egui::RichText::new(short.replace('_', " ")).italics()).on_hover_text(v.replace('_', " "));
                            } else {
                                ui.label(v);
                            }
                        }
                        ui.end_row();
                    }
                });
            });
        });
        self.show_node_report = open;
    }

    /// Find and replace in labels with a regular expression, previewing every change.
    fn regex_dialog(&mut self, ctx: &egui::Context) {
        let Some(Dialog::Regex(r)) = &mut self.dialog else { return };
        let Some(doc) = self.docs.get(self.active) else { return };
        let (mut apply, mut cancel) = (false, false);
        let re = regex::RegexBuilder::new(&r.pattern).case_insensitive(r.ignore_case).build();
        let nodes: Vec<NodeId> = doc
            .tree
            .preorder()
            .into_iter()
            .filter(|&n| if doc.tree.is_tip(n) { r.tips } else { r.internal })
            .filter(|n| !r.selected_only || doc.selection.contains(n))
            .collect();
        let changes = match (&re, r.pattern.is_empty()) {
            (Ok(re), false) => ops::regex_replacements(&doc.tree, re, &r.replace, &nodes),
            _ => Vec::new(),
        };
        egui::Window::new("Find and replace in labels").collapsible(false).default_width(520.0).show(ctx, |ui| {
            egui::Grid::new("regex-fields").num_columns(2).show(ui, |ui| {
                ui.label("Find (regex)");
                ui.add(egui::TextEdit::singleline(&mut r.pattern).font(egui::TextStyle::Monospace).hint_text(r"e.g. ^(\w+_\w+)_.*$").desired_width(320.0));
                ui.end_row();
                ui.label("Replace with");
                ui.add(egui::TextEdit::singleline(&mut r.replace).font(egui::TextStyle::Monospace).hint_text("e.g. $1").desired_width(320.0));
                ui.end_row();
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut r.tips, "Tip labels");
                ui.checkbox(&mut r.internal, "Internal node labels");
                ui.add_enabled(!doc.selection.is_empty(), egui::Checkbox::new(&mut r.selected_only, "Selected nodes only"));
                ui.checkbox(&mut r.ignore_case, "Ignore case");
            });
            ui.small("Rust regex syntax. Use $1, $2 or ${name} in the replacement for captured groups; every match in a label is replaced.");
            if let Err(e) = &re {
                ui.colored_label(Color32::from_rgb(198, 40, 40), format!("Invalid pattern: {}", e.to_string().lines().last().unwrap_or("")));
            } else if !r.pattern.is_empty() {
                ui.label(format!("{} of {} labels will change", changes.len(), nodes.len()));
                egui::ScrollArea::vertical().max_height(260.0).show(ui, |ui| {
                    egui::Grid::new("regex-preview").striped(true).show(ui, |ui| {
                        for (_, old, new) in changes.iter().take(500) {
                            ui.label(old);
                            ui.label("→");
                            ui.label(egui::RichText::new(new).strong());
                            ui.end_row();
                        }
                    });
                });
            }
            ui.horizontal(|ui| {
                if ui.add_enabled(!changes.is_empty(), egui::Button::new(format!("Rename {} labels", changes.len()))).clicked() {
                    apply = true;
                }
                if ui.button("Cancel").clicked() {
                    cancel = true;
                }
            });
        });
        if apply {
            let doc = &mut self.docs[self.active];
            doc.checkpoint();
            for (n, _, new) in &changes {
                doc.set_label(*n, if new.is_empty() { None } else { Some(new.clone()) });
            }
            self.status = format!("Renamed {} labels", changes.len());
        }
        if apply || cancel {
            self.dialog = None;
        }
    }

    /// Choose the support attribute and threshold; shows how many nodes would go.
    fn collapse_dialog(&mut self, ctx: &egui::Context) {
        let Some(Dialog::Collapse { key, threshold }) = &mut self.dialog else { return };
        let Some(doc) = self.docs.get(self.active) else { return };
        let keys = ops::support_keys(&doc.tree);
        let (mut apply, mut cancel) = (false, false);
        egui::Window::new("Collapse weakly supported nodes").collapsible(false).resizable(false).show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Support");
                egui::ComboBox::from_id_salt("collapse-key").selected_text(key.as_str()).show_ui(ui, |ui| {
                    for k in &keys {
                        if ui.selectable_label(key == k, *k).clicked() && key != k {
                            *key = k.to_string();
                            *threshold = default_support_threshold(&doc.tree, k);
                        }
                    }
                });
            });
            let max = if max_support(&doc.tree, key) <= 1.0 { 1.0 } else { 100.0 };
            ui.add(egui::Slider::new(threshold, 0.0..=max).text("collapse below"));
            let mut t = doc.tree.clone();
            let k = ops::collapse_weak(&mut t, 0.0, Some((key.as_str(), *threshold)));
            ui.label(format!("{} node{} would be collapsed into polytomies.", k, if k == 1 { "" } else { "s" }));
            ui.horizontal(|ui| {
                if ui.add_enabled(k > 0, egui::Button::new("Collapse")).clicked() {
                    apply = true;
                }
                if ui.button("Cancel").clicked() {
                    cancel = true;
                }
            });
        });
        if apply {
            let (key, th) = (key.clone(), *threshold);
            let doc = &mut self.docs[self.active];
            doc.checkpoint();
            let k = ops::collapse_weak(&mut doc.tree, 0.0, Some((key.as_str(), th)));
            self.status = format!("Collapsed {} nodes with {} < {}", k, key, crate::tree::format_num(th, 3));
        }
        if apply || cancel {
            self.dialog = None;
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
            ui.label(egui::RichText::new("Reference").strong().small());
            ui.small(kind.reference());
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
                if ui.button("Save tree as…  (Ctrl+S)").on_hover_text("Newick, NEXUS or phyloXML, choosing what to include").clicked() {
                    ui.close_menu();
                    self.dialog = Some(Dialog::SaveTree(self.save_opts.clone()));
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
                if ui.button("Find and replace in labels (regex)…").clicked() {
                    self.dialog = Some(Dialog::Regex(RegexEdit { tips: true, ..Default::default() }));
                    ui.close_menu();
                }
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
                if ui.button("Node report (support per node)…").on_hover_text("Table of every clade / bipartition with its support values, as from RAxML, IQ-TREE, MrBayes or BEAST").clicked() {
                    self.show_node_report = true;
                    ui.close_menu();
                }
                if ui.button("Collapse weakly supported nodes…").clicked() {
                    if let Some(d) = self.docs.get(self.active) {
                        match ops::support_keys(&d.tree).first() {
                            Some(&key) => {
                                let threshold = default_support_threshold(&d.tree, key);
                                self.dialog = Some(Dialog::Collapse { key: key.to_string(), threshold });
                            }
                            None => self.status = "This tree has no support values (posterior, prob, bootstrap or numeric node labels).".into(),
                        }
                    }
                    ui.close_menu();
                }
                if let Some(r) = self.docs.get(self.active).map(|d| d.tree.root) {
                    let mut acts = Vec::new();
                    ui.menu_button("Polytomies", |ui| canvas::polytomy_menu(ui, r, &mut acts));
                    for a in acts {
                        self.apply(a);
                    }
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
                if ui.checkbox(&mut self.settings.show_right_panel, "Show right panel").on_hover_text("Tree summary, node inspector and comparative data").changed() {
                    ui.close_menu();
                }
                ui.separator();
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
            ui.menu_button("Taxa", |ui| {
                if ui.button("Check taxon names (Open Tree of Life)…").on_hover_text("Find misspellings, synonyms and unknown names; optionally rename tips").clicked() {
                    self.start_taxon_check();
                    ui.close_menu();
                }
                ui.separator();
                ui.label(egui::RichText::new("PhyloPic silhouettes").small().weak());
                if ui.checkbox(&mut self.settings.auto_phylopic, "Fetch PhyloPic silhouettes automatically").changed() {
                    self.pics.enabled = self.settings.auto_phylopic;
                }
                if ui.button("Show PhyloPic silhouettes for this tree").clicked() {
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
                if ui.button("Retry failed PhyloPic lookups").clicked() {
                    self.pics.retry_failed();
                    ui.close_menu();
                }
                if ui.button("PhyloPic image credits…").clicked() {
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
                let label = d.name.clone();
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
        if let Some((title, message)) = &self.error {
            let mut close = false;
            egui::Window::new(egui::RichText::new(format!("⚠ {}", title)).color(Color32::from_rgb(198, 40, 40)))
                .collapsible(false)
                .resizable(false)
                .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
                .max_width(520.0)
                .show(ctx, |ui| {
                    ui.label(message);
                    ui.add_space(6.0);
                    ui.horizontal(|ui| {
                        if ui.button("OK").clicked() || ui.input(|i| i.key_pressed(Key::Enter) || i.key_pressed(Key::Escape)) {
                            close = true;
                        }
                        if ui.button("Copy message").clicked() {
                            ui.ctx().copy_text(message.clone());
                        }
                    });
                });
            if close {
                self.error = None;
            }
        }
        self.taxon_dialog(ctx);
        self.node_report(ctx);
        let mut done = false;
        if matches!(self.dialog, Some(Dialog::Transform { .. })) {
            self.transform_dialog(ctx);
        } else if matches!(self.dialog, Some(Dialog::Collapse { .. })) {
            self.collapse_dialog(ctx);
        } else if matches!(self.dialog, Some(Dialog::SaveTree(_))) {
            self.save_tree_dialog(ctx);
        } else if matches!(self.dialog, Some(Dialog::Regex(_))) {
            self.regex_dialog(ctx);
        } else if let Some(dialog) = &mut self.dialog {
            let (title, node, text) = match dialog {
                Dialog::Rename { node, text } => ("Rename", *node, text),
                Dialog::CladeLabel { node, text } => ("Clade label", *node, text),
                Dialog::Transform { .. } | Dialog::Collapse { .. } | Dialog::Regex(_) | Dialog::SaveTree(_) => unreachable!(),
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
                        doc.set_label(node, if value.is_empty() { None } else { Some(value) });
                        doc.show_node_label(node);
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
            let mut hide_right = false;
            if self.settings.show_right_panel {
            // Scrolling both ways keeps wide tables from pinning the panel's width,
            // so it can be dragged narrower as well as wider.
            egui::SidePanel::right("right").resizable(true).default_width(300.0).width_range(200.0..=1000.0).show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.weak("Tree & node details");
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.small_button("✕").on_hover_text("Hide this panel (View › Show right panel brings it back)").clicked() {
                            hide_right = true;
                        }
                    });
                });
                egui::ScrollArea::both().auto_shrink([false, false]).show(ui, |ui| {
                    ui.collapsing("Tree summary", |ui| {
                        let t = &doc.tree;
                        ui.label(format!("{} tips, {} internal nodes", t.num_tips(), t.preorder().len() - t.num_tips()));
                        ui.label(format!("Branch lengths: {}", if t.has_lengths() { "yes" } else { "no" }));
                        ui.label(format!("Rooted: {}", match t.rooted { Some(true) => "yes", Some(false) => "no ([&U])", None => "unspecified" }));
                        let h = t.heights();
                        ui.label(format!("Root height: {}", crate::tree::format_num(h[t.root], 5)));
                        let (c, norm, _) = t.colless(t.root);
                        ui.label(format!("Colless index: {}{}", c, norm.map(|x| format!(" (normalized {})", crate::tree::format_num(x, 3))).unwrap_or_default()))
                            .on_hover_text(panels::COLLESS_HOVER);
                        match t.gamma(t.root) {
                            Ok((g, p)) => ui.label(format!("γ: {} (p = {})", crate::tree::format_num(g, 3), crate::tree::format_num(p, 3))),
                            Err(why) => ui.weak(format!("γ: n/a ({})", why.split('(').next().unwrap_or(&why).trim())),
                        }
                        .on_hover_text(panels::GAMMA_HOVER);
                        let keys = t.attr_keys();
                        if !keys.is_empty() {
                            ui.label(format!("Annotations: {}", keys.join(", ")));
                        }
                    });
                    ui.separator();
                    actions.extend(panels::inspector(ui, doc, &mut self.search));
                    ui.separator();
                    load_data = panels::data_panel(ui, doc, &mut self.status);
                });
            });
            }
            if hide_right {
                self.settings.show_right_panel = false;
            }
        }

        let hovering_files = ctx.input(|i| !i.raw.hovered_files.is_empty());
        egui::CentralPanel::default().frame(egui::Frame::none()).show(ctx, |ui| {
            if let Some(doc) = self.docs.get_mut(self.active) {
                let full = ui.max_rect();
                actions.extend(canvas::show(ui, doc, &mut self.canvas, &mut self.textures, &mut self.pics, self.clip.is_some()));
                if !self.settings.show_right_panel {
                    let r = egui::Rect::from_min_size(egui::pos2(full.right() - 92.0, full.top() + 6.0), egui::vec2(86.0, 22.0));
                    if ui.put(r, egui::Button::new("◀ Details")).on_hover_text("Show the right panel").clicked() {
                        self.settings.show_right_panel = true;
                    }
                }
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
                self.status = match export::export(doc, &self.raster_fonts, &self.fonts, &images, &job) {
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
