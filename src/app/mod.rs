//! The Canopy desktop application (egui/eframe).

pub mod canvas;
pub mod document;
pub mod export;
pub mod panels;
pub mod references;
pub mod table;

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

/// Bundled example (a set of 10 posterior trees of 9 primates), available from
/// Help › Open example trees. Compiled in, so it works without the examples folder.
const EXAMPLE: &str = include_str!("../../examples/Apes_10.nwk");
/// Source of the example data, cited under Help › Open example trees.
const EXAMPLE_CITATION: &str = "Kay, R. F., Organ, C. L., Morse, P. E., Allen, K. L. & Kirk, E. C. (2026) Fossil evidence favors a role for vision in the modular evolution of the primate neocortex. Science 393: 622–627. doi:10.1126/science.aee3541";
/// Reference for the rogue-taxon analysis (Tree › Rogue taxa).
const ROGUE_CITATION: &str = "Klawitter, J., Bouckaert, R. R. & Drummond, A. J. (2024) Skeletons in the forest: using entropy-based rogue detection on Bayesian phylogenetic tree distributions. bioRxiv, doi:10.1101/2024.09.25.615070";

#[derive(serde::Serialize, serde::Deserialize)]
#[serde(default)]
struct Settings {
    /// Add a PhyloPic layer (and so fetch silhouettes) when a tree is opened. Off by
    /// default; stored under a new key so the earlier "on" default doesn't carry over.
    #[serde(rename = "auto_phylopic_v2")]
    auto_phylopic: bool,
    recent: Vec<PathBuf>,
    /// Stored under a new key so themes saved before "Follow system" became the
    /// default don't override it (earlier keys: "theme", "ui_theme").
    #[serde(rename = "ui_theme_v2")]
    theme: Theme,
    /// The right-hand panel (tree summary, inspector, data).
    show_right_panel: bool,
    /// Check GitHub for a newer release at startup.
    check_updates: bool,
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
        Settings { auto_phylopic: false, recent: Vec::new(), theme: Theme::System, show_right_panel: true, check_updates: true }
    }
}

/// Keep the tips of `t` whose labels pass `keep`; false (tree unchanged in
/// meaning) when nothing would go or fewer than two tips would stay.
fn keep_named(t: &mut Tree, keep: impl Fn(&str) -> bool) -> bool {
    let tips = t.tips();
    let kept: Vec<NodeId> = tips.iter().copied().filter(|&x| keep(t.label(x))).collect();
    kept.len() >= 2 && kept.len() < tips.len() && ops::keep_tips(t, &kept).is_ok()
}

/// Status for an edit applied to a whole tree set.
fn set_note(what: &str, done: usize, skipped: usize) -> String {
    let mut s = format!("{} in the shown tree and {} of the set", what, count(done, "tree", "trees"));
    if skipped > 0 {
        s += &format!(" ({} left unchanged)", count(skipped, "tree", "trees"));
    }
    s
}

/// "1 tree", "3 trees": a count with its noun in the singular or plural.
pub(crate) fn count(n: usize, one: &str, many: &str) -> String {
    format!("{} {}", n, if n == 1 { one } else { many })
}

/// How many of `pool` taxa to remove for `percent`: a clade keeps at least one
/// taxon, and the tree at least three of its `total`.
fn thin_count(pool: usize, total: usize, percent: f32, clade: bool) -> usize {
    let k = (pool as f32 * percent / 100.0).round() as usize;
    let k = if clade { k.min(pool.saturating_sub(1)) } else { k };
    k.min(total.saturating_sub(3))
}

/// Add a fossil-range layer to `doc` if it has none.
fn ensure_fossil_layer(doc: &mut Document) {
    if !doc.layers.iter().any(|e| matches!(e.layer, Layer::FossilRanges(_))) {
        doc.add_layer(Layer::default_fossil_ranges());
    }
}

/// Age (Ma) of node `n` in the tree's time frame (`style::time_frame`).
fn fossil_age_of(doc: &Document, n: NodeId) -> f64 {
    let (young, mpu) = time_frame(doc.layers.iter().filter(|e| e.enabled).map(|e| &e.layer), &doc.tree);
    let d = doc.tree.depths(true);
    let max = doc.tree.tips().into_iter().map(|t| d[t]).fold(0.0, f64::max);
    young + (max - d[n]) * mpu
}

/// A consensus or MCC tree being built in the background.
struct SummaryJob {
    label: String,
    doc_name: String,
    layout: crate::layout::LayoutOptions,
    /// The tree set and its burn-in, kept by the new tab for DensiTree.
    sample: (Arc<Vec<Tree>>, usize),
    started: std::time::Instant,
    rx: std::sync::mpsc::Receiver<Result<Tree, String>>,
}

/// State of the rogue-taxon window: the post-burn-in trees, settings and results
/// (entropy, probability of the most probable tree, rogue scores, skeleton rounds).
struct RogueAnalysis {
    name: String,
    /// The tree set (shared) and how many burn-in trees to skip.
    trees: Arc<Vec<Tree>>,
    skip: usize,
    max_size: usize,
    min_prob: f64,
    result: Option<Result<(f64, f64, Vec<crate::rogue::RogueScore>, Vec<crate::rogue::Round>), String>>,
    /// A run in progress (in the background): when it started and its result.
    running: Option<(std::time::Instant, std::sync::mpsc::Receiver<Result<(f64, f64, Vec<crate::rogue::RogueScore>, Vec<crate::rogue::Round>), String>>)>,
}

enum Dialog {
    Rename { node: NodeId, text: String },
    /// Name for a new tip added at `node`.
    AddTaxon { node: NodeId, mode: ops::PasteMode, text: String },
    /// Fossil range of a tip being edited (oldest and youngest ages, Ma).
    FossilRange { node: NodeId, oldest: f64, youngest: f64 },
    /// Keep a random percentage of a tree set's trees (new tab).
    ThinTrees { percent: f32 },
    /// Remove a random percentage of taxa from the tree(s) and data (new tab).
    /// `clade`: only that clade's taxa (at least one of them stays).
    ThinTaxa { percent: f32, clade: Option<NodeId> },
    /// Drop tips whose terminal branch is shorter than `min`.
    DropShort { min: f64 },
    /// Taxon whose PhyloPic silhouette a node shows.
    Silhouette { node: NodeId, text: String },
    CladeLabel { node: NodeId, text: String },
    /// Branch-length transform being previewed; `base` holds the untransformed tree.
    Transform { node: NodeId, kind: ops::BranchTransform, value: f64, base: Tree },
    /// Collapse nodes whose support attribute `key` is below `threshold`.
    Collapse { key: String, threshold: f64 },
    /// Regex find-and-replace on labels.
    Regex(RegexEdit),
    SaveTree(SaveTreeOptions),
    /// Drop species lacking data in the ticked trait columns.
    DropMissing { cols: Vec<(String, bool)> },
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
    /// NEXUS: add a CANOPY block with the display settings.
    canopy: bool,
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
        SaveTreeOptions { format: TreeFormat::Nexus, lengths: true, internal_labels: true, annotations: true, which: SaveWhich::All, canopy: true }
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
    /// Node-density (delta) test result for the active tree, while its window is open.
    delta_report: Option<Result<crate::nodedensity::DeltaTest, String>>,
    /// BayesTraits variable-rates summary for the active tree: (file name, branches, warnings).
    varrates: Option<(String, Vec<crate::varrates::BranchRate>, Vec<String>)>,
    /// MCMC diagnostics of the imported variable-rates runs (with their logs, if found).
    varrates_diag: Option<crate::mcmc::Diagnostics>,
    /// BayesTraits Fabric summary for the active tree: (file name, branches, warnings).
    fabric: Option<(String, Vec<crate::fabric::FabricBranch>, Vec<String>)>,
    fabric_diag: Option<crate::mcmc::Diagnostics>,
    /// The MCMC diagnostics window.
    mcmc_win: Option<McmcWindow>,
    /// The tabs the variable-rates and Fabric results were imported onto.
    varrates_doc: usize,
    fabric_doc: usize,
    /// The variable-rates and Fabric results windows are shown (closing hides them).
    show_varrates: bool,
    show_fabric: bool,
    /// The runs behind the variable-rates and Fabric results.
    varrates_runs: RunFiles<crate::varrates::VarRates>,
    fabric_runs: RunFiles<crate::fabric::Fabric>,
    /// Palettes chosen in the variable-rates and Fabric windows.
    varrates_palette: Palette,
    fabric_palette: Palette,
    /// Plot the median instead of the mean (branch scalars / LS_Beta).
    varrates_median: bool,
    fabric_median: bool,
    /// Rogue-taxon analysis of a tree set (window).
    rogue: Option<RogueAnalysis>,
    /// Consensus / MCC tree being built.
    job: Option<SummaryJob>,
    /// Lineages-through-time window.
    ltt: Option<LttWindow>,
    tangle: Option<TangleWindow>,
    traversal: Option<TraversalWindow>,
    /// Comparative-data editor window open.
    data_editor: bool,
    /// Summary-statistics window when open: (only taxa in the tree, column plotted).
    data_summary: Option<SummaryView>,
    save_opts: SaveTreeOptions,
    /// The startup release check, until it answers.
    update_rx: Option<std::sync::mpsc::Receiver<Option<crate::update::Release>>>,
    /// A newer release on GitHub, if the check found one.
    new_release: Option<crate::update::Release>,
    show_about: bool,
    show_references: bool,
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
        let pics = PhyloPic::new(Arc::new(move || ctx.request_repaint()));
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
            show_references: false,
            ctx: cc.egui_ctx.clone(),
            clip: None,
            clip_text: String::new(),
            logo: None,
            taxon_check: None,
            error: None,
            show_node_report: false,
            delta_report: None,
            varrates: None,
            varrates_diag: None,
            fabric_diag: None,
            mcmc_win: None,
            varrates_doc: 0,
            fabric_doc: 0,
            show_varrates: false,
            show_fabric: false,
            varrates_runs: RunFiles::default(),
            fabric_runs: RunFiles::default(),
            fabric: None,
            varrates_palette: Palette::Viridis,
            fabric_palette: Palette::BlueRed,
            varrates_median: false,
            fabric_median: false,
            rogue: None,
            job: None,
            ltt: None,
            tangle: None,
            traversal: None,
            data_editor: false,
            data_summary: None,
            save_opts: SaveTreeOptions::default(),
            update_rx: None,
            new_release: None,
        };
        if app.settings.check_updates {
            app.start_update_check();
        }
        for f in files {
            app.open_path(&f);
        }
        app
    }

    /// Ask GitHub for the latest release on a background thread (quietly does
    /// nothing when offline).
    fn start_update_check(&mut self) {
        let (tx, rx) = std::sync::mpsc::channel();
        let ctx = self.ctx.clone();
        std::thread::spawn(move || {
            let _ = tx.send(crate::update::check().ok().flatten());
            ctx.request_repaint();
        });
        self.update_rx = Some(rx);
    }

    fn poll_update_check(&mut self) {
        let Some(rx) = &self.update_rx else { return };
        match rx.try_recv() {
            Ok(r) => {
                if let Some(r) = &r {
                    self.status = format!("Canopy {} is available (Help › Download Canopy {}).", r.version, r.version);
                }
                self.new_release = r;
                self.update_rx = None;
            }
            Err(std::sync::mpsc::TryRecvError::Disconnected) => self.update_rx = None,
            Err(std::sync::mpsc::TryRecvError::Empty) => {}
        }
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

    /// Open trees in a new tab: one tree as a plain document, several as a tree set
    /// (posterior sample or bootstrap replicates, guessed from the name).
    fn open_trees(&mut self, name: &str, mut trees: Vec<Tree>, path: Option<&Path>) {
        let n = trees.len();
        if n == 1 {
            let mut d = Document::new(name, trees.remove(0), self.settings.auto_phylopic);
            d.path = path.map(Path::to_path_buf);
            self.add_doc(d);
            self.status = format!("Opened {}", name);
        } else {
            let first = trees[0].clone();
            let mut d = Document::new(name, first, self.settings.auto_phylopic);
            d.path = path.map(Path::to_path_buf);
            let mut post = Posterior::new(trees);
            post.set_kind(SampleKind::guess(name));
            d.posterior = Some(post);
            d.annotate_posterior();
            let kind = d.posterior.as_ref().map(|p| p.kind.name()).unwrap_or_default();
            self.add_doc(d);
            self.status = format!("Opened {} with {} trees ({}; change it in the tree-set panel if wrong)", name, n, kind.to_lowercase());
        }
    }

    pub fn open_path(&mut self, path: &Path) {
        let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        match io::classify(path) {
            FileKind::Tree | FileKind::Unknown => match io::read_trees(path) {
                Ok(trees) => {
                    self.remember(path);
                    self.open_trees(&name, trees, Some(path));
                    // Display settings saved by Canopy in a NEXUS CANOPY block.
                    let text = std::fs::read_to_string(path).unwrap_or_default();
                    // HyPhy results: say which branch values came with the tree.
                    if io::hyphy::is_hyphy_json(&text) {
                        if let Some(d) = self.docs.get(self.active) {
                            let skip = ["tested", "hyphy_node"];
                            let keys: Vec<String> = d.tree.attr_keys().into_iter().filter(|k| !skip.contains(&k.as_str())).collect();
                            let method = d.tree.name.clone().unwrap_or_else(|| "HyPhy".into());
                            self.status = format!("Opened {} results with branch values: {}. Color or label branches by them in the layer settings.", method, keys.join(", "));
                        }
                    }
                    if io::nexus::is_nexus(&text) {
                        if let Some(d) = self.docs.get_mut(self.active) {
                            match d.apply_canopy_block(&text) {
                                Ok(true) => self.status = format!("Opened {} with its Canopy settings", name),
                                Ok(false) => {}
                                Err(e) => self.status = format!("Opened {}; its Canopy settings could not be read: {:#}", name, e),
                            }
                        }
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
            Action::ThinClade(n) => {
                self.dialog = Some(Dialog::ThinTaxa { percent: 10.0, clade: Some(n) });
                return;
            }
            Action::Silhouette(n) => {
                let t = &self.docs[self.active].tree;
                // Start from the node's taxon, else its label (a clade name or the tip).
                let text = match t.nodes[n].attrs.get(crate::scene::PHYLOPIC_KEY) {
                    Some(crate::tree::Attr::Text(s)) => s.clone(),
                    _ => t.label(n).replace('_', " "),
                };
                self.dialog = Some(Dialog::Silhouette { node: n, text });
                return;
            }
            Action::FossilRange(n) => {
                let d = &self.docs[self.active];
                // Start from the existing range, or the tip's own age.
                let (oldest, youngest) = match crate::fossil::range(&d.tree, n) {
                    Some((lo, hi)) => (hi, lo),
                    None => {
                        let age = fossil_age_of(d, n);
                        (age, age)
                    }
                };
                self.dialog = Some(Dialog::FossilRange { node: n, oldest, youngest });
                return;
            }
            Action::AddTaxon(n, mode) => {
                self.dialog = Some(Dialog::AddTaxon { node: n, mode, text: String::new() });
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
            Action::CutClade(n) => {
                if n == self.docs[self.active].tree.root {
                    return;
                }
                self.apply(Action::CopyClade(n));
                let tips = self.docs[self.active].tree.tips_below(n).len();
                self.apply(Action::Drop(vec![n]));
                if self.docs[self.active].tree.reachable()[n] {
                    return; // the drop failed; its status explains why
                }
                self.status = format!("Cut clade of {} tips (also on the clipboard as Newick). Select a node in any tree and paste.", tips);
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
                self.dialog = Some(Dialog::Transform { node: n, kind: ops::BranchTransform::Scale, value: 1.0, base: doc.tree.clone() });
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
        let mut note: Option<String> = None;
        let result: anyhow::Result<()> = (|| {
            match a {
                Action::Extract(_)
                | Action::CladeLabel(_)
                | Action::Rename(_)
                | Action::AssignImage(_)
                | Action::RefetchImage(_)
                | Action::CopyClade(_)
                | Action::CutClade(_)
                | Action::Paste(..)
                | Action::CopyNames(_)
                | Action::CopyData(_)
                | Action::CopyVcv(_)
                | Action::Polytomies(..)
                | Action::AddTaxon(..)
                | Action::FossilRange(_)
                | Action::Silhouette(_)
                | Action::ThinClade(_)
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
                    let names: std::collections::HashSet<String> = nodes.iter().flat_map(|&n| doc.tree.tips_below(n)).map(|t| doc.tree.label(t).to_string()).collect();
                    doc.checkpoint();
                    if let Err(e) = ops::drop_nodes(&mut doc.tree, &nodes) {
                        doc.undo();
                        return Err(e);
                    }
                    doc.selection.clear();
                    // A tree set: the same tips go from every tree.
                    if let Some((done, skipped)) = doc.apply_to_set(&|t| keep_named(t, |l| !names.contains(l))) {
                        doc.refresh_set();
                        note = Some(set_note(&format!("Dropped {}", count(names.len(), "tip", "tips")), done, skipped));
                    }
                }
                Action::KeepOnly(nodes) => {
                    let tips: Vec<NodeId> = nodes.iter().flat_map(|&n| doc.tree.tips_below(n)).collect();
                    let names: std::collections::HashSet<String> = tips.iter().map(|&t| doc.tree.label(t).to_string()).collect();
                    doc.checkpoint();
                    if let Err(e) = ops::keep_tips(&mut doc.tree, &tips) {
                        doc.undo();
                        return Err(e);
                    }
                    doc.selection.clear();
                    if let Some((done, skipped)) = doc.apply_to_set(&|t| keep_named(t, |l| names.contains(l))) {
                        doc.refresh_set();
                        note = Some(set_note(&format!("Kept {}", count(names.len(), "tip", "tips")), done, skipped));
                    }
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
        if let Some(n) = note {
            self.status = n;
        }
        if let Err(e) = result {
            self.status = format!("{:#}", e);
        }
    }

    fn open_dialog(&mut self) {
        let files = rfd::FileDialog::new()
            .add_filter("Trees and data", &["nwk", "newick", "tre", "tree", "trees", "nex", "nexus", "nxs", "t", "con", "treefile", "contree", "xml", "phyloxml", "json", "txt", "csv", "tsv"])
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
            if o.format == TreeFormat::Nexus {
                ui.checkbox(&mut o.canopy, "Canopy settings").on_hover_text("Add a CANOPY block with the layout, layers, colors, comparative data and tree-set settings, restored when the file is opened in Canopy. Other programs ignore it.");
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
            TreeFormat::Nexus => {
                let mut s = nexus::write_nexus(&trees, &opts);
                if o.canopy {
                    match doc.canopy_block() {
                        Ok(b) => s += &b,
                        Err(e) => self.status = format!("Canopy settings not saved: {:#}", e),
                    }
                }
                s
            }
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
        match a {
            panels::PosteriorAction::Rogue => {
                self.rogue = Some(RogueAnalysis { name: doc.name.clone(), trees: p.trees.clone(), skip: p.trees.len() - p.kept().len(), max_size: 1, min_prob: 0.1, result: None, running: None });
            }
            panels::PosteriorAction::OpenSample(i) => {
                let name = format!("{} – sample {}", doc.name, i + 1);
                let mut d = Document::new(&name, p.trees[i].clone(), auto);
                d.view.layout = doc.view.layout.clone();
                self.status = format!("Created {}", name);
                self.add_doc(d);
            }
            panels::PosteriorAction::Consensus | panels::PosteriorAction::Mcc => {
                // Built in the background (reusing the cached clade summary when there
                // is one) so large tree sets don't freeze the window.
                let mcc = matches!(a, panels::PosteriorAction::Mcc);
                let label = if mcc { "MCC".to_string() } else { format!("Consensus {:.0}%", p.threshold * 100.0) };
                let (trees, skip) = (p.trees.clone(), p.trees.len() - p.kept().len());
                let (cached, rooted, key, threshold, heights) = (p.cached_summary(), p.rooted, p.kind.support_key(), p.threshold as f64, p.heights);
                let (tx, rx) = std::sync::mpsc::channel();
                let job_trees = trees.clone();
                std::thread::spawn(move || {
                    let res = (|| -> anyhow::Result<Tree> {
                        let s = match cached {
                            Some(s) => s,
                            None => {
                                let mut s = consensus::summarize(&job_trees[skip..], rooted)?;
                                s.support_key = key;
                                Arc::new(s)
                            }
                        };
                        Ok(if mcc { consensus::mcc_tree(&job_trees[skip..], &s, heights)? } else { consensus::consensus(&s, threshold, heights) })
                    })();
                    let _ = tx.send(res.map_err(|e| format!("{:#}", e)));
                });
                self.status = format!("Building the {} tree from {} trees…", label, trees.len() - skip);
                self.job = Some(SummaryJob { label, doc_name: doc.name.clone(), layout: doc.view.layout.clone(), sample: (trees, skip), started: std::time::Instant::now(), rx });
            }
        }
    }

    /// Collect a finished consensus / MCC tree and open it in a new tab.
    fn poll_job(&mut self, ctx: &egui::Context) {
        let Some(job) = &self.job else { return };
        let res = match job.rx.try_recv() {
            Ok(res) => res,
            Err(std::sync::mpsc::TryRecvError::Empty) => {
                ctx.request_repaint_after(std::time::Duration::from_millis(200));
                return;
            }
            Err(std::sync::mpsc::TryRecvError::Disconnected) => Err("the computation stopped unexpectedly".into()),
        };
        let job = self.job.take().unwrap();
        match res {
            Ok(t) => {
                let name = format!("{} – {}", job.doc_name, job.label);
                let mut d = Document::new(&name, t, self.settings.auto_phylopic);
                d.view.layout = job.layout;
                // Consensus/MCC tabs keep a reference to the post-burn-in sample
                // (shared, not copied) for DensiTree only.
                d.source_sample = Some(job.sample);
                self.status = format!("Created {} in {:.1}s", name, job.started.elapsed().as_secs_f32());
                self.add_doc(d);
            }
            Err(e) => self.status = format!("{}: {}", job.label, e),
        }
        ctx.request_repaint();
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
            // Clipboard: Ctrl+C copies the selected clade, Ctrl+X cuts it; Ctrl+V pastes it
            // (or Newick text from another program) next to the selection.
            let (copy, cut, paste) = ctx.input(|i| {
                let copy = i.events.iter().any(|e| matches!(e, egui::Event::Copy));
                let cut = i.events.iter().any(|e| matches!(e, egui::Event::Cut));
                let paste = i.events.iter().find_map(|e| if let egui::Event::Paste(s) = e { Some(s.clone()) } else { None });
                (copy, cut, paste)
            });
            let selected = self.docs.get(self.active).and_then(|d| d.single_selection());
            // Ctrl+A selects every species.
            if ctx.input_mut(|i| i.consume_shortcut(&KeyboardShortcut::new(Modifiers::COMMAND, Key::A))) {
                if let Some(d) = self.doc() {
                    d.selection = d.tree.tips().into_iter().collect();
                    self.status = format!("Selected all {} species (Ctrl+C copies their names)", d.selection.len());
                }
            }
            // Ctrl+C: a single selected clade is copied for pasting into another
            // tree; otherwise species names (the selection, or all) are copied.
            if copy {
                let single_clade = selected.filter(|&n| self.docs.get(self.active).map(|d| !d.tree.is_tip(n)).unwrap_or(false));
                match single_clade {
                    Some(n) => self.apply(Action::CopyClade(n)),
                    None => self.copy_tip_names(),
                }
            }
            // Ctrl+X: cut the selected clade (copy it, then remove it).
            if cut {
                if let Some(n) = selected {
                    self.apply(Action::CutClade(n));
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
                    let tid = egui::Id::new("taxon-grid");
                    table::headers(ui, tid, &["Tip", "Status", "Open Tree name", "Rename"]);
                    ui.end_row();
                    let keys: Vec<Vec<String>> = tc.results.iter().enumerate().map(|(i, r)| vec![r.query.clone(), format!("{:?}", r.status), r.suggestion().map(|s| s.to_string()).unwrap_or_default(), if tc.rename[i] { "1".into() } else { "0".into() }]).collect();
                    for i in table::order(ui, tid, &keys) {
                        let r = &tc.results[i];
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

    /// Results of the node-density artifact (delta) test for the whole tree, with a
    /// plot of nodes against path length and the fitted curve n = βx^δ.
    fn delta_window(&mut self, ctx: &egui::Context) {
        let Some(result) = &self.delta_report else { return };
        let mut open = true;
        let mut save = false;
        let mut rerun = false;
        egui::Window::new("Node-density artifact (delta test)").open(&mut open).default_width(480.0).show(ctx, |ui| {
            let r = match result {
                Ok(r) => r,
                Err(e) => {
                    ui.colored_label(Color32::from_rgb(198, 40, 40), e);
                    return;
                }
            };
            let col = if r.artifact() { Color32::from_rgb(214, 96, 0) } else { Color32::from_rgb(46, 125, 50) };
            ui.colored_label(col, egui::RichText::new(r.verdict()).strong());
            ui.add_space(4.0);
            let f = crate::tree::format_num;
            egui::Grid::new("delta-grid").num_columns(2).striped(true).show(ui, |ui| {
                ui.label("Species");
                ui.label(r.labels.len().to_string());
                ui.end_row();
                ui.label("Slope β*").on_hover_text("How much root-to-tip path length increases per node along the path");
                ui.label(f(r.beta_star, 4));
                ui.end_row();
                ui.label("Likelihood ratio").on_hover_text("Evidence that path length is related to the number of nodes (directional vs random-walk model)");
                ui.label(format!("{}  (p = {})", f(r.lr, 3), f(r.p_value, 4)));
                ui.end_row();
                ui.label("δ").on_hover_text("Curvature of nodes against path length (n = βx^δ); δ > 1 is the shape the artifact produces");
                ui.label(crate::nodedensity::DeltaTest::delta_text(r.delta));
                ui.end_row();
                ui.label("δ (alternative)").on_hover_text("δ from the equivalent curve x = β*·n^(1/δ); can be the better estimate when δ is very large");
                ui.label(crate::nodedensity::DeltaTest::delta_text(r.delta_alt));
                ui.end_row();
            });
            ui.add_space(4.0);
            // Nodes (y) against path length (x), with the fitted curve.
            let (resp, painter) = ui.allocate_painter(egui::vec2(ui.available_width().min(460.0), 260.0), egui::Sense::hover());
            painter.rect_filled(resp.rect, 0.0, Color32::WHITE);
            let scene = crate::plots::node_density(r, resp.rect.width(), resp.rect.height(), 1.0);
            canvas::paint_scene(&painter, &scene, resp.rect.min, &mut self.textures, &self.pics);
            ui.small("Orange: the fitted curve n = βx^δ. Check that it follows the points; if they form near-vertical stacks it can fit poorly, and the alternative δ may be the better estimate.");
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                if ui.button("Copy report").on_hover_text("Summary and per-species data, tab-separated").clicked() {
                    ui.ctx().copy_text(r.report());
                }
                if ui.button("Save…").clicked() {
                    if let Some(path) = rfd::FileDialog::new().add_filter("Text", &["txt", "tsv"]).set_file_name("delta_test.txt").save_file() {
                        let _ = std::fs::write(path, r.report());
                    }
                }
                if ui.button("Recalculate").on_hover_text("Run again on the tree shown now").clicked() {
                    rerun = true;
                }
                save = ui.button("Save plot…").on_hover_text("PNG (300 dpi), SVG or PDF").clicked();
            });
            ui.add_space(4.0);
            ui.small("Webster, A. J., Payne, R. J. & Pagel, M. (2003) Molecular phylogenies link rates of evolution and speciation. Science 301: 478.");
            ui.small("Venditti, C., Meade, A. & Pagel, M. (2006) Detecting the node-density artifact in phylogeny reconstruction. Systematic Biology 55: 637–643.");
        });
        if save {
            if let Some(Ok(r)) = self.delta_report.clone() {
                self.save_plot(&|w, h, pt| crate::plots::node_density(&r, w, h, pt), "delta_test.png");
            }
        }
        if rerun {
            self.run_delta_test();
        } else if !open {
            self.delta_report = None;
        }
    }

    /// Read a table of fossil ranges onto the tips of the active tree and show them.
    fn import_fossil_ranges(&mut self) {
        if self.docs.get(self.active).is_none() {
            self.status = "Open a tree first.".into();
            return;
        }
        let Some(path) = rfd::FileDialog::new().add_filter("Table", &["csv", "tsv", "txt", "tab"]).pick_file() else { return };
        let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        let text = match std::fs::read_to_string(&path) {
            Ok(t) => t,
            Err(e) => return self.fail(format!("Could not read {}", name), e.to_string()),
        };
        let doc = &mut self.docs[self.active];
        doc.checkpoint();
        match crate::fossil::import(&text, &name, &mut doc.tree) {
            Ok((m, unmatched)) => {
                ensure_fossil_layer(doc);
                self.status = if unmatched.is_empty() {
                    format!("Fossil ranges for {} taxa from {}", m, name)
                } else {
                    let more = if unmatched.len() > 8 { format!(" and {} more", unmatched.len() - 8) } else { String::new() };
                    format!("Fossil ranges for {} taxa from {}; not in the tree: {}{}", m, name, unmatched.iter().take(8).cloned().collect::<Vec<_>>().join(", "), more)
                };
            }
            Err(e) => {
                doc.undo();
                self.fail(format!("Could not read {}", name), format!("{:#}", e));
            }
        }
    }

    /// Choose the silhouette a node shows: a taxon looked up on PhyloPic, or an image file.
    fn silhouette_dialog(&mut self, ctx: &egui::Context) {
        let Some(Dialog::Silhouette { node, text }) = &mut self.dialog else { return };
        let (mut show, mut file, mut remove, mut cancel) = (false, false, false, false);
        egui::Window::new("PhyloPic silhouette").collapsible(false).resizable(false).show(ctx, |ui| {
            ui.label("Taxon to look up on phylopic.org (species, genus or a higher group such as Hominidae)");
            let r = ui.text_edit_singleline(text);
            if r.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                show = true;
            }
            ui.horizontal(|ui| {
                show |= ui.button("Show").clicked();
                file = ui.button("Use image file…").on_hover_text("Your own PNG or JPEG, shown under the name above").clicked();
                remove = ui.button("Remove").clicked();
                cancel = ui.button("Cancel").clicked();
            });
            ui.small("Size or color silhouettes by a trait in the PhyloPic layer's settings.");
        });
        let (n, name) = (*node, text.trim().to_string());
        if show || file || remove || cancel {
            self.dialog = None;
        }
        let Some(d) = self.docs.get_mut(self.active) else { return };
        if remove {
            d.checkpoint();
            d.tree.nodes[n].attrs.remove(crate::scene::PHYLOPIC_KEY);
            return;
        }
        if !(show || file) {
            return;
        }
        let name = if name.is_empty() { format!("node {}", n) } else { name };
        if file {
            let Some(path) = rfd::FileDialog::new().add_filter("Image", &["png", "jpg", "jpeg"]).pick_file() else { return };
            match image::open(&path) {
                Ok(img) => {
                    let key = image_key(&name);
                    self.pics.set_custom(&key, img.to_rgba8());
                    self.textures.invalidate(&key);
                    d.view.tip_images.insert(key, path.display().to_string());
                }
                Err(e) => return self.status = format!("Could not open image: {}", e),
            }
        }
        d.checkpoint();
        d.tree.nodes[n].attrs.insert(crate::scene::PHYLOPIC_KEY.into(), crate::tree::Attr::Text(name.clone()));
        // A PhyloPic layer draws it; one added for a node starts without tip silhouettes.
        if !d.layers.iter().any(|l| matches!(l.layer, Layer::Phylopic(_))) {
            let mut l = Layer::default_phylopic();
            if let Layer::Phylopic(s) = &mut l {
                s.tips = d.tree.is_tip(n);
            }
            d.add_layer(l);
        }
        for e in d.layers.iter_mut() {
            if let Layer::Phylopic(s) = &mut e.layer {
                e.enabled = true;
                if d.tree.is_tip(n) {
                    s.tips = true;
                } else {
                    s.nodes = true;
                }
            }
        }
        self.status = format!("Silhouette of {} at this node{}", name, if file { "" } else { " (looked up on phylopic.org)" });
    }

    /// Edit one tip's fossil range.
    fn fossil_dialog(&mut self, ctx: &egui::Context) {
        let Some(Dialog::FossilRange { node, oldest, youngest }) = &mut self.dialog else { return };
        let (mut ok, mut remove, mut cancel) = (false, false, false);
        let label = self.docs.get(self.active).map(|d| d.tree.label(*node).replace('_', " ")).unwrap_or_default();
        egui::Window::new(format!("Fossil range: {}", label)).collapsible(false).resizable(false).show(ctx, |ui| {
            egui::Grid::new("fossil-range").num_columns(2).show(ui, |ui| {
                ui.label("Oldest (first appearance)");
                ui.add(egui::DragValue::new(oldest).speed(0.1).range(0.0..=5000.0).suffix(" Ma"));
                ui.end_row();
                ui.label("Youngest (last appearance)");
                ui.add(egui::DragValue::new(youngest).speed(0.1).range(0.0..=5000.0).suffix(" Ma"));
                ui.end_row();
            });
            ui.small("Ages in Ma, placed with the geologic timescale layer's settings if there is one, else the time axis's youngest tip age; otherwise branch lengths are read as Ma with the youngest tip at the present.");
            ui.horizontal(|ui| {
                ok = ui.button("OK").clicked();
                remove = ui.button("Remove range").clicked();
                cancel = ui.button("Cancel").clicked();
            });
        });
        let (n, a, b) = (*node, *oldest, *youngest);
        if ok || remove {
            if let Some(d) = self.doc() {
                d.checkpoint();
                if ok {
                    crate::fossil::set_range(&mut d.tree, n, a, b);
                    ensure_fossil_layer(d);
                } else {
                    crate::fossil::clear_range(&mut d.tree, n);
                }
            }
        }
        if ok || remove || cancel {
            self.dialog = None;
        }
    }

    fn import_varrates(&mut self) {
        self.import_bayestraits(crate::bayestraits::RunKind::VariableRates);
    }

    fn import_fabric(&mut self) {
        self.import_bayestraits(crate::bayestraits::RunKind::Fabric);
    }

    /// Read BayesTraits VarRates files and summarize them on the active tree.
    /// Several files are runs of the same analysis: their samples are pooled into
    /// one posterior, with convergence across runs reported (R-hat, ESS, posterior
    /// probabilities per run). The files' entries decide the kind of run
    /// ("Branch" = variable rates, "LS_Beta" = Fabric); files with only Node
    /// scalars are read as the kind `chosen`.
    fn import_bayestraits(&mut self, chosen: crate::bayestraits::RunKind) {
        self.read_bayestraits(chosen, false);
    }

    /// Pick VarRates files and read them (see `import_bayestraits`). From the
    /// MCMC diagnostics window (`from_diagnostics`) their diagnostics tab is
    /// shown, and with no tree open only the diagnostics are built.
    fn read_bayestraits(&mut self, chosen: crate::bayestraits::RunKind, from_diagnostics: bool) {
        use crate::bayestraits::RunKind;
        let diagnostics_only = self.docs.get(self.active).is_none();
        if diagnostics_only && !from_diagnostics {
            self.status = "Open the tree used in the BayesTraits run first.".into();
            return;
        }
        let Some(mut paths) = rfd::FileDialog::new().set_title("BayesTraits VarRates file(s): pick several runs to combine them").add_filter("BayesTraits VarRates", &["txt"]).pick_files() else { return };
        paths.sort();
        // Folder and file name: runs often share a file name (e.g. a full run and a short test run).
        let name_of = |path: &std::path::Path| {
            let file = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
            match path.parent().and_then(|p| p.file_name()) {
                Some(dir) => format!("{} › {}", dir.to_string_lossy(), file),
                None => file,
            }
        };
        let names: Vec<String> = paths.iter().map(|p| name_of(p)).collect();
        let mut texts = Vec::new();
        for (p, nm) in paths.iter().zip(&names) {
            match std::fs::read_to_string(p) {
                Ok(t) => texts.push(t),
                Err(e) => return self.fail(format!("Could not read {}", nm), e.to_string()),
            }
        }
        let mut detected = None;
        for (i, (t, nm)) in texts.iter().zip(&names).enumerate() {
            match crate::bayestraits::detect(t) {
                Ok(d) if i == 0 => detected = d,
                Ok(d) if d != detected => return self.fail("Runs don't match".into(), format!("{} is not the same kind of run as {}", nm, names[0])),
                Ok(_) => {}
                Err(e) => return self.fail(format!("Could not read {}", nm), format!("{:#}", e)),
            }
        }
        let kind = detected.unwrap_or(chosen);
        let how = match detected {
            Some(RunKind::Fabric) => "Fabric run (LS_Beta entries)",
            Some(RunKind::VariableRates) => "variable-rates run (Branch entries)",
            None => "only Node scalars; read as the chosen kind",
        };
        let label = runs_label(&names);
        let doc = self.active;
        match kind {
            RunKind::VariableRates => match texts.iter().map(|t| crate::varrates::parse(t)).collect::<anyhow::Result<Vec<_>>>() {
                Ok(runs) if diagnostics_only => {
                    let mut d = crate::varrates::diagnostics(&runs, &names);
                    d.title = label;
                    d.add_params(sibling_log_params(&paths));
                    self.show_mcmc(d);
                }
                Ok(runs) => {
                    if self.apply_varrates(RunFiles { names, paths, runs }, doc, Some(how)) && from_diagnostics {
                        if let Some(d) = self.varrates_diag.clone() {
                            self.show_mcmc(d);
                        }
                    }
                }
                Err(e) => self.fail(format!("Could not read {}", label), format!("{:#}", e)),
            },
            RunKind::Fabric => match texts.iter().map(|t| crate::fabric::parse(t)).collect::<anyhow::Result<Vec<_>>>() {
                Ok(runs) if diagnostics_only => {
                    let mut d = crate::fabric::diagnostics(&runs, &names);
                    d.title = label;
                    d.add_params(sibling_log_params(&paths));
                    self.show_mcmc(d);
                }
                Ok(runs) => {
                    if self.apply_fabric(RunFiles { names, paths, runs }, doc, Some(how)) && from_diagnostics {
                        if let Some(d) = self.fabric_diag.clone() {
                            self.show_mcmc(d);
                        }
                    }
                }
                Err(e) => self.fail(format!("Could not read {}", label), format!("{:#}", e)),
            },
        }
    }

    /// Pool variable-rates runs into one posterior, summarize it on tab `doc`
    /// and build their diagnostics (replacing an open diagnostics tab of the
    /// earlier runs). False if the runs could not be combined.
    fn apply_varrates(&mut self, files: RunFiles<crate::varrates::VarRates>, doc: usize, how: Option<&str>) -> bool {
        let label = runs_label(&files.names);
        let vr = match crate::varrates::VarRates::combine(&files.runs) {
            Ok(vr) => vr,
            Err(e) => {
                self.fail(format!("Could not combine {}", label), format!("{:#}", e));
                return false;
            }
        };
        let Some(d) = self.docs.get_mut(doc) else { return false };
        let mut diag = crate::varrates::diagnostics(&files.runs, &files.names);
        diag.title = label.clone();
        diag.add_params(sibling_log_params(&files.paths));
        let mut rates = vr.summarize();
        d.checkpoint();
        let warnings = vr.apply(&mut rates, &mut d.tree);
        let matched = rates.iter().filter(|r| r.tree_node.is_some()).count();
        let conv = diag.report.as_ref().map(|c| format!("; {}", runs_verdict(c))).unwrap_or_default();
        let how = how.map(|h| format!(": {}", h)).unwrap_or_default();
        self.status = format!("{}{}; {} iterations; rate scalars on {} of {} branches{}", label, how, vr.iterations.len(), matched, rates.len(), conv);
        self.varrates = Some((format!("{} ({} iterations)", label, vr.iterations.len()), rates, warnings));
        self.varrates_doc = doc;
        self.show_varrates |= !how.is_empty(); // imported, not a run unloaded
        if let Some(old) = self.varrates_diag.replace(diag.clone()) {
            self.replace_mcmc_tab(&old.title, diag);
        }
        self.varrates_runs = files;
        true
    }

    /// As `apply_varrates`, for Fabric runs.
    fn apply_fabric(&mut self, files: RunFiles<crate::fabric::Fabric>, doc: usize, how: Option<&str>) -> bool {
        let label = runs_label(&files.names);
        let fab = match crate::fabric::Fabric::combine(&files.runs) {
            Ok(f) => f,
            Err(e) => {
                self.fail(format!("Could not combine {}", label), format!("{:#}", e));
                return false;
            }
        };
        let Some(d) = self.docs.get_mut(doc) else { return false };
        let mut diag = crate::fabric::diagnostics(&files.runs, &files.names);
        diag.title = label.clone();
        diag.add_params(sibling_log_params(&files.paths));
        let mut rows = fab.summarize();
        d.checkpoint();
        let warnings = fab.apply(&mut rows, &mut d.tree);
        let conv = diag.report.as_ref().map(|c| format!("; {}", runs_verdict(c))).unwrap_or_default();
        let how = how.map(|h| format!(": {}", h)).unwrap_or_default();
        self.status = format!("{}{}; {} iterations summarized{}", label, how, fab.iterations.len(), conv);
        self.fabric = Some((format!("{} ({} iterations)", label, fab.iterations.len()), rows, warnings));
        self.fabric_doc = doc;
        self.show_fabric |= !how.is_empty(); // imported, not a run unloaded
        if let Some(old) = self.fabric_diag.replace(diag.clone()) {
            self.replace_mcmc_tab(&old.title, diag);
        }
        self.fabric_runs = files;
        true
    }

    /// Put `d` in place of the MCMC diagnostics tab titled `old`, if it is open.
    fn replace_mcmc_tab(&mut self, old: &str, d: crate::mcmc::Diagnostics) {
        if let Some(w) = self.mcmc_win.as_mut() {
            if let Some(t) = w.tabs.iter_mut().find(|t| t.title == old) {
                *t = d;
            }
        }
    }

    /// Give the branches colored by `attr` and the node-scalar circles of the
    /// active tree the palette chosen in a BayesTraits window.
    fn set_bayestraits_palette(&mut self, attr: &str, palette: Palette) {
        let Some(d) = self.doc() else { return };
        d.checkpoint();
        for e in d.layers.iter_mut() {
            match &mut e.layer {
                Layer::Tree(s) if s.color_by.as_deref() == Some(attr) => s.palette = palette,
                Layer::NodeScalars(s) => s.palette = palette,
                _ => {}
            }
        }
    }

    /// Switch the branch colors and labels of the active tree from attribute `from`
    /// to `to` (mean ↔ median in the BayesTraits windows).
    fn swap_bayestraits_attr(&mut self, from: &str, to: &str) {
        let Some(d) = self.doc() else { return };
        d.checkpoint();
        for e in d.layers.iter_mut() {
            match &mut e.layer {
                Layer::Tree(s) if s.color_by.as_deref() == Some(from) => s.color_by = Some(to.into()),
                Layer::BranchLengths(s) if s.attr.as_deref() == Some(from) => s.attr = Some(to.into()),
                _ => {}
            }
        }
    }

    /// Per-branch rate scalars from BayesTraits, with actions to color or rescale the tree.
    fn varrates_window(&mut self, ctx: &egui::Context) {
        if !self.show_varrates {
            return;
        }
        let Some((name, rates, warnings)) = &self.varrates else { return };
        let mut open = true;
        let mut open_diag = false;
        let (mut color, mut scaled, mut select) = (false, false, None);
        let (mut branch_labels, mut node_labels) = (false, false);
        let (header, rows) = crate::varrates::table(rates);
        let (mut palette, mut repaint) = (self.varrates_palette, false);
        let mut median = self.varrates_median;
        let median_was = median;
        let attr = if median { "rate_scalar_median" } else { "rate_scalar" };
        let mut legend = self.docs.get(self.active).map(|d| d.view.legend_scale).unwrap_or(1.0);
        let legend_was = legend;
        egui::Window::new(format!("Variable rates: {}", name)).open(&mut open).default_width(560.0).show(ctx, |ui| {
            ui.small("Mean rate scalar per branch across the MCMC iterations (node scalars rescale a clade and its stem branch; branch scalars one branch; overlapping scalars multiply). Branch / Node Scalar PP: posterior probability (proportion of iterations) of a scalar on that branch or at that node. Mean / Median Node Scalar: the value of the node scalar placed at that node, over the iterations with one (1 if none). 95% HPD: the shortest interval holding 95% of the same values (the branch scalar over all iterations, the node scalar over the iterations with one).");
            ui.horizontal_wrapped(|ui| {
                ui.small("Reference: Venditti, C., Meade, A. & Pagel, M. (2011) Multiple routes to mammalian diversity. Nature 479: 393–396.");
                ui.hyperlink_to(egui::RichText::new("doi:10.1038/nature10516").small(), "https://doi.org/10.1038/nature10516");
            });
            if ui.button("MCMC diagnostics…").on_hover_text("Convergence and posterior plots for these runs: R-hat, ESS, density and box plots of every parameter").clicked() {
                open_diag = true;
            }
            for w in warnings {
                ui.colored_label(Color32::from_rgb(198, 40, 40), w);
            }
            ui.horizontal_wrapped(|ui| {
                ui.label("Palette");
                repaint = panels::palette_combo(ui, "vr_palette", &mut palette);
                ui.add(egui::Slider::new(&mut legend, 0.5..=3.0).text("legend size"));
                ui.label("Plot");
                ui.selectable_value(&mut median, false, "Mean");
                ui.selectable_value(&mut median, true, "Median");
                if ui.button(if median { "Color branches by median scalar" } else { "Color branches by mean scalar" }).clicked() {
                    color = true;
                }
                if ui.button(if median { "Show median scalars on branches" } else { "Show mean scalars on branches" }).on_hover_text("Label each rescaled branch with its scalar (unscaled branches, = 1, are left blank)").clicked() {
                    branch_labels = true;
                }
                if ui.button("Show node scalars (% of posterior)").on_hover_text("Circles at nodes where a node scalar was placed, colored by its posterior (proportion of the sample); set a minimum posterior in the layer panel").clicked() {
                    node_labels = true;
                }
                if ui.button("Open mean rate-scaled tree").on_hover_text("New tab: each branch length times its mean scalar").clicked() {
                    scaled = true;
                }
                if ui.button("Copy table").on_hover_text("Tab-separated, 6 decimals").clicked() {
                    ui.ctx().copy_text(io::data::table_to_text(&header, &rows, b'\t'));
                }
                if ui.button("Save CSV…").clicked() {
                    if let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).set_file_name("variable_rates.csv").save_file() {
                        let _ = std::fs::write(path, io::data::table_to_text(&header, &rows, b','));
                    }
                }
            });
            egui::ScrollArea::both().max_height(380.0).show(ui, |ui| {
                egui::Grid::new("varrates").striped(true).show(ui, |ui| {
                    let tid = egui::Id::new("varrates");
                    table::headers(ui, tid, &header);
                    ui.end_row();
                    // Same cells (and precision) as Copy table / Save CSV.
                    for k in table::order(ui, tid, &rows) {
                        let (r, cells) = (&rates[k], &rows[k]);
                        if let Some(n) = r.tree_node {
                            if ui.link(r.id.to_string()).on_hover_text("Select in the tree").clicked() {
                                select = Some(n);
                            }
                        } else {
                            ui.label(r.id.to_string());
                        }
                        for c in &cells[1..] {
                            ui.label(c);
                        }
                        ui.end_row();
                    }
                });
            });
        });
        if let Some(n) = select {
            if let Some(d) = self.doc() {
                d.selection.clear();
                d.selection.insert(n);
            }
        }
        self.varrates_palette = palette;
        self.varrates_median = median;
        if median != median_was {
            let (from, to) = if median { ("rate_scalar", "rate_scalar_median") } else { ("rate_scalar_median", "rate_scalar") };
            self.swap_bayestraits_attr(from, to);
        }
        if legend != legend_was {
            if let Some(d) = self.doc() {
                d.view.legend_scale = legend;
            }
        }
        if repaint {
            self.set_bayestraits_palette(attr, palette);
        }
        if color {
            if let Some(d) = self.doc() {
                d.checkpoint();
                for e in d.layers.iter_mut() {
                    if let Layer::Tree(s) = &mut e.layer {
                        s.color_by = Some(attr.into());
                        s.palette = palette;
                        s.width = s.width.max(1.5);
                    }
                }
            }
        }
        if branch_labels {
            if let Some(d) = self.doc() {
                let exists = d.layers.iter().any(|e| matches!(&e.layer, Layer::BranchLengths(s) if s.attr.as_deref() == Some(attr)));
                if !exists {
                    d.add_layer(Layer::BranchLengths(BranchLengthStyle { size: 7.0, color: Color::rgb(90, 90, 90), digits: 2, attr: Some(attr.into()), hide_value: Some(1.0), min_pp: 0.0, min_value: None, min_abs: false }));
                }
            }
        }
        if node_labels {
            if let Some(d) = self.doc() {
                match d.layers.iter().position(|e| matches!(e.layer, Layer::NodeScalars(_))) {
                    // Already there: select it so its threshold can be adjusted.
                    Some(i) => d.selected_layer = Some(i),
                    None => {
                        let mut l = Layer::default_node_scalars();
                        if let Layer::NodeScalars(s) = &mut l {
                            s.palette = palette;
                        }
                        d.add_layer(l);
                    }
                }
            }
        }
        if scaled {
            if let Some(d) = self.docs.get(self.active) {
                let t = crate::varrates::scaled_tree(&d.tree);
                let name = format!("{} – mean rate-scaled", d.name);
                let mut nd = Document::new(&name, t, self.settings.auto_phylopic);
                nd.view.layout = d.view.layout.clone();
                self.add_doc(nd);
            }
        }
        if open_diag {
            if let Some(d) = self.varrates_diag.clone() {
                self.show_mcmc(d);
            }
        }
        // Closing only hides the window: the results stay loaded (right panel › Loaded files).
        if !open {
            self.show_varrates = false;
        }
    }

    /// The Fabric summary table.
    fn fabric_window(&mut self, ctx: &egui::Context) {
        if !self.show_fabric {
            return;
        }
        let Some((name, rows, warnings)) = &self.fabric else { return };
        let mut open = true;
        let mut open_diag = false;
        let mut select = None;
        let (mut color, mut scaled, mut branch_labels, mut node_labels) = (false, false, false, false);
        let (header, cells) = crate::fabric::table(rows);
        let (mut palette, mut repaint) = (self.fabric_palette, false);
        let mut median = self.fabric_median;
        let median_was = median;
        let attr = if median { "fabric_beta_median" } else { "fabric_beta_mean" };
        let mut legend = self.docs.get(self.active).map(|d| d.view.legend_scale).unwrap_or(1.0);
        let legend_was = legend;
        egui::Window::new(format!("Fabric: {}", name)).open(&mut open).default_width(560.0).show(ctx, |ui| {
            ui.small("Mean and Median (β * BL): LS_Beta on each branch, over the iterations in which the branch carried one (0 if never). β Scalar PP: posterior probability (proportion of iterations) of an LS_Beta on that branch. Node Scalar PP: posterior probability (proportion of iterations) of a node scalar at that node. Mean / Median Node Scalar: the value of that node scalar, over the iterations with one (1 if none). 95% HPD: the shortest interval holding 95% of the same values (LS_Beta and the node scalar over the iterations with one; 0 and 1 if none).");
            ui.horizontal_wrapped(|ui| {
                ui.small("Reference: Pagel, M., O'Donovan, C. & Meade, A. (2022) General statistical model shows that macroevolutionary patterns and processes are consistent with Darwinian gradualism. Nature Communications 13: 1113.");
                ui.hyperlink_to(egui::RichText::new("doi:10.1038/s41467-022-28595-z").small(), "https://doi.org/10.1038/s41467-022-28595-z");
            });
            if ui.button("MCMC diagnostics…").on_hover_text("Convergence and posterior plots for these runs: R-hat, ESS, density and box plots of every parameter").clicked() {
                open_diag = true;
            }
            for w in warnings {
                ui.colored_label(Color32::from_rgb(198, 40, 40), w);
            }
            ui.horizontal_wrapped(|ui| {
                ui.label("Palette");
                repaint = panels::palette_combo(ui, "fab_palette", &mut palette);
                ui.add(egui::Slider::new(&mut legend, 0.5..=3.0).text("legend size"));
                ui.label("Plot");
                ui.selectable_value(&mut median, false, "Mean");
                ui.selectable_value(&mut median, true, "Median");
                if ui.button(if median { "Color branches by median β" } else { "Color branches by mean β" }).clicked() {
                    color = true;
                }
                if ui.button(if median { "Show median β on branches" } else { "Show mean β on branches" }).on_hover_text("Label each branch with its (β * BL) (branches without one, = 0, are left blank)").clicked() {
                    branch_labels = true;
                }
                if ui.button("Show node scalars (% of posterior)").on_hover_text("Circles at nodes where a node scalar was placed, colored by its posterior (proportion of the sample); set a minimum posterior in the layer panel").clicked() {
                    node_labels = true;
                }
                if ui.button("Open mean node-scalar-scaled tree").on_hover_text("New tab: each branch length times the mean product of the node scalars over it (LS_Beta shifts the trait, not the rate, so it is not included)").clicked() {
                    scaled = true;
                }
                if ui.button("Copy table").on_hover_text("Tab-separated").clicked() {
                    ui.ctx().copy_text(io::data::table_to_text(&header, &cells, b'\t'));
                }
                if ui.button("Save CSV…").clicked() {
                    if let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).set_file_name("fabric.csv").save_file() {
                        let _ = std::fs::write(path, io::data::table_to_text(&header, &cells, b','));
                    }
                }
            });
            egui::ScrollArea::both().max_height(380.0).show(ui, |ui| {
                egui::Grid::new("fabric").striped(true).show(ui, |ui| {
                    let tid = egui::Id::new("fabric");
                    table::headers(ui, tid, &header);
                    ui.end_row();
                    // Same cells (and precision) as Copy table / Save CSV.
                    for k in table::order(ui, tid, &cells) {
                        let (r, row) = (&rows[k], &cells[k]);
                        if let Some(n) = r.tree_node {
                            if ui.link(r.id.to_string()).on_hover_text("Select in the tree").clicked() {
                                select = Some(n);
                            }
                        } else {
                            ui.label(r.id.to_string());
                        }
                        for c in &row[1..] {
                            ui.label(c);
                        }
                        ui.end_row();
                    }
                });
            });
        });
        if let Some(n) = select {
            if let Some(d) = self.doc() {
                d.selection.clear();
                d.selection.insert(n);
            }
        }
        self.fabric_palette = palette;
        self.fabric_median = median;
        if median != median_was {
            let (from, to) = if median { ("fabric_beta_mean", "fabric_beta_median") } else { ("fabric_beta_median", "fabric_beta_mean") };
            self.swap_bayestraits_attr(from, to);
        }
        if legend != legend_was {
            if let Some(d) = self.doc() {
                d.view.legend_scale = legend;
            }
        }
        if repaint {
            self.set_bayestraits_palette(attr, palette);
        }
        if color {
            if let Some(d) = self.doc() {
                d.checkpoint();
                for e in d.layers.iter_mut() {
                    if let Layer::Tree(s) = &mut e.layer {
                        s.color_by = Some(attr.into());
                        s.palette = palette;
                        s.width = s.width.max(1.5);
                    }
                }
            }
        }
        if branch_labels {
            if let Some(d) = self.doc() {
                let exists = d.layers.iter().any(|e| matches!(&e.layer, Layer::BranchLengths(s) if s.attr.as_deref() == Some(attr)));
                if !exists {
                    d.add_layer(Layer::BranchLengths(BranchLengthStyle { size: 7.0, color: Color::rgb(90, 90, 90), digits: 2, attr: Some(attr.into()), hide_value: Some(0.0), min_pp: 0.0, min_value: None, min_abs: false }));
                }
            }
        }
        if node_labels {
            if let Some(d) = self.doc() {
                match d.layers.iter().position(|e| matches!(e.layer, Layer::NodeScalars(_))) {
                    Some(i) => d.selected_layer = Some(i),
                    None => {
                        let mut l = Layer::default_node_scalars();
                        if let Layer::NodeScalars(s) = &mut l {
                            s.palette = palette;
                        }
                        d.add_layer(l);
                    }
                }
            }
        }
        if scaled {
            if let Some(d) = self.docs.get(self.active) {
                let t = crate::fabric::scaled_tree(&d.tree);
                let name = format!("{} – mean node-scalar-scaled", d.name);
                let mut nd = Document::new(&name, t, self.settings.auto_phylopic);
                nd.view.layout = d.view.layout.clone();
                self.add_doc(nd);
            }
        }
        if open_diag {
            if let Some(d) = self.fabric_diag.clone() {
                self.show_mcmc(d);
            }
        }
        if !open {
            self.show_fabric = false;
        }
    }

    /// Entropy-based rogue scores and skeleton analysis for the taxa of a tree set
    /// (rogues are flagged in tables; the trees are not changed).
    fn rogue_window(&mut self, ctx: &egui::Context) {
        let Some(ra) = self.rogue.as_mut() else { return };
        let mut open = true;
        let (mut run, mut color) = (false, false);
        let prob = |p: f64| if p < 1e-3 { format!("{:.2e}", p) } else { crate::tree::format_num(p, 3) };
        let table = |scores: &[crate::rogue::RogueScore]| {
            let header: Vec<String> = ["Rogue score", "Clade probability", "Taxa"].iter().map(|s| s.to_string()).collect();
            let rows: Vec<Vec<String>> = scores.iter().map(|s| vec![crate::tree::format_num(s.score, 6), crate::tree::format_num(s.probability, 6), s.taxa.join(" ")]).collect();
            (header, rows)
        };
        egui::Window::new(format!("Rogue taxa: {}", ra.name)).open(&mut open).default_width(600.0).show(ctx, |ui| {
            ui.small("Rogue taxa are taxa (or small clades) whose unstable placement causes most of the uncertainty in a tree set. A clade's rogue score is how much the entropy of the tree set would drop if its taxa were left out; the higher the score, the more of a rogue it is. The skeleton analysis ranks rogues in order: each round it sets aside the top-scoring clade and re-scores the rest. Rogues are only flagged; the trees are not changed.");
            ui.add(egui::Label::new(egui::RichText::new(ROGUE_CITATION).small().weak()).wrap());
            ui.horizontal_wrapped(|ui| {
                ui.label("Score clades of up to");
                ui.add(egui::DragValue::new(&mut ra.max_size).range(1..=10)).on_hover_text("1 = single taxa only");
                ui.label("taxa, with probability ≥");
                ui.add(egui::DragValue::new(&mut ra.min_prob).range(0.0..=1.0).speed(0.01));
                if ui.add_enabled(ra.running.is_none(), egui::Button::new("Run")).clicked() {
                    run = true;
                }
            });
            ui.weak(format!("{} trees after burn-in", ra.trees.len() - ra.skip));
            if let Some((started, _)) = &ra.running {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(format!("Running… {:.0} s", started.elapsed().as_secs_f32()));
                });
            }
            match &ra.result {
                None => {}
                Some(Err(e)) => {
                    ui.colored_label(Color32::from_rgb(198, 40, 40), e);
                }
                Some(Ok((h, map_p, scores, rounds))) => {
                    ui.label(format!("Entropy {} (as uncertain as {} equally likely trees); most probable tree {}", crate::tree::format_num(*h, 3), crate::tree::format_num(h.exp(), 1), prob(*map_p)));
                    ui.horizontal_wrapped(|ui| {
                        if ui.button("Color branches by rogue score").on_hover_text("Black = low, red = high, on the tree in the active tab").clicked() {
                            color = true;
                        }
                        if ui.button("Copy table").on_hover_text("Tab-separated").clicked() {
                            let (h, rows) = table(scores);
                            ui.ctx().copy_text(io::data::table_to_text(&h, &rows, b'\t'));
                        }
                        if ui.button("Save CSV…").clicked() {
                            if let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).set_file_name("rogue_scores.csv").save_file() {
                                let (h, rows) = table(scores);
                                let _ = std::fs::write(path, io::data::table_to_text(&h, &rows, b','));
                            }
                        }
                    });
                    egui::ScrollArea::vertical().id_salt("rogue-scores-scroll").max_height(220.0).show(ui, |ui| {
                        egui::Grid::new("rogue-scores-grid").striped(true).show(ui, |ui| {
                            let tid = egui::Id::new("rogue-scores-grid");
                            table::headers(ui, tid, &["Rank", "Rogue score", "Clade probability", "Taxa"]);
                            ui.end_row();
                            let keys: Vec<Vec<String>> = scores.iter().enumerate().map(|(i, s)| vec![(i + 1).to_string(), s.score.to_string(), s.probability.to_string(), s.taxa.join(", ")]).collect();
                            for i in table::order(ui, tid, &keys) {
                                let s = &scores[i];
                                ui.label((i + 1).to_string());
                                ui.label(crate::tree::format_num(s.score, 3));
                                ui.label(crate::tree::format_num(s.probability, 3));
                                ui.label(egui::RichText::new(s.taxa.join(", ").replace('_', " ")).italics());
                                ui.end_row();
                            }
                        });
                    });
                    ui.add_space(6.0);
                    ui.horizontal(|ui| {
                        ui.strong("Skeleton analysis");
                        if ui.small_button("Copy table").on_hover_text("Tab-separated").clicked() {
                            let header: Vec<String> = ["Round", "Flagged", "Taxa left", "Entropy", "Most probable tree"].iter().map(|s| s.to_string()).collect();
                            let rows: Vec<Vec<String>> = rounds.iter().enumerate().map(|(i, r)| vec![i.to_string(), r.flagged.join(" "), r.taxa_left.to_string(), crate::tree::format_num(r.entropy, 6), format!("{:e}", r.map_probability)]).collect();
                            ui.ctx().copy_text(io::data::table_to_text(&header, &rows, b'\t'));
                        }
                    });
                    ui.small("Entropy and most probable tree are for the taxa not yet flagged.");
                    egui::ScrollArea::vertical().id_salt("rogue-skeleton-scroll").max_height(220.0).show(ui, |ui| {
                        egui::Grid::new("rogue-skeleton-grid").striped(true).show(ui, |ui| {
                            let tid = egui::Id::new("rogue-skeleton-grid");
                            table::headers(ui, tid, &["Round", "Flagged", "Taxa left", "Entropy", "Most probable tree"]);
                            ui.end_row();
                            let keys: Vec<Vec<String>> = rounds.iter().enumerate().map(|(i, r)| vec![i.to_string(), r.flagged.join(", "), r.taxa_left.to_string(), r.entropy.to_string(), r.map_probability.to_string()]).collect();
                            for i in table::order(ui, tid, &keys) {
                                let r = &rounds[i];
                                ui.label(i.to_string());
                                ui.label(egui::RichText::new(if i == 0 { "(all taxa)".to_string() } else { r.flagged.join(", ").replace('_', " ") }).italics());
                                ui.label(r.taxa_left.to_string());
                                ui.label(crate::tree::format_num(r.entropy, 3));
                                ui.label(prob(r.map_probability));
                                ui.end_row();
                            }
                        });
                    });
                }
            }
        });
        if run {
            // In the background, so the window stays usable on large tree sets.
            let (k, t, trees, skip) = (ra.max_size, ra.min_prob, ra.trees.clone(), ra.skip);
            let (tx, rx) = std::sync::mpsc::channel();
            std::thread::spawn(move || {
                let res = crate::rogue::Ccd::from_trees(&trees[skip..]).map(|ccd| (ccd.entropy(), ccd.map_probability(), ccd.rogue_scores(k, t), ccd.skeleton(k, t))).map_err(|e| format!("{:#}", e));
                let _ = tx.send(res);
            });
            ra.running = Some((std::time::Instant::now(), rx));
        }
        if let Some((started, rx)) = &ra.running {
            match rx.try_recv() {
                Ok(res) => {
                    if res.is_ok() {
                        self.status = format!("Rogue scores from {} trees in {:.1}s", ra.trees.len() - ra.skip, started.elapsed().as_secs_f32());
                    }
                    ra.result = Some(res);
                    ra.running = None;
                }
                Err(std::sync::mpsc::TryRecvError::Empty) => ctx.request_repaint_after(std::time::Duration::from_millis(250)),
                Err(std::sync::mpsc::TryRecvError::Disconnected) => {
                    ra.result = Some(Err("the analysis stopped unexpectedly".into()));
                    ra.running = None;
                }
            }
        }
        let scores: Vec<(std::collections::BTreeSet<String>, f64)> = match &ra.result {
            Some(Ok((_, _, s, _))) if color => s.iter().map(|s| (s.taxa.iter().cloned().collect(), s.score)).collect(),
            _ => Vec::new(),
        };
        if !open {
            self.rogue = None;
        }
        if color {
            if let Some(d) = self.doc() {
                d.checkpoint();
                let scores: std::collections::HashMap<std::collections::BTreeSet<String>, f64> = scores.into_iter().collect();
                for n in d.tree.preorder() {
                    let taxa: std::collections::BTreeSet<String> = d.tree.tips_below(n).into_iter().map(|t| d.tree.label(t).to_string()).collect();
                    match scores.get(&taxa) {
                        Some(&s) => d.tree.nodes[n].attrs.insert("rogue_score".into(), crate::tree::Attr::Num(s)),
                        None => d.tree.nodes[n].attrs.remove("rogue_score"),
                    };
                }
                for e in d.layers.iter_mut() {
                    if let Layer::Tree(s) = &mut e.layer {
                        s.color_by = Some("rogue_score".into());
                        s.palette = Palette::BlackRed;
                        s.width = s.width.max(1.5);
                    }
                }
            }
        }
    }

    /// Lineages-through-time plot of the active tree (and, for a tree set, the
    /// other trees in grey behind it).
    /// Tanglegram / mirror tree window: the active tree (left) facing another open tree.
    fn tangle_window(&mut self, ctx: &egui::Context) {
        let Some(w) = self.tangle.as_mut() else { return };
        if self.docs.is_empty() {
            self.tangle = None;
            return;
        }
        w.left = w.left.min(self.docs.len() - 1);
        w.right = w.right.min(self.docs.len() - 1);
        let (li, ri) = (w.left, w.right);
        let (ld, rd) = (&self.docs[li], &self.docs[ri]);
        let names: Vec<String> = self.docs.iter().map(|d| d.name.clone()).collect();
        let matched = crate::tangle::matches(&ld.tree, &rd.tree).len();
        let crossings = crate::tangle::crossings(&ld.tree, &rd.tree);
        let mut open = true;
        let (mut untangle, mut export, mut swap) = (false, false, false);
        let (mut copy_right, mut edit): (bool, Option<usize>) = (false, None);
        let (textures, pics) = (&mut self.textures, &self.pics);
        egui::Window::new("Mirror tree / Tanglegram").open(&mut open).default_width(820.0).default_height(640.0).show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                for (label, idx, salt) in [("Left", &mut w.left, "tg-left"), ("Right", &mut w.right, "tg-right")] {
                    ui.label(label);
                    egui::ComboBox::from_id_salt(salt).selected_text(names[*idx].as_str()).show_ui(ui, |ui| {
                        for (i, n) in names.iter().enumerate() {
                            ui.selectable_value(idx, i, n);
                        }
                    });
                }
                swap = ui.button("⇄").on_hover_text("Swap the two trees").clicked();
                if ui.button("Edit left").on_hover_text("Switch to the left tree's tab: its layers, colors, labels and edits show here as you make them").clicked() {
                    edit = Some(li);
                }
                if ui.button("Edit right").on_hover_text("Switch to the right tree's tab: its layers, colors, labels and edits show here as you make them").clicked() {
                    edit = Some(ri);
                }
            });
            ui.horizontal_wrapped(|ui| {
                ui.checkbox(&mut w.style.lines, "Connecting lines").on_hover_text("Join tips with the same name; off draws a mirror tree");
                if w.style.lines {
                    panels::color_edit(ui, "color", &mut w.style.line_color);
                    ui.add(egui::DragValue::new(&mut w.style.line_width).range(0.1..=6.0).speed(0.05).prefix("width "));
                    ui.checkbox(&mut w.style.clade_colors, "Clade colors").on_hover_text("Lines take the left tree's clade colors where set");
                }
                ui.add(egui::DragValue::new(&mut w.style.gap).range(0.0..=400.0).speed(1.0).prefix("gap ").suffix(" pt"));
            });
            ui.horizontal_wrapped(|ui| {
                ui.label(format!("{} matching tips, {}", matched, count(crossings, "crossing", "crossings")));
                untangle = ui
                    .add_enabled(li != ri && crossings > 0, egui::Button::new("Untangle"))
                    .on_hover_text("Rotate nodes in both trees to reduce crossing lines (undoable in each tab)")
                    .clicked();
                ui.separator();
                ui.add(egui::DragValue::new(&mut w.width_in).range(1.0..=40.0).speed(0.1).suffix(" in"));
                ui.label("×");
                ui.add(egui::DragValue::new(&mut w.height_in).range(1.0..=40.0).speed(0.1).suffix(" in"));
                export = ui.button("Export…").on_hover_text("PNG or TIFF (300 dpi), SVG or PDF").clicked();
            });
            if li == ri {
                ui.horizontal_wrapped(|ui| {
                    ui.colored_label(Color32::from_rgb(198, 120, 0), "Both sides show the same tab, so they look the same.");
                    copy_right = ui.button("Make the right side a separate copy").on_hover_text("Copy this tree, its layers and data into a new tab for the right side, to style or edit on its own").clicked();
                });
            }
            ui.small("Each side is drawn from its own tab: switch to it (Edit left / Edit right) to add or remove tips, label nodes, color branches by a trait or change any layer, and the figure updates here. Tips face the middle; the left tab's background and title are used.");
            let size = ui.available_size().max(egui::vec2(300.0, 200.0));
            // Preview at the export's aspect ratio, as large as fits.
            let aspect = w.width_in / w.height_in;
            let size = if size.x / size.y > aspect { egui::vec2(size.y * aspect, size.y) } else { egui::vec2(size.x, size.x / aspect) };
            let (resp, painter) = ui.allocate_painter(size, egui::Sense::hover());
            painter.rect_filled(resp.rect, 0.0, ld.view.background.map(canvas::c32).unwrap_or(Color32::WHITE));
            let env = canvas::EguiEnv { ctx: ui.ctx(), pics };
            let l = crate::scene::BuildInput { tree: &ld.tree, view: &ld.view, layers: &ld.layers, overlay: &[] };
            let r = crate::scene::BuildInput { tree: &rd.tree, view: &rd.view, layers: &rd.layers, overlay: &[] };
            let scene = crate::tangle::build(&l, &r, &w.style, &env, size.x, size.y, 1.0);
            canvas::paint_scene(&painter, &scene, resp.rect.min, textures, pics);
        });
        if copy_right {
            let copy = self.docs[ri].copy_as(&format!("{} (mirror)", self.docs[ri].name));
            self.add_doc(copy);
            if let Some(w) = self.tangle.as_mut() {
                w.right = self.docs.len() - 1;
            }
            self.status = "The right side is now its own tab; style or edit it there.".into();
        }
        if let Some(i) = edit {
            self.active = i;
        }
        let Some(w) = self.tangle.as_mut() else { return };
        if swap {
            std::mem::swap(&mut w.left, &mut w.right);
        }
        if untangle && li != ri {
            let (a, b) = if li < ri {
                let (x, y) = self.docs.split_at_mut(ri);
                (&mut x[li], &mut y[0])
            } else {
                let (x, y) = self.docs.split_at_mut(li);
                (&mut y[0], &mut x[ri])
            };
            let (mut lt, mut rt) = (a.tree.clone(), b.tree.clone());
            let (before, after) = crate::tangle::untangle(&mut lt, &mut rt);
            a.checkpoint();
            b.checkpoint();
            a.tree = lt;
            b.tree = rt;
            self.status = format!("Untangled: {} → {}", count(before, "crossing", "crossings"), after);
        }
        if export {
            self.export_tangle();
        }
        if !open {
            self.tangle = None;
        }
    }

    /// Save the tanglegram at the window's size: PNG or TIFF at 300 dpi, SVG or PDF.
    fn export_tangle(&mut self) {
        let Some(w) = self.tangle.as_ref() else { return };
        let Some(path) = rfd::FileDialog::new().add_filter("PNG image", &["png"]).add_filter("TIFF image", &["tif", "tiff"]).add_filter("SVG", &["svg"]).add_filter("PDF", &["pdf"]).set_file_name("tanglegram.png").save_file() else { return };
        let path = if path.extension().is_none() { path.with_extension("png") } else { path };
        let format = match path.extension().and_then(|e| e.to_str()).unwrap_or("png").to_ascii_lowercase().as_str() {
            "svg" => export::Format::Svg,
            "pdf" => export::Format::Pdf,
            "tif" | "tiff" => export::Format::Tiff,
            _ => export::Format::Png,
        };
        let dpi = 300.0;
        let job = export::ExportJob { path: path.clone(), format, width_px: (w.width_in * dpi).round() as u32, height_px: (w.height_in * dpi).round() as u32, dpi, transparent: false, tiff_rgb: true, fps: 4.0, every: 1 };
        let (ld, rd) = (&self.docs[w.left], &self.docs[w.right]);
        let l = crate::scene::BuildInput { tree: &ld.tree, view: &ld.view, layers: &ld.layers, overlay: &[] };
        let r = crate::scene::BuildInput { tree: &rd.tree, view: &rd.view, layers: &rd.layers, overlay: &[] };
        let images = self.pics.images();
        let res = export::export_scene(&|env, wd, ht, pt| crate::tangle::build(&l, &r, &w.style, env, wd, ht, pt), ld.view.background, &self.raster_fonts, &self.fonts, &images, &job);
        self.status = match res {
            Ok(()) => format!("Exported {}", path.display()),
            Err(e) => format!("Export failed: {:#}", e),
        };
    }

    /// Tree traversal window: the active tree with pre- or post-order visits
    /// played step by step, and exported as an animated GIF.
    fn traversal_window(&mut self, ctx: &egui::Context) {
        use crate::traversal::{overlay, visits, Order};
        let Some(w) = self.traversal.as_mut() else { return };
        let Some(doc) = self.docs.get(self.active) else {
            self.traversal = None;
            return;
        };
        let mut open = true;
        let mut export_gif = false;
        let (textures, pics) = (&mut self.textures, &self.pics);
        let now = ctx.input(|i| i.time);
        egui::Window::new("Tree traversal").open(&mut open).default_width(760.0).default_height(600.0).show(ctx, |ui| {
            let total = w.total;
            ui.horizontal_wrapped(|ui| {
                for o in Order::ALL {
                    if ui.selectable_label(w.order == o, o.name()).clicked() && w.order != o {
                        w.order = o;
                        w.step = 0;
                    }
                }
                ui.separator();
                if ui.button("⏮").on_hover_text("Start").clicked() {
                    w.step = 0;
                    w.playing = false;
                }
                if ui.button("◀").on_hover_text("Step back").clicked() {
                    w.step = w.step.saturating_sub(1);
                    w.playing = false;
                }
                if ui.button(if w.playing { "⏸" } else { "▶" }).on_hover_text(if w.playing { "Pause" } else { "Play" }).clicked() {
                    w.playing = !w.playing;
                    if w.playing && w.step >= total {
                        w.step = 0;
                    }
                    w.last = now;
                }
                if ui.button("▶|").on_hover_text("Step forward").clicked() {
                    w.step = (w.step + 1).min(total);
                    w.playing = false;
                }
                ui.add(egui::Slider::new(&mut w.speed, 0.25..=8.0).logarithmic(true).text("steps / s"));
                ui.checkbox(&mut w.numbers, "Visit numbers");
            });
            ui.horizontal_wrapped(|ui| {
                ui.label(format!("{}: {}.", w.order.name(), w.order.rule()));
                ui.separator();
                ui.label("GIF");
                ui.add(egui::DragValue::new(&mut w.width_in).range(1.0..=20.0).speed(0.1).suffix(" in"));
                ui.label("×");
                ui.add(egui::DragValue::new(&mut w.height_in).range(1.0..=20.0).speed(0.1).suffix(" in"));
                ui.add(egui::DragValue::new(&mut w.dpi).range(50.0..=300.0).speed(1.0).suffix(" dpi"));
                export_gif = ui.button("Export GIF…").on_hover_text("One frame per visit at the speed set above, starting empty and holding the finished traversal").clicked();
            });
            ui.small("Uses the active tab's tree, layers and styles; a collapsed clade counts as one node.");
            let avail = ui.available_size().max(egui::vec2(300.0, 200.0));
            let aspect = w.width_in / w.height_in;
            let size = if avail.x / avail.y > aspect { egui::vec2(avail.y * aspect, avail.y) } else { egui::vec2(avail.x, avail.x / aspect) };
            let (resp, painter) = ui.allocate_painter(size, egui::Sense::hover());
            painter.rect_filled(resp.rect, 0.0, doc.view.background.map(canvas::c32).unwrap_or(Color32::WHITE));
            let env = canvas::EguiEnv { ctx: ui.ctx(), pics };
            let input = crate::scene::BuildInput { tree: &doc.tree, view: &doc.view, layers: &doc.layers, overlay: &[] };
            let band = crate::traversal::caption_band(1.0);
            let mut scene = crate::scene::build(&input, &env, size.x, (size.y - band).max(50.0), 1.0);
            crate::traversal::make_room(&mut scene, band);
            let v = visits(&doc.tree, &scene, w.order);
            w.total = v.len();
            w.step = w.step.min(w.total);
            overlay(&mut scene, &v, w.step, w.order, w.numbers, 1.0);
            canvas::paint_scene(&painter, &scene, resp.rect.min, textures, pics);
        });
        // Advance while playing, stopping at the end.
        if w.playing {
            if now - w.last >= 1.0 / w.speed as f64 {
                w.last = now;
                w.step += 1;
                if w.step >= w.total {
                    w.step = w.total;
                    w.playing = false;
                }
            }
            ctx.request_repaint_after(std::time::Duration::from_millis(30));
        }
        if export_gif {
            self.export_traversal_gif();
        }
        if !open {
            self.traversal = None;
        }
    }

    /// Save the traversal as an animated GIF: an empty frame, one per visit, and
    /// the finished traversal held for two more.
    fn export_traversal_gif(&mut self) {
        use crate::traversal::{overlay, visits};
        let (Some(w), Some(doc)) = (self.traversal.as_ref(), self.docs.get(self.active)) else { return };
        let stem = format!("{}_{}", doc.name.split('.').next().unwrap_or("tree"), w.order.slug());
        let Some(path) = rfd::FileDialog::new().add_filter("Animated GIF", &["gif"]).set_file_name(format!("{}.gif", stem)).save_file() else { return };
        let path = path.with_extension("gif");
        let job = export::ExportJob { path: path.clone(), format: export::Format::Gif, width_px: (w.width_in * w.dpi).round() as u32, height_px: (w.height_in * w.dpi).round() as u32, dpi: w.dpi, transparent: false, tiff_rgb: true, fps: w.speed, every: 1 };
        let input = crate::scene::BuildInput { tree: &doc.tree, view: &doc.view, layers: &doc.layers, overlay: &[] };
        let (order, numbers, n) = (w.order, w.numbers, w.total);
        let images = self.pics.images();
        let res = export::export_animation(
            n + 3,
            &|env, i, wd, ht, pt| {
                let band = crate::traversal::caption_band(pt);
                let mut sc = crate::scene::build(&input, env, wd, (ht - band).max(50.0), pt);
                crate::traversal::make_room(&mut sc, band);
                let v = visits(&doc.tree, &sc, order);
                overlay(&mut sc, &v, i.min(v.len()), order, numbers, pt);
                sc
            },
            doc.view.background,
            &self.raster_fonts,
            &images,
            &job,
        );
        self.status = match res {
            Ok(()) => format!("Exported {} ({} frames)", path.display(), n + 3),
            Err(e) => format!("Export failed: {:#}", e),
        };
    }

    fn ltt_window(&mut self, ctx: &egui::Context) {
        let Some(w) = self.ltt.as_mut() else { return };
        let Some(doc) = self.docs.get_mut(self.active) else { return };
        // A tree set's clade summary, for its consensus tree (large sets: in the background).
        let summary = doc.posterior.as_mut().and_then(|p| {
            p.poll_summary();
            p.summary()
        });
        let pending = doc.posterior.as_ref().map(|p| p.summary_pending()).unwrap_or(false);
        let doc = &*doc;
        let (young, mpu) = time_frame(doc.layers.iter().filter(|e| e.enabled).map(|e| &e.layer), &doc.tree);
        let mpu = if mpu > 0.0 { mpu } else { 1.0 };
        // Curves as (age, lineages) steps, the oldest first.
        let to_ages = |steps: Vec<(f64, usize)>| -> Vec<(f64, usize)> {
            let deepest = steps.last().map(|s| s.0).unwrap_or(0.0);
            steps.into_iter().map(|(x, k)| (young + (deepest - x) * mpu, k)).collect()
        };
        // The trees behind the solid curve: a tree set after burn-in, or the sample a
        // consensus or MCC tab summarizes.
        let sample: Option<(String, &[Tree])> = match (&doc.posterior, &doc.source_sample) {
            (Some(p), _) => Some((format!("{}|{}|{}", doc.name, p.trees.len(), p.burnin), p.kept())),
            (None, Some((t, skip))) => Some((format!("{}|source|{}|{}", doc.name, t.len(), skip), &t[(*skip).min(t.len())..])),
            _ => None,
        };
        let key = sample.as_ref().map(|s| s.0.clone());
        // Solid: a tree set's majority-rule consensus with mean node heights, else the tree shown.
        let (main, main_label) = match (&doc.posterior, &summary, &key) {
            (Some(p), Some(s), Some(k)) => {
                let k = format!("{}|{}", k, if w.mcc { "mcc" } else { "consensus" });
                if w.main_curve.as_ref().map(|c| &c.0) != Some(&k) {
                    let t = if w.mcc { consensus::mcc_tree(p.kept(), s, consensus::HeightMode::Mean).ok() } else { Some(consensus::consensus(s, 0.5, consensus::HeightMode::Mean)) };
                    w.main_curve = Some((k.clone(), t.unwrap_or_else(|| doc.tree.clone())));
                }
                let label = if w.mcc { "maximum clade credibility (MCC) tree of the set, mean node heights" } else { "majority-rule consensus tree of the set, mean node heights" };
                (to_ages(w.main_curve.as_ref().map(|c| crate::ltt::ltt(&c.1)).unwrap_or_default()), label)
            }
            (Some(_), _, _) => (to_ages(crate::ltt::ltt(&doc.tree)), "the tree shown (the summary tree is being computed)"),
            _ => (to_ages(crate::ltt::ltt(&doc.tree)), "this tree"),
        };
        // All the trees (evenly spaced, at most 1,000), cached.
        if w.show_set {
            if let Some((k, trees)) = &sample {
                if w.set_curves.as_ref().map(|c| &c.0) != Some(k) {
                    let step = trees.len().div_ceil(1000).max(1);
                    let curves = trees.iter().step_by(step).map(|t| to_ages(crate::ltt::ltt(t))).collect();
                    w.set_curves = Some((k.clone(), curves));
                }
            }
        }
        let n_sample = sample.as_ref().map(|s| s.1.len()).unwrap_or(0);
        let others: &[Vec<(f64, usize)>] = match (&w.set_curves, w.show_set, &key) {
            (Some((k, c)), true, Some(cur)) if k == cur => c,
            _ => &[],
        };
        let mut open = true;
        let table = |c: &[(f64, usize)]| {
            let header = vec!["Time before present".to_string(), "Lineages".to_string()];
            let rows: Vec<Vec<String>> = c.iter().map(|(a, k)| vec![crate::tree::format_num(*a, 6), k.to_string()]).collect();
            (header, rows)
        };
        // Timescale rows: those of the tab's timescale layer, else epochs and periods.
        let geo: Vec<crate::geotime::Level> = if w.geoscale {
            use crate::geotime::Level;
            let rows = doc.layers.iter().filter(|e| e.enabled).find_map(|e| if let Layer::Geoscale(g) = &e.layer { Some([(g.ages, Level::Age), (g.epochs, Level::Epoch), (g.periods, Level::Period), (g.eras, Level::Era)].iter().filter(|r| r.0).map(|r| r.1).collect::<Vec<_>>()) } else { None });
            rows.filter(|r| !r.is_empty()).unwrap_or_else(|| vec![Level::Epoch, Level::Period])
        } else {
            Vec::new()
        };
        let mut save_clicked = false;
        #[allow(clippy::type_complexity)]
        let mut save: Option<(String, Vec<(f64, usize)>, Vec<Vec<(f64, usize)>>, bool, bool, Vec<crate::geotime::Level>, Option<(Tree, f64, f64)>)> = None;
        let (textures, pics) = (&mut self.textures, &self.pics);
        let mut toggled = false;
        egui::Window::new(format!("Lineages through time: {}", doc.name)).open(&mut open).default_width(560.0).show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.checkbox(&mut w.log, "Log scale").on_hover_text("Lineages on a log scale: constant-rate pure-birth diversification is a straight line");
                ui.checkbox(&mut w.smooth, "Smooth").on_hover_text("A smooth curve through the steps (shape-preserving, no overshoot) instead of a step plot");
                ui.checkbox(&mut w.geoscale, "Geologic timescale").on_hover_text("Eras, periods, epochs or ages under the time axis (the rows of the tab's timescale layer, else epochs and periods); times are read as Ma");
                ui.checkbox(&mut w.show_tree, "Show tree").on_hover_text("Draw the solid curve's tree (consensus, MCC or this tree) above the plot on the same time axis, with dashed lines from each node down to the curve");
                if doc.posterior.is_some() {
                    ui.label("Solid:");
                    let was = w.mcc;
                    ui.selectable_value(&mut w.mcc, false, "Consensus").on_hover_text("Majority-rule consensus tree of the set, with mean node heights");
                    ui.selectable_value(&mut w.mcc, true, "MCC").on_hover_text("Maximum clade credibility tree: the tree in the set whose clades have the highest product of frequencies, with mean node heights");
                    toggled |= w.mcc != was;
                }
                if n_sample > 0 {
                    toggled |= ui.checkbox(&mut w.show_set, "Trees in the set").on_hover_text("Every tree of the set (after burn-in; at most 1,000, evenly spaced), faint, behind the solid curve").changed();
                }
                if ui.button("Copy data").on_hover_text("Tab-separated: time before present and number of lineages at each step").clicked() {
                    let (h, rows) = table(&main);
                    ui.ctx().copy_text(io::data::table_to_text(&h, &rows, b'\t'));
                }
                save_clicked = ui.button("Save plot…").on_hover_text("PNG (300 dpi), SVG or PDF, 6 × 4 in, titled with the tab's name").clicked();
                if ui.button("Save CSV…").clicked() {
                    if let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).set_file_name("ltt.csv").save_file() {
                        let (h, rows) = table(&main);
                        let _ = std::fs::write(path, io::data::table_to_text(&h, &rows, b','));
                    }
                }
            });
            let title = format!("Lineages through time: {}", doc.name);
            let size = egui::vec2(ui.available_width().max(300.0), if w.show_tree { 520.0 } else { 320.0 });
            let (resp, painter) = ui.allocate_painter(size, egui::Sense::hover());
            painter.rect_filled(resp.rect, 0.0, Color32::WHITE);
            let solid: &Tree = if doc.posterior.is_some() { w.main_curve.as_ref().map(|c| &c.1).unwrap_or(&doc.tree) } else { &doc.tree };
            let opt = crate::plots::LttPlot { log: w.log, smooth: w.smooth, geo: &geo, title: "", tree: w.show_tree.then_some((solid, young, mpu)) };
            let scene = crate::plots::ltt(&main, others, &opt, size.x, size.y, 1.0);
            canvas::paint_scene(&painter, &scene, resp.rect.min, textures, pics);
            if save_clicked {
                save = Some((title, main.clone(), others.to_vec(), w.log, w.smooth, geo.clone(), w.show_tree.then(|| (solid.clone(), young, mpu))));
            }
            if others.is_empty() {
                ui.small(format!("Solid: {}.", main_label));
            } else {
                let shown = if others.len() < n_sample { format!("{} of the {} trees", others.len(), n_sample) } else { format!("all {} trees", n_sample) };
                ui.small(format!("Solid: {}. Faint: {} of the set.", main_label, shown));
            }
            ui.small("Number of lineages at each time, from the root to the youngest tip (Nee et al. 1992). Ages use the timescale or time axis settings, or the tree's own node heights.");
        });
        if let Some((title, main, others, log, smooth, geo, tree)) = save {
            let opt = crate::plots::LttPlot { log, smooth, geo: &geo, title: &title, tree: tree.as_ref().map(|(t, y, m)| (t, *y, *m)) };
            self.save_plot(&|wd, ht, pt| crate::plots::ltt(&main, &others, &opt, wd, ht, pt), "ltt.png");
        }
        if toggled {
            ctx.request_repaint();
        }
        if pending {
            ctx.request_repaint_after(std::time::Duration::from_millis(200));
        }
        if !open {
            self.ltt = None;
        }
    }

    /// Open `trees` (one tree or a set) in a new tab named `name`, keeping the
    /// active tab's layout and comparative data (without rows for `removed`).
    fn open_derived(&mut self, name: &str, trees: Vec<Tree>, removed: &[String]) {
        let (layout, data) = match self.docs.get(self.active) {
            Some(d) => (d.view.layout.clone(), d.data.as_ref().map(|t| t.without(removed))),
            None => return,
        };
        self.open_trees(name, trees, None);
        if let Some(d) = self.docs.get_mut(self.active) {
            d.view.layout = layout;
            if let Some(t) = data {
                t.join_to_tree(&mut d.tree);
                d.data = Some(t);
            }
        }
    }

    /// Thin a tree set, drop random taxa, or drop tips with short branches.
    fn thinning_dialog(&mut self, ctx: &egui::Context) {
        let Some(doc) = self.docs.get(self.active) else {
            self.dialog = None;
            return;
        };
        let n_tips = doc.tree.num_tips();
        let n_trees = doc.posterior.as_ref().map(|p| p.kept().len());
        let (mut ok, mut cancel) = (false, false);
        match &mut self.dialog {
            Some(Dialog::ThinTrees { percent }) => {
                let n = n_trees.unwrap_or(0);
                let k = ((n as f32 * *percent / 100.0).round() as usize).clamp(1, n.max(1));
                egui::Window::new("Thin trees").collapsible(false).resizable(false).show(ctx, |ui| {
                    ui.add(egui::Slider::new(percent, 1.0..=100.0).suffix("%").text("of the trees kept"));
                    ui.label(format!("Keeps {} of the {} trees after burn-in, chosen at random (in their original order), in a new tab.", k, n));
                    ui.horizontal(|ui| {
                        ok = ui.button("Thin").clicked();
                        cancel = ui.button("Cancel").clicked();
                    });
                });
            }
            Some(Dialog::ThinTaxa { percent, clade }) => {
                let pool = clade.map(|c| doc.tree.tips_below(c).len()).unwrap_or(n_tips);
                let k = thin_count(pool, n_tips, *percent, clade.is_some());
                egui::Window::new(if clade.is_some() { "Thin taxa in clade" } else { "Thin taxa" }).collapsible(false).resizable(false).show(ctx, |ui| {
                    ui.add(egui::Slider::new(percent, 1.0..=99.0).suffix("%").text(if clade.is_some() { "of the clade's taxa removed" } else { "of the taxa removed" }));
                    let set = if n_trees.is_some() { " and every tree in the set" } else { "" };
                    let keep = if clade.is_some() { "At least one taxon of the clade stays." } else { "At least three taxa stay." };
                    let of = if clade.is_some() { format!("the clade's {} taxa", pool) } else { format!("the {} taxa", pool) };
                    ui.label(format!("Removes {} of {}, chosen at random, from the tree{} and from the comparative data, in a new tab. {}", k, of, set, keep));
                    ui.horizontal(|ui| {
                        ok = ui.button("Thin").clicked();
                        cancel = ui.button("Cancel").clicked();
                    });
                });
            }
            Some(Dialog::DropShort { min }) => {
                let short = ops::short_tips(&doc.tree, *min);
                egui::Window::new("Drop tips with short branches").collapsible(false).resizable(false).show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        ui.label("Terminal branches shorter than");
                        ui.add(egui::DragValue::new(min).speed(0.01).range(0.0..=f64::MAX));
                    });
                    let names: Vec<String> = short.iter().take(12).map(|&n| doc.tree.label(n).replace('_', " ")).collect();
                    let more = if short.len() > 12 { format!(" and {} more", short.len() - 12) } else { String::new() };
                    ui.label(format!("{} of {} tips{}{}", short.len(), n_tips, if names.is_empty() { String::new() } else { format!(": {}", names.join(", ")) }, more));
                    ui.horizontal(|ui| {
                        ok = ui.add_enabled(!short.is_empty(), egui::Button::new("Drop")).clicked();
                        cancel = ui.button("Cancel").clicked();
                    });
                });
            }
            _ => return,
        }
        if !(ok || cancel) {
            return;
        }
        let dialog = self.dialog.take();
        if cancel {
            return;
        }
        let seed = ops::random_seed();
        let name = doc.name.clone();
        match dialog {
            Some(Dialog::ThinTrees { percent }) => {
                let Some(p) = doc.posterior.as_ref() else { return };
                let kept = p.kept();
                let k = ((kept.len() as f32 * percent / 100.0).round() as usize).clamp(1, kept.len());
                let trees: Vec<Tree> = ops::random_subset(kept.len(), k, seed).into_iter().map(|i| kept[i].clone()).collect();
                self.open_derived(&format!("{} – {}", name, count(k, "tree", "trees")), trees, &[]);
            }
            Some(Dialog::ThinTaxa { percent, clade }) => {
                let tips = match clade {
                    Some(c) => doc.tree.tips_below(c),
                    None => doc.tree.tips(),
                };
                let k = thin_count(tips.len(), n_tips, percent, clade.is_some());
                let removed: Vec<String> = ops::random_subset(tips.len(), k, seed).into_iter().map(|i| doc.tree.label(tips[i]).to_string()).collect();
                let source: Vec<Tree> = match &doc.posterior {
                    Some(p) => p.trees.iter().cloned().collect(),
                    None => vec![doc.tree.clone()],
                };
                match ops::prune_labels(&source, &removed) {
                    Ok(trees) => {
                        self.open_derived(&format!("{} – {} removed", name, count(removed.len(), "taxon", "taxa")), trees, &removed);
                        self.status = format!("Removed {} at random: {}", count(removed.len(), "taxon", "taxa"), removed.iter().take(8).cloned().collect::<Vec<_>>().join(", "));
                    }
                    Err(e) => self.status = format!("{:#}", e),
                }
            }
            Some(Dialog::DropShort { min }) => {
                let short = ops::short_tips(&doc.tree, min);
                let k = short.len();
                self.apply(Action::Drop(short));
                if self.docs.get(self.active).map(|d| d.tree.num_tips()).unwrap_or(0) < n_tips {
                    self.status = format!("Dropped {} with terminal branches shorter than {}", count(k, "tip", "tips"), crate::tree::format_num(min, 4));
                }
            }
            _ => {}
        }
    }

    /// Spreadsheet-style editor for the comparative data. The taxon column is
    /// read-only, so rows keep matching the tree.
    fn data_editor_window(&mut self, ctx: &egui::Context) {
        if !self.data_editor {
            return;
        }
        let Some(doc) = self.docs.get_mut(self.active) else { return };
        let Some(data) = doc.data.as_mut() else {
            self.data_editor = false;
            return;
        };
        let mut open = true;
        let mut edited: Vec<String> = Vec::new();
        egui::Window::new(format!("Edit comparative data: {}", data.name)).open(&mut open).default_size([640.0, 440.0]).show(ctx, |ui| {
            // Each column wide enough for its header (and its sort arrow), at least 96 px.
            let font = egui::TextStyle::Body.resolve(ui.style());
            let widths: Vec<f32> = data.columns.iter().map(|c| ui.fonts(|f| f.layout_no_wrap(format!("{} ▲", c), font.clone(), Color32::BLACK).size().x + 12.0).max(96.0)).collect();
            ui.small("Click a cell to edit it; changes apply to the tree as you type. The taxon column cannot be edited.");
            egui::ScrollArea::horizontal().id_salt("data-edit-h").show(ui, |ui| {
                // Sorting changes only the order shown; edits go to the row itself.
                let tid = egui::Id::new("data-edit");
                ui.horizontal(|ui| {
                    for (i, c) in data.columns.iter().enumerate() {
                        table::header_sized(ui, tid, i, c, Some(widths[i]));
                    }
                });
                let shown = table::order(ui, tid, &data.rows);
                let row_h = ui.spacing().interact_size.y + 2.0;
                let n = data.rows.len();
                egui::ScrollArea::vertical().id_salt("data-edit-v").show_rows(ui, row_h, n, |ui, range| {
                    for r in range.map(|i| shown[i]) {
                        ui.horizontal(|ui| {
                            for c in 0..data.columns.len() {
                                if c == data.key {
                                    ui.add_sized([widths[c], row_h - 2.0], egui::Label::new(egui::RichText::new(data.rows[r][c].replace('_', " ")).italics()).truncate());
                                } else {
                                    let cell = &mut data.rows[r][c];
                                    if ui.add_sized([widths[c], row_h - 2.0], egui::TextEdit::singleline(cell).id_salt(("data-cell", r, c))).changed() {
                                        edited.push(data.columns[c].clone());
                                    }
                                }
                            }
                        });
                    }
                });
            });
        });
        if !edited.is_empty() {
            edited.dedup();
            doc.rejoin_data(&edited);
        }
        if !open {
            self.data_editor = false;
        }
    }

    /// Summary statistics of the comparative data, one row per column, and the
    /// distribution of one column. With nodes selected, only the taxa in those
    /// clades are summarized.
    fn data_summary_window(&mut self, ctx: &egui::Context) {
        let Some(SummaryView { tree_only, plotted, pir_cols, density, shade, opacity }) = self.data_summary.as_mut() else { return };
        let Some(doc) = self.docs.get(self.active) else { return };
        let Some(data) = doc.data.as_ref() else {
            self.data_summary = None;
            return;
        };
        let names = |nodes: Vec<NodeId>| -> std::collections::HashSet<String> { nodes.into_iter().map(|n| io::data::normalize_name(doc.tree.label(n))).collect() };
        // The taxa summarized: the selected clades, else the tree's (or all rows).
        let selected: Vec<NodeId> = doc.selection.iter().flat_map(|&n| doc.tree.tips_below(n)).collect();
        let scope_set = if !selected.is_empty() {
            Some(names(selected.clone()))
        } else if *tree_only {
            Some(names(doc.tree.tips()))
        } else {
            None
        };
        let scope_label = if !selected.is_empty() {
            format!("Selection: {}", count(selected.len(), "taxon", "taxa"))
        } else if *tree_only {
            format!("All {} in the tree", count(doc.tree.num_tips(), "taxon", "taxa"))
        } else {
            "All rows of the table".to_string()
        };
        let (header, rows) = data.summary(scope_set.as_ref());
        let cols: Vec<(usize, String, io::data::ColumnKind)> = (0..data.columns.len()).filter(|&c| c != data.key).map(|c| (c, data.columns[c].clone(), data.kind(c))).collect();
        if !cols.iter().any(|c| c.1 == *plotted) {
            if let Some(c) = cols.first() {
                *plotted = c.1.clone();
            }
        }
        let (numbers, states) = match cols.iter().find(|c| c.1 == *plotted) {
            Some((c, _, io::data::ColumnKind::Numeric)) => (Some(data.numbers(*c, scope_set.as_ref())), Vec::new()),
            Some((c, _, _)) => (None, data.state_counts(*c, scope_set.as_ref())),
            None => (None, Vec::new()),
        };
        // Phylogenetic imbalance ratio of the ticked categorical columns (the first two to start).
        let cat: Vec<(usize, String)> = cols.iter().filter(|c| c.2 != io::data::ColumnKind::Numeric).map(|c| (c.0, c.1.clone())).collect();
        let ticked = pir_cols.get_or_insert_with(|| cat.iter().take(2).map(|c| c.1.clone()).collect());
        let pir_idx: Vec<usize> = cat.iter().filter(|c| ticked.contains(&c.1)).map(|c| c.0).collect();
        let pir = crate::pir::pir(&doc.tree, data, &pir_idx, scope_set.as_ref());
        let title = format!("{} ({})", plotted, scope_label);
        let (show_density, fill) = (*density, if *shade { *opacity } else { 0.0 });
        let build = |w: f32, h: f32, pt: f32| crate::plots::distribution(numbers.as_deref(), &states, show_density, fill, &title, w, h, pt);
        let bar_width = crate::plots::distribution_width(numbers.is_some(), &states);
        let mut open = true;
        let mut save = false;
        let (textures, pics) = (&mut self.textures, &self.pics);
        egui::Window::new(format!("Summary statistics: {}", data.name)).open(&mut open).default_width(760.0).show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.strong(&scope_label);
                ui.add_enabled(selected.is_empty(), egui::Checkbox::new(tree_only, "Only taxa in the tree"))
                    .on_hover_text("Leave out rows for taxa that are not in the tree")
                    .on_disabled_hover_text("Showing the selected clades; clear the selection (Esc) to summarize the whole tree");
                if ui.button("Copy table").on_hover_text("Tab-separated").clicked() {
                    ui.ctx().copy_text(io::data::table_to_text(&header, &rows, b'\t'));
                }
                if ui.button("Save CSV…").clicked() {
                    if let Some(path) = rfd::FileDialog::new().add_filter("CSV", &["csv"]).set_file_name("summary_statistics.csv").save_file() {
                        let _ = std::fs::write(path, io::data::table_to_text(&header, &rows, b','));
                    }
                }
            });
            ui.small("Select nodes in the tree to summarize only their clades.");
            egui::ScrollArea::both().max_height(300.0).show(ui, |ui| {
                egui::Grid::new("data-summary").striped(true).show(ui, |ui| {
                    let tid = egui::Id::new("data-summary");
                    table::headers(ui, tid, &header);
                    ui.end_row();
                    for r in table::order(ui, tid, &rows).into_iter().map(|k| &rows[k]) {
                        for v in r {
                            ui.label(v);
                        }
                        ui.end_row();
                    }
                });
            });
            ui.small("SD is the sample standard deviation. Quantiles interpolate between values (R's default), so 2.5%–97.5% is the central 95% interval. Missing counts empty, NA, ? and - (and, in numeric columns, values that aren't numbers).");
            if !cat.is_empty() {
                ui.separator();
                ui.horizontal_wrapped(|ui| {
                    ui.strong("Phylogenetic imbalance ratio").on_hover_text(PIR_HOVER);
                    ui.label("of");
                    for (_, c) in &cat {
                        let mut on = ticked.contains(c);
                        if ui.checkbox(&mut on, c).changed() {
                            if on {
                                ticked.push(c.clone());
                            } else {
                                ticked.retain(|x| x != c);
                            }
                        }
                    }
                });
                let f = |x: f64| if x.is_finite() { crate::tree::format_num(x, 3) } else { "–".to_string() };
                match &pir {
                    Some(p) if p.steps == 0 => {
                        ui.weak(format!("No state changes on the tree among {} (all taxa share one state), so CI and PIR are undefined.", count(p.taxa, "taxon", "taxa")));
                    }
                    Some(p) => {
                        ui.horizontal_wrapped(|ui| {
                            ui.strong(format!("PIR = {}", f(p.pir)));
                            ui.label(format!("= NIR {} × CI {}", f(p.nir), f(p.ci)));
                            if p.pir < 0.1 {
                                ui.colored_label(Color32::from_rgb(40, 140, 60), "below 0.1, the recommended maximum");
                            } else {
                                ui.colored_label(Color32::from_rgb(198, 120, 0), "0.1 or more: few independent changes or unbalanced classes");
                            }
                        });
                        let classes: Vec<String> = p.classes.iter().map(|c| format!("{{{}}} {}", c.0, c.1)).collect();
                        ui.small(format!("{}; at least {} state changes, {} by parsimony on the tree. Classes: {}.", count(p.taxa, "taxon", "taxa"), p.min_steps, p.steps, classes.join(", ")));
                        ui.small("Normalized imbalance ratio (NIR) = (largest class − smallest class) / taxa");
                        ui.small("Consistency index (CI) = minimum / parsimony changes.");
                        ui.small("Gardner & Organ (2021) Syst. Biol. 70:1061–1075.");
                    }
                    None => {
                        ui.weak("Tick categorical columns that have data for taxa in the tree.");
                    }
                }
            }
            ui.separator();
            ui.horizontal(|ui| {
                ui.label("Distribution of");
                egui::ComboBox::from_id_salt("dist-col").selected_text(plotted.as_str()).show_ui(ui, |ui| {
                    for c in &cols {
                        ui.selectable_value(plotted, c.1.clone(), &c.1);
                    }
                });
                if numbers.is_some() {
                    ui.checkbox(density, "Density").on_hover_text("Show the kernel density curve over the histogram");
                    if *density {
                        ui.checkbox(shade, "Shade").on_hover_text("Fill the area under the curve");
                        if *shade {
                            ui.add(egui::Slider::new(opacity, 0.05..=1.0).text("opacity").fixed_decimals(2));
                        }
                    }
                }
                save = ui.button("Save plot…").on_hover_text("PNG (300 dpi), SVG or PDF, 6 × 4 in").clicked();
            });
            // The same scene as a saved plot, drawn at the window's size.
            // Bar charts keep their natural width; histograms fill the window.
            let wide = ui.available_width().max(300.0);
            let size = egui::vec2(bar_width.map(|b| b.min(wide)).unwrap_or(wide), 260.0);
            let (resp, painter) = ui.allocate_painter(size, egui::Sense::hover());
            painter.rect_filled(resp.rect, 0.0, Color32::WHITE);
            let scene = build(size.x, size.y, 1.0);
            canvas::paint_scene(&painter, &scene, resp.rect.min, textures, pics);
        });
        if save {
            self.save_plot_sized(&build, "distribution.png", bar_width.unwrap_or(432.0));
        }
        if !open {
            self.data_summary = None;
        }
    }

    /// Save a plot as PNG (300 dpi), SVG or PDF, 6 × 4 inches; `build` draws it at
    /// a given width, height (pixels) and points-to-pixels scale.
    fn save_plot(&mut self, build: &dyn Fn(f32, f32, f32) -> crate::scene::Scene, name: &str) {
        self.save_plot_sized(build, name, 432.0);
    }

    /// As `save_plot`, `width` points wide (4 in high).
    fn save_plot_sized(&mut self, build: &dyn Fn(f32, f32, f32) -> crate::scene::Scene, name: &str, width: f32) {
        let Some(mut path) = rfd::FileDialog::new().add_filter("PNG image", &["png"]).add_filter("SVG", &["svg"]).add_filter("PDF", &["pdf"]).set_file_name(name).save_file() else { return };
        if path.extension().is_none() {
            path.set_extension("png");
        }
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("png").to_ascii_lowercase();
        let images = crate::render::ImageMap::new();
        let bg = Some(Color::WHITE);
        let res: anyhow::Result<()> = (|| {
            match ext.as_str() {
                "svg" | "pdf" => {
                    // Vector output in points: `width` × 4 in.
                    let sc = build(width, 288.0, 1.0);
                    if ext == "svg" {
                        std::fs::write(&path, crate::render::svg::render(&sc, &images, bg, &self.fonts.name))?;
                    } else {
                        std::fs::write(&path, crate::render::pdf::render(&sc, &images, bg, &self.fonts)?)?;
                    }
                }
                _ => {
                    let dpi = 300.0;
                    let sc = build(width / 72.0 * dpi, 4.0 * dpi, dpi / 72.0);
                    let img = crate::render::raster::to_rgba(&crate::render::raster::render(&sc, &self.raster_fonts, &images, bg)?);
                    crate::render::raster::save_png(&img, &path, dpi)?;
                }
            }
            Ok(())
        })();
        self.status = match res {
            Ok(()) => format!("Saved {}", path.display()),
            Err(e) => format!("Could not save the plot: {:#}", e),
        };
    }

    fn run_delta_test(&mut self) {
        if let Some(d) = self.docs.get(self.active) {
            self.delta_report = Some(crate::nodedensity::delta_test(&d.tree).map_err(|e| format!("{:#}", e)));
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
                    let tid = egui::Id::new("node-report");
                    table::headers(ui, tid, &header);
                    ui.end_row();
                    for k in table::order(ui, tid, &rows) {
                        let (r, n) = (&rows[k], ids[k]);
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
                        let tid = egui::Id::new("regex-preview");
                        table::header_sized(ui, tid, 0, "Current", None);
                        ui.label("");
                        table::header_sized(ui, tid, 2, "New", None);
                        ui.end_row();
                        let keys: Vec<Vec<String>> = changes.iter().take(500).map(|(_, old, new)| vec![old.clone(), String::new(), new.clone()]).collect();
                        for (_, old, new) in table::order(ui, tid, &keys).into_iter().map(|k| &changes[k]) {
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

    /// Copy tip names (the selected tips, or all tips) to the clipboard, one per line.
    fn copy_tip_names(&mut self) {
        let Some(d) = self.docs.get(self.active) else { return };
        let tips = if d.selection.is_empty() { d.tree.tips() } else { d.selected_tips() };
        let names: Vec<&str> = tips.iter().map(|&t| d.tree.label(t)).collect();
        self.ctx.copy_text(names.join("\n"));
        self.status = format!("Copied {} species names{}", names.len(), if d.selection.is_empty() { "" } else { " (selection)" });
    }

    fn open_drop_missing(&mut self) {
        match self.docs.get(self.active).and_then(|d| d.data.as_ref()) {
            Some(data) => {
                let cols = data.trait_columns().into_iter().map(|(c, _)| (c, true)).collect();
                self.dialog = Some(Dialog::DropMissing { cols });
            }
            None => self.status = "Load a comparative data table first (Comparative data panel or drop a CSV/TSV).".into(),
        }
    }

    /// Pick trait columns; shows which species lack data in any of them, then drops or selects them.
    fn drop_missing_dialog(&mut self, ctx: &egui::Context) {
        let Some(Dialog::DropMissing { cols }) = &mut self.dialog else { return };
        let Some(doc) = self.docs.get(self.active) else { return };
        let chosen: Vec<String> = cols.iter().filter(|c| c.1).map(|c| c.0.clone()).collect();
        let missing = io::data::tips_missing(&doc.tree, &chosen);
        let total = doc.tree.num_tips();
        let (mut drop, mut select, mut cancel) = (false, false, false);
        egui::Window::new("Drop tips with missing data").collapsible(false).resizable(false).show(ctx, |ui| {
            ui.label("A species is dropped if it has no value in any ticked column (including species missing from the table).");
            ui.horizontal(|ui| {
                if ui.small_button("All").clicked() {
                    cols.iter_mut().for_each(|c| c.1 = true);
                }
                if ui.small_button("None").clicked() {
                    cols.iter_mut().for_each(|c| c.1 = false);
                }
            });
            egui::ScrollArea::vertical().id_salt("missing-cols").max_height(180.0).show(ui, |ui| {
                for (name, on) in cols.iter_mut() {
                    let k = io::data::tips_missing(&doc.tree, std::slice::from_ref(name)).len();
                    ui.checkbox(on, format!("{}  ({} missing)", name, k));
                }
            });
            ui.separator();
            ui.strong(format!("{} of {} species would be dropped", missing.len(), total));
            if !missing.is_empty() {
                egui::ScrollArea::vertical().id_salt("missing-tips").max_height(160.0).show(ui, |ui| {
                    for &n in &missing {
                        ui.label(egui::RichText::new(doc.tree.label(n).replace('_', " ")).italics());
                    }
                });
            }
            let keeps_enough = total - missing.len() >= 2;
            if !keeps_enough && !missing.is_empty() {
                ui.colored_label(Color32::from_rgb(198, 40, 40), "At least two species must remain.");
            }
            ui.horizontal(|ui| {
                let label = format!("Drop {} species", missing.len());
                if ui.add_enabled(!missing.is_empty() && keeps_enough, egui::Button::new(label)).clicked() {
                    drop = true;
                }
                if ui.add_enabled(!missing.is_empty(), egui::Button::new("Select them")).clicked() {
                    select = true;
                }
                if ui.button("Cancel").clicked() {
                    cancel = true;
                }
            });
        });
        if select {
            let d = &mut self.docs[self.active];
            d.selection = missing.iter().copied().collect();
            self.status = format!("Selected {} species with missing data", missing.len());
        }
        if drop {
            let k = missing.len();
            self.apply(Action::Drop(missing));
            self.status = format!("Dropped {} species with missing data (Ctrl+Z to undo)", k);
        }
        if drop || select || cancel {
            self.dialog = None;
        }
    }

    /// Choose the support attribute and threshold; shows how many nodes would go.
    fn collapse_dialog(&mut self, ctx: &egui::Context) {
        let Some(Dialog::Collapse { key, threshold }) = &mut self.dialog else { return };
        let Some(doc) = self.docs.get(self.active) else { return };
        let keys = ops::support_keys(&doc.tree);
        let (mut apply, mut cancel) = (false, false);
        egui::Window::new("Drop weakly supported nodes").collapsible(false).resizable(false).show(ctx, |ui| {
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
        let set_n = self.docs.get(self.active).and_then(|d| d.posterior.as_ref()).map(|p| p.trees.len());
        let (mut changed, mut apply, mut cancel) = (false, false, false);
        egui::Window::new("Transform branch lengths").collapsible(false).resizable(false).show(ctx, |ui| {
            ui.label(if whole { "Scope: whole tree".to_string() } else { format!("Scope: clade of {} tips", base.tips_below(*node).len()) });
            if let Some(k) = set_n {
                ui.weak(format!("On Apply, also every tree in the set ({}), on the clade with the same tips.", count(k, "tree", "trees")));
            }
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
                ops::BranchTransform::Scale => egui::Slider::new(value, 0.001..=1000.0).logarithmic(true),
            };
            changed |= ui.add(slider.text(if *kind == ops::BranchTransform::Scale { "factor" } else { "value" })).changed();
            if *kind == ops::BranchTransform::Scale {
                // Or give the height the clade should have; the factor follows.
                ui.horizontal(|ui| {
                    ui.label(format!("Height {} →", crate::tree::format_num(height, 4)));
                    let mut target = height * *value;
                    if ui.add(egui::DragValue::new(&mut target).speed(height * 0.01).range(0.0..=f64::MAX)).changed() && target > 0.0 {
                        *value = target / height;
                        changed = true;
                    }
                });
            }
            ui.small(kind.description());
            if !kind.reference().is_empty() {
                ui.add_space(4.0);
                ui.label(egui::RichText::new("Reference").strong().small());
                ui.small(kind.reference());
            }
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
            if let Some(Dialog::Transform { node, kind, value, base }) = &self.dialog {
                let (n, k, v) = (*node, *kind, *value);
                let what = format!("Applied {} = {}", k.name(), crate::tree::format_num(v, 3));
                self.status = what.clone();
                // A tree set: the same transform on every tree, on the clade with the
                // same tips (trees where it is not a clade are left as they were).
                let whole = n == base.root;
                if let Some(doc) = self.docs.get_mut(self.active) {
                    let base = base.clone();
                    let res = doc.apply_to_set(&|t| {
                        let root = if whole { Some(t.root) } else { ops::map_node(&base, n, t) };
                        root.map(|r| ops::transform_branches(t, r, k, v).is_ok()).unwrap_or(false)
                    });
                    if let Some((done, skipped)) = res {
                        doc.refresh_set();
                        self.status = set_note(&what, done, skipped);
                    }
                }
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
                        self.close_doc(self.active);
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
                if ui.button("Select all tips  (Ctrl+A)").clicked() {
                    if let Some(d) = self.doc() {
                        d.selection = d.tree.tips().into_iter().collect();
                    }
                    ui.close_menu();
                }
                if ui.button("Copy tip names  (Ctrl+C)").on_hover_text("The selected tips, or all tips if nothing is selected; one name per line").clicked() {
                    self.copy_tip_names();
                    ui.close_menu();
                }
                if ui.button("Clear selection (Esc)").clicked() {
                    if let Some(d) = self.doc() {
                        d.selection.clear();
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
            ui.menu_button("Tree", |ui| {
                let has_set_menu = self.docs.get(self.active).map(|d| d.posterior.is_some()).unwrap_or(false);
                if ui.button("Transform branch lengths (scale, λ, κ, δ, OU)…").clicked() {
                    if let Some(r) = self.docs.get(self.active).map(|d| d.tree.root) {
                        self.apply(Action::Transform(r));
                    }
                    ui.close_menu();
                }
                if ui.button("Import BayesTraits variable rates…").on_hover_text("Read a *.VarRates.txt file and summarize the rate scalar on each branch of this tree").clicked() {
                    ui.close_menu();
                    self.import_varrates();
                }
                if ui.button("Import BayesTraits Fabric…").on_hover_text("Read a Fabric-model *.VarRates.txt file and summarize LS_Beta and node scalars on each branch of this tree").clicked() {
                    ui.close_menu();
                    self.import_fabric();
                }
                if ui.button("Import fossil ranges…").on_hover_text("Read Species,minimum_age,maximum_age (comma- or tab-delimited) and draw each taxon's fossil range as a wide bar").clicked() {
                    ui.close_menu();
                    self.import_fossil_ranges();
                }
                ui.separator();
                if ui.button("MCMC diagnostics…").on_hover_text("Convergence and posterior plots of BayesTraits runs: load log files, or use the imported variable-rates or Fabric runs").clicked() {
                    // Opened with a tab for each set of imported runs (none: it offers to load logs).
                    let w = self.mcmc_win.get_or_insert_with(McmcWindow::new);
                    w.open = true;
                    for d in [self.varrates_diag.clone(), self.fabric_diag.clone()].into_iter().flatten() {
                        if !w.tabs.iter().any(|t| t.title == d.title) {
                            w.show(d);
                        }
                    }
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Tree traversal (pre-/post-order)…").on_hover_text("Watch pre-order or post-order traversal step through this tree, and export it as an animated GIF").clicked() {
                    self.traversal = Some(TraversalWindow { order: crate::traversal::Order::Post, step: 0, total: 0, playing: false, speed: 1.5, last: 0.0, numbers: true, width_in: 6.0, height_in: 4.5, dpi: 100.0 });
                    ui.close_menu();
                }
                if ui.add_enabled(has_set_menu, egui::Button::new("Thin trees…")).on_hover_text("Keep a random percentage of the trees in the set, in a new tab").on_disabled_hover_text("Open a set of trees first").clicked() {
                    self.dialog = Some(Dialog::ThinTrees { percent: 10.0 });
                    ui.close_menu();
                }
                if ui.button("Thin taxa…").on_hover_text("Remove a random percentage of taxa from the tree (or tree set) and the comparative data, in a new tab").clicked() {
                    self.dialog = Some(Dialog::ThinTaxa { percent: 10.0, clade: None });
                    ui.close_menu();
                }
                if let Some(r) = self.docs.get(self.active).map(|d| d.tree.root) {
                    let mut acts = Vec::new();
                    let mut collapse = false;
                    ui.menu_button("Polytomies", |ui| {
                        canvas::polytomy_menu(ui, r, &mut acts);
                        if ui.button("Drop weakly supported nodes…").clicked() {
                            collapse = true;
                            ui.close_menu();
                        }
                    });
                    if collapse {
                        if let Some(d) = self.docs.get(self.active) {
                            match ops::support_keys(&d.tree).first() {
                                Some(&key) => {
                                    let threshold = default_support_threshold(&d.tree, key);
                                    self.dialog = Some(Dialog::Collapse { key: key.to_string(), threshold });
                                }
                                None => self.status = "This tree has no support values (posterior, prob, bootstrap or numeric node labels).".into(),
                            }
                        }
                    }
                    for a in acts {
                        self.apply(a);
                    }
                }
                if ui.button("Drop tips with short branches…").on_hover_text("Drop every tip whose terminal branch is shorter than a value you choose").clicked() {
                    self.dialog = Some(Dialog::DropShort { min: 0.001 });
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Node report (support per node)…").on_hover_text("Table of every clade / bipartition with its support values, as from RAxML, IQ-TREE, MrBayes or BEAST").clicked() {
                    self.show_node_report = true;
                    ui.close_menu();
                }
                if ui.button("Lineages through time…").on_hover_text("Plot the number of lineages through time").clicked() {
                    self.ltt = Some(LttWindow { smooth: true, show_set: true, ..Default::default() });
                    ui.close_menu();
                }
                if ui.button("Node-density artifact (delta test)…").on_hover_text("Tests whether root-to-tip path lengths depend on the number of nodes along them, as the node-density artifact predicts").clicked() {
                    self.run_delta_test();
                    ui.close_menu();
                }
                let has_set = self.docs.get(self.active).map(|d| d.posterior.is_some()).unwrap_or(false);
                if ui
                    .add_enabled(has_set, egui::Button::new("Rogue taxa…"))
                    .on_hover_text("Flag the taxa whose unstable placement causes most of the uncertainty in a tree set, in a table of rogue scores (Klawitter, Bouckaert & Drummond 2024)")
                    .on_disabled_hover_text("Open a set of trees (posterior sample or bootstrap replicates) first")
                    .clicked()
                {
                    self.posterior_action(panels::PosteriorAction::Rogue);
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Mirror tree / Tanglegram…").on_hover_text("Draw this tree facing another open tree, with lines joining tips of the same name (or none, for a mirror tree)").clicked() {
                    // With one tab, the right side starts as a copy of it (a mirror tree to style on its own).
                    let left = self.active;
                    let right = if self.docs.len() > 1 {
                        (left + 1) % self.docs.len()
                    } else if let Some(d) = self.docs.get(left) {
                        let copy = d.copy_as(&format!("{} (mirror)", d.name));
                        self.add_doc(copy);
                        self.active = left;
                        self.docs.len() - 1
                    } else {
                        left
                    };
                    self.tangle = Some(TangleWindow { left, right, style: Default::default(), width_in: 9.0, height_in: 6.0 });
                    ui.close_menu();
                }
            });
            ui.menu_button("Taxa", |ui| {
                if ui.button("Select all species  (Ctrl+A)").clicked() {
                    if let Some(d) = self.doc() {
                        d.selection = d.tree.tips().into_iter().collect();
                    }
                    ui.close_menu();
                }
                if ui.button("Copy species names  (Ctrl+C)").on_hover_text("The selected species, or all species if nothing is selected; one name per line").clicked() {
                    self.copy_tip_names();
                    ui.close_menu();
                }
                let has_data = self.docs.get(self.active).map(|d| d.data.is_some()).unwrap_or(false);
                if ui.add_enabled(has_data, egui::Button::new("Drop tips with missing data…")).on_disabled_hover_text("Load a comparative data table first").clicked() {
                    self.open_drop_missing();
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Check taxon names (Open Tree of Life)…").on_hover_text("Find misspellings, synonyms and unknown names; optionally rename tips").clicked() {
                    self.start_taxon_check();
                    ui.close_menu();
                }
                ui.separator();
                ui.label(egui::RichText::new("PhyloPic silhouettes").small().weak());
                ui.checkbox(&mut self.settings.auto_phylopic, "Fetch PhyloPic silhouettes automatically")
                    .on_hover_text("Add a PhyloPic layer to each tree you open. Off: add one from Layers › Add or the button below");
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
            ui.menu_button("Help", |ui| {
                if ui.button("Open example trees").on_hover_text("A set of 10 trees of 8 ape genera").clicked() {
                    match io::read_trees_str(EXAMPLE) {
                        Ok(trees) => self.open_trees("Example: apes", trees, None),
                        Err(e) => self.fail("Could not open the example".into(), format!("{:#}", e)),
                    }
                    ui.close_menu();
                }
                ui.set_max_width(340.0);
                ui.add(egui::Label::new(egui::RichText::new(format!("Example data: {}", EXAMPLE_CITATION)).small()).wrap());
                ui.separator();

                if let Some(r) = &self.new_release {
                    if ui.button(egui::RichText::new(format!("Download Canopy {}…", r.version)).strong()).on_hover_text(r.url.as_str()).clicked() {
                        ui.ctx().open_url(egui::OpenUrl::new_tab(&r.url));
                        ui.close_menu();
                    }
                }
                if ui.checkbox(&mut self.settings.check_updates, "Check for new releases at startup").on_hover_text("Asks GitHub for the latest Canopy release when Canopy starts").changed() && self.settings.check_updates && self.new_release.is_none() {
                    self.start_update_check();
                }
                ui.separator();
                if ui.button("References…").on_hover_text("The software, methods, data sources and papers Canopy builds on").clicked() {
                    self.show_references = true;
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
                self.close_doc(i);
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
        self.delta_window(ctx);
        self.varrates_window(ctx);
        self.fabric_window(ctx);
        self.rogue_window(ctx);
        self.ltt_window(ctx);
        self.mcmc_window(ctx);
        self.poll_update_check();
        self.tangle_window(ctx);
        self.traversal_window(ctx);
        self.data_editor_window(ctx);
        self.data_summary_window(ctx);
        let mut done = false;
        if matches!(self.dialog, Some(Dialog::Transform { .. })) {
            self.transform_dialog(ctx);
        } else if matches!(self.dialog, Some(Dialog::Collapse { .. })) {
            self.collapse_dialog(ctx);
        } else if matches!(self.dialog, Some(Dialog::SaveTree(_))) {
            self.save_tree_dialog(ctx);
        } else if matches!(self.dialog, Some(Dialog::Regex(_))) {
            self.regex_dialog(ctx);
        } else if matches!(self.dialog, Some(Dialog::DropMissing { .. })) {
            self.drop_missing_dialog(ctx);
        } else if matches!(self.dialog, Some(Dialog::ThinTrees { .. } | Dialog::ThinTaxa { .. } | Dialog::DropShort { .. })) {
            self.thinning_dialog(ctx);
        } else if matches!(self.dialog, Some(Dialog::Silhouette { .. })) {
            self.silhouette_dialog(ctx);
        } else if matches!(self.dialog, Some(Dialog::FossilRange { .. })) {
            self.fossil_dialog(ctx);
        } else if let Some(dialog) = &mut self.dialog {
            let (title, node, text) = match dialog {
                Dialog::Rename { node, text } => ("Rename", *node, text),
                Dialog::CladeLabel { node, text } => ("Clade label", *node, text),
                Dialog::AddTaxon { node, mode, text } => (if *mode == ops::PasteMode::Child { "New taxon (child)" } else { "New taxon (sister)" }, *node, text),
                Dialog::Transform { .. } | Dialog::Collapse { .. } | Dialog::Regex(_) | Dialog::SaveTree(_) | Dialog::DropMissing { .. } | Dialog::FossilRange { .. } | Dialog::Silhouette { .. } | Dialog::ThinTrees { .. } | Dialog::ThinTaxa { .. } | Dialog::DropShort { .. } => unreachable!(),
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
                let add = match dialog {
                    Dialog::AddTaxon { mode, .. } => Some(*mode),
                    _ => None,
                };
                let is_label = matches!(dialog, Dialog::CladeLabel { .. });
                if let Some(doc) = self.docs.get_mut(self.active) {
                    if is_label {
                        doc.add_layer(Layer::CladeLabel(CladeLabelStyle::new(node, value)));
                    } else if let Some(mode) = add {
                        let name = if value.trim().is_empty() { "New_taxon".to_string() } else { value.trim().to_string() };
                        doc.checkpoint();
                        match ops::add_taxon(&mut doc.tree, node, &name, mode) {
                            Ok(s) => {
                                doc.selection.clear();
                                doc.selection.insert(s);
                                self.status = format!("Added {}", name);
                            }
                            Err(e) => {
                                doc.undo();
                                self.status = format!("{:#}", e);
                            }
                        }
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
        if self.show_references {
            let mut open = true;
            egui::Window::new("References").open(&mut open).default_width(640.0).default_height(560.0).show(ctx, |ui| {
                ui.horizontal(|ui| {
                    if ui.button("Copy all").on_hover_text("Plain text, to paste elsewhere").clicked() {
                        ui.ctx().copy_text(references::plain_text());
                    }
                });
                egui::ScrollArea::vertical().show(ui, |ui| references::show(ui));
            });
            self.show_references = open;
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
                ui.hyperlink_to(egui::RichText::new(format!("Canopy {}", env!("CARGO_PKG_VERSION"))).heading(), "https://github.com/ChrisOrgan/Canopy").on_hover_text("Canopy on GitHub");
                ui.label("Phylogenetic tree visualization inspired by ggtree, FigTree, Mesquite, and BayesTrees.");
                ui.add_space(6.0);
                ui.strong("Credits");
                ui.label("Created by Chris Organ and Claude Code using cleanroom coding, 2026. Thanks to Jacob Gardner, Jorge Avaria, Kevin Surya for input and advice.");
                ui.add_space(6.0);
                ui.strong("How to cite");
                ui.label(crate::citation::text());
                ui.horizontal(|ui| {
                    ui.hyperlink_to(format!("doi:{}", crate::citation::DOI), format!("https://doi.org/{}", crate::citation::DOI));
                    if ui.small_button("Copy citation").clicked() {
                        ui.ctx().copy_text(crate::citation::text());
                    }
                    if ui.small_button("Copy BibTeX").clicked() {
                        ui.ctx().copy_text(crate::citation::bibtex());
                    }
                });
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
        self.poll_job(ctx);
        // Clade summaries of large tree sets arrive from the background.
        for d in self.docs.iter_mut() {
            if d.poll_summary() {
                ctx.request_repaint();
            }
        }
        if self.docs.iter().any(|d| d.posterior.as_ref().map(|p| p.summary_pending()).unwrap_or(false)) {
            ctx.request_repaint_after(std::time::Duration::from_millis(200));
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

        let files = self.loaded_files();
        let mut unload_file = None;
        let mut actions = Vec::new();
        let mut load_data = false;
        let mut open_missing = false;
        if let Some(doc) = self.docs.get_mut(self.active) {
            egui::SidePanel::left("left").resizable(true).default_width(300.0).show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    actions.extend(panels::layout_panel(ui, doc));
                    ui.separator();
                    panels::zoom_panel(ui, doc);
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
                    ui.heading("Details");
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.small_button("✕").on_hover_text("Hide this panel (View › Show right panel brings it back)").clicked() {
                            hide_right = true;
                        }
                    });
                });
                egui::ScrollArea::both().auto_shrink([false, false]).show(ui, |ui| {
                    egui::CollapsingHeader::new("Tree Details").id_salt("tree-details").default_open(true).show(ui, |ui| {
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
                        ui.add_space(4.0);
                        ui.label(egui::RichText::new(format!("Loaded files ({})", files.len())).strong())
                            .on_hover_text("Support files held in memory. Closing a tab also unloads the BayesTraits results imported onto it.");
                        unload_file = Self::loaded_files_ui(ui, &files);
                    });
                    ui.separator();
                    actions.extend(panels::inspector(ui, doc, &mut self.search));
                    ui.separator();
                    load_data = panels::data_panel(ui, doc, &mut self.status, &mut open_missing, &mut self.data_editor, &mut self.data_summary);
                });
            });
            }
            if hide_right {
                self.settings.show_right_panel = false;
            }
        }
        match unload_file {
            Some((u, true)) => self.show_file(u),
            Some((u, false)) => self.unload(u),
            None => {}
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
        if open_missing {
            self.open_drop_missing();
        }

        self.dialogs(ctx);
        if let Some(doc) = self.docs.get(self.active) {
            let images = self.pics.images();
            if let Some(job) = self.export.show(ctx, doc, &self.raster_fonts, &images) {
                let saved = if job.format == export::Format::Gif {
                    export::export_gif(&mut self.docs[self.active], &self.raster_fonts, &images, &job).map(|n| format!("{} frames, {}×{} px", n, job.width_px, job.height_px))
                } else {
                    export::export(doc, &self.raster_fonts, &self.fonts, &images, &job).map(|_| format!("{}×{} px, {} dpi", job.width_px, job.height_px, job.dpi))
                };
                self.status = match saved {
                    Ok(what) => format!("Exported {} ({})", job.path.display(), what),
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

#[cfg(test)]
mod example_tests {
    #[test]
    fn thinning_counts() {
        // 50% of a 4-taxon clade in a 20-taxon tree: 2. A clade keeps one taxon.
        assert_eq!(super::thin_count(4, 20, 50.0, true), 2);
        assert_eq!(super::thin_count(4, 20, 99.0, true), 3);
        // The whole tree keeps at least three taxa.
        assert_eq!(super::thin_count(5, 5, 99.0, false), 2);
        assert_eq!(super::thin_count(10, 10, 10.0, false), 1);
    }

    #[test]
    fn bundled_example_is_a_tree_set() {
        let trees = crate::io::read_trees_str(super::EXAMPLE).unwrap();
        assert_eq!(trees.len(), 10);
        assert!(trees.iter().all(|t| t.num_tips() == 8 && t.find_label("Homo").is_some()));
    }
}

/// State of the summary-statistics window.
#[derive(Clone)]
pub struct SummaryView {
    /// Only rows of taxa in the tree (when nothing is selected).
    tree_only: bool,
    /// Column plotted.
    plotted: String,
    /// Categorical columns in the phylogenetic imbalance ratio (None: the first two).
    pir_cols: Option<Vec<String>>,
    /// Histograms: the density curve, and the area under it shaded at `opacity`.
    density: bool,
    shade: bool,
    opacity: f32,
}

/// Hover for the phylogenetic imbalance ratio.
const PIR_HOVER: &str = "How suitable categorical data are for phylogenetic comparative models: high when states changed few times on the tree and the state combinations are unevenly represented (Gardner & Organ 2021).";

impl Default for SummaryView {
    fn default() -> Self {
        SummaryView { tree_only: true, plotted: String::new(), pir_cols: None, density: true, shade: true, opacity: crate::plots::SHADE_OPACITY }
    }
}

/// References for the convergence measures: (measure, citation, DOI).
const CONVERGENCE_REFS: [(&str, &str, &str); 5] = [
    ("R-hat", "Gelman, A. & Rubin, D. B. (1992) Inference from iterative simulation using multiple sequences. Statistical Science 7: 457–472.", "10.1214/ss/1177011136"),
    ("R-hat, ESS", "Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. & Bürkner, P.-C. (2021) Rank-normalization, folding, and localization: an improved R-hat for assessing convergence of MCMC. Bayesian Analysis 16: 667–718.", "10.1214/20-BA1221"),
    ("ESS", "Geyer, C. J. (1992) Practical Markov chain Monte Carlo. Statistical Science 7: 473–483.", "10.1214/ss/1177011137"),
    ("ESS ≥ 200", "Rambaut, A., Drummond, A. J., Xie, D., Baele, G. & Suchard, M. A. (2018) Posterior summarization in Bayesian phylogenetics using Tracer 1.7. Systematic Biology 67: 901–904.", "10.1093/sysbio/syy032"),
    ("PP SD (ASDSF)", "Ronquist, F., Teslenko, M., van der Mark, P., et al. (2012) MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. Systematic Biology 61: 539–542.", "10.1093/sysbio/sys029"),
];

/// MCMC diagnostics window state: one tab per set of runs (log files, imported
/// VarRates runs), the active tab, the plotted parameter and plot kind.
struct McmcWindow {
    /// Shown on screen; closing the window hides it and keeps its tabs loaded.
    open: bool,
    tabs: Vec<crate::mcmc::Diagnostics>,
    active: usize,
    param: usize,
    kind: crate::plots::McmcPlot,
    /// Shade the area under density curves, at this opacity.
    shade: bool,
    opacity: f32,
}

impl McmcWindow {
    fn new() -> Self {
        McmcWindow { open: true, tabs: Vec::new(), active: 0, param: 0, kind: crate::plots::McmcPlot::Density, shade: true, opacity: crate::plots::SHADE_OPACITY }
    }

    /// Show `d` in its own tab, replacing a tab of the same name.
    fn show(&mut self, d: crate::mcmc::Diagnostics) {
        match self.tabs.iter().position(|t| t.title == d.title) {
            Some(i) => {
                self.tabs[i] = d;
                self.active = i;
            }
            None => {
                self.tabs.push(d);
                self.active = self.tabs.len() - 1;
            }
        }
        self.param = 0;
        self.open = true;
    }
}

/// Posterior probabilities of scalars compared across variable-rates or Fabric
/// runs: findings against the usual targets, then a per-branch table.
fn report_ui(ui: &mut egui::Ui, r: &crate::mcmc::Report, id: &str) {
    for (ok, line) in r.findings() {
        ui.horizontal_wrapped(|ui| {
            if ok {
                ui.colored_label(Color32::from_rgb(40, 140, 60), "✔");
            } else {
                ui.colored_label(Color32::from_rgb(198, 120, 0), "⚠");
            }
            ui.label(line);
        });
    }
    let (header, rows) = r.table();
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new("Scalars across runs").strong());
        if ui.small_button("Copy table").on_hover_text("Tab-separated").clicked() {
            ui.ctx().copy_text(io::data::table_to_text(&header, &rows, b'\t'));
        }
    });
    egui::ScrollArea::both().id_salt((id, "scroll")).max_height(200.0).show(ui, |ui| {
        egui::Grid::new((id, "grid")).striped(true).show(ui, |ui| {
            let tid = egui::Id::new((id, "table"));
            table::headers(ui, tid, &header);
            ui.end_row();
            for k in table::order(ui, tid, &rows) {
                for v in &rows[k] {
                    ui.label(v);
                }
                ui.end_row();
            }
        });
    });
}

/// A term, its explanation and the references behind it (indices into
/// `CONVERGENCE_REFS`), each with its DOI link.
fn explain_with_refs(ui: &mut egui::Ui, term: &str, text: &str, refs: &[usize]) {
    ui.horizontal_wrapped(|ui| {
        ui.spacing_mut().item_spacing.x = 0.0;
        ui.label(egui::RichText::new(format!("{}  ", term)).small().strong());
        ui.label(egui::RichText::new(text).small());
    });
    for &i in refs {
        let (_, cite, doi) = CONVERGENCE_REFS[i];
        ui.horizontal_wrapped(|ui| {
            ui.spacing_mut().item_spacing.x = 0.0;
            ui.add_space(14.0);
            ui.label(egui::RichText::new(format!("{} ", cite)).small().weak());
            ui.hyperlink_to(egui::RichText::new(format!("doi:{}", doi)).small(), format!("https://doi.org/{}", doi));
        });
    }
}

impl CanopyApp {
    /// Tree › MCMC diagnostics: BayesTraits logs or imported VarRates runs, one
    /// tab per set, with a table of every parameter (mean, median, 95% HPD, R-hat,
    /// ESS), density, box-and-whisker and trace plots run by run, and for
    /// variable-rates or Fabric runs the comparison of posterior probabilities.
    fn mcmc_window(&mut self, ctx: &egui::Context) {
        let Some(w) = self.mcmc_win.as_mut().filter(|w| w.open) else { return };
        let mut open = true;
        let (mut load_logs, mut use_vr, mut use_fab, mut save) = (false, false, false, false);
        let mut close: Option<usize> = None;
        let (textures, pics) = (&mut self.textures, &self.pics);
        egui::Window::new("MCMC diagnostics").id(egui::Id::new("mcmc-diagnostics")).open(&mut open).default_width(720.0).default_height(720.0).show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                load_logs = ui.button("Load BayesTraits log files…").on_hover_text("One or more *.Log.txt files: runs of the same analysis, compared and pooled (opens a new tab)").clicked();
                use_vr = ui.button("Load variable-rates files…").on_hover_text("BayesTraits VarRates files of a variable-rates analysis (several runs are compared and pooled; their logs are added when found next to them). With a tree open, they are also summarized on it, as with Tree › Import BayesTraits variable rates").clicked();
                use_fab = ui.button("Load Fabric files…").on_hover_text("BayesTraits VarRates files of a Fabric analysis (several runs are compared and pooled; their logs are added when found next to them). With a tree open, they are also summarized on it, as with Tree › Import BayesTraits Fabric").clicked();
            });
            if w.tabs.is_empty() {
                ui.weak("Load BayesTraits log files, or the VarRates files of variable-rates or Fabric runs.");
                return;
            }
            // Tabs: one per set of runs.
            w.active = w.active.min(w.tabs.len() - 1);
            ui.horizontal_wrapped(|ui| {
                for (i, d) in w.tabs.iter().enumerate() {
                    if ui.selectable_label(w.active == i, d.title.as_str()).clicked() && w.active != i {
                        w.active = i;
                        w.param = 0;
                    }
                    if ui.small_button("×").on_hover_text("Close this tab").clicked() {
                        close = Some(i);
                    }
                    ui.add_space(6.0);
                }
            });
            ui.separator();
            let d = &w.tabs[w.active];
            if d.params.is_empty() {
                ui.weak("No parameters in these runs.");
                return;
            }
            egui::ScrollArea::vertical().id_salt("mcmc-scroll").show(ui, |ui| {
                // Runs.
                let total: usize = d.runs.iter().map(|r| r.1).sum();
                ui.label(egui::RichText::new(format!("{} ({} samples)", count(d.runs.len(), "run", "runs"), total)).strong());
                for (name, n) in &d.runs {
                    ui.small(format!("{}: {} samples", name, n));
                }
                if d.runs.iter().any(|r| r.1 != d.runs[0].1) {
                    ui.small("Runs differ in length: R-hat and ESS compare them over the shortest.");
                }
                // Parameters: click a name to plot it.
                ui.add_space(6.0);
                let (header, rows) = d.table();
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new("Parameters").strong());
                    if ui.small_button("Copy table").on_hover_text("Tab-separated").clicked() {
                        ui.ctx().copy_text(io::data::table_to_text(&header, &rows, b'\t'));
                    }
                });
                egui::ScrollArea::both().id_salt("mcmc-params").max_height(220.0).show(ui, |ui| {
                    egui::Grid::new("mcmc-params-grid").striped(true).show(ui, |ui| {
                        let tid = egui::Id::new("mcmc-params-table");
                        table::headers(ui, tid, &header);
                        ui.end_row();
                        for k in table::order(ui, tid, &rows) {
                            for (j, v) in rows[k].iter().enumerate() {
                                if j == 0 {
                                    if ui.selectable_label(w.param == k, v).on_hover_text("Plot this parameter").clicked() {
                                        w.param = k;
                                    }
                                } else {
                                    ui.label(v);
                                }
                            }
                            ui.end_row();
                        }
                    });
                });
                // Plot of the chosen parameter.
                ui.add_space(6.0);
                w.param = w.param.min(d.params.len() - 1);
                ui.horizontal_wrapped(|ui| {
                    egui::ComboBox::from_id_salt("mcmc-param").selected_text(d.params[w.param].name.as_str()).width(260.0).show_ui(ui, |ui| {
                        for (i, p) in d.params.iter().enumerate() {
                            ui.selectable_value(&mut w.param, i, p.name.as_str());
                        }
                    });
                    ui.selectable_value(&mut w.kind, crate::plots::McmcPlot::Density, "Density");
                    ui.selectable_value(&mut w.kind, crate::plots::McmcPlot::Box, "Box and whiskers");
                    ui.selectable_value(&mut w.kind, crate::plots::McmcPlot::Trace, "Trace");
                    if w.kind == crate::plots::McmcPlot::Density {
                        ui.checkbox(&mut w.shade, "Shade").on_hover_text("Fill the area under each curve");
                        if w.shade {
                            ui.add(egui::Slider::new(&mut w.opacity, 0.05..=1.0).text("opacity").fixed_decimals(2));
                        }
                    }
                    save = ui.button("Save plot…").on_hover_text("PNG (300 dpi), SVG or PDF, 6 × 4 in").clicked();
                });
                let names: Vec<String> = d.runs.iter().map(|r| r.0.clone()).collect();
                let p = &d.params[w.param];
                // Box plots keep their natural width; the others fill the window.
                let wide = ui.available_width().max(300.0);
                let size = egui::vec2(crate::plots::mcmc_width(w.kind, &names).map(|b| b.min(wide)).unwrap_or(wide), 300.0);
                let (resp, painter) = ui.allocate_painter(size, egui::Sense::hover());
                painter.rect_filled(resp.rect, 0.0, Color32::WHITE);
                let scene = crate::plots::mcmc_param(&p.chains, &names, w.kind, if w.shade { w.opacity } else { 0.0 }, &p.name, size.x, size.y, 1.0);
                canvas::paint_scene(&painter, &scene, resp.rect.min, textures, pics);
                // Scalars across runs (variable rates, Fabric).
                if let Some(r) = &d.report {
                    ui.add_space(6.0);
                    report_ui(ui, r, "mcmc-report");
                }
                // What the measures mean, each with the work it comes from.
                ui.add_space(6.0);
                ui.label(egui::RichText::new("What these mean").strong().small());
                explain_with_refs(ui, "R-hat", "Do the runs agree? It compares the spread between runs with the spread within them (a single run is split in half). About 1.00 means they sample the same distribution. Values greater than 1.1 indicate inadequate convergence.", &[0, 1]);
                explain_with_refs(ui, "ESS", "How much information? The number of independent draws the (autocorrelated) samples are worth. 400 or more is good; Tracer's rule of thumb is 200.", &[2, 1, 3]);
                explain_with_refs(ui, "95% HPD", "The shortest interval holding 95% of the pooled samples.", &[]);
                if d.report.is_some() {
                    explain_with_refs(ui, "PP SD", "Do the runs put scalars in the same places? The spread across runs of each branch's probability of a scalar (or β). Averaged (like MrBayes' ASDSF), below 0.01 is very good and below 0.05 acceptable.", &[4]);
                }
                explain_with_refs(ui, "Plots", "One colour per run. Density curves of converged runs overlap; a trace should look like stationary noise, with no trend and the runs mixed together. Box plots (Tufte's form): lines are the whiskers to 1.5 IQR, the gap between them the quartiles, the dot the median, and the red line beside them the 95% HPD.", &[]);
                ui.small("BayesTraits has already removed burn-in from its logs.");
            });
        });
        if let (Some(i), Some(w)) = (close, self.mcmc_win.as_mut()) {
            w.tabs.remove(i);
            w.active = w.active.min(w.tabs.len().saturating_sub(1));
            w.param = 0;
        }
        if save {
            if let Some(w) = self.mcmc_win.as_ref().filter(|w| !w.tabs.is_empty()) {
                let d = &w.tabs[w.active.min(w.tabs.len() - 1)];
                let names: Vec<String> = d.runs.iter().map(|r| r.0.clone()).collect();
                let p = d.params[w.param.min(d.params.len().saturating_sub(1))].clone();
                let kind = w.kind;
                let shade = if w.shade { w.opacity } else { 0.0 };
                let file = match kind {
                    crate::plots::McmcPlot::Density => "mcmc_density.png",
                    crate::plots::McmcPlot::Box => "mcmc_box.png",
                    crate::plots::McmcPlot::Trace => "mcmc_trace.png",
                };
                let width = crate::plots::mcmc_width(kind, &names).unwrap_or(432.0);
                self.save_plot_sized(&|wd, ht, pt| crate::plots::mcmc_param(&p.chains, &names, kind, shade, &p.name, wd, ht, pt), file, width);
            }
        }
        if use_vr {
            self.read_bayestraits(crate::bayestraits::RunKind::VariableRates, true);
        }
        if use_fab {
            self.read_bayestraits(crate::bayestraits::RunKind::Fabric, true);
        }
        if load_logs {
            self.load_mcmc_logs();
        }
        if !open {
            if let Some(w) = self.mcmc_win.as_mut() {
                w.open = false;
            }
        }
    }

    /// Show diagnostics in the MCMC diagnostics window (opening it), in a tab.
    fn show_mcmc(&mut self, d: crate::mcmc::Diagnostics) {
        self.mcmc_win.get_or_insert_with(McmcWindow::new).show(d);
    }

    /// Pick BayesTraits log files (runs of one analysis) and show their diagnostics.
    fn load_mcmc_logs(&mut self) {
        let Some(mut paths) = rfd::FileDialog::new().set_title("BayesTraits log file(s): pick several runs to compare them").add_filter("BayesTraits log", &["txt", "log"]).pick_files() else { return };
        paths.sort();
        let mut logs = Vec::new();
        for p in &paths {
            let name = p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
            match std::fs::read_to_string(p).map_err(anyhow::Error::from).and_then(|t| crate::bayestraits::parse_log(&t)) {
                Ok(cols) => logs.push((name, cols)),
                Err(e) => return self.fail(format!("Could not read {}", name), format!("{:#}", e)),
            }
        }
        let title = if logs.len() == 1 { logs[0].0.clone() } else { logs.iter().map(|l| l.0.as_str()).collect::<Vec<_>>().join(" + ") };
        let d = crate::mcmc::Diagnostics::from_logs(&title, &logs);
        self.status = format!("MCMC diagnostics: {} parameters from {}", d.params.len(), count(logs.len(), "run", "runs"));
        self.show_mcmc(d);
    }
}

/// A support file held in memory, as listed under Loaded files.
#[derive(Clone, Copy)]
enum Unload {
    Data(usize),
    Fossils(usize),
    VarRates,
    /// One run of several variable-rates runs.
    VarRatesRun(usize),
    Fabric,
    FabricRun(usize),
    /// An MCMC diagnostics tab of log files.
    Mcmc(usize),
    /// One run (log file) of an MCMC diagnostics tab.
    McmcRun(usize, usize),
}

/// A row under Loaded files: what to unload, kind, name, detail, and whether
/// it is one run of the set above it.
type LoadedFile = (Unload, &'static str, String, String, bool);

impl CanopyApp {
    /// Every support file held in memory: comparative data and fossil ranges per
    /// tab, BayesTraits variable-rates and Fabric results (each run of several
    /// on its own row) and MCMC diagnostics tabs of log files (each log too).
    fn loaded_files(&self) -> Vec<LoadedFile> {
        let mut out = Vec::new();
        for (i, d) in self.docs.iter().enumerate() {
            if let Some(t) = &d.data {
                out.push((Unload::Data(i), "Comparative data", t.name.clone(), format!("{} rows, {} columns · tab {}", t.rows.len(), t.columns.len().saturating_sub(1), d.name), false));
            }
            let fossils = d.tree.tips().into_iter().filter(|&t| crate::fossil::range(&d.tree, t).is_some()).count();
            if fossils > 0 {
                out.push((Unload::Fossils(i), "Fossil ranges", format!("{} tips", fossils), format!("tab {}", d.name), false));
            }
        }
        let tab = |i: usize| self.docs.get(i).map(|d| format!(" · tab {}", d.name)).unwrap_or_default();
        let runs = |out: &mut Vec<LoadedFile>, names: &[String], unload: fn(usize) -> Unload| {
            if names.len() > 1 {
                for (i, n) in names.iter().enumerate() {
                    out.push((unload(i), "Run", n.clone(), String::new(), true));
                }
            }
        };
        if let Some((name, rates, _)) = &self.varrates {
            out.push((Unload::VarRates, "Variable rates", name.clone(), format!("{} branches summarized{}", rates.len(), tab(self.varrates_doc)), false));
            runs(&mut out, &self.varrates_runs.names, Unload::VarRatesRun);
        }
        if let Some((name, rows, _)) = &self.fabric {
            out.push((Unload::Fabric, "Fabric", name.clone(), format!("{} branches summarized{}", rows.len(), tab(self.fabric_doc)), false));
            runs(&mut out, &self.fabric_runs.names, Unload::FabricRun);
        }
        if let Some(w) = &self.mcmc_win {
            for (i, d) in w.tabs.iter().enumerate() {
                // Diagnostics of the imported runs go with them (rows above).
                if self.is_imported_diag(&d.title) {
                    continue;
                }
                let samples: usize = d.runs.iter().map(|r| r.1).sum();
                out.push((Unload::Mcmc(i), "MCMC diagnostics", d.title.clone(), format!("{}, {} samples, {} parameters", count(d.runs.len(), "run", "runs"), samples, d.params.len()), false));
                if d.runs.len() > 1 {
                    for (r, (n, s)) in d.runs.iter().enumerate() {
                        out.push((Unload::McmcRun(i, r), "Run", n.clone(), format!("{} samples", s), true));
                    }
                }
            }
        }
        out
    }

    /// Whether an MCMC diagnostics tab belongs to the imported variable-rates or Fabric runs.
    fn is_imported_diag(&self, title: &str) -> bool {
        [&self.varrates_diag, &self.fabric_diag].iter().any(|d| d.as_ref().map(|d| d.title == title).unwrap_or(false))
    }

    /// The Loaded files list: each support file with a button to unload it (and
    /// one to show results in their window); returns the one clicked, true for Show.
    fn loaded_files_ui(ui: &mut egui::Ui, files: &[LoadedFile]) -> Option<(Unload, bool)> {
        let mut unload = None;
        if files.is_empty() {
            ui.weak("No support files loaded.");
        }
        for (u, kind, what, detail, sub) in files {
            ui.horizontal_wrapped(|ui| {
                if *sub {
                    ui.add_space(18.0);
                }
                let hover = if *sub { "Leave this run out and pool the others again" } else { "Remove it from Canopy's memory" };
                if ui.small_button("Unload").on_hover_text(hover).clicked() {
                    unload = Some((*u, false));
                }
                if matches!(u, Unload::VarRates | Unload::Fabric | Unload::Mcmc(_)) && ui.small_button("Show").on_hover_text("Open its window").clicked() {
                    unload = Some((*u, true));
                }
                if *sub {
                    ui.label(what.as_str());
                } else {
                    ui.label(egui::RichText::new(*kind).strong());
                    ui.label(what.as_str());
                }
                if !detail.is_empty() {
                    ui.weak(detail.as_str());
                }
            });
        }
        unload
    }

    /// Drop the variable-rates results and their MCMC diagnostics tab.
    fn unload_varrates(&mut self) {
        self.varrates = None;
        self.varrates_runs = RunFiles::default();
        if let Some(d) = self.varrates_diag.take() {
            self.close_mcmc_tab(&d.title);
        }
    }

    /// Drop the Fabric results and their MCMC diagnostics tab.
    fn unload_fabric(&mut self) {
        self.fabric = None;
        self.fabric_runs = RunFiles::default();
        if let Some(d) = self.fabric_diag.take() {
            self.close_mcmc_tab(&d.title);
        }
    }

    fn close_mcmc_tab(&mut self, title: &str) {
        if let Some(w) = self.mcmc_win.as_mut() {
            w.tabs.retain(|t| t.title != title);
            w.active = w.active.min(w.tabs.len().saturating_sub(1));
            w.param = 0;
        }
    }

    /// Open the window of loaded results.
    fn show_file(&mut self, u: Unload) {
        match u {
            Unload::VarRates => self.show_varrates = true,
            Unload::Fabric => self.show_fabric = true,
            Unload::Mcmc(i) => {
                if let Some(w) = self.mcmc_win.as_mut() {
                    w.open = true;
                    w.active = i.min(w.tabs.len().saturating_sub(1));
                    w.param = 0;
                }
            }
            _ => {}
        }
    }

    fn unload(&mut self, u: Unload) {
        const KEPT: &str = "values already drawn on the tree stay until you remove those layers";
        match u {
            Unload::Data(i) => {
                // The table and the values it put on the tips.
                if let Some(d) = self.docs.get_mut(i) {
                    if let Some(t) = d.data.take() {
                        let cols: Vec<String> = t.columns.iter().enumerate().filter(|(c, _)| *c != t.key).map(|(_, n)| n.clone()).collect();
                        d.rejoin_data(&cols);
                        self.status = format!("Unloaded comparative data {}", t.name);
                    }
                }
            }
            Unload::Fossils(i) => {
                if let Some(d) = self.docs.get_mut(i) {
                    d.checkpoint();
                    for t in d.tree.tips() {
                        crate::fossil::clear_range(&mut d.tree, t);
                    }
                    self.status = "Removed the fossil ranges".into();
                }
            }
            Unload::VarRates => {
                self.unload_varrates();
                self.status = format!("Unloaded the variable-rates results ({})", KEPT);
            }
            Unload::Fabric => {
                self.unload_fabric();
                self.status = format!("Unloaded the Fabric results ({})", KEPT);
            }
            Unload::VarRatesRun(i) => {
                let name = self.varrates_runs.names.get(i).cloned().unwrap_or_default();
                let mut rest = std::mem::take(&mut self.varrates_runs);
                rest.remove(i);
                if self.apply_varrates(rest, self.varrates_doc, None) {
                    self.status = format!("Unloaded {}; {}", name, self.status);
                }
            }
            Unload::FabricRun(i) => {
                let name = self.fabric_runs.names.get(i).cloned().unwrap_or_default();
                let mut rest = std::mem::take(&mut self.fabric_runs);
                rest.remove(i);
                if self.apply_fabric(rest, self.fabric_doc, None) {
                    self.status = format!("Unloaded {}; {}", name, self.status);
                }
            }
            Unload::Mcmc(i) => {
                if let Some(w) = self.mcmc_win.as_mut() {
                    if i < w.tabs.len() {
                        w.tabs.remove(i);
                        w.active = w.active.min(w.tabs.len().saturating_sub(1));
                        w.param = 0;
                    }
                }
                self.status = "Unloaded the MCMC diagnostics".into();
            }
            Unload::McmcRun(t, r) => {
                if let Some(d) = self.mcmc_win.as_mut().and_then(|w| w.tabs.get_mut(t)) {
                    let name = d.runs.get(r).map(|x| x.0.clone()).unwrap_or_default();
                    d.remove_run(r);
                    if let Some(w) = self.mcmc_win.as_mut() {
                        w.param = 0;
                    }
                    self.status = format!("Unloaded {}", name);
                }
            }
        }
    }

    /// Close tab `i` and unload the BayesTraits results imported onto it.
    fn close_doc(&mut self, i: usize) {
        if i >= self.docs.len() {
            return;
        }
        self.docs.remove(i);
        if self.varrates.is_some() {
            if self.varrates_doc == i {
                self.unload_varrates();
            } else if self.varrates_doc > i {
                self.varrates_doc -= 1;
            }
        }
        if self.fabric.is_some() {
            if self.fabric_doc == i {
                self.unload_fabric();
            } else if self.fabric_doc > i {
                self.fabric_doc -= 1;
            }
        }
    }
}

/// BayesTraits runs read from VarRates files (name, path and parsed run each),
/// kept so that one run can be unloaded and the rest pooled again.
struct RunFiles<T> {
    names: Vec<String>,
    paths: Vec<PathBuf>,
    runs: Vec<T>,
}

impl<T> Default for RunFiles<T> {
    fn default() -> Self {
        RunFiles { names: Vec::new(), paths: Vec::new(), runs: Vec::new() }
    }
}

impl<T> RunFiles<T> {
    /// Leave out run `i`.
    fn remove(&mut self, i: usize) {
        if i < self.runs.len() {
            self.names.remove(i);
            self.paths.remove(i);
            self.runs.remove(i);
        }
    }
}

/// The name of a set of runs: the file, or "n runs (a + b)".
fn runs_label(names: &[String]) -> String {
    if names.len() == 1 { names[0].clone() } else { format!("{} runs ({})", names.len(), names.join(" + ")) }
}

/// Short verdict of a convergence report for the status bar.
fn runs_verdict(r: &crate::mcmc::Report) -> &'static str {
    if r.findings().iter().all(|f| f.0) { "runs agree" } else { "check convergence" }
}

/// The BayesTraits logs next to VarRates files (`X.VarRates.txt` → `X.Log.txt`),
/// if every run has one, as diagnostics parameters (the log-likelihood, already
/// in the VarRates, left out).
fn sibling_log_params(paths: &[std::path::PathBuf]) -> Vec<crate::mcmc::Param> {
    let mut logs = Vec::new();
    for p in paths {
        let s = p.to_string_lossy();
        let Some(stem) = s.strip_suffix(".VarRates.txt") else { return Vec::new() };
        let Ok(cols) = std::fs::read_to_string(format!("{}.Log.txt", stem)).map_err(anyhow::Error::from).and_then(|t| crate::bayestraits::parse_log(&t)) else { return Vec::new() };
        logs.push((String::new(), cols));
    }
    crate::mcmc::Diagnostics::from_logs("", &logs).params.into_iter().filter(|p| p.name != "Lh").collect()
}

/// Tree traversal window state.
struct TraversalWindow {
    order: crate::traversal::Order,
    /// Visits shown so far.
    step: usize,
    /// Visits in the current tree (updated as it is drawn).
    total: usize,
    playing: bool,
    /// Steps per second, also the GIF frame rate.
    speed: f32,
    /// Time of the last step while playing.
    last: f64,
    numbers: bool,
    width_in: f32,
    height_in: f32,
    dpi: f32,
}

/// Tanglegram window state: the two tabs (indices), the drawing style and the export size.
struct TangleWindow {
    left: usize,
    right: usize,
    style: crate::tangle::TangleStyle,
    width_in: f32,
    height_in: f32,
}

/// State of the lineages-through-time window.
#[derive(Default)]
struct LttWindow {
    log: bool,
    smooth: bool,
    /// Geologic timescale under the time axis.
    geoscale: bool,
    show_set: bool,
    /// Curves of the trees in the set, cached by (tab, set size, burn-in).
    set_curves: Option<(String, Vec<Vec<(f64, usize)>>)>,
    /// The solid curve of a tree set: its MCC tree (else its consensus).
    mcc: bool,
    /// Draw the solid curve's tree above the plot, its nodes joined to the curve.
    show_tree: bool,
    /// The summary tree (consensus or MCC) of a tree set, cached by the same key.
    main_curve: Option<(String, Tree)>,
}
