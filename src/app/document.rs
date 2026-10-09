//! An open document: a tree, its view state and layers, attached data,
//! an optional tree set (posterior sample or bootstrap replicates), selection
//! and undo history.

use crate::consensus::{self, HeightMode, PosteriorSummary};
use std::sync::Arc;
use crate::io::data::DataTable;
use crate::ops;
use crate::style::*;
use crate::tree::{Attr, NodeId, Tree};
use std::collections::BTreeSet;
use std::path::PathBuf;

#[derive(Clone)]
pub struct Snapshot {
    pub tree: Tree,
    pub view: ViewState,
    pub layers: Vec<LayerEntry>,
}

/// What a tree set is: it decides the wording, the burn-in and the name of the
/// clade-support attribute.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SampleKind {
    /// Bayesian MCMC sample (BEAST, MrBayes): burn-in, posterior probabilities.
    Posterior,
    /// Bootstrap replicates (RAxML, IQ-TREE): no burn-in, bootstrap proportions.
    Bootstrap,
}

impl SampleKind {
    pub fn name(self) -> &'static str {
        match self {
            SampleKind::Posterior => "Bayesian posterior sample",
            SampleKind::Bootstrap => "Bootstrap replicates",
        }
    }

    /// Attribute holding clade frequencies.
    pub fn support_key(self) -> &'static str {
        match self {
            SampleKind::Posterior => "posterior",
            SampleKind::Bootstrap => "bootstrap",
        }
    }

    /// Guess from the file name: bootstrap files usually say so
    /// (RAxML_bootstrap.*, *.boottrees, *.ufboot).
    pub fn guess(file_name: &str) -> SampleKind {
        if file_name.to_ascii_lowercase().contains("boot") {
            SampleKind::Bootstrap
        } else {
            SampleKind::Posterior
        }
    }
}

/// A set of trees: a Bayesian posterior sample or bootstrap replicates.
pub struct Posterior {
    pub kind: SampleKind,
    /// Shared, so consensus/MCC tabs can keep the sample (for DensiTree) cheaply.
    pub trees: Arc<Vec<Tree>>,
    /// Fraction of samples discarded as burn-in.
    pub burnin: f32,
    pub rooted: bool,
    pub threshold: f32,
    pub heights: HeightMode,
    pub browse: usize,
    pub info: String,
    /// The document state before stepping through samples began.
    pub reference: Option<Snapshot>,
    /// Ladderize each sample so successive trees are comparable.
    pub ladderize: bool,
    pub playing: bool,
    /// Samples per second during playback.
    pub speed: f32,
    pub last_step: Option<std::time::Instant>,
    /// Clade frequencies used to label each sample, with the settings they were computed for.
    summary: Option<(f32, bool, Arc<PosteriorSummary>)>,
}

impl Posterior {
    pub fn new(trees: Vec<Tree>) -> Self {
        Self::shared(Arc::new(trees))
    }

    pub fn shared(trees: Arc<Vec<Tree>>) -> Self {
        let rooted = trees.first().map(|t| t.rooted != Some(false)).unwrap_or(true);
        Posterior {
            kind: SampleKind::Posterior,
            trees,
            burnin: 0.0,
            rooted,
            threshold: 0.5,
            heights: HeightMode::Mean,
            browse: 0,
            info: String::new(),
            reference: None,
            ladderize: true,
            playing: false,
            speed: 4.0,
            last_step: None,
            summary: None,
        }
    }

    /// Set the kind of tree set; bootstrap replicates have no burn-in.
    pub fn set_kind(&mut self, kind: SampleKind) {
        self.kind = kind;
        if kind == SampleKind::Bootstrap {
            self.burnin = 0.0;
        }
        self.summary = None;
    }

    /// Clade frequencies over the post-burn-in trees (cached).
    pub fn summary(&mut self) -> Option<Arc<PosteriorSummary>> {
        if let Some((b, r, s)) = &self.summary {
            if *b == self.burnin && *r == self.rooted {
                return Some(s.clone());
            }
        }
        match consensus::summarize(self.kept(), self.rooted) {
            Ok(mut s) => {
                s.support_key = self.kind.support_key();
                let s = Arc::new(s);
                self.summary = Some((self.burnin, self.rooted, s.clone()));
                Some(s)
            }
            Err(e) => {
                self.info = format!("Cannot compute support: {:#}", e);
                None
            }
        }
    }

    pub fn browsing(&self) -> bool {
        self.reference.is_some()
    }

    pub fn kept(&self) -> &[Tree] {
        let skip = ((self.trees.len() as f32) * self.burnin).floor() as usize;
        &self.trees[skip.min(self.trees.len().saturating_sub(1))..]
    }
}

pub struct Document {
    pub name: String,
    pub path: Option<PathBuf>,
    pub tree: Tree,
    pub view: ViewState,
    pub layers: Vec<LayerEntry>,
    pub data: Option<DataTable>,
    pub posterior: Option<Posterior>,
    /// For a consensus or MCC tree: the sample it summarizes (shared) and the
    /// number of burn-in trees to skip. Used only for DensiTree; the document
    /// itself is a single tree.
    pub source_sample: Option<(Arc<Vec<Tree>>, usize)>,
    pub selection: BTreeSet<NodeId>,
    pub selected_layer: Option<usize>,
    pub zoom: [f32; 2],
    pub pan: [f32; 2],
    /// Canvas size in points at the last frame (for zooming about its centre).
    pub view_size: [f32; 2],
    undo: Vec<Snapshot>,
    redo: Vec<Snapshot>,
}

impl Document {
    pub fn new(name: &str, tree: Tree, phylopic: bool) -> Self {
        let layers = smart_layers(&tree, phylopic);
        Document {
            name: name.to_string(),
            path: None,
            tree,
            view: ViewState::default(),
            layers,
            data: None,
            posterior: None,
            source_sample: None,
            selection: BTreeSet::new(),
            selected_layer: None,
            zoom: [1.0, 1.0],
            pan: [0.0, 0.0],
            view_size: [0.0, 0.0],
            undo: Vec::new(),
            redo: Vec::new(),
        }
    }

    /// Multiply the zoom by (fx, fy), keeping the point `about` (canvas
    /// coordinates; None = the canvas centre) fixed on screen.
    pub fn zoom_by(&mut self, fx: f32, fy: f32, about: Option<[f32; 2]>) {
        let c = about.unwrap_or([self.view_size[0] / 2.0, self.view_size[1] / 2.0]);
        let nz = [(self.zoom[0] * fx).clamp(0.2, 60.0), (self.zoom[1] * fy).clamp(0.2, 200.0)];
        for i in 0..2 {
            let f = nz[i] / self.zoom[i];
            self.pan[i] = c[i] - (c[i] - self.pan[i]) * f;
        }
        self.zoom = nz;
    }

    fn snapshot(&self) -> Snapshot {
        Snapshot { tree: self.tree.clone(), view: self.view.clone(), layers: self.layers.clone() }
    }

    /// Record the current state before a change.
    pub fn checkpoint(&mut self) {
        self.undo.push(self.snapshot());
        if self.undo.len() > 100 {
            self.undo.remove(0);
        }
        self.redo.clear();
    }

    fn restore(&mut self, s: Snapshot) {
        self.tree = s.tree;
        self.view = s.view;
        self.layers = s.layers;
        let reach = self.tree.reachable();
        self.selection.retain(|&n| n < reach.len() && reach[n]);
        if self.selected_layer.map(|i| i >= self.layers.len()).unwrap_or(false) {
            self.selected_layer = None;
        }
    }

    pub fn undo(&mut self) -> bool {
        match self.undo.pop() {
            Some(s) => {
                self.redo.push(self.snapshot());
                self.restore(s);
                true
            }
            None => false,
        }
    }

    pub fn redo(&mut self) -> bool {
        match self.redo.pop() {
            Some(s) => {
                self.undo.push(self.snapshot());
                self.restore(s);
                true
            }
            None => false,
        }
    }

    pub fn can_undo(&self) -> bool {
        !self.undo.is_empty()
    }
    pub fn can_redo(&self) -> bool {
        !self.redo.is_empty()
    }

    /// Post-burn-in posterior trees, for DensiTree layers: this document's own
    /// tree set, or the sample a consensus/MCC tree was built from.
    pub fn overlay(&self) -> &[Tree] {
        match (&self.posterior, &self.source_sample) {
            (Some(p), _) => p.kept(),
            (None, Some((trees, skip))) => &trees[(*skip).min(trees.len())..],
            _ => &[],
        }
    }

    /// Tips under the selection (selected tips plus tips of selected clades).
    pub fn selected_tips(&self) -> Vec<NodeId> {
        let mut v: Vec<NodeId> = self.selection.iter().flat_map(|&n| self.tree.tips_below(n)).collect();
        v.sort();
        v.dedup();
        v
    }

    pub fn single_selection(&self) -> Option<NodeId> {
        if self.selection.len() == 1 {
            self.selection.iter().next().copied()
        } else {
            None
        }
    }

    /// Show posterior sample `i` in this window, carrying highlights, clade
    /// labels, clade colors and collapsed clades over by taxon set.
    pub fn show_sample(&mut self, i: usize) {
        let Some(p) = self.posterior.as_mut() else { return };
        if p.trees.is_empty() {
            return;
        }
        let i = i.min(p.trees.len() - 1);
        if p.reference.is_none() {
            p.reference = Some(Snapshot { tree: self.tree.clone(), view: self.view.clone(), layers: self.layers.clone() });
            // Show support values while browsing (once; the user may remove it).
            ensure_support_layer(&mut self.layers, p.kind.support_key());
        }
        let mut t = p.trees[i].clone();
        if p.ladderize {
            let r = t.root;
            ops::ladderize(&mut t, r, true);
        }
        if let Some(d) = &self.data {
            d.join_to_tree(&mut t);
        }
        // Label every clade with its support (posterior or bootstrap) across the set.
        if let Some(s) = p.summary() {
            let _ = consensus::annotate(&mut t, &s, HeightMode::Keep);
        }
        p.browse = i;
        let reference = p.reference.clone().unwrap();
        let prev = std::mem::replace(&mut self.tree, t);
        remap(&prev, &reference, &self.tree, &mut self.view, &mut self.layers);
        self.selection.clear();
    }

    /// Label the displayed tree with clade support from the tree set
    /// (used when a tree set is opened and when burn-in changes).
    pub fn annotate_posterior(&mut self) {
        let Some(p) = self.posterior.as_mut() else { return };
        if let Some(s) = p.summary() {
            // Drop values left from the set's other kind (posterior vs bootstrap).
            let other = if p.kind == SampleKind::Posterior { SampleKind::Bootstrap } else { SampleKind::Posterior };
            for n in self.tree.preorder() {
                self.tree.nodes[n].attrs.remove(other.support_key());
            }
            let _ = consensus::annotate(&mut self.tree, &s, HeightMode::Keep);
            ensure_support_layer(&mut self.layers, p.kind.support_key());
        }
    }

    /// Stop browsing and return to the tree shown before.
    pub fn leave_samples(&mut self) {
        let Some(p) = self.posterior.as_mut() else { return };
        let Some(reference) = p.reference.take() else { return };
        p.playing = false;
        let prev = std::mem::replace(&mut self.tree, reference.tree.clone());
        remap(&prev, &reference, &self.tree, &mut self.view, &mut self.layers);
        self.selection.clear();
        self.annotate_posterior();
    }

    /// Set a node's label. For an internal node a numeric label is its support
    /// (as when a tree is read), so `support` follows it, and a label that
    /// replaces a numeric one drops the support it carried.
    pub fn set_label(&mut self, n: NodeId, label: Option<String>) {
        let numbers = |l: Option<&str>| -> Option<Vec<f64>> { l.and_then(|l| l.split('/').map(|p| p.trim().parse::<f64>().ok()).collect()) };
        let was = numbers(self.tree.nodes[n].label.as_deref());
        let now = numbers(label.as_deref());
        self.tree.nodes[n].label = label;
        if self.tree.is_tip(n) {
            return;
        }
        let attrs = &mut self.tree.nodes[n].attrs;
        if was.is_some() || now.is_some() {
            attrs.retain(|k, _| k != "support" && !k.starts_with("support_"));
        }
        for (i, v) in now.unwrap_or_default().into_iter().enumerate() {
            let key = if i == 0 { "support".to_string() } else { format!("support_{}", i + 1) };
            attrs.insert(key, Attr::Num(v));
        }
    }

    /// After editing internal node `n`'s label, make sure the figure shows it:
    /// a visible node-label layer for its label, or for `support` when the
    /// label is a number and support labels are already shown.
    pub fn show_node_label(&mut self, n: NodeId) {
        if self.tree.is_tip(n) || self.tree.nodes[n].label.is_none() {
            return;
        }
        let numeric = self.tree.nodes[n].attrs.contains_key("support");
        let shown = |attr: &str| self.layers.iter().any(|e| e.enabled && matches!(&e.layer, Layer::NodeLabels(s) if s.attr == attr));
        if shown("label") || (numeric && shown("support")) {
            return;
        }
        let key = if numeric { "support" } else { "label" };
        // Re-enable a hidden layer before adding another.
        if let Some(e) = self.layers.iter_mut().find(|e| matches!(&e.layer, Layer::NodeLabels(s) if s.attr == key)) {
            e.enabled = true;
        } else {
            self.layers.push(LayerEntry::new(Layer::default_node_labels(key)));
        }
    }

    pub fn add_layer(&mut self, layer: Layer) {
        self.checkpoint();
        // Highlights belong underneath everything else.
        if matches!(layer, Layer::Highlight(_)) {
            self.layers.insert(0, LayerEntry::new(layer));
            self.selected_layer = Some(0);
        } else {
            self.layers.push(LayerEntry::new(layer));
            self.selected_layer = Some(self.layers.len() - 1);
        }
    }
}

/// Make sure a node-label layer showing clade support (`key`) exists and is visible.
fn ensure_support_layer(layers: &mut Vec<LayerEntry>, key: &str) {
    let mut found = false;
    for e in layers.iter_mut() {
        if let Layer::NodeLabels(s) = &mut e.layer {
            // A label layer for the other kind of support switches over.
            if s.attr == key || (matches!(s.attr.as_str(), "posterior" | "bootstrap") && !found) {
                s.attr = key.to_string();
                e.enabled = true;
                found = true;
            }
        }
    }
    if !found {
        layers.push(LayerEntry::new(Layer::default_node_labels(key)));
    }
}

fn layer_node(l: &Layer) -> Option<NodeId> {
    match l {
        Layer::Highlight(h) => Some(h.node),
        Layer::CladeLabel(c) => Some(c.node),
        _ => None,
    }
}

/// Re-point node-based annotations from `prev` to `new`. Layers that existed
/// when browsing started are resolved from the reference tree, so a clade
/// that is not monophyletic in one sample reappears in the next one that has it.
fn remap(prev: &Tree, reference: &Snapshot, new: &Tree, view: &mut ViewState, layers: &mut [LayerEntry]) {
    view.collapsed = view.collapsed.iter().filter_map(|&n| ops::map_node(prev, n, new)).collect();
    view.clade_colors = view.clade_colors.iter().filter_map(|(&n, &c)| ops::map_node(prev, n, new).map(|m| (m, c))).collect();
    for (i, e) in layers.iter_mut().enumerate() {
        let Some(cur) = layer_node(&e.layer) else { continue };
        let from_ref = reference
            .layers
            .get(i)
            .filter(|r| std::mem::discriminant(&r.layer) == std::mem::discriminant(&e.layer))
            .and_then(|r| layer_node(&r.layer));
        let mapped = match from_ref {
            Some(rn) => ops::map_node(&reference.tree, rn, new),
            None => ops::map_node(prev, cur, new),
        };
        // An out-of-range id hides the layer until its clade reappears.
        let n = mapped.unwrap_or(usize::MAX);
        match &mut e.layer {
            Layer::Highlight(h) => h.node = n,
            Layer::CladeLabel(c) => c.node = n,
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::newick::parse_newick_multi;

    #[test]
    fn stepped_samples_show_support() {
        let trees = parse_newick_multi("((A:1,B:1):1,(C:1,D:1):1);((A:1,B:1):1,(C:1,D:1):1);((A:1,C:1):1,(B:1,D:1):1);").unwrap();
        let mut d = Document::new("t", trees[0].clone(), false);
        let mut p = Posterior::new(trees);
        p.burnin = 0.0;
        d.posterior = Some(p);
        d.show_sample(1);
        let a = d.tree.find_label("A").unwrap();
        let ab = d.tree.nodes[a].parent.unwrap();
        let post = d.tree.nodes[ab].attrs.get("posterior").and_then(|v| v.as_f64()).unwrap();
        assert!((post - 2.0 / 3.0).abs() < 1e-9);
        assert!(d.layers.iter().any(|e| e.enabled && matches!(&e.layer, Layer::NodeLabels(s) if s.attr == "posterior")));
        d.leave_samples();
        assert!(d.posterior.as_ref().unwrap().reference.is_none());
    }

    #[test]
    fn opened_tree_set_shows_support() {
        let trees = parse_newick_multi("((A:1,B:1):1,(C:1,D:1):1);((A:1,B:1):1,(C:1,D:1):1);((A:1,C:1):1,(B:1,D:1):1);").unwrap();
        let mut d = Document::new("t", trees[0].clone(), false);
        d.posterior = Some(Posterior::new(trees));
        d.posterior.as_mut().unwrap().burnin = 0.0;
        d.annotate_posterior();
        let a = d.tree.find_label("A").unwrap();
        let ab = d.tree.nodes[a].parent.unwrap();
        assert!(d.tree.nodes[ab].attrs.contains_key("posterior"));
        assert!(d.layers.iter().any(|e| e.enabled && matches!(&e.layer, Layer::NodeLabels(s) if s.attr == "posterior")));
    }

    #[test]
    fn bootstrap_set_uses_all_trees_and_bootstrap_labels() {
        assert_eq!(SampleKind::guess("RAxML_bootstrap.run1"), SampleKind::Bootstrap);
        assert_eq!(SampleKind::guess("primates.trees"), SampleKind::Posterior);
        let trees = parse_newick_multi("((A:1,B:1):1,(C:1,D:1):1);((A:1,B:1):1,(C:1,D:1):1);((A:1,C:1):1,(B:1,D:1):1);((A:1,B:1):1,(C:1,D:1):1);").unwrap();
        let mut d = Document::new("t", trees[0].clone(), false);
        let mut p = Posterior::new(trees);
        p.set_kind(SampleKind::Bootstrap);
        assert_eq!(p.kept().len(), 4);
        d.posterior = Some(p);
        d.annotate_posterior();
        let ab = d.tree.nodes[d.tree.find_label("A").unwrap()].parent.unwrap();
        assert_eq!(d.tree.nodes[ab].attrs.get("bootstrap").and_then(|v| v.as_f64()), Some(0.75));
        assert!(!d.tree.nodes[ab].attrs.contains_key("posterior"));
        assert!(d.layers.iter().any(|e| e.enabled && matches!(&e.layer, Layer::NodeLabels(s) if s.attr == "bootstrap")));
    }

    #[test]
    fn edited_node_labels_reach_the_figure() {
        // Bootstrap-style tree: the node-label layer shows `support`.
        let t = crate::io::newick::parse_newick("((A:1,B:1)100:1,(C:1,D:1)75:1);").unwrap();
        let mut d = Document::new("t", t, false);
        let ab = d.tree.nodes[d.tree.find_label("A").unwrap()].parent.unwrap();
        d.set_label(ab, Some("95".into()));
        d.show_node_label(ab);
        assert_eq!(d.tree.value(ab, "support"), Some(Attr::Num(95.0)));
        // A name replaces the number, and a label layer appears to show it.
        d.set_label(ab, Some("Hominidae".into()));
        d.show_node_label(ab);
        assert!(!d.tree.nodes[ab].attrs.contains_key("support"));
        assert!(d.layers.iter().any(|e| e.enabled && matches!(&e.layer, Layer::NodeLabels(s) if s.attr == "label")));
        // No duplicate layer on further edits.
        d.set_label(ab, Some("Hominids".into()));
        d.show_node_label(ab);
        assert_eq!(d.layers.iter().filter(|e| matches!(&e.layer, Layer::NodeLabels(s) if s.attr == "label")).count(), 1);
    }

    #[test]
    fn summary_tree_is_single_but_keeps_densitree_sample() {
        let trees = Arc::new(parse_newick_multi("((A:1,B:1):1,C:2);((A:1,C:1):1,B:2);((A:1,B:1):1,C:2);((B:1,C:1):1,A:2);").unwrap());
        let mut d = Document::new("MCC", trees[0].clone(), false);
        d.source_sample = Some((trees.clone(), 1));
        assert!(d.posterior.is_none());
        assert_eq!(d.overlay().len(), 3);
    }
}

/// Default layers chosen from what the tree contains (like a sensible ggtree recipe).
pub fn smart_layers(tree: &Tree, phylopic: bool) -> Vec<LayerEntry> {
    let keys = tree.attr_keys();
    let has = |k: &str| keys.iter().any(|x| x == k);
    let mut layers = vec![LayerEntry::new(Layer::default_tree())];
    if has("height_95%_HPD") {
        layers.insert(0, LayerEntry::new(Layer::default_node_bars()));
    }
    layers.push(LayerEntry::new(Layer::default_tip_labels()));
    if has("posterior") {
        layers.push(LayerEntry::new(Layer::default_node_labels("posterior")));
    } else if has("support") {
        layers.push(LayerEntry::new(Layer::default_node_labels("support")));
    } else if has("bootstrap") {
        layers.push(LayerEntry::new(Layer::default_node_labels("bootstrap")));
    }
    if phylopic {
        layers.push(LayerEntry::new(Layer::default_phylopic()));
    }
    if has("height_95%_HPD") || has("height") {
        layers.push(LayerEntry::new(Layer::default_axis()));
    } else {
        layers.push(LayerEntry::new(Layer::default_scale_bar()));
    }
    layers
}
