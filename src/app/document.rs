//! An open document: a tree, its view state and layers, attached data,
//! an optional posterior tree sample, selection and undo history.

use crate::consensus::{self, HeightMode, PosteriorSummary};
use std::sync::Arc;
use crate::io::data::DataTable;
use crate::ops;
use crate::style::*;
use crate::tree::{NodeId, Tree};
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

#[derive(Clone, Serialize, Deserialize)]
pub struct Snapshot {
    pub tree: Tree,
    pub view: ViewState,
    pub layers: Vec<LayerEntry>,
}

pub struct Posterior {
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
            trees,
            burnin: 0.1,
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

    /// Clade frequencies over the post-burn-in trees (cached).
    pub fn summary(&mut self) -> Option<Arc<PosteriorSummary>> {
        if let Some((b, r, s)) = &self.summary {
            if *b == self.burnin && *r == self.rooted {
                return Some(s.clone());
            }
        }
        match consensus::summarize(self.kept(), self.rooted) {
            Ok(s) => {
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
    pub selection: BTreeSet<NodeId>,
    pub selected_layer: Option<usize>,
    pub zoom: [f32; 2],
    pub pan: [f32; 2],
    pub dirty: bool,
    undo: Vec<Snapshot>,
    redo: Vec<Snapshot>,
}

#[derive(Serialize, Deserialize)]
struct ProjectFile {
    format: String,
    version: u32,
    name: String,
    snapshot: Snapshot,
    data: Option<DataTable>,
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
            selection: BTreeSet::new(),
            selected_layer: None,
            zoom: [1.0, 1.0],
            pan: [0.0, 0.0],
            dirty: false,
            undo: Vec::new(),
            redo: Vec::new(),
        }
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
        self.dirty = true;
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

    pub fn save_project(&self, path: &Path) -> Result<()> {
        let pf = ProjectFile {
            format: "canopy".into(),
            version: 1,
            name: self.name.clone(),
            snapshot: self.snapshot(),
            data: self.data.clone(),
        };
        std::fs::write(path, serde_json::to_string(&pf)?).with_context(|| format!("writing {}", path.display()))
    }

    pub fn load_project(path: &Path) -> Result<Document> {
        let pf: ProjectFile = serde_json::from_str(&std::fs::read_to_string(path)?).context("not a Canopy project file")?;
        let mut d = Document::new(&pf.name, pf.snapshot.tree, false);
        d.view = pf.snapshot.view;
        d.layers = pf.snapshot.layers;
        d.data = pf.data;
        d.path = Some(path.to_path_buf());
        Ok(d)
    }

    /// Post-burn-in posterior trees, for DensiTree layers.
    pub fn overlay(&self) -> &[Tree] {
        self.posterior.as_ref().map(|p| p.kept()).unwrap_or(&[])
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
            ensure_support_layer(&mut self.layers);
        }
        let mut t = p.trees[i].clone();
        if p.ladderize {
            let r = t.root;
            ops::ladderize(&mut t, r, true);
        }
        if let Some(d) = &self.data {
            d.join_to_tree(&mut t);
        }
        // Label every clade with its posterior probability across the sample.
        if let Some(s) = p.summary() {
            let _ = consensus::annotate(&mut t, &s, HeightMode::Keep);
        }
        p.browse = i;
        let reference = p.reference.clone().unwrap();
        let prev = std::mem::replace(&mut self.tree, t);
        remap(&prev, &reference, &self.tree, &mut self.view, &mut self.layers);
        self.selection.clear();
    }

    /// Label the displayed tree with posterior clade support from the sample
    /// (used when a tree set is opened and when burn-in changes).
    pub fn annotate_posterior(&mut self) {
        let Some(p) = self.posterior.as_mut() else { return };
        if let Some(s) = p.summary() {
            let _ = consensus::annotate(&mut self.tree, &s, HeightMode::Keep);
            ensure_support_layer(&mut self.layers);
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

/// Make sure a node-label layer showing "posterior" exists and is visible.
fn ensure_support_layer(layers: &mut Vec<LayerEntry>) {
    let mut found = false;
    for e in layers.iter_mut() {
        if matches!(&e.layer, Layer::NodeLabels(s) if s.attr == "posterior") {
            e.enabled = true;
            found = true;
        }
    }
    if !found {
        layers.push(LayerEntry::new(Layer::default_node_labels("posterior")));
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
