//! Canopy: phylogenetic tree visualization and editing.
//!
//! The library is split into a GUI-independent core (tree model, I/O,
//! tree operations, posterior summaries, layouts, scene building and
//! rendering) and the `app` module containing the egui front-end.

pub mod ancestral;
pub mod app;
pub mod bayestraits;
pub mod cli;
pub mod consensus;
pub mod fabric;
pub mod fossil;
pub mod geotime;
pub mod io;
pub mod layout;
pub mod logo;
pub mod ltt;
pub mod nodedensity;
pub mod ops;
pub mod phylopic;
pub mod pir;
pub mod plots;
pub mod rogue;
pub mod render;
pub mod scene;
pub mod style;
pub mod taxonomy;
pub mod tree;
pub mod varrates;

pub use tree::{Attr, Node, NodeId, Tree};
