//! Renderers for scenes: high-resolution raster (PNG/TIFF), SVG and PDF.
//! The interactive egui painter lives in `app::canvas`.

pub mod fonts;
pub mod pdf;
pub mod raster;
pub mod svg;

use image::RgbaImage;
use std::collections::HashMap;
use std::sync::Arc;

/// Loaded silhouette images by key.
pub type ImageMap = HashMap<String, Arc<RgbaImage>>;
