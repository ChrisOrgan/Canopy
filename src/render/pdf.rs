//! Vector PDF export: the scene is written as SVG (in points), parsed with
//! usvg using Canopy's plot font, and converted by svg2pdf. Text stays text
//! (the font is embedded as a subset) and silhouettes stay embedded images.

use super::fonts::FontBytes;
use super::{svg, ImageMap};
use crate::scene::Scene;
use crate::style::Color;
use anyhow::{anyhow, Result};
use svg2pdf::usvg;

/// Render a scene built in points (72 per inch) to PDF bytes.
pub fn render(scene: &Scene, images: &ImageMap, background: Option<Color>, fonts: &FontBytes) -> Result<Vec<u8>> {
    let text = svg::render(scene, images, background, &fonts.name);
    let mut opt = usvg::Options::default();
    {
        let db = opt.fontdb_mut();
        db.load_font_data(fonts.regular.clone());
        for f in [&fonts.italic, &fonts.bold].into_iter().flatten() {
            db.load_font_data(f.clone());
        }
        db.set_sans_serif_family(fonts.name.clone());
    }
    let tree = usvg::Tree::from_str(&text, &opt).map_err(|e| anyhow!("could not lay out the figure for PDF: {}", e))?;
    svg2pdf::to_pdf(&tree, svg2pdf::ConversionOptions::default(), svg2pdf::PageOptions { dpi: 72.0 }).map_err(|e| anyhow!("PDF conversion failed: {:?}", e))
}
