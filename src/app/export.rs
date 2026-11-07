//! Figure export dialog: PNG / TIFF (with DPI metadata), SVG and PDF.

use super::document::Document;
use crate::render::fonts::FontBytes;
use crate::render::raster::{self, RasterEnv, RasterFonts};
use crate::render::{pdf, svg, ImageMap};
use crate::scene::{self, BuildInput, Scene, SceneEnv};
use crate::style::Color;
use anyhow::{bail, Result};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Format {
    Png,
    Tiff,
    Svg,
    Pdf,
    /// Animated GIF of the trees in a tree set, as stepped through in the GUI.
    Gif,
}

impl Format {
    fn ext(self) -> &'static str {
        match self {
            Format::Png => "png",
            Format::Tiff => "tif",
            Format::Svg => "svg",
            Format::Pdf => "pdf",
            Format::Gif => "gif",
        }
    }

    fn vector(self) -> bool {
        matches!(self, Format::Svg | Format::Pdf)
    }
}

pub struct ExportDialog {
    pub open: bool,
    pub format: Format,
    pub width_in: f32,
    pub height_in: f32,
    pub dpi: f32,
    pub transparent: bool,
    pub tiff_rgb: bool,
    pub mm: bool,
    /// Animated GIF: screen resolution, frames per second, and every nth tree.
    pub gif_dpi: f32,
    pub gif_fps: f32,
    pub gif_every: usize,
    /// Sizes in points scaled so text and lines take the same share of the
    /// figure as in the tree window (else true point sizes).
    pub match_view: bool,
    /// Size (points) of the tree as drawn in the window: canvas times zoom.
    view: [f32; 2],
    preview: Option<egui::TextureHandle>,
    last_preview: Option<Instant>,
}

impl Default for ExportDialog {
    fn default() -> Self {
        ExportDialog {
            open: false,
            format: Format::Tiff,
            width_in: 6.85,
            height_in: 8.0,
            dpi: 600.0,
            transparent: false,
            tiff_rgb: true,
            mm: false,
            gif_dpi: 100.0,
            gif_fps: 4.0,
            gif_every: 1,
            match_view: true,
            view: [0.0, 0.0],
            preview: None,
            last_preview: None,
        }
    }
}

pub struct ExportJob {
    pub path: PathBuf,
    pub format: Format,
    pub width_px: u32,
    pub height_px: u32,
    pub dpi: f32,
    pub transparent: bool,
    pub tiff_rgb: bool,
    /// Animated GIF: frames per second and every nth tree.
    pub fps: f32,
    pub every: usize,
    /// Multiplier on every size in points (text, lines, points): 1 gives true
    /// point sizes; the dialog's "as in the tree window" makes sizes the same
    /// share of the figure as on screen.
    pub scale: f32,
}

/// Indices of the trees an animated GIF shows: every `every`th tree of the set.
pub fn gif_frames(n: usize, every: usize) -> Vec<usize> {
    (0..n).step_by(every.max(1)).collect()
}

/// Save the trees of the document's tree set as an animated GIF, one frame
/// per tree as "Step through trees" shows them (same layers, highlights and
/// ladderizing), looping forever. Returns the number of frames.
pub fn export_gif(doc: &mut Document, fonts: &RasterFonts, images: &ImageMap, job: &ExportJob) -> Result<usize> {
    use image::codecs::gif::{GifEncoder, Repeat};
    use image::{Delay, Frame};
    if job.width_px == 0 || job.height_px == 0 || job.width_px > 4000 || job.height_px > 4000 {
        bail!("an animated GIF must be between 1 and 4000 pixels per side; lower the resolution");
    }
    let Some(p) = doc.posterior.as_mut() else { bail!("open a set of trees (posterior sample or bootstrap replicates) to export an animation") };
    // Every frame is labelled with clade support: have the summary ready first.
    p.summary_blocking();
    let frames = gif_frames(p.trees.len(), job.every);
    let (was_browsing, was) = (p.browsing(), p.browse);
    let mut enc = GifEncoder::new_with_speed(std::io::BufWriter::new(std::fs::File::create(&job.path)?), 10);
    enc.set_repeat(Repeat::Infinite)?;
    let delay = Delay::from_numer_denom_ms((1000.0 / job.fps.max(0.1)).round() as u32, 1);
    let pt = job.dpi / 72.0 * job.scale;
    let bg = Some(doc.view.background.unwrap_or(Color::WHITE));
    let mut result = Ok(());
    for &i in &frames {
        doc.show_sample(i);
        let env = RasterEnv { fonts, images };
        let input = BuildInput { tree: &doc.tree, view: &doc.view, layers: &doc.layers, overlay: doc.overlay() };
        let sc = scene::build(&input, &env, job.width_px as f32, job.height_px as f32, pt);
        result = raster::render(&sc, fonts, images, bg).and_then(|pm| Ok(enc.encode_frame(Frame::from_parts(raster::to_rgba(&pm), 0, 0, delay))?));
        if result.is_err() {
            break;
        }
    }
    // Back to what was shown before.
    if was_browsing {
        doc.show_sample(was);
    } else {
        doc.leave_samples();
    }
    result.map(|_| frames.len())
}

/// Indices (into the whole set) of the trees a trait-mapping GIF shows: every
/// `every`th tree after the burn-in.
pub fn posterior_frames(n: usize, burnin: f32, every: usize) -> Vec<usize> {
    let skip = ((n as f32) * burnin).floor() as usize;
    (skip.min(n.saturating_sub(1))..n).step_by(every.max(1)).collect()
}

/// Map Trait to Posterior: an animated GIF of the tree set after burn-in, each
/// tree with trait `key` reconstructed on it (ML under Brownian motion when
/// `numeric`, Fitch parsimony otherwise) and its branches colored by it, with
/// the tab's other layers as shown. Returns the number of frames.
pub fn export_trait_gif(doc: &mut Document, key: &str, numeric: bool, fonts: &RasterFonts, images: &ImageMap, job: &ExportJob) -> Result<usize> {
    use image::codecs::gif::{GifEncoder, Repeat};
    use image::{Delay, Frame};
    if job.width_px == 0 || job.height_px == 0 || job.width_px > 4000 || job.height_px > 4000 {
        bail!("an animated GIF must be between 1 and 4000 pixels per side; lower the resolution");
    }
    if doc.data.is_none() {
        bail!("load trait data first");
    }
    let Some(p) = doc.posterior.as_mut() else { bail!("open a set of trees (posterior sample or bootstrap replicates) to map a trait onto it") };
    p.summary_blocking();
    let frames = posterior_frames(p.trees.len(), p.burnin, job.every);
    let (was_browsing, was) = (p.browsing(), p.browse);
    let mut enc = GifEncoder::new_with_speed(std::io::BufWriter::new(std::fs::File::create(&job.path)?), 10);
    enc.set_repeat(Repeat::Infinite)?;
    let delay = Delay::from_numer_denom_ms((1000.0 / job.fps.max(0.1)).round() as u32, 1);
    let pt = job.dpi / 72.0 * job.scale;
    let bg = Some(doc.view.background.unwrap_or(Color::WHITE));
    let mut result = Ok(());
    for &i in &frames {
        // The sampled tree (trait data joined to its tips), the trait estimated
        // at its internal nodes, and the branches colored by it.
        doc.show_sample(i);
        if numeric {
            crate::ancestral::reconstruct_continuous(&mut doc.tree, key);
        } else {
            crate::ancestral::reconstruct_discrete(&mut doc.tree, key);
        }
        let mut layers = doc.layers.clone();
        for e in layers.iter_mut() {
            if let crate::style::Layer::Tree(s) = &mut e.layer {
                if s.color_by.as_deref() != Some(key) {
                    s.color_by = Some(key.to_string());
                    s.palette = if numeric { crate::style::Palette::Viridis } else { crate::style::Palette::OkabeIto };
                }
                s.width = s.width.max(1.5);
            }
        }
        let env = RasterEnv { fonts, images };
        let input = BuildInput { tree: &doc.tree, view: &doc.view, layers: &layers, overlay: doc.overlay() };
        let sc = scene::build(&input, &env, job.width_px as f32, job.height_px as f32, pt);
        result = raster::render(&sc, fonts, images, bg).and_then(|pm| Ok(enc.encode_frame(Frame::from_parts(raster::to_rgba(&pm), 0, 0, delay))?));
        if result.is_err() {
            break;
        }
    }
    if was_browsing {
        doc.show_sample(was);
    } else {
        doc.leave_samples();
    }
    result.map(|_| frames.len())
}

/// Save an animated GIF of `frames` scenes: `frame` draws frame `i` with a
/// text-measuring environment at a width and height (pixels) and points-to-pixels
/// scale. Size, resolution and frames per second come from `job`.
pub fn export_animation(frames: usize, frame: &dyn Fn(&dyn SceneEnv, usize, f32, f32, f32) -> Scene, background: Option<Color>, fonts: &RasterFonts, images: &ImageMap, job: &ExportJob) -> Result<()> {
    use image::codecs::gif::{GifEncoder, Repeat};
    use image::{Delay, Frame};
    if job.width_px == 0 || job.height_px == 0 || job.width_px > 4000 || job.height_px > 4000 {
        bail!("an animated GIF must be between 1 and 4000 pixels per side; lower the resolution");
    }
    let mut enc = GifEncoder::new_with_speed(std::io::BufWriter::new(std::fs::File::create(&job.path)?), 10);
    enc.set_repeat(Repeat::Infinite)?;
    let delay = Delay::from_numer_denom_ms((1000.0 / job.fps.max(0.1)).round() as u32, 1);
    let pt = job.dpi / 72.0 * job.scale;
    let bg = Some(background.unwrap_or(Color::WHITE));
    let env = RasterEnv { fonts, images };
    for i in 0..frames {
        let sc = frame(&env, i, job.width_px as f32, job.height_px as f32, pt);
        let pm = raster::render(&sc, fonts, images, bg)?;
        enc.encode_frame(Frame::from_parts(raster::to_rgba(&pm), 0, 0, delay))?;
    }
    Ok(())
}

/// Render and save a figure.
pub fn export(doc: &Document, fonts: &RasterFonts, font: &FontBytes, images: &ImageMap, job: &ExportJob) -> Result<()> {
    let input = BuildInput { tree: &doc.tree, view: &doc.view, layers: &doc.layers, overlay: doc.overlay() };
    export_scene(&|env, w, h, pt| scene::build(&input, env, w, h, pt), doc.view.background, fonts, font, images, job)
}

/// Render and save any figure: `build` draws it with a text-measuring
/// environment at a width and height (pixels) and points-to-pixels scale.
pub fn export_scene(build: &dyn Fn(&dyn SceneEnv, f32, f32, f32) -> Scene, background: Option<Color>, fonts: &RasterFonts, font: &FontBytes, images: &ImageMap, job: &ExportJob) -> Result<()> {
    if job.width_px == 0 || job.height_px == 0 || job.width_px > 30000 || job.height_px > 30000 {
        bail!("figure size must be between 1 and 30000 pixels per side");
    }
    let pt = job.dpi / 72.0 * job.scale;
    let env = RasterEnv { fonts, images };
    let bg = if job.transparent { None } else { Some(background.unwrap_or(Color::WHITE)) };
    match job.format {
        Format::Svg | Format::Pdf => {
            // Vector output is in points: 72 per inch.
            let w = job.width_px as f32 / job.dpi * 72.0;
            let h = job.height_px as f32 / job.dpi * 72.0;
            let sc = build(&env, w, h, job.scale);
            if job.format == Format::Svg {
                std::fs::write(&job.path, svg::render(&sc, images, bg, &font.name))?;
            } else {
                std::fs::write(&job.path, pdf::render(&sc, images, bg, font)?)?;
            }
        }
        Format::Png | Format::Tiff => {
            let sc = build(&env, job.width_px as f32, job.height_px as f32, pt);
            let pm = raster::render(&sc, fonts, images, bg)?;
            let img = raster::to_rgba(&pm);
            if job.format == Format::Png {
                raster::save_png(&img, &job.path, job.dpi)?;
            } else {
                raster::save_tiff(&img, &job.path, job.dpi, job.tiff_rgb && !job.transparent)?;
            }
        }
        Format::Gif => bail!("animated GIFs are saved with export_gif"),
    }
    Ok(())
}

/// (width, height) in inches with height / width = `ratio`, keeping `width`
/// unless the height would leave 0.5–40 in; then the height is clamped and
/// the width follows.
fn fit_aspect(width: f32, ratio: f32) -> (f32, f32) {
    let h = width * ratio;
    if (0.5..=40.0).contains(&h) {
        (width, h)
    } else {
        let h = h.clamp(0.5, 40.0);
        ((h / ratio).clamp(0.5, 40.0), h)
    }
}

impl ExportDialog {
    /// Resolution for the current format (GIFs have their own, for screens).
    fn dpi(&self) -> f32 {
        if self.format == Format::Gif { self.gif_dpi } else { self.dpi }
    }

    /// Open the dialog with the height set so the figure has the tree
    /// window's proportions (`view` = canvas size in points), keeping the width.
    pub fn open_for(&mut self, view: [f32; 2]) {
        self.open = true;
        self.view = view;
        self.match_window(view);
    }

    /// Height (or, past the 40 in limit, width) set to the tree window's proportions.
    fn match_window(&mut self, [w, h]: [f32; 2]) {
        if w < 1.0 || h < 1.0 {
            return;
        }
        let (width, height) = fit_aspect(self.width_in, h / w);
        self.width_in = width;
        self.height_in = height;
        self.preview = None;
    }

    /// Multiplier on sizes in points: with `match_view`, the figure's size in
    /// points over the window's (the smaller of width and height), so text
    /// and lines look as they do on screen; else 1 (true point sizes).
    fn scale(&self) -> f32 {
        let [w, h] = self.view;
        if !self.match_view || w < 1.0 || h < 1.0 {
            return 1.0;
        }
        (self.width_in * 72.0 / w).min(self.height_in * 72.0 / h)
    }

    fn px(&self) -> (u32, u32) {
        ((self.width_in * self.dpi()).round() as u32, (self.height_in * self.dpi()).round() as u32)
    }

    /// Returns a job when the user confirms.
    pub fn show(&mut self, ctx: &egui::Context, doc: &Document, fonts: &RasterFonts, images: &ImageMap) -> Option<ExportJob> {
        if !self.open {
            return None;
        }
        // Follow the window as it is resized or zoomed.
        self.view = doc.drawn_size();
        let mut job = None;
        let mut open = self.open;
        egui::Window::new("Export figure").open(&mut open).resizable(true).default_width(560.0).show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Format");
                ui.selectable_value(&mut self.format, Format::Tiff, "TIFF");
                ui.selectable_value(&mut self.format, Format::Png, "PNG");
                ui.selectable_value(&mut self.format, Format::Svg, "SVG (vector)");
                ui.selectable_value(&mut self.format, Format::Pdf, "PDF (vector)");
                if doc.posterior.is_some() {
                    ui.selectable_value(&mut self.format, Format::Gif, "GIF (animated)").on_hover_text("One frame per tree of the set, as Step through trees shows them");
                } else if self.format == Format::Gif {
                    self.format = Format::Png;
                }
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut self.mm, "Millimetres");
                let (f, unit) = if self.mm { (25.4, "mm") } else { (1.0, "in") };
                let mut w = self.width_in * f;
                let mut h = self.height_in * f;
                ui.add(egui::DragValue::new(&mut w).speed(0.05 * f).range(0.5 * f..=40.0 * f).suffix(unit).prefix("W "));
                ui.add(egui::DragValue::new(&mut h).speed(0.05 * f).range(0.5 * f..=40.0 * f).suffix(unit).prefix("H "));
                self.width_in = w / f;
                self.height_in = h / f;
            });
            ui.horizontal_wrapped(|ui| {
                ui.label("Presets:");
                // Column widths keep the current proportions.
                for (name, w) in [("1 column (85 mm)", 85.0 / 25.4), ("1.5 col (114 mm)", 114.0 / 25.4), ("2 columns (174 mm)", 174.0 / 25.4)] {
                    if ui.small_button(name).clicked() {
                        let ratio = self.height_in / self.width_in;
                        (self.width_in, self.height_in) = fit_aspect(w, ratio);
                    }
                }
            });
            if self.format == Format::Gif {
                let n = doc.posterior.as_ref().map(|p| p.trees.len()).unwrap_or(0);
                ui.horizontal(|ui| {
                    ui.label("Resolution");
                    for d in [72.0, 100.0, 150.0] {
                        ui.selectable_value(&mut self.gif_dpi, d, format!("{} dpi", d));
                    }
                    ui.add(egui::DragValue::new(&mut self.gif_dpi).range(36.0..=300.0).suffix(" dpi"));
                });
                ui.horizontal(|ui| {
                    ui.add(egui::DragValue::new(&mut self.gif_fps).range(0.5..=30.0).speed(0.1).suffix(" trees/s"));
                    ui.label("Every");
                    ui.add(egui::DragValue::new(&mut self.gif_every).range(1..=n.max(1)));
                    ui.label(format!("th tree: {} frames", gif_frames(n, self.gif_every).len()));
                });
            } else if !self.format.vector() {
                ui.horizontal(|ui| {
                    ui.label("Resolution");
                    for d in [300.0, 600.0, 1200.0] {
                        ui.selectable_value(&mut self.dpi, d, format!("{} dpi", d));
                    }
                    ui.add(egui::DragValue::new(&mut self.dpi).range(72.0..=2400.0).suffix(" dpi"));
                });
            }
            if self.format != Format::Gif {
                ui.checkbox(&mut self.transparent, "Transparent background");
            }
            if self.format == Format::Tiff {
                ui.checkbox(&mut self.tiff_rgb, "RGB without alpha (journal-safe)");
            }
            ui.checkbox(&mut self.match_view, "Sizes as in the tree window")
                .on_hover_text("Text, lines and points take the same share of the figure as they do on screen. Untick for true point sizes (a 9 pt label prints at 9 pt).");
            let scale = self.scale();
            let sizes = if self.match_view { format!("sizes at {:.0}% of their points, as in the window", scale * 100.0) } else { "text sizes are in points".to_string() };
            let (pw, ph) = self.px();
            if self.format.vector() {
                ui.label(format!("{:.0} × {:.0} mm page  ·  vector: lines and text stay sharp at any zoom  ·  {}", self.width_in * 25.4, self.height_in * 25.4, sizes));
            } else {
                ui.label(format!(
                    "{} × {} px  ·  {:.1} MP  ·  {} ({} px per pt)",
                    pw,
                    ph,
                    pw as f32 * ph as f32 / 1e6,
                    sizes,
                    format_args!("{:.2}", self.dpi() / 72.0 * scale)
                ));
            }
            ui.separator();

            // Live preview at low resolution with the same physical layout.
            let refresh = self.last_preview.map(|t| t.elapsed().as_millis() > 400).unwrap_or(true);
            if refresh {
                let pdpi = (420.0 / self.width_in.max(self.height_in)).min(110.0);
                let env = RasterEnv { fonts, images };
                let input = BuildInput { tree: &doc.tree, view: &doc.view, layers: &doc.layers, overlay: doc.overlay() };
                let (w, h) = (self.width_in * pdpi, self.height_in * pdpi);
                let sc = scene::build(&input, &env, w, h, pdpi / 72.0 * scale);
                if let Ok(pm) = raster::render(&sc, fonts, images, Some(Color::WHITE)) {
                    let img = raster::to_rgba(&pm);
                    let ci = egui::ColorImage::from_rgba_unmultiplied([img.width() as usize, img.height() as usize], img.as_raw());
                    self.preview = Some(ctx.load_texture("export-preview", ci, egui::TextureOptions::LINEAR));
                }
                self.last_preview = Some(Instant::now());
                ctx.request_repaint_after(std::time::Duration::from_millis(450));
            }
            if let Some(t) = &self.preview {
                let size = t.size_vec2();
                let fit = (520.0 / size.x).min(420.0 / size.y).min(1.0);
                ui.add(egui::Image::new((t.id(), size * fit)).bg_fill(egui::Color32::WHITE));
            }
            ui.separator();
            if ui.button("Export…").clicked() {
                let default = format!("{}.{}", sanitize(&doc.name), self.format.ext());
                let mut dlg = rfd::FileDialog::new().set_file_name(&default);
                dlg = match self.format {
                    Format::Png => dlg.add_filter("PNG image", &["png"]),
                    Format::Tiff => dlg.add_filter("TIFF image", &["tif", "tiff"]),
                    Format::Svg => dlg.add_filter("SVG", &["svg"]),
                    Format::Pdf => dlg.add_filter("PDF", &["pdf"]),
                    Format::Gif => dlg.add_filter("Animated GIF", &["gif"]),
                };
                if let Some(mut path) = dlg.save_file() {
                    if path.extension().is_none() {
                        path.set_extension(self.format.ext());
                    }
                    job = Some(ExportJob {
                        path,
                        format: self.format,
                        width_px: pw,
                        height_px: ph,
                        dpi: self.dpi(),
                        transparent: self.transparent && self.format != Format::Gif,
                        tiff_rgb: self.tiff_rgb,
                        fps: self.gif_fps,
                        every: self.gif_every,
                        scale,
                    });
                }
            }
        });
        self.open = open;
        job
    }
}

fn sanitize(s: &str) -> String {
    let stem = Path::new(s).file_stem().map(|x| x.to_string_lossy().into_owned()).unwrap_or_else(|| s.to_string());
    stem.chars().map(|c| if c.is_alphanumeric() || c == '-' || c == '_' { c } else { '_' }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn export_size_follows_window_proportions() {
        // A window twice as tall as wide: 6 in wide becomes 12 in tall.
        let mut d = ExportDialog { width_in: 6.0, ..Default::default() };
        d.open_for([400.0, 800.0]);
        assert!(d.open && (d.width_in, d.height_in) == (6.0, 12.0));
        // Too tall for 40 in: the height stops there and the width shrinks.
        assert_eq!(fit_aspect(6.0, 10.0), (4.0, 40.0));
        // Sizes as in the window: 6 in = 432 pt across a 400-pt window scales
        // points by 1.08, so the export is the window scene scaled up.
        assert!((d.scale() - 432.0 / 400.0).abs() < 1e-5);
        d.match_view = false;
        assert_eq!(d.scale(), 1.0);
        // No canvas yet: the size is left alone.
        let mut d = ExportDialog::default();
        d.open_for([0.0, 0.0]);
        assert_eq!((d.width_in, d.height_in), (6.85, 8.0));
    }

    #[test]
    fn animated_gif_has_a_frame_per_tree() {
        use image::AnimationDecoder;
        let trees = crate::io::read_trees_str("((A:1,B:1):1,C:2);\n((A:1,C:1):1,B:2);\n((B:1,C:1):1,A:2);\n").unwrap();
        let mut doc = Document::new("set", trees[0].clone(), false);
        doc.posterior = Some(super::super::document::Posterior::new(trees));
        let font = FontBytes::load();
        let fonts = RasterFonts::new(&font).unwrap();
        let path = std::env::temp_dir().join(format!("canopy-test-{}.gif", std::process::id()));
        let job = ExportJob { path: path.clone(), format: Format::Gif, width_px: 200, height_px: 150, dpi: 72.0, transparent: false, tiff_rgb: true, fps: 2.0, every: 1, scale: 1.0 };
        assert_eq!(export_gif(&mut doc, &fonts, &ImageMap::new(), &job).unwrap(), 3);
        let frames = image::codecs::gif::GifDecoder::new(std::io::BufReader::new(std::fs::File::open(&path).unwrap())).unwrap().into_frames().collect_frames().unwrap();
        let _ = std::fs::remove_file(&path);
        assert_eq!(frames.len(), 3);
        assert_ne!(frames[0].buffer(), frames[1].buffer(), "successive trees should look different");
        // The document goes back to the tree it showed before.
        assert!(!doc.posterior.as_ref().unwrap().browsing());
        assert_eq!(gif_frames(10, 3), [0, 3, 6, 9]);
    }

    #[test]
    fn trait_mapped_on_each_posterior_tree() {
        use image::AnimationDecoder;
        let trees = crate::io::read_trees_str("((A:1,B:1):1,(C:1,D:1):1);\n((A:1,C:1):1,(B:1,D:1):1);\n((A:1,D:1):1,(B:1,C:1):1);\n((A:1,B:1):1,(C:1,D:1):1);\n").unwrap();
        let mut doc = Document::new("set", trees[0].clone(), false);
        let mut p = super::super::document::Posterior::new(trees);
        p.burnin = 0.25;
        doc.posterior = Some(p);
        doc.data = Some(crate::io::data::DataTable::from_str("taxon,mass,diet\nA,1,x\nB,2,x\nC,8,y\nD,9,y\n", "t", None).unwrap());
        let font = FontBytes::load();
        let fonts = RasterFonts::new(&font).unwrap();
        for (key, numeric) in [("mass", true), ("diet", false)] {
            let path = std::env::temp_dir().join(format!("canopy-trait-{}-{}.gif", key, std::process::id()));
            let job = ExportJob { path: path.clone(), format: Format::Gif, width_px: 200, height_px: 150, dpi: 72.0, transparent: false, tiff_rgb: true, fps: 2.0, every: 1, scale: 1.0 };
            // Burn-in leaves three of the four trees.
            assert_eq!(export_trait_gif(&mut doc, key, numeric, &fonts, &ImageMap::new(), &job).unwrap(), 3);
            let frames = image::codecs::gif::GifDecoder::new(std::io::BufReader::new(std::fs::File::open(&path).unwrap())).unwrap().into_frames().collect_frames().unwrap();
            let _ = std::fs::remove_file(&path);
            assert_eq!(frames.len(), 3);
            assert_ne!(frames[0].buffer(), frames[1].buffer(), "different trees, different reconstructions");
        }
        // Back to the tab's own tree, its layers unchanged (no trait coloring).
        assert!(!doc.posterior.as_ref().unwrap().browsing());
        assert!(doc.layers.iter().all(|e| !matches!(&e.layer, crate::style::Layer::Tree(s) if s.color_by.is_some())));
        assert_eq!(posterior_frames(10, 0.25, 3), [2, 5, 8]);
    }
}
