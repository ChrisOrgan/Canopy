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
    let pt = job.dpi / 72.0;
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
    let pt = job.dpi / 72.0;
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
    let pt = job.dpi / 72.0;
    let env = RasterEnv { fonts, images };
    let bg = if job.transparent { None } else { Some(background.unwrap_or(Color::WHITE)) };
    match job.format {
        Format::Svg | Format::Pdf => {
            // Vector output is in points: 72 per inch.
            let w = job.width_px as f32 / job.dpi * 72.0;
            let h = job.height_px as f32 / job.dpi * 72.0;
            let sc = build(&env, w, h, 1.0);
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

impl ExportDialog {
    /// Resolution for the current format (GIFs have their own, for screens).
    fn dpi(&self) -> f32 {
        if self.format == Format::Gif { self.gif_dpi } else { self.dpi }
    }

    fn px(&self) -> (u32, u32) {
        ((self.width_in * self.dpi()).round() as u32, (self.height_in * self.dpi()).round() as u32)
    }

    /// Returns a job when the user confirms.
    pub fn show(&mut self, ctx: &egui::Context, doc: &Document, fonts: &RasterFonts, images: &ImageMap) -> Option<ExportJob> {
        if !self.open {
            return None;
        }
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
                for (name, w) in [("1 column (85 mm)", 85.0 / 25.4), ("1.5 col (114 mm)", 114.0 / 25.4), ("2 columns (174 mm)", 174.0 / 25.4)] {
                    if ui.small_button(name).clicked() {
                        self.width_in = w;
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
            let (pw, ph) = self.px();
            if self.format.vector() {
                ui.label(format!(
                    "{:.0} × {:.0} mm page  ·  vector: lines and text stay sharp at any zoom  ·  text sizes are in points",
                    self.width_in * 25.4,
                    self.height_in * 25.4
                ));
            } else {
                ui.label(format!(
                    "{} × {} px  ·  {:.1} MP  ·  text sizes are in points ({} px per pt)",
                    pw,
                    ph,
                    pw as f32 * ph as f32 / 1e6,
                    format_args!("{:.2}", self.dpi() / 72.0)
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
                let sc = scene::build(&input, &env, w, h, pdpi / 72.0);
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
                let scale = (520.0 / size.x).min(420.0 / size.y).min(1.0);
                ui.add(egui::Image::new((t.id(), size * scale)).bg_fill(egui::Color32::WHITE));
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
    fn animated_gif_has_a_frame_per_tree() {
        use image::AnimationDecoder;
        let trees = crate::io::read_trees_str("((A:1,B:1):1,C:2);\n((A:1,C:1):1,B:2);\n((B:1,C:1):1,A:2);\n").unwrap();
        let mut doc = Document::new("set", trees[0].clone(), false);
        doc.posterior = Some(super::super::document::Posterior::new(trees));
        let font = FontBytes::load();
        let fonts = RasterFonts::new(&font).unwrap();
        let path = std::env::temp_dir().join(format!("canopy-test-{}.gif", std::process::id()));
        let job = ExportJob { path: path.clone(), format: Format::Gif, width_px: 200, height_px: 150, dpi: 72.0, transparent: false, tiff_rgb: true, fps: 2.0, every: 1 };
        assert_eq!(export_gif(&mut doc, &fonts, &ImageMap::new(), &job).unwrap(), 3);
        let frames = image::codecs::gif::GifDecoder::new(std::io::BufReader::new(std::fs::File::open(&path).unwrap())).unwrap().into_frames().collect_frames().unwrap();
        let _ = std::fs::remove_file(&path);
        assert_eq!(frames.len(), 3);
        assert_ne!(frames[0].buffer(), frames[1].buffer(), "successive trees should look different");
        // The document goes back to the tree it showed before.
        assert!(!doc.posterior.as_ref().unwrap().browsing());
        assert_eq!(gif_frames(10, 3), [0, 3, 6, 9]);
    }
}
