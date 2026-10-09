//! Figure export dialog: PNG / TIFF (with DPI metadata), SVG and PDF.

use super::document::Document;
use crate::render::fonts::FontBytes;
use crate::render::raster::{self, RasterEnv, RasterFonts};
use crate::render::{pdf, svg, ImageMap};
use crate::scene::{self, BuildInput};
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
}

impl Format {
    fn ext(self) -> &'static str {
        match self {
            Format::Png => "png",
            Format::Tiff => "tif",
            Format::Svg => "svg",
            Format::Pdf => "pdf",
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
}

/// Render and save a figure.
pub fn export(doc: &Document, fonts: &RasterFonts, font: &FontBytes, images: &ImageMap, job: &ExportJob) -> Result<()> {
    if job.width_px == 0 || job.height_px == 0 || job.width_px > 30000 || job.height_px > 30000 {
        bail!("figure size must be between 1 and 30000 pixels per side");
    }
    let pt = job.dpi / 72.0;
    let env = RasterEnv { fonts, images };
    let input = BuildInput { tree: &doc.tree, view: &doc.view, layers: &doc.layers, overlay: doc.overlay() };
    let bg = if job.transparent { None } else { Some(doc.view.background.unwrap_or(Color::WHITE)) };
    match job.format {
        Format::Svg | Format::Pdf => {
            // Vector output is in points: 72 per inch.
            let w = job.width_px as f32 / job.dpi * 72.0;
            let h = job.height_px as f32 / job.dpi * 72.0;
            let sc = scene::build(&input, &env, w, h, 1.0);
            if job.format == Format::Svg {
                std::fs::write(&job.path, svg::render(&sc, images, bg, &font.name))?;
            } else {
                std::fs::write(&job.path, pdf::render(&sc, images, bg, font)?)?;
            }
        }
        Format::Png | Format::Tiff => {
            let sc = scene::build(&input, &env, job.width_px as f32, job.height_px as f32, pt);
            let pm = raster::render(&sc, fonts, images, bg)?;
            let img = raster::to_rgba(&pm);
            if job.format == Format::Png {
                raster::save_png(&img, &job.path, job.dpi)?;
            } else {
                raster::save_tiff(&img, &job.path, job.dpi, job.tiff_rgb && !job.transparent)?;
            }
        }
    }
    Ok(())
}

impl ExportDialog {
    fn px(&self) -> (u32, u32) {
        ((self.width_in * self.dpi).round() as u32, (self.height_in * self.dpi).round() as u32)
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
            if !self.format.vector() {
                ui.horizontal(|ui| {
                    ui.label("Resolution");
                    for d in [300.0, 600.0, 1200.0] {
                        ui.selectable_value(&mut self.dpi, d, format!("{} dpi", d));
                    }
                    ui.add(egui::DragValue::new(&mut self.dpi).range(72.0..=2400.0).suffix(" dpi"));
                });
            }
            ui.checkbox(&mut self.transparent, "Transparent background");
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
                    format_args!("{:.2}", self.dpi / 72.0)
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
                        dpi: self.dpi,
                        transparent: self.transparent,
                        tiff_rgb: self.tiff_rgb,
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
