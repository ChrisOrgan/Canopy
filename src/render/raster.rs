//! Anti-aliased raster rendering with tiny-skia, and PNG/TIFF writers that
//! embed DPI metadata (pHYs chunk / TIFF resolution tags).

use super::fonts::FontBytes;
use super::ImageMap;
use crate::scene::{Prim, Scene, SceneEnv, P};
use crate::style::Color;
use ab_glyph::{Font, FontVec, PxScale, ScaleFont};
use anyhow::{Context, Result};
use image::RgbaImage;
use std::path::Path;
use tiny_skia::{FillRule, Paint, PathBuilder, Pixmap, PixmapPaint, Stroke, StrokeDash, Transform};

pub struct RasterFonts {
    regular: FontVec,
    italic: Option<FontVec>,
    bold: Option<FontVec>,
}

impl RasterFonts {
    pub fn new(bytes: &FontBytes) -> Result<Self> {
        let regular = FontVec::try_from_vec(bytes.regular.clone()).context("loading font")?;
        let italic = bytes.italic.clone().and_then(|b| FontVec::try_from_vec(b).ok());
        let bold = bytes.bold.clone().and_then(|b| FontVec::try_from_vec(b).ok());
        Ok(RasterFonts { regular, italic, bold })
    }

    fn pick(&self, italic: bool, bold: bool) -> (&FontVec, bool) {
        if bold {
            if let Some(b) = &self.bold {
                return (b, false);
            }
        }
        if italic {
            if let Some(i) = &self.italic {
                return (i, false);
            }
            return (&self.regular, true);
        }
        (&self.regular, false)
    }

    pub fn width(&self, text: &str, size: f32, italic: bool, bold: bool) -> f32 {
        let (font, _) = self.pick(italic, bold);
        let sf = font.as_scaled(PxScale::from(size));
        let mut w = 0.0;
        let mut prev = None;
        for ch in text.chars() {
            let id = sf.glyph_id(ch);
            if let Some(p) = prev {
                w += sf.kern(p, id);
            }
            w += sf.h_advance(id);
            prev = Some(id);
        }
        w
    }
}

pub struct RasterEnv<'a> {
    pub fonts: &'a RasterFonts,
    pub images: &'a ImageMap,
}

impl SceneEnv for RasterEnv<'_> {
    fn text_width(&self, text: &str, size: f32, italic: bool) -> f32 {
        self.fonts.width(text, size, italic, false)
    }
    fn image_aspect(&self, key: &str) -> Option<f32> {
        self.images.get(key).map(|i| i.width() as f32 / i.height().max(1) as f32)
    }
}

fn paint(c: Color) -> Paint<'static> {
    let mut p = Paint::default();
    p.set_color_rgba8(c.r, c.g, c.b, c.a);
    p.anti_alias = true;
    p
}

fn polyline(pts: &[P], close: bool) -> Option<tiny_skia::Path> {
    let mut pb = PathBuilder::new();
    let first = pts.first()?;
    pb.move_to(first[0], first[1]);
    for p in &pts[1..] {
        pb.line_to(p[0], p[1]);
    }
    if close {
        pb.close();
    }
    pb.finish()
}

fn stroke(width: f32, dashed: bool) -> Stroke {
    let mut s = Stroke { width: width.max(0.1), line_cap: tiny_skia::LineCap::Butt, line_join: tiny_skia::LineJoin::Round, ..Default::default() };
    if dashed {
        s.dash = StrokeDash::new(vec![width * 2.0 + 1.0, width * 2.0 + 1.0], 0.0);
        s.line_cap = tiny_skia::LineCap::Butt;
    }
    s
}

#[allow(clippy::too_many_arguments)]
fn draw_text(pm: &mut Pixmap, fonts: &RasterFonts, pos: P, text: &str, size: f32, color: Color, angle: f32, halign: f32, valign: f32, italic: bool, bold: bool) {
    if text.is_empty() || size < 0.5 {
        return;
    }
    let (font, synth_italic) = fonts.pick(italic, bold);
    let scale = PxScale::from(size);
    let sf = font.as_scaled(scale);
    let ascent = sf.ascent();
    let descent = sf.descent();
    let h = ascent - descent;
    let pad = 2.0f32 + if synth_italic { size * 0.25 } else { 0.0 };
    let w = fonts.width(text, size, italic, bold);
    let bw = (w + 2.0 * pad).ceil() as u32 + 1;
    let bh = (h + 2.0 * pad).ceil() as u32 + 1;
    let Some(mut buf) = Pixmap::new(bw, bh) else { return };
    let mut cov = vec![0f32; (bw * bh) as usize];
    let mut caret = 0.0;
    let mut prev = None;
    for ch in text.chars() {
        let id = sf.glyph_id(ch);
        if let Some(p) = prev {
            caret += sf.kern(p, id);
        }
        let g = id.with_scale_and_position(scale, ab_glyph::point(pad + caret, pad + ascent));
        caret += sf.h_advance(id);
        prev = Some(id);
        if let Some(og) = font.outline_glyph(g) {
            let b = og.px_bounds();
            og.draw(|x, y, c| {
                let px = b.min.x as i32 + x as i32;
                let py = b.min.y as i32 + y as i32;
                if px >= 0 && py >= 0 && (px as u32) < bw && (py as u32) < bh {
                    let i = (py as u32 * bw + px as u32) as usize;
                    cov[i] = (cov[i] + c).min(1.0);
                }
            });
        }
    }
    let data = buf.data_mut();
    let a = color.a as f32 / 255.0;
    for (i, &c) in cov.iter().enumerate() {
        if c > 0.0 {
            let al = c * a;
            data[i * 4] = (color.r as f32 * al).round() as u8;
            data[i * 4 + 1] = (color.g as f32 * al).round() as u8;
            data[i * 4 + 2] = (color.b as f32 * al).round() as u8;
            data[i * 4 + 3] = (al * 255.0).round() as u8;
        }
    }
    let mut t = Transform::from_translate(pos[0], pos[1]).pre_rotate(angle.to_degrees()).pre_translate(-halign * w - pad, -valign * h - pad);
    if synth_italic {
        // Shear about the baseline.
        let base = pad + ascent;
        t = t.pre_translate(0.0, base).pre_concat(Transform::from_row(1.0, 0.0, -0.2, 1.0, 0.0, 0.0)).pre_translate(0.0, -base);
    }
    let pp = PixmapPaint { quality: tiny_skia::FilterQuality::Bilinear, ..Default::default() };
    pm.draw_pixmap(0, 0, buf.as_ref(), &pp, t, None);
}

fn draw_image(pm: &mut Pixmap, img: &RgbaImage, center: P, height: f32, tint: Option<Color>) {
    let (w, h) = img.dimensions();
    let Some(mut src) = Pixmap::new(w, h) else { return };
    let data = src.data_mut();
    for (i, p) in img.pixels().enumerate() {
        let [r, g, b, a] = p.0;
        let (r, g, b) = match tint {
            Some(c) => (c.r, c.g, c.b),
            None => (r, g, b),
        };
        let af = a as f32 / 255.0;
        data[i * 4] = (r as f32 * af).round() as u8;
        data[i * 4 + 1] = (g as f32 * af).round() as u8;
        data[i * 4 + 2] = (b as f32 * af).round() as u8;
        data[i * 4 + 3] = a;
    }
    let s = height / h as f32;
    let dw = w as f32 * s;
    let t = Transform::from_translate(center[0] - dw / 2.0, center[1] - height / 2.0).pre_scale(s, s);
    let pp = PixmapPaint { quality: tiny_skia::FilterQuality::Bicubic, ..Default::default() };
    pm.draw_pixmap(0, 0, src.as_ref(), &pp, t, None);
}

/// Render a scene. `background = None` gives a transparent image.
pub fn render(scene: &Scene, fonts: &RasterFonts, images: &ImageMap, background: Option<Color>) -> Result<Pixmap> {
    let w = scene.width.ceil().max(1.0) as u32;
    let h = scene.height.ceil().max(1.0) as u32;
    let mut pm = Pixmap::new(w, h).context("image too large")?;
    if let Some(bg) = background {
        pm.fill(tiny_skia::Color::from_rgba8(bg.r, bg.g, bg.b, bg.a));
    }
    let id = Transform::identity();
    for prim in &scene.prims {
        match prim {
            Prim::Path { pts, color, width, dashed } => {
                if let Some(p) = polyline(pts, false) {
                    pm.stroke_path(&p, &paint(*color), &stroke(*width, *dashed), id, None);
                }
            }
            Prim::Poly { pts, fill, stroke: st } => {
                if let Some(p) = polyline(pts, true) {
                    pm.fill_path(&p, &paint(*fill), FillRule::Winding, id, None);
                    if let Some((c, w)) = st {
                        pm.stroke_path(&p, &paint(*c), &stroke(*w, false), id, None);
                    }
                }
            }
            Prim::Strip { a, b, fill, stroke: st } => {
                let mut pts = a.clone();
                pts.extend(b.iter().rev());
                if let Some(p) = polyline(&pts, true) {
                    pm.fill_path(&p, &paint(*fill), FillRule::Winding, id, None);
                    if let Some((c, w)) = st {
                        pm.stroke_path(&p, &paint(*c), &stroke(*w, false), id, None);
                    }
                }
            }
            Prim::Circle { c, r, fill, stroke: st } => {
                if let Some(p) = PathBuilder::from_circle(c[0], c[1], r.max(0.1)) {
                    pm.fill_path(&p, &paint(*fill), FillRule::Winding, id, None);
                    if let Some((sc, w)) = st {
                        pm.stroke_path(&p, &paint(*sc), &stroke(*w, false), id, None);
                    }
                }
            }
            Prim::Text { pos, text, size, color, angle, halign, valign, italic, bold } => {
                draw_text(&mut pm, fonts, *pos, text, *size, *color, *angle, *halign, *valign, *italic, *bold)
            }
            Prim::Image { key, center, height, tint } => {
                if let Some(img) = images.get(key) {
                    draw_image(&mut pm, img, *center, *height, *tint);
                }
            }
        }
    }
    Ok(pm)
}

/// Convert a premultiplied pixmap into straight-alpha RGBA.
pub fn to_rgba(pm: &Pixmap) -> RgbaImage {
    let mut img = RgbaImage::new(pm.width(), pm.height());
    for (dst, src) in img.pixels_mut().zip(pm.pixels()) {
        let c = src.demultiply();
        *dst = image::Rgba([c.red(), c.green(), c.blue(), c.alpha()]);
    }
    img
}

pub fn save_png(img: &RgbaImage, path: &Path, dpi: f32) -> Result<()> {
    let file = std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let w = std::io::BufWriter::new(file);
    let mut enc = png::Encoder::new(w, img.width(), img.height());
    enc.set_color(png::ColorType::Rgba);
    enc.set_depth(png::BitDepth::Eight);
    let ppm = (dpi / 0.0254).round() as u32;
    enc.set_pixel_dims(Some(png::PixelDimensions { xppu: ppm, yppu: ppm, unit: png::Unit::Meter }));
    let mut writer = enc.write_header()?;
    writer.write_image_data(img.as_raw())?;
    Ok(())
}

/// Write a TIFF with LZW compression and resolution tags. `rgb` drops the
/// alpha channel (some journals reject RGBA TIFFs).
pub fn save_tiff(img: &RgbaImage, path: &Path, dpi: f32, rgb: bool) -> Result<()> {
    use tiff::encoder::{colortype, compression::Lzw, Rational, TiffEncoder};
    use tiff::tags::ResolutionUnit;
    let file = std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut enc = TiffEncoder::new(std::io::BufWriter::new(file))?;
    let res = Rational { n: (dpi * 100.0).round() as u32, d: 100 };
    if rgb {
        let data: Vec<u8> = img
            .pixels()
            .flat_map(|p| {
                // Composite onto white.
                let a = p.0[3] as f32 / 255.0;
                let f = |c: u8| (c as f32 * a + 255.0 * (1.0 - a)).round() as u8;
                [f(p.0[0]), f(p.0[1]), f(p.0[2])]
            })
            .collect();
        let mut im = enc.new_image_with_compression::<colortype::RGB8, _>(img.width(), img.height(), Lzw)?;
        im.resolution(ResolutionUnit::Inch, res);
        im.write_data(&data)?;
    } else {
        let mut im = enc.new_image_with_compression::<colortype::RGBA8, _>(img.width(), img.height(), Lzw)?;
        im.resolution(ResolutionUnit::Inch, res);
        im.write_data(img.as_raw())?;
    }
    Ok(())
}
