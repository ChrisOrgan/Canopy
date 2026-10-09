//! SVG export (vector; editable in Illustrator/Inkscape).

use super::ImageMap;
use crate::scene::{Prim, Scene, P};
use crate::style::Color;
use base64::Engine;
use std::fmt::Write;

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

fn rgb(c: Color) -> String {
    format!("rgb({},{},{})", c.r, c.g, c.b)
}

fn op(c: Color) -> String {
    if c.a == 255 {
        String::new()
    } else {
        format!(" {}-opacity=\"{:.3}\"", "OPKIND", c.a as f32 / 255.0)
    }
}

fn fill_attr(c: Color) -> String {
    format!("fill=\"{}\"{}", rgb(c), op(c).replace("OPKIND", "fill"))
}

fn stroke_attr(c: Color, w: f32) -> String {
    format!("stroke=\"{}\" stroke-width=\"{:.2}\"{}", rgb(c), w, op(c).replace("OPKIND", "stroke"))
}

fn pts_str(pts: &[P]) -> String {
    pts.iter().map(|p| format!("{:.2},{:.2}", p[0], p[1])).collect::<Vec<_>>().join(" ")
}

/// `font_family` is written into the SVG (e.g. "Arial").
pub fn render(scene: &Scene, images: &ImageMap, background: Option<Color>, font_family: &str) -> String {
    let mut s = String::new();
    let _ = writeln!(
        s,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w:.1}\" height=\"{h:.1}\" viewBox=\"0 0 {w:.1} {h:.1}\">",
        w = scene.width,
        h = scene.height
    );
    if let Some(bg) = background {
        let _ = writeln!(s, "<rect width=\"100%\" height=\"100%\" {}/>", fill_attr(bg));
    }
    let mut encoded: std::collections::HashMap<&str, String> = std::collections::HashMap::new();
    for prim in &scene.prims {
        match prim {
            Prim::Path { pts, color, width, dashed } => {
                let dash = if *dashed { format!(" stroke-dasharray=\"{:.1},{:.1}\"", width * 2.0 + 1.0, width * 2.0 + 1.0) } else { String::new() };
                let _ = writeln!(s, "<polyline points=\"{}\" fill=\"none\" {} stroke-linejoin=\"round\"{}/>", pts_str(pts), stroke_attr(*color, *width), dash);
            }
            Prim::Poly { pts, fill, stroke } => {
                let st = stroke.map(|(c, w)| format!(" {}", stroke_attr(c, w))).unwrap_or_default();
                let _ = writeln!(s, "<polygon points=\"{}\" {}{}/>", pts_str(pts), fill_attr(*fill), st);
            }
            Prim::Strip { a, b, fill, stroke } => {
                let mut pts = a.clone();
                pts.extend(b.iter().rev());
                let st = stroke.map(|(c, w)| format!(" {}", stroke_attr(c, w))).unwrap_or_default();
                let _ = writeln!(s, "<polygon points=\"{}\" {}{}/>", pts_str(&pts), fill_attr(*fill), st);
            }
            Prim::Circle { c, r, fill, stroke } => {
                let st = stroke.map(|(c, w)| format!(" {}", stroke_attr(c, w))).unwrap_or_default();
                let _ = writeln!(s, "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" {}{}/>", c[0], c[1], r, fill_attr(*fill), st);
            }
            Prim::Text { pos, text, size, color, angle, halign, valign, italic, bold } => {
                let anchor = if *halign < 0.25 {
                    "start"
                } else if *halign > 0.75 {
                    "end"
                } else {
                    "middle"
                };
                // Approximate vertical alignment with a baseline shift in ems.
                let dy = 0.75 - 0.95 * valign;
                let _ = writeln!(
                    s,
                    "<text x=\"0\" y=\"0\" dy=\"{:.2}em\" transform=\"translate({:.2},{:.2}) rotate({:.2})\" font-family=\"{}\" font-size=\"{:.2}\" text-anchor=\"{}\"{}{} {}>{}</text>",
                    dy,
                    pos[0],
                    pos[1],
                    angle.to_degrees(),
                    esc(font_family),
                    size,
                    anchor,
                    if *italic { " font-style=\"italic\"" } else { "" },
                    if *bold { " font-weight=\"bold\"" } else { "" },
                    fill_attr(*color),
                    esc(text)
                );
            }
            Prim::Image { key, center, height, tint } => {
                let Some(img) = images.get(key) else { continue };
                let data = encoded.entry(key.as_str()).or_insert_with(|| {
                    let mut img = (**img).clone();
                    if let Some(c) = tint {
                        for p in img.pixels_mut() {
                            p.0[0] = c.r;
                            p.0[1] = c.g;
                            p.0[2] = c.b;
                        }
                    }
                    let mut buf = std::io::Cursor::new(Vec::new());
                    let _ = img.write_to(&mut buf, image::ImageFormat::Png);
                    base64::engine::general_purpose::STANDARD.encode(buf.into_inner())
                });
                let w = height * img.width() as f32 / img.height().max(1) as f32;
                let _ = writeln!(
                    s,
                    "<image x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" href=\"data:image/png;base64,{}\"/>",
                    center[0] - w / 2.0,
                    center[1] - height / 2.0,
                    w,
                    height,
                    data
                );
            }
        }
    }
    s.push_str("</svg>\n");
    s
}
