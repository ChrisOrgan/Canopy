//! Rasterizes the logo SVG. Shared by the app (`logo.rs`) and `build.rs`
//! (which turns it into the Windows .exe icon), so it depends only on resvg.

use resvg::{tiny_skia, usvg};

/// Render `svg` into a `size`×`size` straight-alpha RGBA buffer, centered and
/// scaled to fit. With `tile`, the mark sits on a white rounded square.
pub fn render_rgba(svg: &str, size: u32, tile: bool) -> Vec<u8> {
    let tree = usvg::Tree::from_str(svg, &usvg::Options::default()).expect("assets/logo.svg is not valid SVG");
    let mut pm = tiny_skia::Pixmap::new(size, size).expect("icon size");
    let s = size as f32;
    if tile {
        let r = s * 0.22;
        let mut pb = tiny_skia::PathBuilder::new();
        pb.move_to(r, 0.0);
        pb.line_to(s - r, 0.0);
        pb.quad_to(s, 0.0, s, r);
        pb.line_to(s, s - r);
        pb.quad_to(s, s, s - r, s);
        pb.line_to(r, s);
        pb.quad_to(0.0, s, 0.0, s - r);
        pb.line_to(0.0, r);
        pb.quad_to(0.0, 0.0, r, 0.0);
        pb.close();
        if let Some(p) = pb.finish() {
            let mut paint = tiny_skia::Paint::default();
            paint.set_color_rgba8(255, 255, 255, 255);
            paint.anti_alias = true;
            pm.fill_path(&p, &paint, tiny_skia::FillRule::Winding, tiny_skia::Transform::identity(), None);
        }
    }
    let avail = s * if tile { 0.8 } else { 1.0 };
    let (w, h) = (tree.size().width(), tree.size().height());
    let k = avail / w.max(h);
    let t = tiny_skia::Transform::from_row(k, 0.0, 0.0, k, (s - w * k) / 2.0, (s - h * k) / 2.0);
    resvg::render(&tree, t, &mut pm.as_mut());
    let mut out = Vec::with_capacity((size * size * 4) as usize);
    for px in pm.pixels() {
        let c = px.demultiply();
        out.extend_from_slice(&[c.red(), c.green(), c.blue(), c.alpha()]);
    }
    out
}
