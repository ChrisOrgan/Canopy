// Builds the Windows .exe icon from assets/logo.svg (the single logo source)
// and embeds it. Nothing to run by hand: edit the SVG and rebuild.

#[path = "src/logo_raster.rs"]
mod logo_raster;

use std::path::PathBuf;

/// Write an .ico containing PNG-compressed frames.
fn write_ico(svg: &str, path: &std::path::Path) -> std::io::Result<()> {
    let sizes = [16u32, 24, 32, 48, 64, 128, 256];
    let mut pngs = Vec::new();
    for &s in &sizes {
        let rgba = logo_raster::render_rgba(svg, s, true);
        let mut buf = Vec::new();
        {
            let mut enc = png::Encoder::new(&mut buf, s, s);
            enc.set_color(png::ColorType::Rgba);
            enc.set_depth(png::BitDepth::Eight);
            let mut w = enc.write_header()?;
            w.write_image_data(&rgba)?;
        }
        pngs.push(buf);
    }
    let mut ico = Vec::new();
    ico.extend_from_slice(&[0, 0, 1, 0]);
    ico.extend_from_slice(&(sizes.len() as u16).to_le_bytes());
    let mut offset = 6 + 16 * sizes.len() as u32;
    for (s, data) in sizes.iter().zip(&pngs) {
        let dim = if *s >= 256 { 0 } else { *s as u8 };
        ico.extend_from_slice(&[dim, dim, 0, 0]);
        ico.extend_from_slice(&1u16.to_le_bytes()); // color planes
        ico.extend_from_slice(&32u16.to_le_bytes()); // bits per pixel
        ico.extend_from_slice(&(data.len() as u32).to_le_bytes());
        ico.extend_from_slice(&offset.to_le_bytes());
        offset += data.len() as u32;
    }
    for d in &pngs {
        ico.extend_from_slice(d);
    }
    std::fs::write(path, ico)
}

fn main() {
    println!("cargo:rerun-if-changed=assets/logo.svg");
    println!("cargo:rerun-if-changed=src/logo_raster.rs");
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows") {
        return;
    }
    let svg = match std::fs::read_to_string("assets/logo.svg") {
        Ok(s) => s,
        Err(e) => {
            println!("cargo:warning=assets/logo.svg not readable ({e}); building without an .exe icon");
            return;
        }
    };
    let ico = PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("canopy.ico");
    if let Err(e) = write_ico(&svg, &ico) {
        println!("cargo:warning=could not build the icon: {e}");
        return;
    }
    let mut res = winresource::WindowsResource::new();
    res.set_icon(ico.to_str().unwrap());
    if let Err(e) = res.compile() {
        println!("cargo:warning=could not embed the application icon: {e}");
    }
}
