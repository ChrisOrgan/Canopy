//! The Canopy logo. `assets/logo.svg` is the single source: it is compiled
//! into the program and rendered at whatever size is needed (window icon,
//! About dialog). `build.rs` turns the same file into the .exe icon.

use image::RgbaImage;

#[path = "logo_raster.rs"]
mod logo_raster;

/// The logo SVG, as shipped in `assets/logo.svg`.
pub const SVG: &str = include_str!("../assets/logo.svg");

/// Render the logo at `size`×`size` pixels. With `tile`, on a white rounded square.
pub fn render(size: u32, tile: bool) -> RgbaImage {
    RgbaImage::from_raw(size, size, logo_raster::render_rgba(SVG, size, tile)).expect("logo buffer size")
}

/// The logo SVG; with `tile`, wrapped on a white rounded square (app-icon style).
pub fn svg(tile: bool) -> String {
    if !tile {
        return SVG.to_string();
    }
    // Nest the original drawing (without its XML prolog) inside a square tile.
    let inner = SVG.find("<svg").map(|i| &SVG[i..]).unwrap_or(SVG);
    let inner = inner.replacen("<svg", "<svg x=\"20\" y=\"20\" width=\"160\" height=\"160\"", 1);
    let inner = inner.replacen("width=\"100%\" height=\"100%\" ", "", 1);
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 200 200\" width=\"512\" height=\"512\">\n<title>Canopy</title>\n<rect width=\"200\" height=\"200\" rx=\"44\" fill=\"#FFFFFF\"/>\n{}\n</svg>\n",
        inner.trim()
    )
}

/// Window icon for eframe.
pub fn icon_data(size: u32) -> egui::IconData {
    let img = render(size, true);
    egui::IconData { width: img.width(), height: img.height(), rgba: img.into_raw() }
}

#[cfg(test)]
mod tests {
    #[test]
    fn renders_logo() {
        let img = super::render(64, false);
        // The mark should cover a reasonable part of the square.
        let opaque = img.pixels().filter(|p| p.0[3] > 128).count();
        assert!(opaque > 64 * 64 / 20, "only {} opaque pixels", opaque);
        let tile = super::svg(true);
        assert!(tile.contains("<rect") && tile.contains("viewBox=\"0 0 1741 1993\""));
    }
}
