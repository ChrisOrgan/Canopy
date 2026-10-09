//! Export the logo in other formats, from `assets/logo.svg`:
//! `cargo run --example make_icon`
//! Writes into exports/: canopy-icon.svg (on a white rounded tile),
//! canopy-256.png (tile), canopy-logo-1024.png (mark only), and
//! canopy.iconset/ (the PNG sizes macOS `iconutil` turns into an .icns).
//! The .exe and window icons need nothing: they are built from the SVG.

fn main() -> anyhow::Result<()> {
    std::fs::create_dir_all("exports/canopy.iconset")?;
    std::fs::write("exports/canopy-icon.svg", canopy::logo::svg(true))?;
    canopy::logo::render(256, true).save("exports/canopy-256.png")?;
    canopy::logo::render(1024, false).save("exports/canopy-logo-1024.png")?;
    for size in [16u32, 32, 128, 256, 512] {
        canopy::logo::render(size, true).save(format!("exports/canopy.iconset/icon_{0}x{0}.png", size))?;
        canopy::logo::render(size * 2, true).save(format!("exports/canopy.iconset/icon_{0}x{0}@2x.png", size))?;
    }
    println!("wrote exports/canopy-icon.svg, canopy-256.png, canopy-logo-1024.png, canopy.iconset/");
    Ok(())
}
