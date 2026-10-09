//! Export the logo in other formats, from `assets/logo.svg`:
//! `cargo run --example make_icon`
//! Writes into exports/: canopy-icon.svg (on a white rounded tile),
//! canopy-256.png (tile) and canopy-logo-1024.png (mark only).
//! The .exe and window icons need nothing: they are built from the SVG.

fn main() -> anyhow::Result<()> {
    std::fs::create_dir_all("exports")?;
    std::fs::write("exports/canopy-icon.svg", canopy::logo::svg(true))?;
    canopy::logo::render(256, true).save("exports/canopy-256.png")?;
    canopy::logo::render(1024, false).save("exports/canopy-logo-1024.png")?;
    println!("wrote exports/canopy-icon.svg, canopy-256.png, canopy-logo-1024.png");
    Ok(())
}
