//! Font discovery shared by the GUI and the exporters, so on-screen text
//! metrics match exported figures. Prefers Arial/Helvetica-like system fonts
//! (the usual journal requirement) and falls back to egui's bundled font.

use std::path::Path;

pub struct FontBytes {
    pub regular: Vec<u8>,
    pub italic: Option<Vec<u8>>,
    pub bold: Option<Vec<u8>>,
    pub name: String,
}

const CANDIDATES: &[(&str, &str, &str, &str)] = &[
    ("Arial", "C:\\Windows\\Fonts\\arial.ttf", "C:\\Windows\\Fonts\\ariali.ttf", "C:\\Windows\\Fonts\\arialbd.ttf"),
    ("Arial", "/Library/Fonts/Arial.ttf", "/Library/Fonts/Arial Italic.ttf", "/Library/Fonts/Arial Bold.ttf"),
    (
        "Arial",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    ),
    (
        "Liberation Sans",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    ),
    (
        "DejaVu Sans",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    ),
];

impl FontBytes {
    pub fn load() -> FontBytes {
        for (name, r, i, b) in CANDIDATES {
            if let Ok(regular) = std::fs::read(r) {
                let opt = |p: &str| if Path::new(p).exists() { std::fs::read(p).ok() } else { None };
                return FontBytes { regular, italic: opt(i), bold: opt(b), name: name.to_string() };
            }
        }
        let defs = egui::FontDefinitions::default();
        let regular = defs
            .font_data
            .get("Ubuntu-Light")
            .or_else(|| defs.font_data.values().next())
            .map(|f| f.font.to_vec())
            .unwrap_or_default();
        FontBytes { regular, italic: None, bold: None, name: "Ubuntu".into() }
    }
}
