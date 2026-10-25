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


/// System fonts with wide symbol coverage (arrows, geometric shapes, dingbats),
/// used by the interface as fallbacks for glyphs its own fonts lack.
const SYMBOL_FONTS: &[&str] = &["C:\\Windows\\Fonts\\seguisym.ttf", "/System/Library/Fonts/Apple Symbols.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"];

/// The first symbol font found on this system.
pub fn symbol_font() -> Option<Vec<u8>> {
    SYMBOL_FONTS.iter().find_map(|p| std::fs::read(p).ok())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ab_glyph::Font;

    /// Every non-ASCII character in the interface code has a glyph in the fonts
    /// bundled with egui that the interface uses (its proportional fonts, then
    /// Hack), so none shows as a box on any system. The system fonts installed
    /// as further fallbacks are not counted: they differ between platforms.
    #[test]
    fn interface_symbols_have_glyphs() {
        let defs = egui::FontDefinitions::default();
        let mut names = defs.families[&egui::FontFamily::Proportional].clone();
        names.push("Hack".into());
        let chain: Vec<Vec<u8>> = names.iter().map(|n| defs.font_data[n].font.to_vec()).collect();
        let has = |c: char| chain.iter().any(|f| ab_glyph::FontRef::try_from_slice(f).map(|f| f.glyph_id(c).0 != 0).unwrap_or(false));
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
        let mut missing = std::collections::BTreeSet::new();
        for dir in [root.join("app"), root.clone()] {
            for e in std::fs::read_dir(dir).unwrap() {
                let p = e.unwrap().path();
                if p.extension().map(|x| x == "rs").unwrap_or(false) {
                    for line in std::fs::read_to_string(&p).unwrap().lines() {
                        // Code only: comments may use any character.
                        let code = line.split("//").next().unwrap_or("");
                        missing.extend(code.chars().filter(|&c| c as u32 > 127 && !('\u{300}'..='\u{36f}').contains(&c) && !has(c)));
                    }
                }
            }
        }
        assert!(missing.is_empty(), "no glyph for {:?}", missing);
    }
}

