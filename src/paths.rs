//! Where Canopy keeps its own files (settings, the PhyloPic cache): a
//! CanopyData folder beside the program, not the system's per-user folders.

use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// Name of the folder beside the program.
const FOLDER: &str = "CanopyData";

/// The CanopyData folder beside the executable (beside the .app bundle on
/// macOS, whose inside must not change), created if needed. Only if that
/// location can't be written (a read-only install) does Canopy fall back to
/// the per-user data folder, so settings can still be saved.
pub fn data_dir() -> PathBuf {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    DIR.get_or_init(|| {
        beside_program().unwrap_or_else(|| {
            let d = dirs::data_local_dir().unwrap_or_else(std::env::temp_dir).join("Canopy");
            let _ = std::fs::create_dir_all(&d);
            d
        })
    })
    .clone()
}

/// The folder for `exe`: CanopyData next to it, or next to its .app bundle.
fn folder_for(exe: &Path) -> Option<PathBuf> {
    let mut dir = exe.parent()?.to_path_buf();
    if dir.ends_with("Contents/MacOS") {
        dir = dir.parent()?.parent()?.parent()?.to_path_buf();
    }
    Some(dir.join(FOLDER))
}

fn beside_program() -> Option<PathBuf> {
    let d = folder_for(&std::env::current_exe().ok()?)?;
    std::fs::create_dir_all(&d).ok()?;
    let probe = d.join(".write-test");
    std::fs::write(&probe, b"").ok()?;
    let _ = std::fs::remove_file(&probe);
    Some(d)
}

/// The settings file (recent files, theme, panels), written by eframe.
pub fn settings_file() -> PathBuf {
    let file = data_dir().join("settings.ron");
    // Settings saved by earlier versions in the per-user folder carry over once.
    if !file.exists() {
        if let Some(old) = dirs::data_dir().map(|d| [d.join("Canopy").join("data").join("app.ron"), d.join("Canopy").join("app.ron"), d.join("canopy").join("app.ron")]) {
            if let Some(o) = old.iter().find(|o| o.exists()) {
                let _ = std::fs::copy(o, &file);
            }
        }
    }
    file
}

/// The PhyloPic silhouette cache, created if needed.
pub fn phylopic_dir() -> Option<PathBuf> {
    let d = data_dir().join("phylopic");
    std::fs::create_dir_all(&d).ok()?;
    Some(d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn folder_beside_the_program() {
        assert_eq!(folder_for(Path::new("/opt/canopy/canopy")), Some(PathBuf::from("/opt/canopy/CanopyData")));
        // macOS: beside the app bundle, not inside it.
        assert_eq!(folder_for(Path::new("/Applications/Canopy.app/Contents/MacOS/canopy")), Some(PathBuf::from("/Applications/CanopyData")));
    }
}
