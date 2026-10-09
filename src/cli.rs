//! Headless command line: render a figure without opening the GUI.
//!
//! ```text
//! canopy export TREE OUT.(png|tif|svg) [--layout rectangular|slanted|roundrect|ellipse|
//!        dendrogram|circular|fan|inward|radial] [--width IN] [--height IN] [--dpi N]
//!        [--data TRAITS.csv] [--color-by COLUMN] [--heatmap COL1,COL2] [--phylopic]
//!        [--cladogram] [--align] [--title TEXT] [--mcc | --consensus P] [--burnin F]
//! ```

use crate::app::document::{smart_layers, Document};
use crate::app::export::{export, ExportJob, Format};
use crate::consensus::{self, HeightMode};
use crate::io::{self, data::DataTable};
use crate::layout::LayoutKind;
use crate::phylopic;
use crate::render::fonts::FontBytes;
use crate::render::raster::RasterFonts;
use crate::render::ImageMap;
use crate::scene::image_key;
use crate::style::*;
use anyhow::{anyhow, bail, Result};
use std::path::PathBuf;
use std::sync::Arc;

pub const USAGE: &str = "usage: canopy export TREE OUT.(png|tif|svg) [--layout NAME] [--width IN] [--height IN] [--dpi N] \
[--data TABLE] [--color-by COL] [--heatmap COLS] [--phylopic] [--cladogram] [--align] [--title TEXT] \
[--mcc | --consensus P] [--burnin F] [--densitree] [--densitree-n N] [--geoscale]";

pub fn run(args: &[String]) -> Result<()> {
    let mut pos = Vec::new();
    let mut opt = std::collections::HashMap::new();
    let mut flags = std::collections::HashSet::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if let Some(k) = a.strip_prefix("--") {
            match k {
                "phylopic" | "cladogram" | "align" | "mcc" | "densitree" | "geoscale" => {
                    flags.insert(k.to_string());
                }
                _ => {
                    let v = args.get(i + 1).ok_or_else(|| anyhow!("--{} needs a value", k))?;
                    opt.insert(k.to_string(), v.clone());
                    i += 1;
                }
            }
        } else {
            pos.push(a.clone());
        }
        i += 1;
    }
    let [tree_path, out] = pos.as_slice() else { bail!(USAGE) };
    let out = PathBuf::from(out);
    let mut trees = io::read_trees(tree_path.as_ref())?;
    let sample = if flags.contains("densitree") { trees.clone() } else { Vec::new() };
    let burnin: f32 = opt.get("burnin").map(|s| s.parse()).transpose()?.unwrap_or(0.1);
    let tree = if flags.contains("mcc") || opt.contains_key("consensus") {
        let skip = (trees.len() as f32 * burnin) as usize;
        let kept = &trees[skip.min(trees.len() - 1)..];
        let rooted = kept[0].rooted != Some(false);
        let s = consensus::summarize(kept, rooted)?;
        if flags.contains("mcc") {
            consensus::mcc_tree(kept, &s, HeightMode::Mean)?
        } else {
            let p: f64 = opt["consensus"].parse()?;
            consensus::consensus(&s, p, HeightMode::Mean)
        }
    } else {
        trees.remove(0)
    };

    let mut doc = Document::new("figure", tree, false);
    doc.layers = smart_layers(&doc.tree, flags.contains("phylopic"));
    if flags.contains("densitree") && sample.len() > 1 {
        let mut p = crate::app::document::Posterior::new(sample);
        p.burnin = burnin;
        doc.posterior = Some(p);
        let mut l = Layer::default_densitree();
        if let (Layer::DensiTree(d), Some(n)) = (&mut l, opt.get("densitree-n")) {
            d.max_trees = n.parse()?;
        }
        doc.layers.insert(0, LayerEntry::new(l));
    }
    if flags.contains("geoscale") {
        let mut l = Layer::default_geoscale();
        if let Layer::Geoscale(g) = &mut l {
            g.boundaries = true;
        }
        doc.layers.push(LayerEntry::new(l));
    }
    if let Some(l) = opt.get("layout") {
        doc.view.layout.kind = match l.to_ascii_lowercase().as_str() {
            "rectangular" | "rect" => LayoutKind::Rectangular,
            "slanted" => LayoutKind::Slanted,
            "roundrect" => LayoutKind::Roundrect,
            "ellipse" => LayoutKind::Ellipse,
            "dendrogram" => LayoutKind::Dendrogram,
            "circular" => LayoutKind::Circular,
            "fan" => LayoutKind::Fan,
            "inward" | "inward_circular" => LayoutKind::InwardCircular,
            "radial" | "unrooted" | "equal_angle" => LayoutKind::Radial,
            other => bail!("unknown layout '{}'", other),
        };
    }
    doc.view.layout.use_lengths = !flags.contains("cladogram");
    if let Some(t) = opt.get("title") {
        doc.view.title = t.clone();
    }
    for e in doc.layers.iter_mut() {
        if let Layer::TipLabels(s) = &mut e.layer {
            s.align = flags.contains("align");
        }
    }
    if let Some(d) = opt.get("data") {
        let table = DataTable::from_path(d.as_ref())?;
        let (m, _) = table.join_to_tree(&mut doc.tree);
        eprintln!("joined {} rows to tips", m);
    }
    if let Some(col) = opt.get("color-by") {
        let mut l = Layer::default_tip_points();
        if let Layer::TipPoints(s) = &mut l {
            s.color_by = Some(col.clone());
        }
        doc.layers.push(LayerEntry::new(l));
    }
    if let Some(cols) = opt.get("heatmap") {
        doc.layers.push(LayerEntry::new(Layer::Heatmap(HeatmapStyle {
            columns: cols.split(',').map(|s| s.trim().to_string()).collect(),
            offset: 4.0,
            cell_width: 12.0,
            palette: Palette::Viridis,
            show_names: true,
            name_size: 8.0,
        })));
    }

    let mut images: ImageMap = ImageMap::new();
    if flags.contains("phylopic") {
        for t in doc.tree.tips() {
            let key = image_key(doc.tree.label(t));
            match phylopic::fetch_blocking(&key) {
                Ok(Some((img, meta))) => {
                    eprintln!("phylopic: {} -> {} ({}, {})", key, meta.matched_name, meta.attribution.unwrap_or_default(), meta.license.unwrap_or_default());
                    images.insert(key, Arc::new(img));
                }
                Ok(None) => eprintln!("phylopic: {} not found", key),
                Err(e) => eprintln!("phylopic: {} failed: {:#}", key, e),
            }
        }
    }

    let format = match out.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase().as_str() {
        "png" => Format::Png,
        "tif" | "tiff" => Format::Tiff,
        "svg" => Format::Svg,
        other => bail!("unsupported output format '{}'", other),
    };
    let dpi: f32 = opt.get("dpi").map(|s| s.parse()).transpose()?.unwrap_or(300.0);
    let w: f32 = opt.get("width").map(|s| s.parse()).transpose()?.unwrap_or(6.85);
    let h: f32 = opt.get("height").map(|s| s.parse()).transpose()?.unwrap_or(w * 1.1);
    let fonts = FontBytes::load();
    let rf = RasterFonts::new(&fonts)?;
    let job = ExportJob {
        path: out.clone(),
        format,
        width_px: (w * dpi).round() as u32,
        height_px: (h * dpi).round() as u32,
        dpi,
        transparent: false,
        tiff_rgb: true,
    };
    export(&doc, &rf, &fonts.name, &images, &job)?;
    eprintln!("wrote {} ({}x{} px at {} dpi)", out.display(), job.width_px, job.height_px, dpi);
    Ok(())
}
