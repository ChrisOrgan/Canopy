//! Headless command line: render a figure without opening the GUI.
//!
//! ```text
//! canopy export TREE OUT.(png|tif|svg|pdf|gif) [--layout rectangular|slanted|roundrect|ellipse|
//!        dendrogram|circular|fan|inward|radial] [--width IN] [--height IN] [--dpi N]
//!        [--data TRAITS.csv] [--color-by COLUMN] [--point-offset PT] [--heatmap COL1,COL2] [--bars COL] [--phylopic]
//!        [--cladogram] [--align] [--title TEXT] [--title-size PT] [--background #RRGGBB] [--tips SIDE] [--open-angle DEG] [--mcc | --consensus P] [--burnin F]
//! ```

use crate::app::document::{smart_layers, Document};
use crate::app::export::{export, ExportJob, Format};
use crate::consensus::{self, HeightMode};
use crate::io::{self, data::DataTable};
use crate::layout::{LayoutKind, TipSide};
use crate::phylopic;
use crate::render::fonts::FontBytes;
use crate::render::raster::RasterFonts;
use crate::render::ImageMap;
use crate::scene::image_key;
use crate::style::*;
use anyhow::{anyhow, bail, Result};
use std::path::PathBuf;
use std::sync::Arc;

pub const USAGE: &str = "usage: canopy export TREE OUT.(png|tif|svg|pdf|gif) [--layout NAME] [--width IN] [--height IN] [--dpi N] \
[--data TABLE] [--color-by COL] [--point-offset PT] [--heatmap COLS] [--bars COL] [--phylopic] [--cladogram] [--align] [--title TEXT] [--title-size PT] [--fossils RANGES.csv] [--clade-label TIP1,TIP2=TEXT] [--clade-gradient HEX] [--background HEX] [--tips SIDE] [--open-angle DEG] \
[--mcc | --consensus P] [--burnin F] [--densitree] [--densitree-n N] [--geoscale] [--geoscale-ages] [--fps N] [--every N] (a .gif steps through every tree of a tree set)";

/// `canopy check-names TREE`: check tip names against the Open Tree of Life
/// taxonomy and print a tab-separated report (label, status, accepted name, OTT id).
pub fn check_names(args: &[String]) -> Result<()> {
    let [tree_path] = args else { bail!("usage: canopy check-names TREE") };
    let tree = io::read_trees(tree_path.as_ref())?.remove(0);
    let tips = tree.tips();
    let names: Vec<String> = tips.iter().map(|&t| crate::taxonomy::query_name(tree.label(t))).collect();
    let checks = crate::taxonomy::check_names(&names)?;
    println!("label\tstatus\taccepted_name\tott_id\tflags");
    for (&t, c) in tips.iter().zip(&checks) {
        println!(
            "{}\t{}\t{}\t{}\t{}",
            tree.label(t),
            c.status.label(),
            c.accepted.as_deref().unwrap_or(""),
            c.ott_id.map(|i| i.to_string()).unwrap_or_default(),
            c.flags.join(",")
        );
    }
    Ok(())
}

pub fn run(args: &[String]) -> Result<()> {
    let mut pos = Vec::new();
    let mut opt = std::collections::HashMap::new();
    let mut flags = std::collections::HashSet::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if let Some(k) = a.strip_prefix("--") {
            match k {
                "phylopic" | "cladogram" | "align" | "mcc" | "densitree" | "geoscale" | "geoscale-ages" => {
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
    let gif = out.extension().map(|e| e.eq_ignore_ascii_case("gif")).unwrap_or(false);
    let all = if gif { trees.clone() } else { Vec::new() };
    let sample = if flags.contains("densitree") { trees.clone() } else { Vec::new() };
    let burnin: f32 = opt.get("burnin").map(|s| s.parse()).transpose()?.unwrap_or(0.0);
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
    if gif {
        if all.len() < 2 {
            bail!("an animated GIF needs a file with several trees");
        }
        doc.posterior = Some(crate::app::document::Posterior::new(all));
    }
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
            g.youngest_age = doc.tree.youngest_tip_age().unwrap_or(0.0);
            g.ages = flags.contains("geoscale-ages");
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
    if let Some(a) = opt.get("open-angle") {
        doc.view.layout.open_angle = a.parse()?;
    }
    if let Some(t) = opt.get("tips") {
        doc.view.layout.tips = match t.to_ascii_lowercase().as_str() {
            "right" => TipSide::Right,
            "left" => TipSide::Left,
            "top" => TipSide::Top,
            "bottom" => TipSide::Bottom,
            other => bail!("unknown --tips side '{}' (right, left, top or bottom)", other),
        };
    }
    doc.view.layout.use_lengths = !flags.contains("cladogram");
    if let Some(t) = opt.get("title") {
        doc.view.title = t.clone();
    }
    if let Some(s) = opt.get("title-size") {
        doc.view.title_size = s.parse()?;
    }
    if let Some(c) = opt.get("background") {
        doc.view.background = Some(Color::hex(c));
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
    // --clade-label "TipA,TipB=Text": a label on the two tips' MRCA (gradient with --clade-gradient HEX[AA]).
    if let Some(spec) = opt.get("clade-label") {
        let (tips, text) = spec.split_once('=').ok_or_else(|| anyhow!("--clade-label wants TIP1,TIP2=TEXT"))?;
        let ids: Vec<crate::tree::NodeId> = tips.split(',').map(|t| doc.tree.find_label(t.trim()).ok_or_else(|| anyhow!("no tip \"{}\"", t.trim()))).collect::<Result<_>>()?;
        let node = doc.tree.mrca(&ids).ok_or_else(|| anyhow!("no common ancestor for {}", tips))?;
        let mut s = CladeLabelStyle::new(node, text.to_string());
        if let Some(c) = opt.get("clade-gradient") {
            s.gradient = true;
            s.gradient_color = Color::hex(c);
        }
        doc.layers.push(LayerEntry::new(Layer::CladeLabel(s)));
    }
    if let Some(f) = opt.get("fossils") {
        let text = std::fs::read_to_string(f)?;
        let (m, un) = crate::fossil::import(&text, f, &mut doc.tree)?;
        eprintln!("fossil ranges for {} taxa{}", m, if un.is_empty() { String::new() } else { format!("; not in the tree: {}", un.join(", ")) });
        doc.layers.push(LayerEntry::new(Layer::default_fossil_ranges()));
    }
    if let Some(col) = opt.get("color-by") {
        let mut l = Layer::default_tip_points();
        if let Layer::TipPoints(s) = &mut l {
            s.color_by = Some(col.clone());
            s.offset = opt.get("point-offset").map(|v| v.parse()).transpose()?.unwrap_or(0.0);
        }
        doc.layers.push(LayerEntry::new(l));
    }
    if let Some(col) = opt.get("bars") {
        doc.layers.push(LayerEntry::new(Layer::default_bars(col)));
    }
    if let Some(cols) = opt.get("heatmap") {
        doc.layers.push(LayerEntry::new(Layer::Heatmap(HeatmapStyle {
            columns: cols.split(',').map(|s| s.trim().to_string()).collect(),
            offset: 4.0,
            cell_width: 12.0,
            palette: Palette::Viridis,
            reverse_palette: false,
            show_names: true,
            name_size: 8.0,
        })));
    }

    let mut images: ImageMap = ImageMap::new();
    if flags.contains("phylopic") {
        // Tips by name, and nodes given a taxon (`phylopic` attribute).
        let mut keys: Vec<String> = doc.tree.tips().into_iter().map(|t| image_key(doc.tree.label(t))).collect();
        for n in doc.tree.preorder() {
            if let Some(crate::tree::Attr::Text(name)) = doc.tree.nodes[n].attrs.get(crate::scene::PHYLOPIC_KEY) {
                keys.push(image_key(name));
            }
        }
        keys.dedup();
        for key in keys {
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
        "pdf" => Format::Pdf,
        "gif" => Format::Gif,
        other => bail!("unsupported output format '{}'", other),
    };
    let dpi: f32 = opt.get("dpi").map(|s| s.parse()).transpose()?.unwrap_or(if gif { 100.0 } else { 300.0 });
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
        fps: opt.get("fps").map(|s| s.parse()).transpose()?.unwrap_or(4.0),
        every: opt.get("every").map(|s| s.parse()).transpose()?.unwrap_or(1),
    };
    if gif {
        let n = crate::app::export::export_gif(&mut doc, &rf, &images, &job)?;
        eprintln!("wrote {} ({} frames, {}x{} px)", out.display(), n, job.width_px, job.height_px);
        return Ok(());
    }
    export(&doc, &rf, &fonts, &images, &job)?;
    if matches!(format, Format::Svg | Format::Pdf) {
        eprintln!("wrote {} ({:.2} x {:.2} in, vector)", out.display(), w, h);
    } else {
        eprintln!("wrote {} ({}x{} px at {} dpi)", out.display(), job.width_px, job.height_px, dpi);
    }
    Ok(())
}

/// `canopy icons [DIR]`: write the logo in other formats (default DIR: exports):
/// canopy-icon.svg (white rounded tile), canopy-256.png (tile), canopy-logo-1024.png
/// (mark only) and canopy.iconset/ (the PNG sizes macOS `iconutil` turns into an .icns).
/// Used by the release workflow; the .exe and window icons are built from the SVG directly.
pub fn export_icons(args: &[String]) -> Result<()> {
    let dir = std::path::PathBuf::from(args.first().map(String::as_str).unwrap_or("exports"));
    let iconset = dir.join("canopy.iconset");
    std::fs::create_dir_all(&iconset)?;
    std::fs::write(dir.join("canopy-icon.svg"), crate::logo::svg(true))?;
    crate::logo::render(256, true).save(dir.join("canopy-256.png"))?;
    crate::logo::render(1024, false).save(dir.join("canopy-logo-1024.png"))?;
    for size in [16u32, 32, 128, 256, 512] {
        crate::logo::render(size, true).save(iconset.join(format!("icon_{0}x{0}.png", size)))?;
        crate::logo::render(size * 2, true).save(iconset.join(format!("icon_{0}x{0}@2x.png", size)))?;
    }
    eprintln!("wrote logo files to {}", dir.display());
    Ok(())
}

/// `canopy delta-test TREE`: the node-density artifact (delta) test on the first tree.
pub fn delta_test(args: &[String]) -> Result<()> {
    let [tree_path] = args else { bail!("usage: canopy delta-test TREE") };
    let trees = io::read_trees(tree_path.as_ref())?;
    let r = crate::nodedensity::delta_test(&trees[0])?;
    print!("{}", r.report());
    Ok(())
}

/// `canopy varrates TREE FILE` / `canopy fabric TREE FILE`: per-branch summary
/// of a BayesTraits run as a tab-separated table. The file's entries decide the
/// summary ("Branch" = variable rates, "LS_Beta" = Fabric); a file with only
/// Node scalars is read as the command's kind.
pub fn bayestraits(args: &[String], chosen: crate::bayestraits::RunKind) -> Result<()> {
    use crate::bayestraits::RunKind;
    let [tree_path, path] = args else { bail!("usage: canopy {} TREE VARRATES.txt", if chosen == RunKind::Fabric { "fabric" } else { "varrates" }) };
    let mut tree = io::read_trees(tree_path.as_ref())?.remove(0);
    let text = std::fs::read_to_string(path)?;
    let detected = crate::bayestraits::detect(&text)?;
    let kind = detected.unwrap_or(chosen);
    if detected.is_some() && kind != chosen {
        eprintln!("note: this is {} output; showing the {} summary", kind.name(), kind.name());
    }
    let (header, cells) = match kind {
        RunKind::VariableRates => {
            let vr = crate::varrates::parse(&text)?;
            let mut rates = vr.summarize();
            for w in vr.apply(&mut rates, &mut tree) {
                eprintln!("warning: {}", w);
            }
            crate::varrates::table(&rates)
        }
        RunKind::Fabric => {
            let fab = crate::fabric::parse(&text)?;
            let mut rows = fab.summarize();
            for w in fab.apply(&mut rows, &mut tree) {
                eprintln!("warning: {}", w);
            }
            crate::fabric::table(&rows)
        }
    };
    print!("{}", io::data::table_to_text(&header, &cells, b'\t'));
    Ok(())
}

/// `canopy rogue TREES [--burnin F] [--max-size K] [--min-prob T]`: entropy-based
/// rogue scores and skeleton analysis of a tree set (Klawitter, Bouckaert &
/// Drummond 2024). Rogues are only flagged; the trees are not changed.
pub fn rogue(args: &[String]) -> Result<()> {
    let mut opt = std::collections::HashMap::new();
    let mut pos = Vec::new();
    let mut i = 0;
    while i < args.len() {
        if let Some(k) = args[i].strip_prefix("--") {
            opt.insert(k.to_string(), args.get(i + 1).cloned().ok_or_else(|| anyhow!("--{} needs a value", k))?);
            i += 1;
        } else {
            pos.push(args[i].clone());
        }
        i += 1;
    }
    let [path] = pos.as_slice() else { bail!("usage: canopy rogue TREES [--burnin F] [--max-size K] [--min-prob T]") };
    let trees = io::read_trees(path.as_ref())?;
    let burnin: f64 = opt.get("burnin").map(|s| s.parse()).transpose()?.unwrap_or(0.1);
    let k: usize = opt.get("max-size").map(|s| s.parse()).transpose()?.unwrap_or(1);
    let t: f64 = opt.get("min-prob").map(|s| s.parse()).transpose()?.unwrap_or(0.1);
    let skip = ((trees.len() as f64) * burnin) as usize;
    let started = std::time::Instant::now();
    let ccd = crate::rogue::Ccd::from_trees(&trees[skip.min(trees.len() - 1)..])?;
    println!("# {} trees after burn-in, {} taxa; entropy {:.4}", trees.len() - skip, ccd.num_taxa(), ccd.entropy());
    println!("# rogue scores (clades of at most {} taxa with probability >= {})", k, t);
    println!("score\tprobability\ttaxa");
    for s in ccd.rogue_scores(k, t).iter() {
        println!("{:.4}\t{:.3}\t{}", s.score, s.probability, s.taxa.join(" "));
    }
    println!("# skeleton analysis (rogues flagged in order; entropy and MAP probability for the taxa not yet flagged)");
    println!("round\tflagged\ttaxa_left\tentropy\tMAP_probability");
    for (r, x) in ccd.skeleton(k, t).iter().enumerate() {
        println!("{}\t{}\t{}\t{:.4}\t{:.4e}", r, x.flagged.join(" "), x.taxa_left, x.entropy, x.map_probability);
    }
    println!("# Klawitter, J., Bouckaert, R. R. & Drummond, A. J. (2024) Skeletons in the forest: using entropy-based rogue detection on Bayesian phylogenetic tree distributions. bioRxiv, doi:10.1101/2024.09.25.615070");
    eprintln!("({:.1}s)", started.elapsed().as_secs_f32());
    Ok(())
}

/// `canopy ltt TREE [OUT]`: lineages through time of the first tree, as
/// tab-separated (depth from the root, lineages) steps; with OUT (PNG, TIFF, SVG
/// or PDF), the plot instead, with every tree of a multi-tree file faint behind
/// the first.
pub fn ltt(args: &[String]) -> Result<()> {
    let geoscale = args.iter().any(|a| a == "--geoscale");
    let show_tree = args.iter().any(|a| a == "--tree");
    let pos: Vec<&String> = args.iter().filter(|a| *a != "--geoscale" && *a != "--tree").collect();
    let (path, out) = match pos.as_slice() {
        [p] => (*p, None),
        [p, o] => (*p, Some(PathBuf::from(o))),
        _ => bail!("usage: canopy ltt TREE [OUT.png|svg|pdf] [--geoscale] [--tree]"),
    };
    let trees = io::read_trees(path.as_ref())?;
    if let Some(out) = out {
        let ages = |t: &crate::tree::Tree| {
            let s = crate::ltt::ltt(t);
            let deepest = s.last().map(|x| x.0).unwrap_or(0.0);
            s.into_iter().map(|(x, k)| (deepest - x, k)).collect::<Vec<_>>()
        };
        let main = ages(&trees[0]);
        let others: Vec<Vec<(f64, usize)>> = if trees.len() > 1 { trees.iter().map(ages).collect() } else { Vec::new() };
        let format = match out.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase().as_str() {
            "png" => Format::Png,
            "tif" | "tiff" => Format::Tiff,
            "svg" => Format::Svg,
            "pdf" => Format::Pdf,
            other => bail!("unsupported output format '{}'", other),
        };
        let fonts = FontBytes::load();
        let rf = RasterFonts::new(&fonts)?;
        let job = ExportJob { path: out.clone(), format, width_px: 1800, height_px: if show_tree { 1800 } else { 1200 }, dpi: 300.0, transparent: false, tiff_rgb: true, fps: 4.0, every: 1 };
        let geo = if geoscale { vec![crate::geotime::Level::Epoch, crate::geotime::Level::Period] } else { Vec::new() };
        let title = format!("Lineages through time: {}", std::path::Path::new(path).file_name().and_then(|f| f.to_str()).unwrap_or(""));
        crate::app::export::export_scene(&|_, w, h, pt| crate::plots::ltt(&main, &others, &crate::plots::LttPlot { log: false, smooth: true, geo: &geo, title: &title, tree: show_tree.then_some((&trees[0], 0.0, 1.0)) }, w, h, pt), None, &rf, &fonts, &ImageMap::new(), &job)?;
        eprintln!("wrote {}", out.display());
        return Ok(());
    }
    println!("depth\tlineages");
    for (x, k) in crate::ltt::ltt(&trees[0]) {
        println!("{}\t{}", crate::tree::format_num(x, 6), k);
    }
    Ok(())
}

/// `canopy pir TREE DATA COLUMN [COLUMN ..]`: phylogenetic imbalance ratio
/// (Gardner & Organ 2021) of categorical columns on the first tree.
pub fn pir(args: &[String]) -> Result<()> {
    let [tree, data, cols @ ..] = args else { bail!("usage: canopy pir TREE DATA COLUMN [COLUMN ..]") };
    if cols.is_empty() {
        bail!("usage: canopy pir TREE DATA COLUMN [COLUMN ..]");
    }
    let trees = io::read_trees(tree.as_ref())?;
    let table = DataTable::from_path(data.as_ref())?;
    let idx: Vec<usize> = cols.iter().map(|c| table.columns.iter().position(|x| x == c).ok_or_else(|| anyhow!("no column \"{}\"", c))).collect::<Result<_>>()?;
    let p = crate::pir::pir(&trees[0], &table, &idx, None).ok_or_else(|| anyhow!("no taxa in the tree have values in every column"))?;
    let f = |x: f64| crate::tree::format_num(x, 3);
    println!("taxa\t{}\nmin_steps\t{}\nsteps\t{}\nCI\t{}\nNIR\t{}\nPIR\t{}", p.taxa, p.min_steps, p.steps, f(p.ci), f(p.nir), f(p.pir));
    for (c, n) in &p.classes {
        println!("class {{{}}}\t{}", c, n);
    }
    Ok(())
}

/// `canopy tanglegram LEFT RIGHT OUT [--untangle] [--mirror] [--title T]
/// [--width IN] [--height IN] [--dpi N]`: the first tree of each file face to
/// face, tips joined by name (no lines with --mirror). PNG, TIFF, SVG or PDF.
pub fn tanglegram(args: &[String]) -> Result<()> {
    let usage = "usage: canopy tanglegram LEFT RIGHT OUT [--untangle] [--mirror] [--title T] [--width IN] [--height IN] [--dpi N]";
    let (mut pos, mut opt, mut flags) = (Vec::new(), std::collections::HashMap::new(), std::collections::HashSet::new());
    let mut i = 0;
    while i < args.len() {
        match args[i].strip_prefix("--") {
            Some(k @ ("untangle" | "mirror")) => {
                flags.insert(k.to_string());
            }
            Some(k) => {
                opt.insert(k.to_string(), args.get(i + 1).ok_or_else(|| anyhow!("--{} needs a value", k))?.clone());
                i += 1;
            }
            None => pos.push(args[i].clone()),
        }
        i += 1;
    }
    let [left, right, out] = pos.as_slice() else { bail!(usage) };
    let mut l = io::read_trees(left.as_ref())?.remove(0);
    let mut r = io::read_trees(right.as_ref())?.remove(0);
    if flags.contains("untangle") {
        let (b, a) = crate::tangle::untangle(&mut l, &mut r);
        eprintln!("crossings: {} -> {}", b, a);
    }
    let (ll, rl) = (smart_layers(&l, false), smart_layers(&r, false));
    let mut view = crate::style::ViewState::default();
    view.title = opt.get("title").cloned().unwrap_or_default();
    let rview = crate::style::ViewState::default();
    let style = crate::tangle::TangleStyle { lines: !flags.contains("mirror"), ..Default::default() };
    let out = PathBuf::from(out);
    let format = match out.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase().as_str() {
        "png" => Format::Png,
        "tif" | "tiff" => Format::Tiff,
        "svg" => Format::Svg,
        "pdf" => Format::Pdf,
        other => bail!("unsupported output format '{}'", other),
    };
    let dpi: f32 = opt.get("dpi").map(|s| s.parse()).transpose()?.unwrap_or(300.0);
    let w: f32 = opt.get("width").map(|s| s.parse()).transpose()?.unwrap_or(9.0);
    let h: f32 = opt.get("height").map(|s| s.parse()).transpose()?.unwrap_or(6.0);
    let fonts = FontBytes::load();
    let rf = RasterFonts::new(&fonts)?;
    let job = ExportJob { path: out.clone(), format, width_px: (w * dpi).round() as u32, height_px: (h * dpi).round() as u32, dpi, transparent: false, tiff_rgb: true, fps: 4.0, every: 1 };
    let li = crate::scene::BuildInput { tree: &l, view: &view, layers: &ll, overlay: &[] };
    let ri = crate::scene::BuildInput { tree: &r, view: &rview, layers: &rl, overlay: &[] };
    crate::app::export::export_scene(&|env, w, h, pt| crate::tangle::build(&li, &ri, &style, env, w, h, pt), None, &rf, &fonts, &ImageMap::new(), &job)?;
    eprintln!("wrote {} ({} matching tips)", out.display(), crate::tangle::matches(&l, &r).len());
    Ok(())
}

/// `canopy traversal TREE OUT [--post] [--step K] [--fps N] [--no-numbers]
/// [--width IN] [--height IN] [--dpi N]`: pre-order (or post-order) traversal of
/// the first tree, as an animated GIF (OUT.gif) or one frame after K visits
/// (PNG, TIFF, SVG or PDF; all visits by default).
pub fn traversal(args: &[String]) -> Result<()> {
    use crate::traversal::{overlay, visits, Order};
    let usage = "usage: canopy traversal TREE OUT [--post] [--step K] [--fps N] [--no-numbers] [--width IN] [--height IN] [--dpi N]";
    let (mut pos, mut opt, mut flags) = (Vec::new(), std::collections::HashMap::new(), std::collections::HashSet::new());
    let mut i = 0;
    while i < args.len() {
        match args[i].strip_prefix("--") {
            Some(k @ ("post" | "no-numbers")) => {
                flags.insert(k.to_string());
            }
            Some(k) => {
                opt.insert(k.to_string(), args.get(i + 1).ok_or_else(|| anyhow!("--{} needs a value", k))?.clone());
                i += 1;
            }
            None => pos.push(args[i].clone()),
        }
        i += 1;
    }
    let [tree, out] = pos.as_slice() else { bail!(usage) };
    let t = io::read_trees(tree.as_ref())?.remove(0);
    let layers = smart_layers(&t, false);
    let view = crate::style::ViewState::default();
    let order = if flags.contains("post") { Order::Post } else { Order::Pre };
    let numbers = !flags.contains("no-numbers");
    let out = PathBuf::from(out);
    let ext = out.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
    let format = match ext.as_str() {
        "gif" => Format::Gif,
        "png" => Format::Png,
        "tif" | "tiff" => Format::Tiff,
        "svg" => Format::Svg,
        "pdf" => Format::Pdf,
        other => bail!("unsupported output format '{}'", other),
    };
    let gif = format == Format::Gif;
    let dpi: f32 = opt.get("dpi").map(|s| s.parse()).transpose()?.unwrap_or(if gif { 100.0 } else { 300.0 });
    let w: f32 = opt.get("width").map(|s| s.parse()).transpose()?.unwrap_or(6.0);
    let h: f32 = opt.get("height").map(|s| s.parse()).transpose()?.unwrap_or(4.5);
    let fonts = FontBytes::load();
    let rf = RasterFonts::new(&fonts)?;
    let fps: f32 = opt.get("fps").map(|s| s.parse()).transpose()?.unwrap_or(1.5);
    let job = ExportJob { path: out.clone(), format, width_px: (w * dpi).round() as u32, height_px: (h * dpi).round() as u32, dpi, transparent: false, tiff_rgb: true, fps, every: 1 };
    let input = crate::scene::BuildInput { tree: &t, view: &view, layers: &layers, overlay: &[] };
    let frame = |env: &dyn crate::scene::SceneEnv, step: usize, w: f32, h: f32, pt: f32| {
        let mut sc = crate::scene::build(&input, env, w, h, pt);
        let v = visits(&t, &sc, order);
        overlay(&mut sc, &v, step, order, numbers, pt);
        sc
    };
    let n = visits(&t, &crate::scene::build(&input, &crate::scene::ApproxEnv, 600.0, 400.0, 1.0), order).len();
    if gif {
        // An empty first frame, one per visit, and the finished traversal held for two more.
        crate::app::export::export_animation(n + 3, &|env, i, w, h, pt| frame(env, i.min(n), w, h, pt), None, &rf, &ImageMap::new(), &job)?;
        eprintln!("wrote {} ({} frames)", out.display(), n + 3);
    } else {
        let step: usize = opt.get("step").map(|s| s.parse()).transpose()?.unwrap_or(n);
        crate::app::export::export_scene(&|env, w, h, pt| frame(env, step, w, h, pt), None, &rf, &fonts, &ImageMap::new(), &job)?;
        eprintln!("wrote {} (step {} of {})", out.display(), step.min(n), n);
    }
    Ok(())
}
