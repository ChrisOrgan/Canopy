# Canopy

Interactive phylogenetic tree visualization in Rust, modelled on R's **ggtree**: a layer-based
figure grammar with a desktop GUI, tree editing, comparative data, Bayesian posterior summaries,
PhyloPic silhouettes, and publication-quality PNG / TIFF / SVG export.

## Build & run

Needs Rust (stable) and, on Windows, the MSVC C++ build tools.

```bash
cargo run --release                                  # GUI with an example tree
cargo run --release -- examples/primates_posterior.trees
cargo test                                           # core test suite
cargo run --example make_posterior                   # regenerate the example posterior sample
```

### The logo

`assets/logo.svg` is the single source for the logo. The window icon and About box render it at runtime, and `build.rs` turns it into the Windows `.exe` icon. To change the logo, replace that file and rebuild:

```bash
cargo build --release
```

To export PNG and tile-SVG versions (into `exports/`): `cargo run --example make_icon`.
### Headless export (scripting / `ggsave`-style)

```bash
cargo run --release -- export examples/primates.nwk fig.tif --data examples/primate_traits.csv --heatmap body_mass_kg,brain_g --color-by diet --phylopic --dpi 600
cargo run --release -- export examples/primates_posterior.trees mcc.png --mcc --burnin 0.1 --layout fan
```

Options: `--layout rectangular|slanted|roundrect|ellipse|dendrogram|circular|fan|inward|radial`,
`--width/--height` (inches), `--dpi`, `--data`, `--color-by`, `--heatmap`, `--phylopic`,
`--cladogram`, `--align`, `--title`, `--mcc` or `--consensus 0.5`, `--burnin`.
(Release builds on Windows are GUI-subsystem binaries, so CLI messages are only visible from debug builds.)

## Features

| Requirement | What's there |
|---|---|
| Newick & NEXUS | Quoted labels, comments, BEAST/MrBayes `[&k=v,{a,b}]` and NHX annotations, TRANSLATE tables, `[&R]/[&U]`, multi-tree files; Newick/NEXUS writers (annotations preserved) |
| ggtree functionality | Layouts: rectangular, slanted, roundrect, ellipse, dendrogram, circular, fan, inward circular, unrooted (equal angle); cladogram mode; root edge; flip. Layers: `geom_tree`, `geom_tiplab` (align + leaders), `geom_nodelab`, `geom_tippoint`, `geom_nodepoint` (with thresholds), `geom_range` (HPD bars), `geom_hilight`, `geom_cladelab`, `geom_treescale`, `theme_tree2` axis, `gheatmap`, `geom_facet` bars, `geom_phylopic`; `aes(color=…)` mapping with discrete/continuous palettes and legends; `groupClade` colors; `collapse`, `rotate`, `flip`. See [ROADMAP.md](ROADMAP.md) |
| GUI | egui desktop app with tabs, layer stack, property editors, inspector, undo/redo, project files (`.canopy`) |
| Interactive editing & comparative data | Click/shift-click selection, search, context menu: reroot, midpoint root, rotate, ladderize, collapse, drop/keep tips, extract clade, rename, color clades; CSV/TSV trait tables joined by tip label → tip colors, points, heatmaps, bars, and branch coloring by ML (Brownian motion) or Fitch ancestral reconstruction |
| Bayesian posterior samples | Burn-in, clade frequencies (rooted clades or unrooted splits), majority-rule / extended-majority consensus, MCC tree; posterior, mean/median heights and 95% HPD annotations (TreeAnnotator-style) |
| Publication figures | TIFF (LZW, RGB, DPI tags) and PNG (pHYs DPI), plus SVG; sizes in in/mm with journal column presets; text sized in points; Arial used when available; live preview |
| Drag & drop | Drop trees, data tables or projects onto the window; drop an image onto a tip to use it as its silhouette; drag a clade onto another branch to move it (SPR); drag layers to reorder |
| PhyloPic | Automatic background lookup by binomial with genus fallback, disk cache (`%LOCALAPPDATA%/Canopy/phylopic`), recoloring, and a credits dialog listing contributors and licenses |

### Also

- **Step through posterior samples in one window:** the posterior panel's ⏮ ◀ ▶ ▶| buttons, sample slider, ← / → keys and Play (adjustable trees per second). Optional ladderizing keeps successive trees comparable; highlights, clade labels, colors and collapsed clades follow their taxa from tree to tree. Use "Back to original tree" to return.
- **Species & data from a node:** clicking a node lists its species with all their data in the inspector, with Copy names, Copy table (tab-separated for Excel) and Save CSV. The same commands are in the right-click menu.
- **Branch-length transforms:** Tree › Transform branch lengths, or right-click a clade › Transform. λ, κ, δ and Ornstein–Uhlenbeck (α) preview live as you drag the slider; Apply keeps the result (undoable) and Cancel restores the tree. References are shown in the dialog and listed below.
- **Copy & paste clades between trees:** select a clade and press Ctrl+C (or right-click › Copy clade), then select a node in another tab and press Ctrl+V to paste as sister, or right-click › Paste clade for sister / child / replace. Copied clades also go on the system clipboard as Newick, and Newick copied from other programs can be pasted (with nothing selected it opens as a new tab).

- **DensiTree:** with a posterior set open (or in a consensus/MCC tab made from one), Layers › Add › DensiTree overlays up to N post-burn-in trees, translucent and aligned at the tips, under the main tree. Optionally colors the three most frequent topologies blue, red and green.
- **Geologic timescale:** Layers › Add › Geologic timescale adds ICS 2023 epoch/period/era bars (standard colors) under time-scaled trees, with optional boundary lines across the plot. Set "Ma per branch-length unit" and "Youngest tip age" if the tree isn't in Ma with tips at the present. Circular layouts show the finest level as background rings.
- Headless: `canopy export trees.nex fig.png --mcc --densitree --geoscale`

## Mouse & keys

- Drag empty space to pan · scroll to pan · **Ctrl+scroll** to zoom (add **Alt** to zoom vertically only)
- Click to select · Shift/Ctrl-click for multi-select · double-click to rename · right-click for the menu
- Drag a node onto another branch to regraft it there
- Ctrl+O open · Ctrl+S save project · Ctrl+E export · Ctrl+Z/Ctrl+Y undo/redo · Del drop selection · Esc clear selection

## Code layout

```
src/tree.rs        arena tree with stable node ids and annotations
src/io/            newick.rs, nexus.rs, data.rs (trait tables)
src/ops.rs         reroot, midpoint, ladderize, drop, regraft, collapse
src/consensus.rs   clade counting, consensus, MCC, HPD
src/ancestral.rs   continuous ML and Fitch reconstruction
src/layout.rs      logical layouts; src/scene.rs projects them and builds drawing primitives
src/render/        tiny-skia raster (PNG/TIFF), SVG, font discovery
src/phylopic.rs    PhyloPic API client and cache
src/app/           egui GUI: canvas, panels, export dialog, documents
src/cli.rs         headless `canopy export`
```

The same `Scene` drives the screen and every exporter, so exported figures match the canvas layout.

## References for branch-length transforms

- **λ and δ:** Pagel, M. (1997) Inferring evolutionary processes from phylogenies. *Zoologica Scripta* 26: 331–348. Pagel, M. (1999) Inferring the historical patterns of biological evolution. *Nature* 401: 877–884.
- **κ:** Pagel, M. (1994) Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters. *Proc. R. Soc. Lond. B* 255: 37–45.
- **Ornstein–Uhlenbeck:** Hansen, T. F. (1997) Stabilizing selection and the comparative analysis of adaptation. *Evolution* 51: 1341–1351. Butler, M. A. & King, A. A. (2004) Phylogenetic comparative analysis: a modeling approach for adaptive evolution. *Am. Nat.* 164: 683–695.
- **Implementation** follows geiger: Pennell, M. W. et al. (2014) geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. *Bioinformatics* 30: 2216–2218.