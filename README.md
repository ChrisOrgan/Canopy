# <img alt="Canopy" src="/assets/Canopy.png">

Canopy is an interactive phylogenetic tree visualization solution built in Rust with the help of Claude Code. It uses a layer-based
figure approach in a desktop GUI, tree editing, comparative data, Bayesian posterior summaries,
PhyloPic silhouettes, and publication-quality PNG / TIFF / SVG / PDF export.

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
`--cladogram`, `--align`, `--title`, `--background #RRGGBB`, `--tips right|left|top|bottom`, `--open-angle DEG` (fan), `--mcc` or `--consensus 0.5`, `--burnin`.
(Release builds on Windows are GUI-subsystem binaries, so CLI messages are only visible from debug builds.)

## Features

| Requirement | What's there |
|---|---|
| Newick, NEXUS & phyloXML | Quoted labels, comments, BEAST/MrBayes `[&k=v,{a,b}]` and NHX annotations, RAxML/IQ-TREE support labels, TRANSLATE tables, `[&R]/[&U]`, multi-tree files; phyloXML (confidence, taxonomy, properties); Newick/NEXUS/phyloXML writers (annotations preserved) |
| ggtree functionality | Layouts: rectangular, slanted, roundrect, ellipse, dendrogram, circular, fan, inward circular, unrooted (equal angle); cladogram mode; root edge; flip. Layers: `geom_tree`, `geom_tiplab` (align + leaders), `geom_nodelab`, `geom_tippoint`, `geom_nodepoint` (with thresholds), `geom_range` (HPD bars), `geom_hilight`, `geom_cladelab`, `geom_treescale`, `theme_tree2` axis, `gheatmap`, `geom_facet` bars, `geom_phylopic`; `aes(color=…)` mapping with discrete/continuous palettes and legends; `groupClade` colors; `collapse`, `rotate`, `flip`. See [ROADMAP.md](ROADMAP.md) |
| GUI | egui desktop app with tabs, layer stack, property editors, inspector, undo/redo |
| Interactive editing & comparative data | Click/shift-click selection, search, context menu: reroot, midpoint root, rotate, ladderize, collapse, drop/keep tips, extract clade, rename, color clades; CSV/TSV trait tables joined by tip label → tip colors, points, heatmaps, bars, and branch coloring by ML (Brownian motion) or Fitch ancestral reconstruction |
| Tree sets (posterior or bootstrap) | Bayesian posterior samples (burn-in) or bootstrap replicates (all trees); clade frequencies (rooted clades or unrooted splits), majority-rule / extended-majority consensus, MCC tree; posterior or bootstrap support, mean/median heights and 95% HPD annotations (TreeAnnotator-style) |
| Publication figures | TIFF (LZW, RGB, DPI tags) and PNG (pHYs DPI), plus vector SVG and PDF (text kept as text with the font embedded); sizes in in/mm with journal column presets; text sized in points; Arial used when available; live preview |
| Drag & drop | Drop trees or data tables onto the window; drop an image onto a tip to use it as its silhouette; drag a clade onto another branch to move it (SPR); drag layers to reorder |
| PhyloPic | Automatic background lookup by binomial with genus fallback, disk cache (`%LOCALAPPDATA%/Canopy/phylopic`), recoloring, and a credits dialog listing contributors and licenses |

### Also

- **Step through a tree set in one window:** the tree-set panel's ⏮ ◀ ▶ ▶| buttons, sample slider, ← / → keys and Play (adjustable trees per second). Optional ladderizing keeps successive trees comparable; highlights, clade labels, colors and collapsed clades follow their taxa from tree to tree. Use "Back to original tree" to return.
- **Species & data from a node:** clicking a node lists its species with all their data in the inspector, with Copy names, Copy table (tab-separated for Excel) and Save CSV. The same commands are in the right-click menu.
- **Variance–covariance matrix:** for an internal node, the inspector's "Variance–covariance matrix" section shows the clade's phylogenetic VCV (as ape's `vcv`, measured from that node: shared branch length for each pair of tips, distance from the node on the diagonal). Copy matrix and Save CSV export it at full precision; right-click › Copy variance–covariance matrix does the same.
- **Open clade in new tab:** with a posterior tree set, the new tab gets the clade from each post-burn-in sample that contains exactly that clade, i.e. the trees behind its posterior support (20% support in 100 trees gives 20 trees). Support, consensus, MCC and DensiTree then work within the clade. The status bar gives the count.
- **Tree orientation:** Layout panel › Tips at Right / Left / Top / Bottom (rectangular, slanted and round-rect layouts). Top and Bottom rotate the tree so tip labels run vertically; the scale bar turns with it. Headless: `--tips top`.
- **Fan layout:** Layout › Fan draws the tree over part of a circle (default 180°, a half circle above the root), scaled and centred to fill the canvas. "open angle°" sets the gap and "rotate°" turns it. A timescale becomes half-rings with the age axis along the fan's edge. Headless: `--layout fan --open-angle 120`.
- **Zoom:** the Layout panel's − / + buttons zoom about the canvas centre; Fit resets.
- **Check taxon names:** Taxa › Check taxon names (Open Tree of Life)… sends the tip names to the [Open Tree of Life TNRS](https://github.com/OpenTreeOfLife/germinator/wiki/Open-Tree-of-Life-Web-APIs) and lists each as ok, synonym (with the accepted name), possible misspelling (with the closest match) or not found. Tick suggestions to rename tips (undoable); Copy report gives a table. Headless: `canopy check-names tree.nwk`.
- **Hard and soft polytomies:** Tree › Polytomies (or right-click a clade › Polytomies in clade): *Make soft* resolves each polytomy into a binary ladder of zero-length branches (ape `multi2di`); *Make hard* collapses zero-length internal branches into a single node (ape `di2multi`).
- **Bootstrap or posterior tree sets:** a multi-tree file opens as a tree set whose type is guessed from the file name ("boot" → bootstrap replicates) and can be changed in the tree-set panel. Bootstrap sets use every tree (no burn-in) and label clades `bootstrap`; posterior samples label clades `posterior` and have a burn-in slider (0% by default).
- **RAxML / IQ-TREE support and the node report:** support from RAxML bipartitions files (node labels, or `[100]` after branch lengths in `bipartitionsBranchLabels`) and IQ-TREE `SH-aLRT/UFBoot` labels (as `support` and `support_2`) is read automatically. Tree › Node report lists every clade (bipartition) with its support values, branch length, height and taxa; click a node number to select it, Copy table or Save CSV.
- **Regex label editing:** Edit › Find and replace in labels (regex)… with capture groups (`$1`), for tips, internal labels or the selection only, with a preview of every change before renaming (undoable).
- **phyloXML:** open `.xml`/`.phyloxml` files (confidence values, taxonomy names, colors, dates and properties become annotations) and save them with File › Save tree as… (support as `<confidence>`, other annotations as `canopy:` properties so they survive a round trip).
- **Save tree as… (Ctrl+S):** writes Newick, NEXUS or phyloXML with checkboxes for branch lengths, internal node labels and annotations, e.g. a plain NEXUS tree with only topology, tip names and branch lengths. In a tab holding a tree set it saves all trees in the set (or only those after burn-in, or only the tree shown). The last choice is remembered. Canopy has no project files: save trees with this, figures with Export, and reload data tables from their CSV.
- **Tree shape: Colless and γ:** the inspector shows the Colless index and Pybus & Harvey's γ (with a two-tailed p-value against a constant-rate pure-birth model) for any clade, and Tree summary (top of the right panel) shows them for the whole tree. Hover over either name for a one-line reading guide. References:
  - Colless, D. H. (1982) Review of *Phylogenetics: the theory and practice of phylogenetic systematics*, by E. O. Wiley. *Systematic Zoology* 31: 100–104.
  - Pybus, O. G. & Harvey, P. H. (2000) Testing macro-evolutionary models using incomplete molecular phylogenies. *Proc. R. Soc. Lond. B* 267: 2267–2272.
- **Collapse weakly supported nodes:** Tree › Collapse weakly supported nodes… collapses nodes below a support threshold into polytomies. It uses whichever support the tree has (`posterior` from BEAST or Canopy, `prob` from MrBayes, numeric node labels or `bootstrap`), suggests 0.5 or 50, and shows how many nodes will go before you apply.
- **One row per taxon:** a trait table with the same taxon in more than one row is rejected with the duplicated names and row numbers (names compared ignoring case, spaces = underscores).
- **Right panel:** drag its left edge to resize; hide it with ✕ or View › Show right panel, and bring it back with the "◀ Details" button at the top right of the canvas (remembered between sessions).
- **Background color:** Layout panel › Background; exports use it unless "Transparent background" is ticked.
- **Branch-length transforms:** Tree › Transform branch lengths, or right-click a clade › Transform. λ, κ, δ and Ornstein–Uhlenbeck (α) preview live as you drag the slider; Apply keeps the result (undoable) and Cancel restores the tree. References are shown in the dialog and listed below.
- **Copy & paste clades between trees:** select a clade and press Ctrl+C (or right-click › Copy clade), then select a node in another tab and press Ctrl+V to paste as sister, or right-click › Paste clade for sister / child / replace. Copied clades also go on the system clipboard as Newick, and Newick copied from other programs can be pasted (with nothing selected it opens as a new tab).

- **DensiTree:** with a posterior set open (or in a consensus/MCC tab made from one), Layers › Add › DensiTree overlays up to N post-burn-in trees, translucent and aligned at the tips, under the main tree. Optionally colors the three most frequent topologies blue, red and green.
- **Geologic timescale:** Layers › Add › Geologic timescale adds ICS 2023 epoch/period/era bars (standard colors) under time-scaled trees (to the right of the tree when tips are at the top or bottom), with optional boundary lines across the plot. Set "Ma per branch-length unit" and "Youngest tip age" if the tree isn't in Ma with tips at the present. Circular layouts show the finest level as background rings.
- **Tip points away from the tips:** in the Tip points layer, "distance from tip" moves the points out along each tip (in any layout); tip labels shift out past them. Headless: `--color-by COL --point-offset 12`.
- **Bar charts:** Layers › Add › Bar chart (or "bars" next to a numeric column in the data panel) draws a value per tip. Bars start at zero, so negative values (e.g. residuals or contrasts) extend the other way in their own color, with a zero line; a scale with round tick values and the column name sits under the bars (beside them when tips face up or down). Headless: `--bars COLUMN`.
- Headless: `canopy export trees.nex fig.png --mcc --densitree --geoscale`

## Mouse & keys

- Drag empty space to pan · scroll to pan · **+** / **−** to zoom in and out (add **Alt** to zoom vertically only)
- Click to select · Shift/Ctrl-click for multi-select · double-click to rename · right-click for the menu
- Drag a node onto another branch to regraft it there
- Ctrl+O open · Ctrl+S save tree(s) · Ctrl+E export · Ctrl+Z/Ctrl+Y undo/redo · Del drop selection · Esc clear selection

## Code layout

```
src/tree.rs        arena tree with stable node ids and annotations
src/io/            newick.rs, nexus.rs, data.rs (trait tables)
src/ops.rs         reroot, midpoint, ladderize, drop, regraft, collapse
src/consensus.rs   clade counting, consensus, MCC, HPD
src/ancestral.rs   continuous ML and Fitch reconstruction
src/layout.rs      logical layouts; src/scene.rs projects them and builds drawing primitives
src/render/        tiny-skia raster (PNG/TIFF), SVG, PDF (via svg2pdf), font discovery
src/phylopic.rs    PhyloPic API client and cache
src/app/           egui GUI: canvas, panels, export dialog, documents
src/cli.rs         headless `canopy export`
```

The same `Scene` drives the screen and every exporter, so exported figures match the canvas layout.

## References for branch-length transforms

- **λ:** Pagel, M. (1999) Inferring the historical patterns of biological evolution. *Nature* 401: 877–884.
- **κ:** Pagel, M. (1994) Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters. *Proc. R. Soc. Lond. B* 255: 37–45.
- **δ:** Pagel, M. (1997) Inferring evolutionary processes from phylogenies. *Zoologica Scripta* 26: 331–348.
- **Ornstein–Uhlenbeck:** Hansen, T. F. (1997) Stabilizing selection and the comparative analysis of adaptation. *Evolution* 51: 1341–1351.
