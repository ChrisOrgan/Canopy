# Canopy roadmap

What Canopy covers from ggtree/treeio/tidytree and related R packages (phangorn, geiger, deeptime), and what's left.
✅ done · 🟡 partial · ⬜ not yet.

## Layouts
- ✅ rectangular, slanted, roundrect, ellipse, dendrogram, circular, fan (`open.angle`), inward_circular, equal_angle
- ✅ `branch.length = "none"`, root edge, `scale_x_reverse`-style flip, rotation
- ⬜ daylight layout, `ape` unrooted variants, `layout_rectangular(branch.length = <attr>)` (scale by any attribute)
- ⬜ `scaleClade` (stretch a clade vertically), `open_tree`/`rotate_tree` animations

## Geoms & annotation
- ✅ geom_tree, geom_tiplab (align with dotted leaders), geom_nodelab, geom_tippoint, geom_nodepoint, geom_range, geom_hilight, geom_cladelab (bar + text), geom_treescale, theme_tree2, gheatmap, geom_facet (bars), geom_phylopic (resizable, optional fit to tip spacing)
- ✅ Branch-length labels on branches; node labels colored by a threshold (posterior support red < 0.5, green ≥ 0.5 by default; colors and cut-off editable)
- ✅ Geologic timescale (deeptime `coord_geo`): ICS 2023 eras, periods and epochs with CGMW colors, boundary lines, Ma conversion; circular layouts draw solid concentric interval rings with names and an age axis
- 🟡 aes mapping: color mapped for tree/labels/points/heatmap; ⬜ size, shape, alpha and linetype mappings
- 🟡 geom_hilight: rectangle/sector/hull; ⬜ gradient fill, `type = "encircle"` smoothing
- 🟡 Geologic timescale: rectangular-type and circular layouts; ⬜ unrooted layout, stages/ages level, several levels at once in circular layouts
- ⬜ geom_strip (label spanning two arbitrary tips), geom_taxalink (curved links between taxa)
- ⬜ geom_inset / nodepie / nodebar (pie or bar charts at nodes, e.g. ancestral state probabilities)
- ⬜ geom_tiplab2 / geom_text2 parse expressions (plotmath, mixed italic/roman labels)
- ⬜ geom_cladelab with images, geom_balance, geom_label background boxes
- ⬜ facet_plot / geom_facet with boxplots, dots, density, alignments (msaplot)
- ⬜ multiple heatmaps with independent scales, `colnames_angle`, `legend` positioning options
- ⬜ geom_rootpoint, geom_point2 subset expressions (only threshold filters today)

## Data & treeio
- ✅ Newick, NEXUS, BEAST/MrBayes annotations, NHX, CSV/TSV joins (`%<+%`)
- ✅ Per-clade species & data table (root distance, nodes to root, height, traits; 3 decimals) with copy/CSV export
- ⬜ read.jplace, PAML (rst/mlc), HyPhy, CODEML, r8s, RAxML bipartitions, IQ-TREE `.iqtree` reports, phyloXML, Jtree/JSON
- ⬜ Sequence alignments (FASTA/PHYLIP) for `msaplot`
- ⬜ Pattern/regex label editing, label lookup tables (rename tips from a two-column file)

## Analysis
- ✅ reroot, midpoint root, ladderize, rotate/flip, drop/keep tips, extract clade, SPR regraft, collapse weak nodes
- ✅ Copy/paste clades between trees (sister, child or replace; Newick via the system clipboard)
- ✅ Branch-length transforms λ, κ, δ and Ornstein–Uhlenbeck α (whole tree or clade, live preview, with references)
- ✅ Majority/extended-majority consensus, MCC, posterior, height HPD; ML continuous and Fitch discrete ancestral states
- ✅ Posterior browsing: step/play through samples in one window with per-sample clade support; DensiTree overlay (optionally colored by topology)
- ✅ Clade shape statistics: Colless index (raw and normalized), cherries
- 🟡 Tree shape: ⬜ Sackin index, gamma statistic, γ/LTT plots, Colless for polytomies
- ⬜ Common-ancestor heights (TreeAnnotator `-heights ca`), strict consensus option, tree-to-tree distances (RF)
- ⬜ Stochastic character mapping, Mk-model marginal ancestral states (for nodepie)
- ⬜ Tanglegrams / cophylo comparison of two trees
- ⬜ Fitting λ/κ/δ/OU by maximum likelihood to trait data (transforms are currently set by hand)

## Export & UX
- ✅ PNG/TIFF with DPI metadata, SVG, live preview, journal size presets
- ✅ Headless `canopy export` (layouts, data, heatmap, PhyloPic, MCC/consensus, DensiTree, timescale)
- ✅ App icon and logo (window, executable, About), empty start window with Help › Open example tree
- ⬜ PDF/EPS export, CMYK TIFF, embedded font selection UI (Arial auto-detected today)
- ⬜ Undo for style edits (structural edits are undoable), multi-select drag of layers
- ⬜ Large-tree performance work (scene caching and culling for more than 10k tips; caching DensiTree geometry)
- ⬜ Clade-level PhyloPic images, choosing among alternative PhyloPic images
- ⬜ Console output from release builds for `canopy export` (currently a GUI-subsystem binary on Windows)
