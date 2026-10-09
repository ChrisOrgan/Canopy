# Canopy

Desktop app (Rust, egui/eframe 0.29) for viewing, editing and publishing phylogenetic trees, modelled on R's ggtree. Owner: Chris Organ. See README.md for features and ROADMAP.md for what's done and what's next; keep both up to date when features change.

## Build, run, test

```bash
cargo build                 # debug build
cargo run --release         # run the GUI (opens empty; Help > Open example tree)
cargo test --lib            # unit tests; keep them passing
cargo run --example make_posterior   # regenerate examples/primates_posterior.trees
cargo run --example make_icon        # export logo PNG/SVG variants into exports/
```

On the original Windows PC, Rust is installed via rustup but not on PATH in new shells. Prefix commands with `$env:Path = "$env:USERPROFILE\.cargo\bin;$env:Path"` (PowerShell). There's no Python on that machine, so write helper scripts in Rust (as `examples/`).

## Verifying changes

- Check visual changes with the headless exporter and look at the PNG. It uses the same scene builder as the GUI:
  `target/debug/canopy export examples/primates.nwk out.png --layout circular --phylopic --dpi 100`
  Options include `--data examples/primate_traits.csv --heatmap col1,col2 --color-by col --mcc --consensus 0.5 --burnin 0.1 --densitree --geoscale --align --cladogram --title T`.
- Don't drive the GUI with simulated keystrokes or clicks while the user may have Canopy open: input can land in their window. Ask first, or verify headlessly.
- Add a unit test next to the code for new logic (scene tests inspect `Scene::prims`).

## Architecture (src/)

- `tree.rs`: arena tree. Structural edits never remove nodes from the arena; detached nodes become unreachable, so `NodeId`s stay stable (selections and layers keep working). Use `compact()` only for output.
- `io/`: `newick.rs` (BEAST/NHX annotations), `nexus.rs` (TRANSLATE, multi-tree), `data.rs` (CSV/TSV join to tips; `tips_table` for the species/data table).
- `ops.rs`: reroot, midpoint, ladderize, drop/keep, regraft (drag-and-drop SPR), graft/paste, `map_node` (match clades between trees by tip set), branch-length transforms (λ, κ, δ, OU).
- `consensus.rs`: clade counting, majority/extended consensus, MCC, HPD. `ancestral.rs`: ML continuous, Fitch discrete.
- `layout.rs`: logical x (depth) / y (slot) per node; `LayoutKind::MENU` is what the UI offers (other variants stay for old projects/CLI). The slanted layout defaults to a V-shaped cladogram.
- `scene.rs`: turns tree + `ViewState` + layers into renderer-independent `Prim`s in pixels. All layouts go through `Geo` (Lin / Polar / Free). Sizes are in points times `pt` (GUI pt = 1; export pt = dpi/72).
- `render/`: tiny-skia raster (PNG/TIFF with DPI), SVG, font discovery (Arial if present). `app/canvas.rs` paints scenes with egui.
- `style.rs`: colors, palettes and the `Layer` enum (ggtree-style layers). `geotime.rs`: ICS 2023 intervals.
- `app/`: `mod.rs` (menus, dialogs, actions, settings), `document.rs` (undo, projects, posterior browsing with annotation remapping), `panels.rs`, `canvas.rs`, `export.rs`.
- `phylopic.rs`: PhyloPic API v2 client with a disk cache. `cli.rs`: `canopy export`.
- Logo: `assets/logo.svg` is the single source. `logo.rs` renders it via resvg; `build.rs` builds the .exe icon from it.

## Adding a layer type

1. `style.rs`: style struct, `Layer` variant, `name()` / `description()` arms, a `default_*()` constructor. Use `#[serde(default)]` for new fields so old `.canopy` projects still load.
2. `scene.rs`: a `draw_*` function plus a dispatch arm in `build()`. Reserve margins or columns there if it needs space.
3. `app/panels.rs`: editor arm in `layer_editor` and an entry in `add_layer_menu`.
4. Optionally a CLI flag in `cli.rs`, a test, and README/ROADMAP lines.

## Conventions

- egui 0.29 API (`ComboBox::from_id_salt`, `DragValue::range`, `Stroke::new(1.0f32, ..)` needs the f32 suffix).
- Match the existing style: compact functions, doc comments on public items, no new dependencies without reason.
- Numbers shown to users are rounded to 3 decimals; heights below 1e-4 display as 0 (`io::data::table_value`).
- Commit and push only when asked.
