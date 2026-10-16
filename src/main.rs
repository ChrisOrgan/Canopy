#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use std::path::PathBuf;

fn main() -> eframe::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    // Release builds on Windows have no console of their own: with arguments
    // (command-line use), write to the console Canopy was started from.
    #[cfg(all(windows, not(debug_assertions)))]
    if !args.is_empty() {
        extern "system" {
            fn AttachConsole(process_id: u32) -> i32;
        }
        const ATTACH_PARENT_PROCESS: u32 = u32::MAX;
        unsafe {
            AttachConsole(ATTACH_PARENT_PROCESS);
        }
    }
    if args.first().map(|a| a == "--cite" || a == "cite").unwrap_or(false) {
        println!("{}\n\nBibTeX:\n{}", canopy::citation::text(), canopy::citation::bibtex());
        return Ok(());
    }
    if args.first().map(|a| a == "export").unwrap_or(false) {
        if let Err(e) = canopy::cli::run(&args[1..]) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "check-names").unwrap_or(false) {
        if let Err(e) = canopy::cli::check_names(&args[1..]) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "delta-test").unwrap_or(false) {
        if let Err(e) = canopy::cli::delta_test(&args[1..]) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "fabric").unwrap_or(false) {
        if let Err(e) = canopy::cli::bayestraits(&args[1..], canopy::bayestraits::RunKind::Fabric) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "varrates").unwrap_or(false) {
        if let Err(e) = canopy::cli::bayestraits(&args[1..], canopy::bayestraits::RunKind::VariableRates) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "ltt").unwrap_or(false) {
        if let Err(e) = canopy::cli::ltt(&args[1..]) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "traversal").unwrap_or(false) {
        if let Err(e) = canopy::cli::traversal(&args[1..]) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "tanglegram").unwrap_or(false) {
        if let Err(e) = canopy::cli::tanglegram(&args[1..]) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "pir").unwrap_or(false) {
        if let Err(e) = canopy::cli::pir(&args[1..]) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "rogue").unwrap_or(false) {
        if let Err(e) = canopy::cli::rogue(&args[1..]) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "icons").unwrap_or(false) {
        if let Err(e) = canopy::cli::export_icons(&args[1..]) {
            eprintln!("error: {:#}", e);
            std::process::exit(1);
        }
        return Ok(());
    }
    if args.first().map(|a| a == "--help" || a == "-h").unwrap_or(false) {
        println!("canopy [FILES...]   open the GUI\n{}", canopy::cli::USAGE);
        return Ok(());
    }
    let files: Vec<PathBuf> = args.into_iter().map(PathBuf::from).collect();
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Canopy")
            .with_inner_size([1440.0, 900.0])
            .with_min_inner_size([800.0, 500.0])
            .with_drag_and_drop(true)
            .with_icon(canopy::logo::icon_data(256)),
        ..Default::default()
    };
    eframe::run_native("Canopy", options, Box::new(move |cc| Ok(Box::new(canopy::app::CanopyApp::new(cc, files)))))
}
