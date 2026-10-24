//! Help › References: the works Canopy builds on, as light Markdown (section
//! lines in **bold**, items as "- " lines with **bold** labels and *italic*
//! titles), drawn with clickable URLs and DOIs.

use egui::RichText;

pub const REFERENCES: &str = r#"**Overall design & function space**
- **ggtree:** Yu, G., Smith, D. K., Zhu, H., Guan, Y. & Lam, T. T.-Y. (2017) ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. *Methods in Ecology and Evolution* 8: 28–36. https://doi.org/10.1111/2041-210X.12628
- **FigTree:** Rambaut, A. FigTree: a graphical viewer of phylogenetic trees. https://github.com/rambaut/figtree
- **Mesquite:** Maddison, W. P. & Maddison, D. R. Mesquite: a modular system for evolutionary analysis. https://www.mesquiteproject.org
- **BayesTrees:** Meade, A. & Pagel, M. BayesTrees. http://www.evolution.reading.ac.uk/BayesTrees.html
- **TreeAnnotator (BEAST):** Drummond, A. J., & Rambaut, A. (2007). BEAST: Bayesian evolutionary analysis by sampling trees. *BMC Evolutionary Biology*, 7, 214. https://doi.org/10.1186/1471-2148-7-214

**Tree format**
- **phyloXML format:** Han, M. V. & Zmasek, C. M. (2009) phyloXML: XML for evolutionary biology and comparative genomics. *BMC Bioinformatics* 10: 356.
- **Newick:** https://phylipweb.github.io/phylip/newicktree.html
- **Nexus:** Maddison, D. R., D. L. Swofford, and W. P. Maddison. Nexus: An Extensible File Format for Systematic Information. *Systematic Biology* 46, no. 4 (1997): 590–621. https://doi.org/10.1093/sysbio/46.4.590

**Tree inference support**
- **BEAST 2:** Bouckaert, Remco, Joseph Heled, Denise Kühnert, et al. 2014. “BEAST 2: A Software Platform for Bayesian Evolutionary Analysis.” *PLoS Comput Biol* 10 (4): e1003537.
- **HyPhy:** Kosakovsky Pond, Sergei L., Art F. Y. Poon, Ryan Velazquez, et al. 2020. “HyPhy 2.5—A Customizable Platform for Evolutionary Hypothesis Testing Using Phylogenies.” *Molecular Biology and Evolution* 37 (1): 295–99. https://doi.org/10.1093/molbev/msz197
- **MrBayes:** Ronquist, F., M. Teslenko, P. van der Mark, et al. 2012. “MrBayes 3.2: Efficient Bayesian Phylogenetic Inference and Model Choice across a Large Model Space.” *Systematic Biology*.
- **RAxML:** Stamatakis, Alexandros. 2014. “RAxML Version 8: A Tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies.” *Bioinformatics*, January 21, btu033. https://doi.org/10.1093/bioinformatics/btu033
- **IQ-TREE 3:** Wong, Thomas K. F., Nhan Ly-Trong, Huaiyan Ren, et al. 2026. “IQ-TREE 3: Phylogenomic Inference Software Using Complex Evolutionary Models.” *Molecular Biology and Evolution* 43 (5): msag117. https://doi.org/10.1093/molbev/msag117

**Tree methods**
- **Fitch parsimony:** Fitch, W. M. (1971) Toward defining the course of evolution: minimum change for a specific tree topology. *Systematic Zoology* 20: 406–416.
- **Parsimony with polytomies:** Hartigan, J. A. (1973) Minimum mutation fits to a given tree. *Biometrics* 29: 53–65.
- **Equal-angle unrooted layout:** Felsenstein, J. (2004) *Inferring Phylogenies*. Sinauer.
- **Ancestral states (ML, Brownian motion):** equivalent to phytools fastAnc. Revell, L. J. (2012) phytools: an R package for phylogenetic comparative biology (and other things). *Methods Ecol. Evol.* 3: 217–223.
- **ape equivalents (di2multi, multi2di, vcv):** Paradis, E. & Schliep, K. (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35: 526–528.
- **DensiTree layer:** Bouckaert, R. R. (2010) DensiTree: making sense of sets of phylogenetic trees. *Bioinformatics* 26: 1372–1373.

**Example data**
- Kay, R. F., Organ, C. L., Morse, P. E., Allen, K. L. & Kirk, E. C. (2026) Fossil evidence favors a role for vision in the modular evolution of the primate neocortex. *Science* 393: 622–627. doi:10.1126/science.aee3541

**Data sources and services**
- PhyloPic (phylopic.org): silhouettes, with each image's contributor and license in Help › PhyloPic image credits.
- Open Tree of Life taxonomy and its name-resolution service: the taxon name check.
- Macrostrat (macrostrat.org, CC-BY 4.0): geologic timescale names, boundaries and colors, from the International Chronostratigraphic Chart v2024/12.

**Branch-length transforms**
- **λ:** Pagel, M. (1999) Inferring the historical patterns of biological evolution. *Nature* 401: 877–884.
- **κ:** Pagel, M. (1994) Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters. *Proc. R. Soc. Lond. B* 255: 37–45.
- **δ:** Pagel, M. (1997) Inferring evolutionary processes from phylogenies. *Zoologica Scripta* 26: 331–348.
- **Ornstein–Uhlenbeck:** Hansen, T. F. (1997) Stabilizing selection and the comparative analysis of adaptation. *Evolution* 51: 1341–1351.

**BayesTraits imports**
- Venditti, C., Meade, A. & Pagel, M. (2011) Multiple routes to mammalian diversity. *Nature* 479: 393–396. doi:10.1038/nature10516
- Pagel, M., O'Donovan, C. & Meade, A. (2022) General statistical model shows that macroevolutionary patterns and processes are consistent with Darwinian gradualism. *Nature Communications* 13: 1113. doi:10.1038/s41467-022-28595-z

**Node-density artifact**
- Webster, A. J., Payne, R. J. & Pagel, M. (2003) Molecular phylogenies link rates of evolution and speciation. *Science* 301: 478.
- Venditti, C., Meade, A. & Pagel, M. (2006) Detecting the node-density artifact in phylogeny reconstruction. *Systematic Biology* 55: 637–643.

**MCMC convergence (combined runs)**
- Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. & Bürkner, P.-C. (2021) Rank-normalization, folding, and localization: an improved R-hat for assessing convergence of MCMC. *Bayesian Analysis* 16: 667–718. https://doi.org/10.1214/20-BA1221
- Gelman, A. & Rubin, D. B. (1992) Inference from iterative simulation using multiple sequences. *Statistical Science* 7: 457–472. https://doi.org/10.1214/ss/1177011136
- Geyer, C. J. (1992) Practical Markov chain Monte Carlo. *Statistical Science* 7: 473–483. https://doi.org/10.1214/ss/1177011137
- Rambaut, A., Drummond, A. J., Xie, D., Baele, G. & Suchard, M. A. (2018) Posterior summarization in Bayesian phylogenetics using Tracer 1.7. *Systematic Biology* 67: 901–904. https://doi.org/10.1093/sysbio/syy032
- Ronquist, F., Teslenko, M., van der Mark, P., et al. (2012) MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. *Systematic Biology* 61: 539–542. https://doi.org/10.1093/sysbio/sys029

**Rogue taxa**
- Klawitter, J., Bouckaert, R. R. & Drummond, A. J. (2024) Skeletons in the forest: using entropy-based rogue detection on Bayesian phylogenetic tree distributions. *bioRxiv*. doi:10.1101/2024.09.25.615070

**Tree & data metrics**
- **Lineages through time:** Nee, S., Mooers, A. Ø. & Harvey, P. H. (1992) Tempo and mode of evolution revealed from molecular phylogenies. *PNAS* 89: 8322–8326.
- **γ statistic:** Pybus, O. G. & Harvey, P. H. (2000) Testing macro-evolutionary models using incomplete molecular phylogenies. *Proc. R. Soc. Lond. B* 267: 2267–2272.
- **Colless index:** Colless, D. H. (1982) Review of *Phylogenetics: the theory and practice of phylogenetic systematics*. *Systematic Zoology* 31: 100–104.
- **Consistency index:** Kluge, A. G. & Farris, J. S. (1969) Quantitative phyletics and the evolution of anurans. *Systematic Zoology* 18: 1–32.
- **Phylogenetic imbalance ratio:** Gardner, J. D. & Organ, C. L. (2021) Evolutionary sample size and consilience in phylogenetic comparative analysis. *Systematic Biology* 70: 1061–1075.

**Plotting**
- **Plot design:** Tufte, E. R. (2001) *The Visual Display of Quantitative Information*, 2nd ed. Graphics Press, Cheshire, CT.
- **Okabe–Ito palette:** Okabe, M. & Ito, K. (2008) Color Universal Design (CUD): how to make figures and presentations that are friendly to colorblind people. https://jfly.uni-koeln.de/color/
- **Palette families:** ColorBrewer palettes (Brewer, C. A., colorbrewer2.org), and viridis and magma (van der Walt, S. & Smith, N., matplotlib).
"#;

/// A piece of an item: text (bold, italic) or a link (shown text, target).
#[derive(Debug, PartialEq)]
enum Piece {
    Text(String, bool, bool),
    Link(String, String),
}

/// Split an item into text with **bold** / *italic* runs and links (http(s)
/// URLs, and "doi:" names linked through doi.org); a trailing period stays text.
fn pieces(line: &str) -> Vec<Piece> {
    let mut out = Vec::new();
    let (mut bold, mut italic) = (false, false);
    let mut buf = String::new();
    let flush = |buf: &mut String, out: &mut Vec<Piece>, bold: bool, italic: bool| {
        if !buf.is_empty() {
            out.push(Piece::Text(std::mem::take(buf), bold, italic));
        }
    };
    let mut rest = line;
    while !rest.is_empty() {
        if let Some(r) = rest.strip_prefix("**") {
            flush(&mut buf, &mut out, bold, italic);
            bold = !bold;
            rest = r;
        } else if let Some(r) = rest.strip_prefix('*') {
            flush(&mut buf, &mut out, bold, italic);
            italic = !italic;
            rest = r;
        } else if rest.starts_with("https://") || rest.starts_with("http://") || rest.starts_with("doi:") {
            flush(&mut buf, &mut out, bold, italic);
            let end = rest.find(char::is_whitespace).unwrap_or(rest.len());
            let mut link = &rest[..end];
            while link.ends_with(['.', ',', ')', ';']) {
                link = &link[..link.len() - 1];
            }
            let target = match link.strip_prefix("doi:") {
                Some(d) => format!("https://doi.org/{}", d),
                None => link.to_string(),
            };
            out.push(Piece::Link(link.to_string(), target));
            rest = &rest[link.len()..];
        } else {
            let c = rest.chars().next().unwrap();
            buf.push(c);
            rest = &rest[c.len_utf8()..];
        }
    }
    flush(&mut buf, &mut out, bold, italic);
    out
}

/// The references as plain text (no Markdown marks), for copying.
pub fn plain_text() -> String {
    REFERENCES.replace("**", "").replace('*', "")
}

/// Draw the references: section headings, then their items.
pub fn show(ui: &mut egui::Ui) {
    for line in REFERENCES.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        match line.strip_prefix("- ") {
            None => {
                ui.add_space(8.0);
                ui.label(RichText::new(line.trim_matches('*')).strong().size(15.0));
            }
            Some(item) => {
                ui.horizontal_wrapped(|ui| {
                    ui.spacing_mut().item_spacing.x = 0.0;
                    ui.label("•  ");
                    for p in pieces(item) {
                        match p {
                            Piece::Text(t, b, i) => {
                                let mut r = RichText::new(t);
                                if b {
                                    r = r.strong();
                                }
                                if i {
                                    r = r.italics();
                                }
                                ui.label(r);
                            }
                            Piece::Link(shown, target) => {
                                ui.hyperlink_to(shown, target);
                            }
                        }
                    }
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn items_split_into_styles_and_links() {
        let p = pieces("**λ:** Pagel, M. (1999) *Nature* 401. doi:10.1038/x.");
        assert_eq!(p[0], Piece::Text("λ:".into(), true, false));
        assert_eq!(p[2], Piece::Text("Nature".into(), false, true));
        assert_eq!(p[4], Piece::Link("doi:10.1038/x".into(), "https://doi.org/10.1038/x".into()));
        assert_eq!(p[5], Piece::Text(".".into(), false, false));
        assert!(!plain_text().contains('*'));
    }
}
