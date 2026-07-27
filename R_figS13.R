# FigS13 — Exploratory sequence- and structure-derived targetability descriptors (R)
# Purpose: panel a (structural similarity + domain conservation) +
#          panel b (composite prioritization score ranking)
#          individual 90x95mm panels -> 180x95mm composite
# Note:    content moved from former Fig. 4c/d when the structural analysis was
#          demoted to exploratory descriptor status (see Methods / Limitations)
# Usage:   Rscript R_figS13.R   (run from the project root containing paralog_sl_predictor/)

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)

# ── Constants ──
BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7
PANEL_W <- 90; PANEL_H <- 95  # mm
OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Colors ──
BLUE  <- "#2171B5"; RED   <- "#CB181D"; ORANGE <- "#E6550D"
TEAL  <- "#0D7377"; GRAY  <- "#636363"

# ── Theme ──
theme_sci <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.line    = element_line(linewidth = 0.4),
  axis.ticks   = element_line(linewidth = 0.3),
  axis.text    = element_text(size = TICK_FS),
  axis.title   = element_text(size = BASE_FS),
  legend.text  = element_text(size = LEGEND_FS),
  legend.title = element_blank(),
  legend.background = element_blank(),
  legend.key        = element_blank(),
  legend.key.size   = unit(3, "mm"),
  plot.margin  = margin(4, 4, 4, 4, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

# ═══════════════════════════════════════════════════════════════
# PANEL A — Structural similarity and domain conservation
# ═══════════════════════════════════════════════════════════════
panel_a <- function() {
  struct_path <- "paralog_sl_predictor/output/alphafold_structural_analysis.csv"
  if (file.exists(struct_path)) {
    st <- read_csv(struct_path, show_col_types = FALSE) %>%
      filter(domain_similarity > 0) %>%  # exclude zero-domain pairs (e.g. BRCA1/BRCA2)
      arrange(desc(structural_similarity)) %>% head(8)
  } else {
    stop("paralog_sl_predictor/output/alphafold_structural_analysis.csv not found — ",
         "run the pipeline first; simulated fallbacks are forbidden")
  }
  st$pair <- factor(paste0(st$gene_a, "/", st$gene_b),
                    levels = rev(paste0(st$gene_a, "/", st$gene_b)))
  df <- st %>% select(pair, structural_similarity, domain_similarity) %>%
    pivot_longer(-pair) %>%
    mutate(name = recode(name, structural_similarity = "Structural",
                         domain_similarity = "Domain"))

  ggplot(df, aes(value, pair, fill = name)) +
    geom_col(position = position_dodge(0.7), width = 0.55) +
    scale_fill_manual(values = c(Structural = TEAL, Domain = ORANGE)) +
    labs(x = "Score", y = NULL) +
    # legend below the panel: inside placement overlapped the short bottom bars;
    # reverse legend so entry order matches top-to-bottom bar order
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_sci + theme(legend.position = "bottom",
                      legend.direction = "horizontal",
                      legend.margin = margin(0, 0, 0, 0),
                      plot.margin = margin(4, 4, 4, 20, "pt"))
}

# ═══════════════════════════════════════════════════════════════
# PANEL B — Composite prioritization score ranking
# ═══════════════════════════════════════════════════════════════
panel_b <- function() {
  struct_path <- "paralog_sl_predictor/output/alphafold_structural_analysis.csv"
  if (file.exists(struct_path)) {
    st <- read_csv(struct_path, show_col_types = FALSE)
    if ("clinical_targetability" %in% names(st)) {
      cand <- st %>% arrange(desc(clinical_targetability)) %>% head(10) %>%
        mutate(label = paste0(driver, "->", paralog), score = clinical_targetability)
    } else { cand <- NULL }
  } else { cand <- NULL }
  if (is.null(cand) || nrow(cand) == 0) {
    stop("clinical_targetability column missing in ",
         "alphafold_structural_analysis.csv — run the structural/targetability ",
         "pipeline first; simulated fallbacks are forbidden")
  }
  cand$label <- factor(cand$label, levels = rev(cand$label))
  # Highlight matches the text narrative: ARID1A->ARID1B is the leading
  # SELECTIVE candidate; NF1->RASA2 ranks first on score only through a
  # near-zero pan-essential denominator (selectivity ~ 0)
  cand$cat <- "Others"
  cand$cat[cand$driver == "NF1"    & cand$paralog == "RASA2"]  <- "Rank 1 (non-selective)"
  cand$cat[cand$driver == "ARID1A" & cand$paralog == "ARID1B"] <- "Leading selective candidate"
  cand$cat <- factor(cand$cat, levels = c("Leading selective candidate",
                                          "Rank 1 (non-selective)", "Others"))
  cat_colors <- c("Leading selective candidate" = RED,
                  "Rank 1 (non-selective)" = ORANGE, "Others" = BLUE)
  cand$txt_col <- cat_colors[as.character(cand$cat)]

  ggplot(cand, aes(score, label, fill = cat)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = sprintf("%.3f", score), color = txt_col),
              hjust = -0.1, size = 2.5, fontface = "bold", show.legend = FALSE) +
    # direct in-bar category labels instead of a legend: a 3-entry legend
    # overlapped the bottom bars and their value labels in this 90mm panel
    geom_text(data = cand %>% filter(cat != "Others"),
              aes(x = score - 0.02, label = as.character(cat)),
              hjust = 1, size = 2.5, color = "white", fontface = "bold",
              family = "Arial", show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_manual(values = cat_colors, guide = "none") +
    geom_vline(xintercept = 0.5, linewidth = 0.3, color = GRAY,
               linetype = "dashed", alpha = 0.3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Composite prioritization score", y = NULL) +
    theme_sci
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== FigS13 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b()

ggsave(file.path(OUT_DIR, "FigS13_panel_a.pdf"), pa,
       width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS13_panel_b.pdf"), pb,
       width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)

p <- cowplot::plot_grid(pa, pb, ncol = 2, rel_widths = c(0.5, 0.5),
                        labels = c("a","b"),
                        label_size = 9, label_fontface = "bold",
                        label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "FigS13_Structure_Targetability.pdf"), p,
       width = 180, height = 95, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS13_Structure_Targetability.svg"), p,
       width = 180, height = 95, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS13_Structure_Targetability.tiff"), p,
       width = 180, height = 95, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS13_Structure_Targetability.pdf (180x95mm) ✓")
