# Fig4 — Dependency-Window Classification (R)
# Purpose: single panel (DWS bubble plot, former panel b) → 120×95mm, no panel letter
# Note:    former panel a (PRISM assay-validity anchors) moved to Supplementary
#          Fig. S7b (R_figS7_prism.R) after manuscript restructuring; former
#          panels c/d (structure-derived descriptors) moved to R_figS13.R
# Usage:   Rscript R_fig4.R

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)

# ── Constants ──
BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7; ANNOT_FS <- 5.5
PANEL_W   <- 90; PANEL_H <- 95  # mm (single-panel PDF source)
FIG_W     <- 120; FIG_H   <- 95 # mm (final single-panel figure)
OUT_DIR   <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Colors ──
BLUE  <- "#2171B5"; RED   <- "#CB181D"; GREEN <- "#238B45"
ORANGE <- "#E6550D"; GRAY  <- "#636363"; TEAL  <- "#0D7377"

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
  plot.margin  = margin(4, 4, 4, 4, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

save_panel <- function(p, name) {
  ggsave(file.path(OUT_DIR, paste0("Fig4_panel_", name, ".pdf")), p,
         width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
  message(sprintf("  panel %s ✓", name))
}

# ═══════════════════════════════════════════════════════════════
# PANEL — Therapeutic Window (DWS bubble plot)
# ═══════════════════════════════════════════════════════════════
panel_b <- function() {
  tw_path <- "paralog_sl_predictor/output/therapeutic_window_paralog_classification.csv"
  if (file.exists(tw_path)) {
    tw <- read_csv(tw_path, show_col_types = FALSE)
  } else {
    stop("paralog_sl_predictor/output/therapeutic_window_paralog_classification.csv not found — ",
         "run the pipeline first; simulated fallbacks are forbidden")
  }
  tw <- tw %>% mutate(
    class_label = case_when(
      classification == "HIGH_SELECTIVITY" ~ "HIGH",
      classification == "MODERATE"         ~ "MODERATE",
      classification == "LOW_SELECTIVITY"  ~ "LOW",
      classification == "PAN_ESSENTIAL"    ~ "PAN",
      TRUE                                 ~ classification),
    class_label = factor(class_label, c("HIGH","MODERATE","LOW","PAN")),
    size = 2 + abs(mean_selectivity) * 15)

  ggplot(tw, aes(mean_ti, mean_selectivity, color = class_label, size = size)) +
    geom_point(alpha = 0.75) +
    ggrepel::geom_text_repel(
              data = filter(tw, mean_ti > 2 | class_label %in% c("HIGH","PAN")),
              aes(label = paste0(driver, "->", paralog)),
              size = 2.5, show.legend = FALSE, max.overlaps = 20,
              box.padding = 0.5, point.padding = 0.8, force = 2,
              min.segment.length = 0.2, segment.size = 0.3, seed = 42) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = GRAY, alpha = 0.4) +
    geom_vline(xintercept = 1, linewidth = 0.3, color = GRAY, alpha = 0.3, linetype = "dashed") +
    scale_color_manual(values = c(HIGH = RED, MODERATE = ORANGE, LOW = BLUE, PAN = GRAY),
                      limits = c("HIGH","MODERATE","LOW","PAN"), drop = FALSE) +
    scale_size(range = c(1.5, 6), guide = "none") +
    labs(x = "Dependency Window Score (DWS)", y = "Selectivity") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme_sci + theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== Fig4 Panel Generation (R) ===")
pb <- panel_b()

save_panel(pb, "b")

# Single-panel figure (former panel b only) — no panel letter
p <- pb

ggsave(file.path(OUT_DIR, "Fig4_Translational.pdf"), p,
       width = FIG_W, height = FIG_H, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig4_Translational.svg"), p,
       width = FIG_W, height = FIG_H, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "Fig4_Translational.tiff"), p,
       width = FIG_W, height = FIG_H, units = "mm", device = ragg::agg_tiff, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig4_Translational.png"), p,
       width = FIG_W, height = FIG_H, units = "mm", device = ragg::agg_png, dpi = 300)
REVIEW_DIR <- "figure_review"
dir.create(REVIEW_DIR, showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(OUT_DIR, "Fig4_Translational.png"),
          file.path(REVIEW_DIR, "Fig4_Translational.png"), overwrite = TRUE)
message("Fig4_Translational.pdf (120x95mm) ✓")
