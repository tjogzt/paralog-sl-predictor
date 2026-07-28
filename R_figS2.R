# FigS2 — Cross-Cancer DD AUROC Detail (R)
# Purpose: 3×60mm square panels → 180×60mm composite
# Usage:   Rscript R_figS2.R
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)

OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7
PANEL_W <- 60; PANEL_H <- 60

BLUE  <- "#2171B5"; RED   <- "#CB181D"
GRAY  <- "#636363"; DARK <- "#252525"

theme_sci <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.line    = element_line(linewidth = 0.35),
  axis.ticks   = element_line(linewidth = 0.3),
  axis.text    = element_text(size = TICK_FS),
  axis.title   = element_text(size = BASE_FS),
  legend.text  = element_text(size = LEGEND_FS),
  legend.title = element_blank(),
  legend.background = element_blank(),
  legend.key        = element_blank(),
  plot.margin  = margin(3, 3, 3, 3, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

# Load real data
solid_path <- "paralog_sl_predictor/output/solid_tumor_summary.csv"
if (file.exists(solid_path)) {
  df_raw <- read_csv(solid_path, show_col_types = FALSE) %>% filter(!is.na(dd_auroc))
} else {
  stop("paralog_sl_predictor/output/solid_tumor_summary.csv not found — ",
       "run the pipeline first; simulated fallbacks are forbidden")
}

df_raw <- df_raw %>% arrange(dd_auroc)
df_raw$cancer_short <- df_raw$cancer  # full names; abbreviate() produced unreadable labels

# ── Panel A: High AUROC (>= 0.7) ──
panel_a <- function() {
  df <- df_raw %>% filter(dd_auroc >= 0.7)
  df$cancer_short <- factor(df$cancer_short, levels = df$cancer_short)
  ggplot(df, aes(dd_auroc, cancer_short)) +
    geom_col(fill = RED, width = 0.6) +
    geom_text(aes(label = sprintf("%.3f", dd_auroc)), hjust = -0.1, size = 2.5, color = RED) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(x = "DD AUROC", y = NULL, title = "Strong signal") +
    theme_sci + theme(plot.title = element_text(size = 7, face = "bold"))
}

# ── Panel B: Moderate (0.5-0.7) plus below-chance (< 0.5) AUROC ──
# Under signed-DD scoring Breast falls below 0.5 on the primary frame; it is
# shown in gray so no evaluable lineage is silently dropped.
panel_b <- function() {
  df <- df_raw %>% filter(dd_auroc < 0.7)
  df$cancer_short <- factor(df$cancer_short, levels = df$cancer_short)
  df$clr <- ifelse(df$dd_auroc >= 0.5, BLUE, GRAY)
  ggplot(df, aes(dd_auroc, cancer_short)) +
    geom_col(aes(fill = clr), width = 0.6) +
    scale_fill_identity() +
    geom_text(aes(label = sprintf("%.3f", dd_auroc), color = clr), hjust = -0.1, size = 2.5) +
    scale_color_identity() +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = GRAY, alpha = 0.6) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(x = "DD AUROC", y = NULL, title = "Moderate / below-chance") +
    theme_sci + theme(plot.title = element_text(size = 7, face = "bold"))
}

# ── Panel C is not used: below-chance lineages are folded into panel b. ──

# ── MAIN ──
message("=== FigS2 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b()

ggsave(file.path(OUT_DIR, "FigS2_panel_a.pdf"), pa, width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS2_panel_b.pdf"), pb, width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
message("  panels saved")

# plot_grid allocates widths properly — fixed draw_plot cells clipped the
# right-most value labels when full lineage names widened panel b's gtable.
p <- cowplot::plot_grid(pa, pb, nrow = 1,
                        labels = c("a","b"),
                        label_size = 8, label_fontface = "bold", label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "FigS2_CrossCancer_AUROC.pdf"), p,
       width = 120, height = 60, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS2_CrossCancer_AUROC.svg"), p,
       width = 120, height = 60, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS2_CrossCancer_AUROC.tiff"), p,
       width = 120, height = 60, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS2_CrossCancer_AUROC.pdf (120x60mm, 2 panels; below-chance lineages shown in panel b) ✓")
