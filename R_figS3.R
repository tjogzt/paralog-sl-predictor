# FigS3 — CPTAC Per-Cohort Scatter Plots (R, real UCEC data)
# Purpose: 3×60mm panels → 180×60mm, with regression lines
library(ggplot2)
library(cowplot)
library(dplyr)
library(jsonlite)

OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FS <- 7; TICK_FS <- 6; LEGEND_FS <- 6
PANEL_W <- 60; PANEL_H <- 60

BLUE  <- "#2171B5"; RED   <- "#CB181D"; GRAY  <- "#636363"

theme_sci <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.line    = element_line(linewidth = 0.35),
  axis.ticks   = element_line(linewidth = 0.3),
  axis.text    = element_text(size = TICK_FS),
  axis.title   = element_text(size = BASE_FS),
  plot.margin  = margin(3, 3, 3, 3, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

# Load real UCEC CPTAC data
d <- fromJSON("paralog_sl_predictor/data/cptac_cache/UCEC_protein_data.json")

make_real_scatter <- function(gene_a, gene_b) {
  x <- unlist(d[[gene_a]])
  y <- unlist(d[[gene_b]])
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)

  r_val <- cor(x, y)
  t_stat <- r_val * sqrt((n - 2) / (1 - r_val^2))
  p_val  <- 2 * pt(abs(t_stat), n - 2, lower.tail = FALSE)
  sig_label <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "(ns)"
  title_lab <- sprintf("%s/%s | r=%.3f %s", gene_a, gene_b, r_val, sig_label)

  df <- tibble(x = x, y = y)
  ggplot(df, aes(x, y)) +
    geom_point(size = 1, alpha = 0.4, color = BLUE) +
    geom_smooth(method = "lm", se = FALSE, color = RED, linewidth = 0.5) +
    labs(x = gene_a, y = gene_b, title = title_lab) +
    theme_sci + theme(plot.title = element_text(size = 7, face = "bold"))
}

# ── MAIN ──
message("=== FigS3 (Real UCEC CPTAC Data) ===")
# Use pairs that exist in the JSON
pa <- make_real_scatter("EP300", "CREBBP")
pb <- make_real_scatter("PIK3CA", "PIK3CB")
pc <- make_real_scatter("PIK3R1", "CRKL")

ggsave(file.path(OUT_DIR, "FigS3_panel_a.pdf"), pa, width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS3_panel_b.pdf"), pb, width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS3_panel_c.pdf"), pc, width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
message("  panels saved")

p <- ggdraw() +
  draw_plot(pa, x = 0,     y = 0, width = 1/3, height = 1) +
  draw_plot(pb, x = 1/3,   y = 0, width = 1/3, height = 1) +
  draw_plot(pc, x = 2/3,   y = 0, width = 1/3, height = 1) +
  draw_plot_label(c("a","b","c"), x = c(0, 1/3, 2/3), y = c(1, 1, 1),
                  size = 8, fontface = "bold")

ggsave(file.path(OUT_DIR, "FigS3_CPTAC_PerCohort.pdf"), p,
       width = 180, height = 60, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS3_CPTAC_PerCohort.svg"), p,
       width = 180, height = 60, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS3_CPTAC_PerCohort.tiff"), p,
       width = 180, height = 60, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS3_CPTAC_PerCohort.pdf (180×60mm) ✓")
