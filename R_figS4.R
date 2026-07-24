# FigS4 — CNV Independence Analysis (R, real DepMap data)
# Purpose: Scatter CNV vs CERES for each paralog gene, show R² < 0.10
# Usage:   Rscript R_figS4.R
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)

OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FS <- 7; TICK_FS <- 6

# Colors
BLUE  <- "#2171B5"; RED   <- "#CB181D"; GRAY  <- "#636363"

theme_sci <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.line    = element_line(linewidth = 0.3),
  axis.ticks   = element_line(linewidth = 0.2),
  axis.text    = element_text(size = TICK_FS),
  axis.title   = element_text(size = TICK_FS),
  plot.margin  = margin(2, 2, 2, 2, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

# ── Load precomputed artifacts (regression_controls.py) ──
# R² values: per-gene lm(CERES ~ CNV) on all CNV∩dep cell lines (n = 1192).
# Scatter points: 300 sampled cell lines per gene (same data, subsampled for plotting).
r2_tab <- read_csv("paralog_sl_predictor/output/cnv_independence.csv", show_col_types = FALSE)
scatter <- read_csv("paralog_sl_predictor/output/cnv_scatter_sample.csv", show_col_types = FALSE)
message(sprintf("Loaded R² for %d genes; %d scatter points",
                nrow(r2_tab), nrow(scatter)))

# Paralog genes to analyze (panel order)
paralog_genes <- c("ARID1A","ARID1B","PIK3CA","PIK3CB","PIK3R1","CRKL",
                   "EP300","CREBBP","KRAS","HRAS","PTEN","TNS2",
                   "SMARCA4","SMARCA2","PPP2R1A","PPP2R1B",
                   "KMT2D","KMT2C","TP53","TP63","FBXW7","FBXW2",
                   "STK11","SIK1")

make_cnv_plot <- function(gene) {
  r2 <- r2_tab$r2[r2_tab$gene == gene]
  df <- scatter %>% filter(gene == !!gene)
  if (length(r2) == 0 || nrow(df) < 10) {
    return(NULL)
  }
  ggplot(df, aes(cnv, dep)) +
    geom_point(size = 0.3, alpha = 0.3, color = BLUE) +
    geom_smooth(method = "lm", se = FALSE, color = RED, linewidth = 0.3) +
    labs(x = "Relative copy number", y = "Gene effect score",
         title = sprintf("%s  R²=%.3f", gene, r2)) +
    theme_sci + theme(plot.title = element_text(size = 6, face = "bold"))
}

# ── MAIN ──
message("Generating CNV plots...")
plots <- lapply(paralog_genes, make_cnv_plot)
plots <- plots[!sapply(plots, is.null)]  # remove genes with no data
n_plots <- length(plots)
message(sprintf("  %d genes plotted (skipped %d with no data)", n_plots, length(paralog_genes) - n_plots))

# Grid: 4 columns, dynamic rows
ncol <- 4
nrow <- ceiling(n_plots / ncol)
h_mm <- nrow * 45  # ~45mm per row
message(sprintf("  Layout: %d×%d, canvas 180×%dmm", ncol, nrow, h_mm))
p <- plot_grid(plotlist = plots, ncol = ncol, align = "hv",
               labels = letters[1:n_plots],
               label_size = 6, label_fontface = "bold")

ggsave(file.path(OUT_DIR, "FigS4_CNV_Independence.pdf"), p,
       width = 180, height = h_mm, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS4_CNV_Independence.svg"), p,
       width = 180, height = 180, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS4_CNV_Independence.tiff"), p,
       width = 180, height = 180, units = "mm", device = ragg::agg_tiff, dpi = 600)
message(sprintf("FigS4_CNV_Independence.pdf (180×%dmm) ✓", h_mm))
