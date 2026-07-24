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

BASE_FS <- 7; TICK_FS <- 6; LEGEND_FS <- 5.5
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
  df_raw <- tibble(
    cancer = c("Biliary Tract","Mesothelioma","Colorectal","Esophagogastric","SCLC",
               "Pancreatic","HNSCC","Breast","Neuroblastoma","Bladder Urothelial",
               "Hepatocellular","Ovarian","Cervical","Glioma","Endometrial",
               "NSCLC","Melanoma"),
    dd_auroc = c(.960,.917,.812,.805,.750,.750,.746,.734,.722,.698,
                 .683,.674,.643,.600,.576,.572,.462))
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
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "DD AUROC", y = NULL, title = "Strong signal") +
    theme_sci + theme(plot.title = element_text(size = 7, face = "bold"))
}

# ── Panel B: Moderate AUROC (0.5-0.7) ──
panel_b <- function() {
  df <- df_raw %>% filter(dd_auroc >= 0.5 & dd_auroc < 0.7)
  df$cancer_short <- factor(df$cancer_short, levels = df$cancer_short)
  ggplot(df, aes(dd_auroc, cancer_short)) +
    geom_col(fill = BLUE, width = 0.6) +
    geom_text(aes(label = sprintf("%.3f", dd_auroc)), hjust = -0.1, size = 2.5, color = BLUE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "DD AUROC", y = NULL, title = "Moderate signal") +
    theme_sci + theme(plot.title = element_text(size = 7, face = "bold"))
}

# ── Panel C: Weak AUROC (<0.5) ──
panel_c <- function() {
  df <- df_raw %>% filter(dd_auroc < 0.5)
  df$cancer_short <- factor(df$cancer_short, levels = df$cancer_short)
  ggplot(df, aes(dd_auroc, cancer_short)) +
    geom_col(fill = GRAY, width = 0.6) +
    geom_text(aes(label = sprintf("%.3f", dd_auroc)), hjust = -0.1, size = 2.5, color = GRAY) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "DD AUROC", y = NULL, title = "Weak signal") +
    theme_sci + theme(plot.title = element_text(size = 7, face = "bold"))
}

# ── MAIN ──
message("=== FigS2 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b(); pc <- panel_c()

ggsave(file.path(OUT_DIR, "FigS2_panel_a.pdf"), pa, width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS2_panel_b.pdf"), pb, width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS2_panel_c.pdf"), pc, width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
message("  panels saved")

# plot_grid allocates widths properly — fixed draw_plot cells clipped the
# right-most value labels when full lineage names widened panel b's gtable.
p <- cowplot::plot_grid(pa, pb, pc, nrow = 1,
                        labels = c("a","b","c"),
                        label_size = 8, label_fontface = "bold")

ggsave(file.path(OUT_DIR, "FigS2_CrossCancer_AUROC.pdf"), p,
       width = 200, height = 60, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS2_CrossCancer_AUROC.svg"), p,
       width = 200, height = 60, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS2_CrossCancer_AUROC.tiff"), p,
       width = 200, height = 60, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS2_CrossCancer_AUROC.pdf (200×60mm) ✓")
