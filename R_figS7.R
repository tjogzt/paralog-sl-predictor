# FigS7 — Therapeutic Window Analysis (R)
# Purpose: 4×90mm panels → 180×180mm composite
# Usage:   Rscript R_figS7.R
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)

OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7
PANEL_W <- 90; PANEL_H <- 90

RED   <- "#CB181D"; BLUE  <- "#2171B5"; GRAY  <- "#636363"
GREEN <- "#238B45"; ORANGE <- "#E6550D"; DARK <- "#252525"

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
  ggsave(file.path(OUT_DIR, paste0("FigS7_panel_", name, ".pdf")), p,
         width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
  message(sprintf("  panel %s ✓", name))
}

# ── Data ──
tw_path <- "paralog_sl_predictor/output/therapeutic_window_paralog_classification.csv"
tw_all_path <- "paralog_sl_predictor/output/therapeutic_window_all_results.csv"

if (file.exists(tw_path)) {
  tw <- read_csv(tw_path, show_col_types = FALSE)
} else { tw <- NULL }
if (file.exists(tw_all_path)) {
  tw_all <- read_csv(tw_all_path, show_col_types = FALSE)
} else { tw_all <- NULL }

class_colors <- c(HIGH = RED, MODERATE = ORANGE, LOW = BLUE, PAN = GRAY)
class_labels <- c(HIGH_SELECTIVITY = "HIGH", MODERATE = "MODERATE",
                  LOW_SELECTIVITY = "LOW", PAN_ESSENTIAL = "PAN")

# ═══════════════════════════════════════════════════════════════
# PANEL A — TI Ranking
# ═══════════════════════════════════════════════════════════════
panel_a <- function() {
  if (is.null(tw)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=5) + theme_void())
  }
  tw_sorted <- tw %>% arrange(mean_ti)
  # two-line labels (house style, same as FigS6c); margin(l) guarantees headroom
  # against the y-label column's slight underestimate of rendered text width,
  # which clipped leading characters at the figure's left edge
  tw_sorted$label <- factor(paste0(tw_sorted$driver, "->\n", tw_sorted$paralog),
                            levels = paste0(tw_sorted$driver, "->\n", tw_sorted$paralog))
  tw_sorted$class_short <- factor(tw_sorted$classification,
                                  levels = names(class_labels),
                                  labels = class_labels)

  ggplot(tw_sorted, aes(mean_ti, label, fill = class_short)) +
    geom_col(width = 0.55) +
    scale_fill_manual(values = class_colors, drop = FALSE) +
    labs(x = "Mean Therapeutic Index (TI)", y = NULL) +
    theme_sci + theme(legend.position = "bottom",
                      plot.margin = margin(10, 4, 4, 4, "pt"),
                      axis.text.y = element_text(size = 7, lineheight = 0.85,
                                                 margin = margin(l = 10, r = 2)),
                      legend.key.size = unit(3.5, "mm"),
                      legend.spacing.x = unit(3, "mm"),
                      legend.text = element_text(size = 7, margin = margin(l = 1, r = 4, unit = "pt")))
}

# ═══════════════════════════════════════════════════════════════
# PANEL B — Safety Classification Pie
# ═══════════════════════════════════════════════════════════════
panel_b <- function() {
  if (is.null(tw)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=5) + theme_void())
  }
  counts <- tw %>% count(classification)
  counts$class_short <- factor(counts$classification,
                               levels = names(class_labels),
                               labels = class_labels)
  counts$pct <- counts$n / sum(counts$n) * 100
  counts$ymax <- cumsum(counts$pct)
  counts$ymin <- c(0, head(counts$ymax, -1))
  counts$label_pos <- (counts$ymax + counts$ymin) / 2

  ggplot(counts, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = class_short)) +
    geom_rect(color = "white", linewidth = 0.3) +
    geom_text(aes(x = 4.5, y = label_pos, label = sprintf("%s\n%d (%.0f%%)", class_short, n, pct)),
              size = 2.5, hjust = 0.5) +
    scale_fill_manual(values = class_colors, drop = FALSE) +
    coord_polar(theta = "y") +
    xlim(c(1, 4.8)) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(0, 0, 0, 0, "pt"),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
}

# ═══════════════════════════════════════════════════════════════
# PANEL C — TI by Cancer Context
# ═══════════════════════════════════════════════════════════════
panel_c <- function() {
  if (is.null(tw_all)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=5) + theme_void())
  }
  ctx <- tw_all %>%
    mutate(context_short = substr(context, 1, 12))

  ggplot(ctx, aes(context_short, therapeutic_index)) +
    geom_boxplot(fill = BLUE, alpha = 0.3, outlier.size = 0.8, linewidth = 0.3) +
    geom_hline(aes(yintercept = 1.0, linetype = "TI = 1"), linewidth = 0.4, color = RED, alpha = 0.5) +
    scale_linetype_manual(values = c("TI = 1" = "dashed")) +
    labs(x = NULL, y = "Therapeutic Index") +
    theme_sci + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                      legend.position = c(0.98, 0.98), legend.justification = c(1, 1))
}

# ═══════════════════════════════════════════════════════════════
# PANEL D — Selectivity vs TI Scatter
# ═══════════════════════════════════════════════════════════════
panel_d <- function() {
  if (is.null(tw)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=5) + theme_void())
  }
  tw_d <- tw %>% filter(!is.na(mean_selectivity), !is.na(mean_ti))
  tw_d$class_short <- factor(tw_d$classification, levels = names(class_labels), labels = class_labels)

  ggplot(tw_d, aes(mean_selectivity, mean_ti, color = class_short)) +
    geom_point(aes(size = abs(mean_dd) * 50), alpha = 0.75) +
    ggrepel::geom_text_repel(
      data = tw_d %>% group_by(class_short) %>% slice_max(abs(mean_dd), n = 1),
      aes(label = paste0(driver, "->", paralog)),
      size = 2.5, show.legend = FALSE, max.overlaps = 20, box.padding = 0.3,
      point.padding = 0.5) +
    geom_hline(yintercept = 1, linewidth = 0.3, color = GRAY, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0.15, linewidth = 0.3, color = GRAY, linetype = "dashed", alpha = 0.3) +
    scale_color_manual(values = class_colors, drop = FALSE) +
    scale_size(range = c(1.5, 6), guide = "none") +
    labs(x = "Selectivity", y = "Therapeutic Index (TI)") +
    theme_sci + theme(legend.position = "bottom")
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== FigS7 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b(); pc <- panel_c(); pd <- panel_d()

save_panel(pa, "a"); save_panel(pb, "b"); save_panel(pc, "c"); save_panel(pd, "d")

# Left column: TI ranking (full height, two-line labels); right column: b/c/d
# stacked. Nested plot_grid (not ggdraw/draw_plot): draw_plot anchored the
# panel region to the viewport edge and the y-label column was clipped at the
# figure's left edge regardless of cell width or label margins.
right_col <- cowplot::plot_grid(pb, pc, pd, ncol = 1,
                                labels = c("b","c","d"),
                                label_size = 9, label_fontface = "bold",
                                label_fontfamily = "Arial")
p <- cowplot::plot_grid(pa, right_col, ncol = 2, rel_widths = c(0.42, 0.58),
                        labels = c("a",""),
                        label_size = 9, label_fontface = "bold",
                        label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "FigS7_TherapeuticWindow.pdf"), p,
       width = 180, height = 180, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS7_TherapeuticWindow.svg"), p,
       width = 180, height = 180, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS7_TherapeuticWindow.tiff"), p,
       width = 180, height = 180, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS7_TherapeuticWindow.pdf (180×180mm) ✓")
