# FigS8 — Therapeutic Window Analysis (R)
# Purpose: 3 panels (a DWS ranking / b tier donut / c by-context) → 180×130mm composite
# Note:    former panel d (selectivity-vs-DWS scatter) removed 2026-08-03 — it
#          duplicated main-text Fig. 4 (same 21 pairs, same two variables,
#          axes swapped); its unique elements (|mean DD| bubble size, both
#          HIGH thresholds) were absorbed into Fig. 4
# Usage:   Rscript R_figS8_dws.R
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
  ggsave(file.path(OUT_DIR, paste0("FigS8_panel_", name, ".pdf")), p,
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
# PANEL A — DWS Ranking
# ═══════════════════════════════════════════════════════════════
panel_a <- function() {
  if (is.null(tw)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=5) + theme_void())
  }
  tw_sorted <- tw %>% arrange(mean_ti)
  # 2026-08-03: single-line labels. With panel d removed the composite shrank
  # to 180x130mm; 12 of 21 pairs have mean DWS ~ 0 (legitimate "no window"
  # negatives kept for completeness with Table S5), so two-line labels wasted
  # vertical space. margin(l) still guards the left-edge clipping fix.
  tw_sorted$label <- factor(paste0(tw_sorted$driver, "->", tw_sorted$paralog),
                            levels = paste0(tw_sorted$driver, "->", tw_sorted$paralog))
  tw_sorted$class_short <- factor(tw_sorted$classification,
                                  levels = names(class_labels),
                                  labels = class_labels)

  ggplot(tw_sorted, aes(mean_ti, label, fill = class_short)) +
    geom_col(width = 0.55) +
    scale_fill_manual(values = class_colors, drop = FALSE) +
    labs(x = "Mean Dependency Window Score (DWS)", y = NULL) +
    theme_sci + theme(legend.position = "bottom",
                      # left margin 8pt + text margin(l=16): single-line labels
                      # are ~2x wider than the old two-line ones; the gutter
                      # underestimates rendered width and clipped the leading
                      # character ("SMARCA4" -> "MARCA4") at the cell edge
                      plot.margin = margin(10, 4, 4, 8, "pt"),
                      axis.text.y = element_text(size = 7,
                                                 margin = margin(l = 16, r = 2)),
                      legend.key.size = unit(3, "mm"),
                      legend.spacing.x = unit(0.5, "mm"),
                      # gap AFTER each label lives in the label's right margin,
                      # so every label sits tight against its own key
                      legend.text = element_text(size = 7, margin = margin(l = 2, r = 5, unit = "pt")))
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
  # Nudge the small top segments' labels sideways (HIGH -> upper right,
  # PAN -> upper left) so the top of the cell holds only the ring; the
  # xlim can then be tightened and the donut grows without stealing
  # height from panel c
  counts$label_y <- counts$label_pos
  counts$label_y[counts$class_short == "HIGH"] <- 7
  counts$label_y[counts$class_short == "PAN"]  <- 86
  # outward horizontal justification: labels extend AWAY from the ring
  counts$hj <- ifelse(counts$class_short %in% c("HIGH", "LOW"), 0, 1)

  ggplot(counts, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = class_short)) +
    geom_rect(color = "white", linewidth = 0.3) +
    geom_text(aes(x = 4.28, y = label_y, hjust = hj,
                  label = sprintf("%s\n%d (%.0f%%)", class_short, n, pct)),
              size = 2.5) +
    scale_fill_manual(values = class_colors, drop = FALSE) +
    coord_polar(theta = "y", clip = "off") +
    # labels sit just beyond the scale range; oob_keep retains them (plain
    # xlim() would drop them as NA) so the ring gets ~98% of the cell
    scale_x_continuous(limits = c(1.2, 4.1), oob = scales::oob_keep) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(2, 10, 2, 10, "pt"),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
}

# ═══════════════════════════════════════════════════════════════
# PANEL C — DWS by Cancer Context
# ═══════════════════════════════════════════════════════════════
panel_c <- function() {
  if (is.null(tw_all)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=5) + theme_void())
  }
  ctx <- tw_all %>%
    mutate(context_short = substr(context, 1, 12))

  ggplot(ctx, aes(context_short, therapeutic_index)) +
    geom_boxplot(fill = BLUE, alpha = 0.3, outlier.size = 0.8, linewidth = 0.3) +
    geom_hline(aes(yintercept = 1.0, linetype = "DWS = 1"), linewidth = 0.4, color = RED, alpha = 0.5) +
    scale_linetype_manual(values = c("DWS = 1" = "dashed")) +
    labs(x = NULL, y = "Dependency Window Score") +
    theme_sci + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                      # 2026-07-30 review: the inside top-right legend overlapped
                      # the PanCancer outlier point (~9.3); move it above the
                      # panel (own allocated strip, no data collision)
                      legend.position = "top",
                      legend.justification = "right",
                      legend.direction = "horizontal",
                      legend.margin = margin(0, 0, 0, 0),
                      legend.box.margin = margin(0, 0, -3, 0))
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== FigS8 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b(); pc <- panel_c()

save_panel(pa, "a"); save_panel(pb, "b"); save_panel(pc, "c")

# Left column: DWS ranking (full height, single-line labels); right column: b/c
# stacked. Nested plot_grid (not ggdraw/draw_plot): draw_plot anchored the
# panel region to the viewport edge and the y-label column was clipped at the
# figure's left edge regardless of cell width or label margins.
right_col <- cowplot::plot_grid(pb, pc, ncol = 1,
                                labels = c("b","c"),
                                label_size = 9, label_fontface = "bold",
                                label_fontfamily = "Arial")
p <- cowplot::plot_grid(pa, right_col, ncol = 2, rel_widths = c(0.42, 0.58),
                        labels = c("a",""),
                        label_size = 9, label_fontface = "bold",
                        label_fontfamily = "Arial")

# 180x130mm: 21 single-line rows need ~5.4mm each; right cells stay ~65mm
# (donut still larger than in the 4-panel 180x180 era)
ggsave(file.path(OUT_DIR, "FigS8_TherapeuticWindow.pdf"), p,
       width = 180, height = 130, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS8_TherapeuticWindow.svg"), p,
       width = 180, height = 130, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS8_TherapeuticWindow.tiff"), p,
       width = 180, height = 130, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS8_TherapeuticWindow.pdf (180×130mm) ✓")
