# FigS6 — Mutation Type Stratification (R)
# Purpose: 4×90mm panels → 180×180mm composite
# Usage:   Rscript R_figS6.R
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
GREEN <- "#238B45"; DARK <- "#252525"

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
  ggsave(file.path(OUT_DIR, paste0("FigS6_panel_", name, ".pdf")), p,
         width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
  message(sprintf("  panel %s ✓", name))
}

# ── Data ──
mut_path <- "paralog_sl_predictor/output/muttype_all_results.csv"
if (file.exists(mut_path)) {
  mut <- read_csv(mut_path, show_col_types = FALSE)
} else {
  mut <- NULL
}

# ═══════════════════════════════════════════════════════════════
# PANEL A — Per-Cancer AUC Comparison
# ═══════════════════════════════════════════════════════════════
panel_a <- function() {
  if (is.null(mut)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=5) + theme_void() + theme(plot.background = element_rect(fill="white", color=NA)))
  }
  cancers <- c("Ovarian","Endometrial","Colorectal","Breast")
  plot_data <- data.frame()
  for (canc in cancers) {
    sub <- mut[mut$cancer == canc, ]
    if (nrow(sub) < 5) next
    yt <- as.integer(sub$is_known_sl)
    if (sum(yt) < 2) next
    for (col in c("dd_all","dd_trunc","dd_miss")) {
      scores <- abs(as.numeric(sub[[col]]))
      scores[is.na(scores)] <- 0
      if (length(unique(yt)) < 2) next
      auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(yt, scores))),
                      error = function(e) NA)
      plot_data <- rbind(plot_data, data.frame(cancer = canc, type = col, auc = auc))
    }
  }
  if (nrow(plot_data) == 0) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Insufficient data", size=5) + theme_void() + theme(plot.background = element_rect(fill="white", color=NA)))
  }
  plot_data$type <- factor(plot_data$type, levels = c("dd_all","dd_trunc","dd_miss"),
                           labels = c("All","Truncating","Missense"))
  plot_data$cancer <- factor(plot_data$cancer, levels = cancers)
  type_colors <- c(All = GRAY, Truncating = RED, Missense = BLUE)

  ggplot(plot_data, aes(cancer, auc, fill = type)) +
    geom_col(position = position_dodge(0.7), width = 0.55) +
    scale_fill_manual(values = type_colors) +
    geom_hline(yintercept = 0.5, linewidth = 0.3, color = GRAY, linetype = "dashed", alpha = 0.3) +
    labs(x = NULL, y = "DD AUROC") +
    theme_sci + theme(legend.position = "bottom")
}

# (panel a — legend below)

# ═══════════════════════════════════════════════════════════════
# PANEL B — |DD| Magnitude: Truncating − Missense
# ═══════════════════════════════════════════════════════════════
panel_b <- function() {
  if (is.null(mut)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=5) + theme_void() + theme(plot.background = element_rect(fill="white", color=NA)))
  }
  sub <- mut[!is.na(mut$dd_trunc) & !is.na(mut$dd_miss), ]
  if (nrow(sub) < 10) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="n<10", size=5) + theme_void() + theme(plot.background = element_rect(fill="white", color=NA)))
  }
  diff_vals <- abs(sub$dd_trunc) - abs(sub$dd_miss)
  df <- data.frame(diff = diff_vals)
  mean_diff <- mean(diff_vals, na.rm = TRUE)

  ggplot(df, aes(diff)) +
    geom_histogram(bins = 25, fill = "#8E44AD", alpha = 0.7, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = DARK) +
    geom_vline(xintercept = mean_diff, linewidth = 0.5, color = RED, linetype = "dashed") +
    annotate("text", x = mean_diff + 0.08, y = Inf, label = sprintf("Mean = %+.3f", mean_diff),
             size = 2.8, color = RED, hjust = 0, vjust = 1.5) +
    # right-side headroom keeps the outermost tick label inside the panel
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.14))) +
    labs(x = "|DD_trunc| − |DD_miss|", y = "Frequency") +
    theme_sci
}

# ═══════════════════════════════════════════════════════════════
# PANEL C — Known SL Pairs: Truncating vs Missense DD
# ═══════════════════════════════════════════════════════════════
panel_c <- function() {
  if (is.null(mut)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=5) + theme_void() + theme(plot.background = element_rect(fill="white", color=NA)))
  }
  known <- mut[mut$is_known_sl == TRUE & !is.na(mut$dd_trunc) & !is.na(mut$dd_miss), ]
  if (nrow(known) < 2) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Known SL pairs insufficient", size=4) + theme_void() + theme(plot.background = element_rect(fill="white", color=NA)))
  }
  # Deduplicate by driver-paralog pair, take max |DD_trunc|
  known <- known %>%
    group_by(driver, paralog) %>%
    summarise(dd_trunc = max(abs(dd_trunc), na.rm = TRUE),
              dd_miss  = max(abs(dd_miss), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(dd_trunc)) %>% head(8)
  # two-line labels keep the combined gtable within the 180mm device width
  # (one-line labels overflowed and were clipped at the figure's left edge)
  known$pair_label <- paste0(known$driver, "->\n", known$paralog)
  known$pair_label <- factor(known$pair_label, levels = rev(known$pair_label))

  df <- bind_rows(
    data.frame(pair = known$pair_label, type = "Truncating", dd = abs(known$dd_trunc)),
    data.frame(pair = known$pair_label, type = "Missense",   dd = abs(known$dd_miss)))
  df$type <- factor(df$type, levels = c("Truncating","Missense"))
  df$pair <- factor(df$pair, levels = known$pair_label)

  ggplot(df, aes(dd, pair, fill = type)) +
    geom_col(position = position_dodge(0.7), width = 0.55) +
    # Pair labels drawn INSIDE the plot above each bar group: y-axis text
    # was clipped at the composite figure's left edge (same issue as Fig3c)
    geom_text(data = known,
              aes(x = 0, y = as.numeric(pair_label) + 0.45,
                  label = paste0(driver, "->", paralog)),
              hjust = 0, size = 2.3, family = "Arial", color = DARK,
              inherit.aes = FALSE) +
    scale_fill_manual(values = c(Truncating = RED, Missense = BLUE)) +
    labs(x = "|DD|", y = NULL) +
    theme_sci + theme(legend.position = "bottom",
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.line.y = element_blank())
}

# (panel c — legend below)

# ═══════════════════════════════════════════════════════════════
# PANEL D — Breast Cancer ROC Exception
# ═══════════════════════════════════════════════════════════════
panel_d <- function() {
  if (is.null(mut)) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No data", size=5) + theme_void() + theme(plot.background = element_rect(fill="white", color=NA)))
  }
  breast <- mut[mut$cancer == "Breast", ]
  if (nrow(breast) < 5) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Breast: insufficient data", size=4) + theme_void() + theme(plot.background = element_rect(fill="white", color=NA)))
  }
  yt <- as.integer(breast$is_known_sl)
  if (sum(yt) < 2) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="Breast: known SL < 2", size=4) + theme_void() + theme(plot.background = element_rect(fill="white", color=NA)))
  }

  roc_data <- data.frame()
  auc_labels <- c()
  auc_colors <- c()
  for (col in c("dd_all","dd_trunc","dd_miss")) {
    scores <- abs(as.numeric(breast[[col]]))
    scores[is.na(scores)] <- 0
    if (length(unique(yt)) < 2) next
    roc_obj <- tryCatch(pROC::roc(yt, scores), error = function(e) NULL)
    if (is.null(roc_obj)) next
    auc_val <- as.numeric(pROC::auc(roc_obj))
    label <- switch(col, dd_all = "All", dd_trunc = "Truncating", dd_miss = "Missense")
    label <- sprintf("%s (AUC=%.3f)", label, auc_val)
    df_roc <- data.frame(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities,
                         type = label)
    roc_data <- rbind(roc_data, df_roc)
  }
  if (nrow(roc_data) == 0) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="ROC computation failed", size=5) + theme_void() + theme(plot.background = element_rect(fill="white", color=NA)))
  }

  roc_colors <- c()
  for (lvl in unique(roc_data$type)) {
    if (grepl("All", lvl)) roc_colors[lvl] <- GRAY
    else if (grepl("Truncating", lvl)) roc_colors[lvl] <- RED
    else roc_colors[lvl] <- BLUE
  }

  ggplot(roc_data, aes(fpr, tpr, color = type)) +
    geom_line(linewidth = 0.6) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, color = GRAY, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = roc_colors) +
    labs(x = "False Positive Rate", y = "True Positive Rate",
         title = "Breast Cancer") +
    theme_sci + theme(plot.title = element_text(size = 7, face = "bold"),
                      legend.position = c(0.98, 0.02), legend.justification = c(1, 0))
}

# (panel d — legend bottom-right with AUC)

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== FigS6 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b(); pc <- panel_c(); pd <- panel_d()

save_panel(pa, "a"); save_panel(pb, "b"); save_panel(pc, "c"); save_panel(pd, "d")

# plot_grid measures each gtable and allocates widths — long y labels in
# panel c are no longer clipped (draw_plot with fixed cells clipped them).
p <- cowplot::plot_grid(pa, pb, pc, pd, ncol = 2,
                        labels = c("a","b","c","d"),
                        label_size = 9, label_fontface = "bold", label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "FigS6_MutationType.pdf"), p,
       width = 180, height = 180, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS6_MutationType.svg"), p,
       width = 180, height = 180, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS6_MutationType.tiff"), p,
       width = 180, height = 180, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS6_MutationType.pdf (180×180mm) ✓")
