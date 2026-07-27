# Fig4 — Drug + Dependency-Window Classification (R)
# Purpose: Generate individual 90×95mm panels → review → 180×95mm composite
# Note:    former panels c/d (structure-derived descriptors) moved to R_figS13.R
# Usage:   Rscript R_fig4.R

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)

# ── Constants ──
BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7; ANNOT_FS <- 5.5
PANEL_W   <- 90; PANEL_H <- 95  # mm
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
# PANEL A — PRISM Drug Selectivity
# ═══════════════════════════════════════════════════════════════
panel_a <- function() {
  prism_path <- "paralog_sl_predictor/output/prism_top_hits.csv"
  if (file.exists(prism_path)) {
    pr_raw <- read_csv(prism_path, show_col_types = FALSE) %>%
      mutate(
        abs_delta = abs(delta_auc),
        drug_upper = toupper(drug),
        drug_class = case_when(
          grepl("MEK|AZD8330|TRAMETINIB|RO-4987655", drug_upper) ~ "MEKi",
          grepl("MTOR|EVEROLIMUS|TEMSIROLIMUS|AKT|IPATASERTIB|GSK-2141795", drug_upper) ~ "mTOR/AKTi",
          grepl("HDAC|PANOBINOSTAT", drug_upper) ~ "HDACi",
          TRUE ~ "Other"))
    # Pick top 1-2 per class by |ΔAUC|
    pr <- bind_rows(
      pr_raw %>% filter(drug_class == "MEKi")      %>% slice_max(abs_delta, n = 2),
      pr_raw %>% filter(drug_class == "mTOR/AKTi")  %>% slice_max(abs_delta, n = 2),
      pr_raw %>% filter(drug_class == "HDACi")      %>% slice_max(abs_delta, n = 2),
      pr_raw %>% filter(drug_class == "Other")      %>% slice_max(abs_delta, n = 3)
    ) %>% distinct(drug, driver, paralog, .keep_all = TRUE) %>%
      slice_max(abs_delta, n = 8)
  } else {
    stop("paralog_sl_predictor/output/prism_top_hits.csv not found — ",
         "run the pipeline first; simulated fallbacks are forbidden")
  }
  pr <- pr %>% mutate(
    label = paste0(drug, "\n", driver, "->", paralog),
    drug_class = factor(drug_class, levels = c("MEKi","mTOR/AKTi","HDACi","Other")))
  pr$label <- factor(pr$label, levels = rev(unique(pr$label)))

  ggplot(pr, aes(abs_delta, label, fill = drug_class)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c(MEKi = RED, `mTOR/AKTi` = BLUE, HDACi = ORANGE, Other = GREEN),
                      drop = FALSE) +
    labs(x = "|ΔAUC|", y = NULL) +
    theme_sci + theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
                      plot.margin = margin(4, 4, 8, 20, "pt"))
}

# ═══════════════════════════════════════════════════════════════
# PANEL B — Therapeutic Window
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
pa <- panel_a(); pb <- panel_b()

save_panel(pa, "a"); save_panel(pb, "b")

p <- cowplot::plot_grid(pa, pb, ncol = 2, rel_widths = c(0.5, 0.5),
                        labels = c("a","b"),
                        label_size = 9, label_fontface = "bold",
                        label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "Fig4_Translational.pdf"), p,
       width = 180, height = 95, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig4_Translational.svg"), p,
       width = 180, height = 95, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "Fig4_Translational.tiff"), p,
       width = 180, height = 95, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("Fig4_Translational.pdf (180x95mm) ✓")
