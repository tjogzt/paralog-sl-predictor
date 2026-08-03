# FigS7 — PRISM Drug Selectivity Heatmap + Assay-Validity Anchors (R)
# Purpose: panel a = heatmap of top drug × paralog-SL pair ΔAUC values
#          panel b = targeted-agent assay-validity anchors (|ΔAUC| bar chart,
#                    moved from main-text Fig. 4a after restructuring)
# Canvas:  180×120mm composite, panel letters a/b
# Usage:   Rscript R_figS7_prism.R
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyr)
library(readr)
library(reshape2)

OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FS <- 7; TICK_FS <- 7

# ── Colors ──
BLUE  <- "#2171B5"; RED   <- "#CB181D"; GREEN <- "#238B45"
ORANGE <- "#E6550D"; GRAY  <- "#636363"

theme_sci <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.line    = element_blank(),
  axis.ticks   = element_line(linewidth = 0.3),
  axis.text.x  = element_text(size = TICK_FS, angle = 45, hjust = 1, vjust = 1),
  axis.text.y  = element_text(size = TICK_FS, hjust = 1),
  axis.title   = element_blank(),
  legend.position = "right",
  legend.key.height = unit(0.5, "cm"),
  legend.key.width  = unit(0.3, "cm"),
  plot.margin  = margin(4, 4, 4, 4, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

# Bar-chart theme (matches the former R_fig4.R style: visible axis lines)
theme_bar <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.line    = element_line(linewidth = 0.4),
  axis.ticks   = element_line(linewidth = 0.3),
  axis.text    = element_text(size = TICK_FS),
  axis.title   = element_text(size = BASE_FS),
  legend.text  = element_text(size = BASE_FS),
  legend.title = element_blank(),
  legend.background = element_blank(),
  legend.key        = element_blank(),
  plot.margin  = margin(4, 4, 8, 20, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

# ═══════════════════════════════════════════════════════════════
# PANEL A — PRISM ΔAUC heatmap (unchanged)
# ═══════════════════════════════════════════════════════════════
panel_a <- function() {
  prism_path <- "paralog_sl_predictor/output/prism_top_hits.csv"
  if (file.exists(prism_path)) {
    pr <- read_csv(prism_path, show_col_types = FALSE)
  } else {
    stop("prism_top_hits.csv not found")
  }

  # Select top drugs and top pairs by effect size
  top_drugs <- pr %>%
    group_by(drug) %>%
    summarise(min_delta = min(delta_auc, na.rm = TRUE)) %>%
    arrange(min_delta) %>%
    head(20) %>%
    pull(drug)

  top_pairs <- pr %>%
    mutate(pair_label = paste0(driver, "->", paralog)) %>%
    group_by(pair_label) %>%
    summarise(min_delta = min(delta_auc, na.rm = TRUE)) %>%
    arrange(min_delta) %>%
    head(15) %>%
    pull(pair_label)

  # Filter to intersection
  pr_sub <- pr %>%
    mutate(pair_label = paste0(driver, "->", paralog)) %>%
    filter(drug %in% top_drugs, pair_label %in% top_pairs) %>%
    group_by(drug, pair_label) %>%
    summarise(delta_auc = min(delta_auc, na.rm = TRUE), .groups = "drop")

  # Shorten drug names (18 chars keeps the y labels inside the composite
  # canvas when panel a shares 180 mm with panel b)
  pr_sub$drug_short <- substr(pr_sub$drug, 1, 18)
  pr_sub$drug_short <- factor(pr_sub$drug_short, levels = rev(sort(unique(pr_sub$drug_short))))
  pr_sub$pair_label <- factor(pr_sub$pair_label, levels = sort(unique(pr_sub$pair_label)))

  # Cap delta_auc for visualization
  pr_sub$delta_capped <- pmax(pr_sub$delta_auc, -0.7)

  ggplot(pr_sub, aes(pair_label, drug_short, fill = delta_capped)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-0.7, 0.1),
                         name = "ΔAUC\n(MUT-WT)") +
    theme_sci
}

# ═══════════════════════════════════════════════════════════════
# PANEL B — Targeted-agent assay-validity anchors (|ΔAUC|)
# ═══════════════════════════════════════════════════════════════
panel_b <- function() {
  # Identical data path to the former main-text Fig. 4a (R_fig4.R prior to
  # the restructuring): top hits per drug class from prism_top_hits.csv plus
  # the caption-named biology-validation rows pulled exactly from
  # prism_full_results.csv. Never hardcoded literals.
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
  # Biology-validation drugs: the caption names trametinib / everolimus /
  # ipatasertib, whose significant PanCancer associations sit below the
  # top-30 display cutoff in prism_top_hits.csv. Pull those exact rows from
  # the full results so the panel shows what the caption claims (no literals:
  # values come straight from prism_full_results.csv).
  full_path <- "paralog_sl_predictor/output/prism_full_results.csv"
  if (!file.exists(full_path))
    stop("paralog_sl_predictor/output/prism_full_results.csv not found")
  bio <- read_csv(full_path, show_col_types = FALSE) %>%
    filter(context == "PanCancer",
           (drug == "TRAMETINIB"  & driver == "KRAS" & paralog == "HRAS") |
           (drug == "EVEROLIMUS"  & driver == "PTEN" & paralog == "TNS2") |
           (drug == "IPATASERTIB" & driver == "PTEN" & paralog == "TNS2")) %>%
    mutate(
      abs_delta = abs(delta_auc),
      drug_upper = toupper(drug),
      drug_class = case_when(
        grepl("MEK|AZD8330|TRAMETINIB|RO-4987655", drug_upper) ~ "MEKi",
        grepl("MTOR|EVEROLIMUS|TEMSIROLIMUS|AKT|IPATASERTIB|GSK-2141795", drug_upper) ~ "mTOR/AKTi",
        grepl("HDAC|PANOBINOSTAT", drug_upper) ~ "HDACi",
        TRUE ~ "Other"))
  if (nrow(bio) != 3)
    stop("expected 3 biology-validation rows in prism_full_results.csv, got ", nrow(bio))
  pr <- bind_rows(pr, bio) %>%
    distinct(drug, driver, paralog, .keep_all = TRUE) %>%
    arrange(desc(abs_delta))
  pr <- pr %>% mutate(
    label = paste0(drug, "\n", driver, "->", paralog),
    drug_class = factor(drug_class, levels = c("MEKi","mTOR/AKTi","HDACi","Other")))
  pr$label <- factor(pr$label, levels = rev(unique(pr$label)))

  ggplot(pr, aes(abs_delta, label, fill = drug_class)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c(MEKi = RED, `mTOR/AKTi` = BLUE, HDACi = ORANGE, Other = GREEN),
                      drop = TRUE) +  # hide classes absent from the plotted hits
                                 # (avoids an empty swatch in the legend)
    labs(x = "|ΔAUC|", y = NULL) +
    theme_bar + theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== FigS7 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b()

# patchwork composition (cowplot::plot_grid clipped the heatmap's long
# y-axis drug labels at the canvas edge)
p <- pa + pb +
  patchwork::plot_layout(ncol = 2, widths = c(0.62, 0.38)) +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(size = 9, face = "bold",
                                          family = "Arial")))

ggsave(file.path(OUT_DIR, "FigS7_PRISM_Selectivity.pdf"), p,
       width = 180, height = 120, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS7_PRISM_Selectivity.svg"), p,
       width = 180, height = 120, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS7_PRISM_Selectivity.tiff"), p,
       width = 180, height = 120, units = "mm", device = ragg::agg_tiff, dpi = 300)
ggsave(file.path(OUT_DIR, "FigS7_PRISM_Selectivity.png"), p,
       width = 180, height = 120, units = "mm", device = ragg::agg_png, dpi = 300)
REVIEW_DIR <- "figure_review"
dir.create(REVIEW_DIR, showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(OUT_DIR, "FigS7_PRISM_Selectivity.png"),
          file.path(REVIEW_DIR, "FigS7_PRISM_Selectivity.png"), overwrite = TRUE)
message("FigS7_PRISM_Selectivity.pdf (180×120mm) ✓")
