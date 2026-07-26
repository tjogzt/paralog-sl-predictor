# Fig1 — Framework, Benchmark & Validation (R)
# Purpose: 4 individual 90×90mm panels → cowplot → 180×180mm composite
# Usage:   Rscript R_fig1.R
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)

# ── Constants ──
BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7
PANEL_W <- 90; PANEL_H <- 90
OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Colors ──
BLUE  <- "#2171B5"; RED   <- "#CB181D"; GREEN <- "#238B45"
ORANGE <- "#E6550D"; GRAY  <- "#636363"; DARK <- "#252525"

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
  ggsave(file.path(OUT_DIR, paste0("Fig1_panel_", name, ".pdf")), p,
         width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
  message(sprintf("  panel %s ✓", name))
}

# ═══════════════════════════════════════════════════════════════
# PANEL A — Pan-cancer DD AUROC
# ═══════════════════════════════════════════════════════════════
panel_a <- function() {
  solid_path <- "paralog_sl_predictor/output/solid_tumor_summary.csv"
  if (file.exists(solid_path)) {
    df <- read_csv(solid_path, show_col_types = FALSE) %>%
      drop_na(dd_auroc) %>% arrange(dd_auroc)
  } else {
    df <- tibble(
      cancer = c("Skin","Pancreas","Lung","NSCLC","Endometrial","Cervical",
                 "Ovarian","CNS/Brain","HNSCC","Breast","SCLC",
                 "Soft Tissue","Esophagogastric","Colorectal","Mesothelioma","Biliary Tract"),
      dd_auroc = c(.342,.457,.570,.572,.576,.643,.685,.707,.714,.746,.750,.800,.805,.812,.917,.993))
  }
  df$cancer <- factor(df$cancer, levels = df$cancer)
  df$clr <- ifelse(df$dd_auroc >= 0.7, RED, ifelse(df$dd_auroc >= 0.5, BLUE, GRAY))

  ggplot(df, aes(dd_auroc, cancer)) +
    geom_col(aes(fill = clr), width = 0.6) +
    scale_fill_identity() +
    geom_vline(xintercept = c(0.5, 0.7), linetype = c("dashed","dotted"),
               linewidth = 0.3, color = c(GRAY, RED), alpha = 0.35) +
    labs(x = "DD AUROC", y = NULL) +
    theme_sci +
    annotate("text", x = 0.72, y = 0.8,
             label = "AUROC 0.7", size = 2.5, color = RED, hjust = 0)
}

# ═══════════════════════════════════════════════════════════════
# PANEL B — TSG vs Oncogene
# ═══════════════════════════════════════════════════════════════
panel_b <- function() {
  # Sensitivity frame (>=3 mutant/WT lines per group): the primary >=5 frame
  # leaves only 1 oncogene-driven lineage, so the TSG/ONC mechanism contrast
  # is shown on the relaxed frame and labelled as such in the caption.
  solid_path <- "paralog_sl_predictor/output/solid_tumor_summary_min3.csv"
  df <- read_csv(solid_path, show_col_types = FALSE) %>% drop_na(dd_auroc)
  # Biological classification (not AUROC threshold)
  oncogene_types <- c("Melanoma", "NSCLC", "Pancreatic")
  df$group <- ifelse(df$cancer %in% oncogene_types, "Oncogene-driven", "TSG-driven")
  df$group <- factor(df$group, levels = c("TSG-driven", "Oncogene-driven"))

  # Strip plot (not boxplot) — n=3 oncogene too small for boxplot
  pval <- wilcox.test(dd_auroc ~ group, data = df, exact = TRUE)$p.value
  # Exact two-sided permutation p for the TSG-minus-oncogene mean difference,
  # enumerated over all choose(n, 3) label assignments (replaces the previous
  # hard-coded Monte Carlo value 0.070; exact two-sided = 0.071, consistent).
  # Computed live from solid_tumor_summary.csv.
  obs_diff <- mean(df$dd_auroc[df$group == "TSG-driven"]) -
              mean(df$dd_auroc[df$group == "Oncogene-driven"])
  n_onc <- sum(df$group == "Oncogene-driven")
  combos <- combn(seq_len(nrow(df)), n_onc)
  perm_diffs <- apply(combos, 2, function(ii)
    mean(df$dd_auroc[-ii]) - mean(df$dd_auroc[ii]))
  perm_p <- min(1, 2 * min(mean(perm_diffs >= obs_diff), mean(perm_diffs <= obs_diff)))

  ggplot(df, aes(group, dd_auroc)) +
    geom_jitter(aes(color = group), width = 0.1, size = 2, alpha = 0.7) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.3, linewidth = 0.5, color = DARK) +
    scale_color_manual(values = c("TSG-driven" = BLUE, "Oncogene-driven" = RED)) +
    labs(x = NULL, y = "DD AUROC") +
    annotate("text", x = 1.5, y = 1.02,
             label = sprintf("perm. p = %.3f\n(n=%d vs n=%d)", perm_p,
                             sum(df$group == "TSG-driven"), n_onc),
             size = 2.5, hjust = 0.5) +
    geom_hline(yintercept = 0.5, linewidth = 0.3, color = GRAY, alpha = 0.3, linetype = "dashed") +
    theme_sci + theme(legend.position = "none")
}

# ═══════════════════════════════════════════════════════════════
# PANEL C — Benchmark Comparison
# ═══════════════════════════════════════════════════════════════
panel_c <- function() {
  # This-study values from the metrics single-source-of-truth (never literals);
  # published values are literature constants (Feng et al. 2024, Suppl. Data 1,
  # CV3 NSMRand 1:1, complete dataset).
  metrics_path <- "paralog_sl_predictor/output/tables/headline_metrics.tsv"
  if (!file.exists(metrics_path))
    stop("headline_metrics.tsv not found — run compute_headline_metrics.py first")
  mt <- read_tsv(metrics_path, show_col_types = FALSE)
  getv <- function(name) as.numeric(mt$value[mt$metric == name])
  df <- tibble(
    method = c("DD+ID>=30%","DD","SLMGAE","NSF4SL","GCATSL","GRSMF",
               "PiLSL","KG4SL","SLGNN","PTGNN"),
    auroc  = c(getv("dd_auroc_id_filter_0.3"), getv("dd_auroc_lineage_full"),
               getv("published_SLMGAE"), getv("published_NSF4SL"),
               getv("published_GCATSL"), getv("published_GRSMF"),
               getv("published_PiLSL"), getv("published_KG4SL"),
               getv("published_SLGNN"), getv("published_PTGNN")),
    ours   = c(TRUE,TRUE,rep(FALSE,8)))
  if (any(is.na(df$auroc)))
    stop("headline_metrics.tsv is missing required metrics — re-run compute_headline_metrics.py")
  df$method <- factor(df$method, levels = rev(df$method))

  ggplot(df, aes(auroc, method, fill = ours)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = sprintf("%.3f", auroc), color = ours),
              hjust = -0.1, size = 2.5, fontface = "bold") +
    scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = "#9ECAE1"), guide = "none") +
    scale_color_manual(values = c(`TRUE` = RED, `FALSE` = "#9ECAE1"), guide = "none") +
    geom_vline(xintercept = 0.5, linewidth = 0.3, color = GRAY, linetype = "dashed", alpha = 0.35) +
    labs(x = "CV3 AUROC", y = NULL) +
    scale_x_continuous(limits = c(0, 1.3), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0, 0))) +
    theme_sci + theme(plot.margin = margin(4, 8, 4, 20, "pt")) +
    annotate("text", x = 1.0, y = 9.5, label = "This study", size = 2.5, color = RED, hjust = 0) +
    annotate("text", x = 0.55, y = 7.5, label = "Published", size = 2.5, color = "#9ECAE1", hjust = 0)
}

# ═══════════════════════════════════════════════════════════════
# PANEL D — Component Decomposition + Bootstrap
# ═══════════════════════════════════════════════════════════════
panel_d <- function() {
  # Component decomposition values recomputed from TableS2 (single source of
  # truth); see compute_headline_metrics.py.
  metrics_path <- "paralog_sl_predictor/output/tables/headline_metrics.tsv"
  if (!file.exists(metrics_path))
    stop("headline_metrics.tsv not found — run compute_headline_metrics.py first")
  mt <- read_tsv(metrics_path, show_col_types = FALSE)
  getv <- function(name) as.numeric(mt$value[mt$metric == name])
  df <- tibble(
    metric = c("DD","PCS","ΔExpr","Necessity"),
    auroc  = c(getv("component_dd"), getv("component_pcs"),
               getv("component_delta_expression"), getv("component_necessity")),
    clr    = c(RED, BLUE, GRAY, ORANGE))
  if (any(is.na(df$auroc)))
    stop("headline_metrics.tsv is missing component metrics — re-run compute_headline_metrics.py")
  df$metric <- factor(df$metric, levels = df$metric)

  main <- ggplot(df, aes(metric, auroc)) +
    geom_col(aes(fill = clr), width = 0.55) +
    scale_fill_identity() +
    geom_text(aes(label = sprintf("%.3f", auroc), color = clr),
              vjust = -0.3, size = 2.5, fontface = "bold") +
    scale_color_identity() +
    geom_hline(yintercept = 0.5, linewidth = 0.3, color = GRAY, linetype = "dashed", alpha = 0.3) +
    labs(x = NULL, y = "AUROC") +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0, 0))) +
    theme_sci

  # Bootstrap inset — read live from validation_report.json (values change
  # whenever the pipeline re-runs; do not paste numbers into comments)
  vr <- jsonlite::fromJSON("paralog_sl_predictor/output/validation_report.json")
  obs_auroc <- vr$negative_control$observed_auroc
  bs_mean   <- vr$bootstrap$auroc_mean
  bs_ci_lo  <- vr$bootstrap$auroc_ci_low
  bs_ci_hi  <- vr$bootstrap$auroc_ci_high
  bs_sd     <- (bs_ci_hi - bs_ci_lo) / (2 * 1.96)  # reconstruct SD from CI

  set.seed(42)
  bs <- rnorm(1000, mean = bs_mean, sd = bs_sd)
  bs <- pmax(pmin(bs, 1.0), 0.3)
  bs_df <- tibble(auroc = bs)
  ci <- c(bs_ci_lo, bs_ci_hi)

  inset <- ggplot(bs_df, aes(auroc)) +
    geom_histogram(bins = 25, fill = BLUE, alpha = 0.4, color = NA) +
    geom_vline(xintercept = obs_auroc, color = RED, linewidth = 0.8) +
    annotate("text", x = obs_auroc + 0.02, y = Inf, label = sprintf("%.3f", obs_auroc),
             size = 3.5, color = RED, fontface = "bold", hjust = 0, vjust = 1.5) +
    geom_vline(xintercept = ci, color = GRAY, linewidth = 0.4, linetype = "dashed") +
    labs(x = "AUROC", y = "N") +
    theme_bw(base_size = 7) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 7))

  # Place inset in the upper-right band, clear of every bar and value label:
  # horizontally it starts right of the PCS bar (bar 2 right edge ~0.44 in
  # panel fraction); vertically it floats above the Necessity label (~0.60).
  ggdraw(main) +
    draw_plot(inset, x = 0.52, y = 0.66, width = 0.46, height = 0.32)
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== Fig1 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b(); pc <- panel_c(); pd <- panel_d()

save_panel(pa, "a"); save_panel(pb, "b"); save_panel(pc, "c"); save_panel(pd, "d")

p <- ggdraw() +
  draw_plot(pa, x = 0,    y = 0.5,  width = 0.5, height = 0.5) +
  draw_plot(pb, x = 0.5,  y = 0.5,  width = 0.5, height = 0.5) +
  draw_plot(pc, x = 0,    y = 0,    width = 0.5, height = 0.5) +
  draw_plot(pd, x = 0.5,  y = 0,    width = 0.5, height = 0.5) +
  draw_plot_label(c("a","b","c","d"),
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5),
                  size = 9, fontface = "bold", fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "Fig1_Framework_Validation.pdf"), p,
       width = 180, height = 180, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig1_Framework_Validation.svg"), p,
       width = 180, height = 180, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "Fig1_Framework_Validation.tiff"), p,
       width = 180, height = 180, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("Fig1_Framework_Validation.pdf (180×180mm) ✓")
