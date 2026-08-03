# Fig1 — Framework, Benchmark & Validation (R)
# Purpose: 4 individual 90×90mm panels → cowplot → 180×180mm composite
# Note:    panel c is the benchmark-selection flowchart (replaced the former
#          cross-study CV3 bar comparison after manuscript restructuring).
#          All flowchart numbers are read from output/headline_metrics.json
#          and output/tables/TableS3_GoldStandard.tsv — never literals.
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
DARKRED <- "#7F0000"

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
    stop("paralog_sl_predictor/output/solid_tumor_summary.csv not found — ",
         "run the pipeline first; simulated fallbacks are forbidden")
  }
  df$cancer <- factor(df$cancer, levels = df$cancer)
  df$clr <- ifelse(df$dd_auroc >= 0.7, RED, ifelse(df$dd_auroc >= 0.5, BLUE, GRAY))

  ggplot(df, aes(dd_auroc, cancer)) +
    geom_col(aes(fill = clr), width = 0.6) +
    scale_fill_identity() +
    geom_vline(xintercept = 0.5, linetype = "dashed",
               linewidth = 0.3, color = GRAY, alpha = 0.35) +
    geom_vline(xintercept = 0.7, linetype = "dotdash",
               linewidth = 0.5, color = DARKRED) +
    labs(x = "DD AUROC", y = NULL) +
    theme_sci +
    annotate("text", x = 0.69, y = 0.8,
             label = "AUROC 0.7", size = 2.5, color = DARKRED, hjust = 1,
             fontface = "bold")
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
             label = sprintf("perm. p = %.3f; exact MW p = %.3f\n(n=%d vs n=%d)",
                             perm_p, pval,
                             sum(df$group == "TSG-driven"), n_onc),
             size = 2.5, hjust = 0.5) +
    geom_hline(yintercept = 0.5, linewidth = 0.3, color = GRAY, alpha = 0.3, linetype = "dashed") +
    theme_sci + theme(legend.position = "none")
}

# ═══════════════════════════════════════════════════════════════
# PANEL C — Benchmark selection flowchart
# ═══════════════════════════════════════════════════════════════
panel_c <- function() {
  # Every number in this flowchart is read from the single-source-of-truth
  # artifacts (simulated/hardcoded values are forbidden):
  #   pair & tier counts  <- output/tables/TableS3_GoldStandard.tsv
  #   entries/positives   <- output/headline_metrics.json (lineage_full)
  #   evaluation lineages <- headline_metrics.json (leave_one_lineage_out)
  #   min-sample rule     <- headline_metrics.json (min_samples string)
  #   framework count     <- headline_metrics.json (framework result blocks)
  gs_path <- "paralog_sl_predictor/output/tables/TableS3_GoldStandard.tsv"
  hm_path <- "paralog_sl_predictor/output/headline_metrics.json"
  if (!file.exists(gs_path))
    stop("TableS3_GoldStandard.tsv not found — run tables.py first")
  if (!file.exists(hm_path))
    stop("headline_metrics.json not found — run compute_headline_metrics.py first")
  gs <- read_tsv(gs_path, show_col_types = FALSE)
  hm <- jsonlite::fromJSON(hm_path)

  n_pairs   <- nrow(gs)                                # 12 curated pairs
  n_tier_a  <- sum(gs$Tier == "A")                     # dual-gene perturbation
  n_tier_b  <- sum(gs$Tier == "B")                     # natural-genotype dep.
  n_tier_c  <- sum(gs$Tier == "C")                     # indirect evidence
  n_comp    <- sum(gs$Tier == "Comparator")            # specificity references
  n_primary <- n_tier_a + n_tier_b                     # primary benchmark pairs
  n_entries <- hm$lineage_full$n_entries               # 110
  n_pos     <- hm$lineage_full$n_positives             # 8
  n_ctrl    <- n_entries - n_pos                       # 102 unlabeled controls
  frame_lineages <- sub("^without_", "",
                        names(hm$leave_one_lineage_out$values))
  # manuscript display order (validated against the artifact keys above)
  lineage_order <- c("Ovarian", "Endometrial", "Cervical")
  frame_lineages <- c(intersect(lineage_order, frame_lineages),
                      setdiff(frame_lineages, lineage_order))
  ms_rule <- gsub(">=", "\u2265",
                  regmatches(hm$min_samples,
                             regexpr(">=[0-9]+ mutant and >=[0-9]+ WT",
                                     hm$min_samples)))
  ms_rule <- gsub(" and ", " + ", ms_rule)
  n_fw <- sum(c("lineage_full", "per_pair_max_from_tables2",
                "per_pair_mean_from_tables2") %in% names(hm))

  # ── Drawing primitives (canvas 0–100 × 0–100 on the square panel) ──
  st_main <- list(fill = "#E8F1FA", border = BLUE, linetype = "solid", text = DARK)
  st_side <- list(fill = "#F7F7F7", border = GRAY, linetype = "dashed", text = GRAY)

  draw_box <- function(xc, yc, w, h, label, st, fs = 2.5) {
    list(
      annotation_custom(
        grid::roundrectGrob(r = grid::unit(0.14, "snpc"),
                            gp = grid::gpar(fill = st$fill, col = st$border,
                                            lwd = 0.6, lty = st$linetype)),
        xmin = xc - w / 2, xmax = xc + w / 2, ymin = yc - h / 2, ymax = yc + h / 2),
      annotate("text", x = xc, y = yc, label = label, size = fs,
               family = "Arial", color = st$text, lineheight = 0.95))
  }
  arw <- grid::arrow(length = grid::unit(1.4, "mm"), type = "closed")
  seg <- function(x1, y1, x2, y2, st = st_main)
    annotate("segment", x = x1, y = y1, xend = x2, yend = y2,
             arrow = arw, color = st$border, linewidth = 0.5,
             linetype = st$linetype)

  layers <- c(
    # main path (solid blue)
    draw_box(37, 93, 62, 9,
             sprintf("%d curated paralog pairs\n(citation-verified)", n_pairs),
             st_main),
    draw_box(18.5, 75.5, 33, 10,
             sprintf("Tier A (%d): dual-gene\nperturbation", n_tier_a), st_main),
    draw_box(53, 75.5, 35, 12,
             sprintf("Tier B (%d):\nnatural-genotype\ndependency", n_tier_b), st_main),
    draw_box(37, 58.5, 54, 8,
             sprintf("%d pairs: primary literature-\nderived benchmark", n_primary), st_main),
    draw_box(37, 43, 70, 10,
             sprintf("Evaluation frame: %s\n%s lines",
                     paste(frame_lineages, collapse = " / "), ms_rule), st_main),
    draw_box(37, 28, 54, 9,
             sprintf("%d entries: %d positive +\n%d unlabeled controls",
                     n_entries, n_pos, n_ctrl), st_main),
    draw_box(37, 12.5, 74, 10,
             sprintf("%d aggregation frameworks: lineage-level (primary)\n\u00b7 per-pair max \u00b7 per-pair mean",
                     n_fw), st_main),
    # side branches (dashed gray — not part of the primary evaluation)
    draw_box(84.5, 78, 29, 10,
             sprintf("Tier C (%d):\nindirect evidence", n_tier_c), st_side),
    draw_box(84.5, 64, 29, 10,
             sprintf("Comparators (%d):\nspecificity references", n_comp), st_side),
    list(
      annotate("text", x = 84.5, y = 55.5,
               label = "excluded from primary\nevaluation",
               size = 2.5, family = "Arial", color = GRAY, fontface = "italic",
               lineheight = 0.95),
      # arrows: fan-out from the curated set
      seg(28, 88.4, 21, 80.7),
      seg(46, 88.4, 50, 81.9),
      seg(68, 93.8, 84.5, 83.2, st_side),
      seg(84.5, 72.8, 84.5, 69.4, st_side),
      # arrows: Tier A + Tier B converge on the primary benchmark
      seg(22, 70.3, 30, 62.7),
      seg(48, 69.3, 44, 62.7),
      # arrows: down the evaluation chain
      seg(37, 54.3, 37, 48.3),
      seg(37, 37.7, 37, 32.8),
      seg(37, 23.3, 37, 17.8)))

  ggplot() + layers +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
    theme_void(base_family = "Arial") +
    theme(plot.margin = margin(2, 2, 2, 2, "pt"),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
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

  # Inset: REAL paired-bootstrap distribution of PCS − DD (10,000 resamples,
  # dumped by compute_headline_metrics.py). Answers the question this panel
  # raises — is the PCS bar really higher than DD? (mean +0.150, 95% CI
  # −0.110 to +0.456 → no significant difference.) The pair-level negative
  # control formerly shown here lives in Supplementary Fig. S8.
  hm <- jsonlite::fromJSON("paralog_sl_predictor/output/headline_metrics.json")
  pb <- hm$component_paired_bootstrap$pcs_minus_dd
  deltas <- readr::read_csv("paralog_sl_predictor/output/component_paired_bootstrap_deltas.csv",
                            show_col_types = FALSE)
  dd_ <- deltas[deltas$component == "pcs_minus_dd", ]
  d_mean  <- pb$mean_delta
  d_ci    <- unlist(pb$ci95)

  inset <- ggplot(dd_, aes(delta)) +
    geom_histogram(bins = 30, fill = BLUE, alpha = 0.4, color = NA) +
    geom_vline(xintercept = 0, color = GRAY, linewidth = 0.4, linetype = "dotted") +
    geom_vline(xintercept = d_mean, color = RED, linewidth = 0.7) +
    geom_vline(xintercept = d_ci, color = GRAY, linewidth = 0.4, linetype = "dashed") +
    annotate("text", x = d_mean, y = Inf,
             label = sprintf("%+.3f", d_mean),
             size = 3.2, color = RED, fontface = "bold", hjust = -0.1, vjust = 1.5) +
    labs(x = "AUROC (PCS \u2212 DD)", y = "N") +
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
       width = 180, height = 180, units = "mm", device = ragg::agg_tiff, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig1_Framework_Validation.png"), p,
       width = 180, height = 180, units = "mm", device = ragg::agg_png, dpi = 300)
REVIEW_DIR <- "figure_review"
dir.create(REVIEW_DIR, showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(OUT_DIR, "Fig1_Framework_Validation.png"),
          file.path(REVIEW_DIR, "Fig1_Framework_Validation.png"), overwrite = TRUE)
message("Fig1_Framework_Validation.pdf (180×180mm) ✓")
