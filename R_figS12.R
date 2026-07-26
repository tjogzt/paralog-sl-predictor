# FigS12 — Mutational co-occurrence by lineage, with and without TMB adjustment (R)
# Purpose: panel a (pan-cancer unadjusted vs TMB-adjusted OR, 90×95mm) +
#          panel b (per-lineage TMB-adjusted OR heatmap, 90×95mm)
#          → 180×95mm composite
# Usage:   Rscript R_figS12.R
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)

# ── Constants ──
BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7
OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Colors ──
BLUE  <- "#2171B5"; RED   <- "#CB181D"; GREEN <- "#238B45"
ORANGE <- "#E6550D"; PURPLE <- "#6A51A3"
GRAY  <- "#636363"; DARK  <- "#252525"

# ── Theme ──
theme_sci <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.line    = element_line(linewidth = 0.4),
  axis.ticks   = element_line(linewidth = 0.3),
  axis.text    = element_text(size = TICK_FS),
  axis.title   = element_text(size = BASE_FS),
  legend.text  = element_text(size = LEGEND_FS),
  legend.background = element_blank(),
  legend.key        = element_blank(),
  plot.margin  = margin(4, 4, 4, 4, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

# ── Data (single source of truth) ──
# output/cooccurrence_by_lineage.json, written by cooccurrence_by_lineage.py:
# Fisher's exact test per pair x lineage on DepMap 26Q1 driver-rule mutation
# status, plus a logistic model logit(paralog_mut ~ driver_mut + log1p(TMB))
# whose driver coefficient gives the TMB-adjusted OR. Never hardcoded literals.
json_path <- "paralog_sl_predictor/output/cooccurrence_by_lineage.json"
if (!file.exists(json_path))
  stop("cooccurrence_by_lineage.json not found — run cooccurrence_by_lineage.py ",
       "first; simulated fallbacks are forbidden")
cl <- jsonlite::fromJSON(json_path)
res <- cl$results %>%
  mutate(ci_lo = vapply(ci, `[`, numeric(1), 1),
         ci_hi = vapply(ci, `[`, numeric(1), 2),
         ci_adj_lo = vapply(ci_adj, function(x) if (length(x) < 2 || is.null(x[[1]])) NA_real_ else x[[1]], numeric(1)),
         ci_adj_hi = vapply(ci_adj, function(x) if (length(x) < 2 || is.null(x[[2]])) NA_real_ else x[[2]], numeric(1)))

pairs <- c("ARID1A/ARID1B","EP300/CREBBP","BRCA1/BRCA2","PIK3CA/PIK3CB","SMARCA4/SMARCA2")
pair_colors <- setNames(c(RED, "#0D7377", ORANGE, BLUE, PURPLE), pairs)

# ═══════════════════════════════════════════════════════════════
# PANEL A — PAN-CANCER: unadjusted vs TMB-adjusted OR
# ═══════════════════════════════════════════════════════════════
panel_a <- function() {
  pan <- res %>% filter(lineage == "PAN-CANCER") %>%
    mutate(pair = factor(pair, levels = rev(pairs)))
  dfl <- bind_rows(
    pan %>% transmute(pair, est = "Unadjusted (Fisher)",
                      or = or, lo = ci_lo, hi = ci_hi, star = ""),
    pan %>% transmute(pair, est = "TMB-adjusted (logistic)",
                      or = tmb_adjusted_or, lo = ci_adj_lo, hi = ci_adj_hi,
                      star = ifelse(p_adj < 0.05, "*", ""))
  ) %>% mutate(est = factor(est, levels = c("Unadjusted (Fisher)",
                                            "TMB-adjusted (logistic)")),
               yn = as.numeric(pair) + ifelse(est == "Unadjusted (Fisher)", 0.17, -0.17))

  ggplot(dfl, aes(y = yn)) +
    geom_vline(xintercept = 1, linewidth = 0.4, color = DARK, alpha = 0.5,
               linetype = "dashed") +
    geom_errorbarh(aes(xmin = lo, xmax = hi, color = pair,
                       alpha = est, linewidth = est), height = 0.12) +
    geom_point(aes(x = or, color = pair, shape = est, alpha = est),
               size = 2) +
    geom_text(data = filter(dfl, star != ""),
              aes(x = hi * 1.18, label = star),
              size = 3.5, color = DARK, family = "Arial") +
    scale_color_manual(values = pair_colors, guide = "none") +
    scale_shape_manual(values = c("Unadjusted (Fisher)" = 1,
                                  "TMB-adjusted (logistic)" = 16), name = NULL) +
    scale_alpha_manual(values = c(0.55, 1), guide = "none") +
    scale_linewidth_manual(values = c(0.6, 0.8), guide = "none") +
    scale_x_continuous(trans = "log10", breaks = c(0.5, 1, 2, 5, 10)) +
    coord_cartesian(xlim = c(0.45, 14)) +
    scale_y_continuous(breaks = seq_along(pairs), labels = rev(pairs),
                       expand = expansion(add = c(0.5, 0.4))) +
    labs(x = "Co-occurrence odds ratio (log scale)", y = NULL) +
    theme_sci +
    theme(legend.position = c(0.03, 0.97), legend.justification = c(0, 1),
          legend.key.size = unit(3, "mm")) +
    guides(shape = guide_legend(override.aes = list(color = DARK, alpha = 1)))
}

# ═══════════════════════════════════════════════════════════════
# PANEL B — Per-lineage TMB-adjusted OR heatmap
# ═══════════════════════════════════════════════════════════════
panel_b <- function() {
  lin_order <- c("Bladder/Urinary Tract","Bowel","Breast","CNS/Brain",
                 "Esophagus/Stomach","Lung","Lymphoid","Ovary/Fallopian Tube",
                 "Skin","Uterus")
  lin_short <- c("Bladder","Bowel","Breast","CNS","Esoph/Stom","Lung",
                 "Lymphoid","Ovary","Skin","Uterus")

  grid <- expand.grid(pair = pairs, lineage = lin_order, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    left_join(res %>% select(pair, lineage, tmb_adjusted_or, p_adj),
              by = c("pair","lineage")) %>%
    mutate(
      skipped = paste(pair, lineage) %in% paste(cl$skipped_lineages$pair, cl$skipped_lineages$lineage),
      in_results = paste(pair, lineage) %in% paste(res$pair, res$lineage),
      cell = case_when(
        !in_results ~ "\u2013",                       # not evaluable
        is.na(tmb_adjusted_or) ~ "n.e.",              # complete separation
        TRUE ~ sprintf("%.3g%s", tmb_adjusted_or,
                       ifelse(!is.na(p_adj) & p_adj < 0.05, "*", ""))),
      pair = factor(pair, levels = rev(pairs)),
      lineage = factor(lineage, levels = lin_order))

  ggplot(grid, aes(lineage, pair, fill = tmb_adjusted_or)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = cell), size = 2.2, family = "Arial", color = DARK) +
    scale_fill_gradient2(low = BLUE, mid = "white", high = RED, midpoint = 0,
                         trans = "log10", na.value = "grey88", guide = "none") +
    scale_x_discrete(labels = setNames(lin_short, lin_order),
                     expand = expansion(add = 0)) +
    scale_y_discrete(expand = expansion(add = 0)) +
    labs(x = NULL, y = NULL) +
    theme_sci +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1, size = TICK_FS),
          axis.text.y = element_text(size = TICK_FS))
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== FigS12 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b()

ggsave(file.path(OUT_DIR, "FigS12_panel_a.pdf"), pa,
       width = 90, height = 95, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS12_panel_b.pdf"), pb,
       width = 90, height = 95, units = "mm", device = cairo_pdf)

p <- cowplot::plot_grid(pa, pb, ncol = 2, rel_widths = c(0.46, 0.54),
                        labels = c("a","b"),
                        label_size = 9, label_fontface = "bold",
                        label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "FigS12_Cooccurrence_TMB.pdf"), p,
       width = 180, height = 95, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS12_Cooccurrence_TMB.svg"), p,
       width = 180, height = 95, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS12_Cooccurrence_TMB.tiff"), p,
       width = 180, height = 95, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS12_Cooccurrence_TMB.pdf (180×95mm) ✓")
