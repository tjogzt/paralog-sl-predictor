# Fig3 — Clinical Stratification (R)
# Purpose: 4 individual 90×90mm panels → cowplot → 180×180mm composite
# Usage:   Rscript R_fig3.R
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)

# ── Constants ──
BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7
PANEL_W <- 90; PANEL_H <- 90
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
  legend.title = element_blank(),
  legend.background = element_blank(),
  legend.key        = element_blank(),
  plot.margin  = margin(4, 4, 4, 4, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

save_panel <- function(p, name) {
  ggsave(file.path(OUT_DIR, paste0("Fig3_panel_", name, ".pdf")), p,
         width = PANEL_W, height = PANEL_H, units = "mm", device = cairo_pdf)
  message(sprintf("  panel %s ✓", name))
}

# ═══════════════════════════════════════════════════════════════
# PANEL A — MSI Stratification
# ═══════════════════════════════════════════════════════════════
panel_a <- function() {
  # Single source of truth: output/msi_key_numbers_min3.json (sensitivity
  # frame, >=3 mutant/WT per group; official DepMap 26Q1 MSIsensor2
  # annotation, MSI-H = MSIscore > 20), written by msi_analysis.py.
  # On the primary >=5 frame the endometrial subgroups are not evaluable
  # and colorectal shows no difference; see manuscript text.
  json_path <- "paralog_sl_predictor/output/msi_key_numbers_min3.json"
  if (!file.exists(json_path))
    stop("msi_key_numbers_min3.json not found — run msi_analysis.py first; ",
         "simulated fallbacks are forbidden")
  kj <- jsonlite::fromJSON(json_path)
  sg <- kj$subgroups
  pick <- function(key) {
    # keys are e.g. "Colorectal_MSI_H" / "Colorectal_MSS"
    status <- if (grepl("_MSI_H$", key)) "MSI-H" else "MSS"
    cancer <- sub("_(MSI_H|MSS)$", "", key)
    tibble(cancer = cancer, msi_status = status,
           auroc = sg[[key]]$dd_auroc, n = sg[[key]]$n_lines)
  }
  df <- bind_rows(pick("Colorectal_MSI_H"), pick("Colorectal_MSS"),
                  pick("Endometrial_MSI_H"), pick("Endometrial_MSS"))
  df$cancer <- factor(df$cancer, levels = c("Colorectal","Endometrial"))
  df$msi_status <- factor(df$msi_status, levels = c("MSI-H","MSS"))

  # For hatched bar, use a dummy value with patterned fill
  df$bar_val <- ifelse(is.na(df$auroc), 0.25, df$auroc)

  ggplot(df, aes(cancer, bar_val, fill = msi_status)) +
    geom_col(position = position_dodge(0.7), width = 0.55) +
    geom_text(aes(y = bar_val + 0.03,
                  label = ifelse(!is.na(auroc), paste0("n=", n), paste0("n=", n, "*\n(insuff.)"))),
              position = position_dodge(0.7), size = 2.5, color = GRAY, vjust = 0) +
    scale_fill_manual(values = c("MSI-H" = "#F4A582", "MSS" = BLUE)) +
    geom_hline(yintercept = 0.5, linewidth = 0.3, color = GRAY, linetype = "dashed", alpha = 0.3) +
    labs(x = NULL, y = "DD AUROC") +
    theme_sci + theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))
}

# ═══════════════════════════════════════════════════════════════
# PANEL B — Mutation Type ΔDD
# ═══════════════════════════════════════════════════════════════
panel_b <- function() {
  # Read directly from output/muttype_{cancer}_results.csv (WT − MUT,
  # manuscript Eq. 1; positive = stronger dependency in the mutant
  # subgroup) — never hardcoded literals.
  specs <- tribble(
    ~file,          ~cancer,       ~driver,   ~paralog,  ~label,
    "ovarian",      "Ovarian",     "ARID1A",  "ARID1B",  "ARID1A->ARID1B\n(Ovarian)",
    "colorectal",   "Colorectal",  "EP300",   "CREBBP",  "EP300->CREBBP\n(Colorectal)",
    "ovarian",      "Ovarian",     "BRCA1",   "BRCA2",   "BRCA1->BRCA2\n(Ovarian)",
    "breast",       "Breast",      "BRCA1",   "BRCA2",   "BRCA1->BRCA2\n(Breast)")
  df <- purrr::pmap_dfr(specs, function(file, cancer, driver, paralog, label) {
    path <- sprintf("paralog_sl_predictor/output/muttype_%s_results.csv", file)
    if (!file.exists(path))
      stop(path, " not found — run the mutation-type analysis first; ",
           "simulated fallbacks are forbidden")
    row <- read_csv(path, show_col_types = FALSE) %>%
      filter(.data$cancer == .env$cancer, .data$driver == .env$driver,
             .data$paralog == .env$paralog)
    if (nrow(row) != 1) stop("expected exactly 1 row for ", label, " in ", path)
    tibble(pair = label,
           Truncating = row$dd_trunc, Missense = row$dd_miss)
  }) %>% pivot_longer(c(Truncating, Missense), names_to = "type", values_to = "dd")
  df$pair <- factor(df$pair, levels = unique(df$pair))
  df$type <- factor(df$type, levels = c("Truncating","Missense"))

  ggplot(df, aes(dd, pair, fill = type)) +
    geom_col(position = position_dodge(0.7), width = 0.55) +
    scale_fill_manual(values = c(Truncating = "#CB181D", Missense = "#9ECAE1")) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = DARK) +
    labs(x = "DD (WT - MUT)", y = NULL) +
    theme_sci + theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))
}

# (panel b — legend upper right)

# ═══════════════════════════════════════════════════════════════
# PANEL C — TCGA Survival Forest Plot
# ═══════════════════════════════════════════════════════════════
panel_c <- function() {
  # Single source of truth: output/tcga_survival_associations.csv,
  # written by tcga_survival.py (Cox PH, median-split high vs low paralog
  # expression, TCGA PanCan Atlas BRCA). Never hardcoded literals.
  tcga_path <- "paralog_sl_predictor/output/tcga_survival_associations.csv"
  if (!file.exists(tcga_path))
    stop("tcga_survival_associations.csv not found — run tcga_survival.py ",
         "first; simulated fallbacks are forbidden")
  genes8 <- c("ARID1B","BRCA2","PIK3CB","CRKL","CREBBP","ATR","SMARCA2","HRAS")
  df <- read_csv(tcga_path, show_col_types = FALSE) %>%
    filter(gene %in% genes8) %>%
    mutate(gene = factor(gene, levels = rev(genes8)),
           sig = ifelse(p_value < 0.05,
                        sprintf("p=%s", formatC(p_value, format = "f", digits = 3)), ""),
           clr = case_when(hr > 1 & sig != "" ~ RED,
                           hr > 1              ~ BLUE,
                           TRUE                ~ GREEN),
           ypos = as.numeric(gene) + 0.18)

  ggplot(df, aes(hr, gene)) +
    geom_vline(xintercept = 1, linewidth = 0.4, color = DARK, alpha = 0.5) +
    # Gene labels drawn INSIDE the plot panel, above each errorbar (left
    # side): the composite figure's left edge clipped the widest y-axis
    # label ("SMARCA2"), so the y-axis text column is removed entirely;
    # labels sit above the bars so whiskers never cross the text.
    geom_text(aes(x = 0.50, y = ypos, label = gene), hjust = 0, size = 2.5,
              family = "Arial", color = DARK) +
    geom_point(aes(color = clr), size = 2) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, color = clr),
                   height = 0.15, linewidth = 0.8) +
    geom_text(aes(label = sig, color = clr, y = ypos),
              size = 2.5, fontface = "bold", hjust = 0.5, vjust = 0) +
    scale_color_identity() +
    scale_x_continuous(limits = c(0.48, 2.4), breaks = c(0.5, 1.0, 1.5, 2.0)) +
    scale_y_discrete(expand = expansion(add = c(0.3, 0.7))) +
    labs(x = "Hazard Ratio\n(high vs low paralog expression)", y = NULL) +
    theme_sci +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank())
}

# ═══════════════════════════════════════════════════════════════
# PANEL D — Mutational Co-occurrence
# ═══════════════════════════════════════════════════════════════
panel_d <- function() {
  # Single source of truth: output/cooccurrence_analysis.csv, written by
  # cooccurrence_analysis.py (Fisher's exact test on DepMap 26Q1 driver-
  # rule mutation status across the 1,208 dependency-profiled cell lines).
  co_path <- "paralog_sl_predictor/output/cooccurrence_analysis.csv"
  if (!file.exists(co_path))
    stop("cooccurrence_analysis.csv not found — run cooccurrence_analysis.py ",
         "first; simulated fallbacks are forbidden")
  pair_colors <- c("ARID1A/ARID1B" = RED, "PIK3CA/PIK3CB" = BLUE,
                   "BRCA1/BRCA2" = ORANGE, "EP300/CREBBP" = "#0D7377",
                   "SMARCA4/SMARCA2" = PURPLE)
  df <- read_csv(co_path, show_col_types = FALSE) %>%
    arrange(desc(odds_ratio)) %>%
    mutate(pair = factor(pair, levels = pair),
           clr = pair_colors[as.character(pair)],
           star = ifelse(p_value < 0.05, "*", ""))

  ggplot(df, aes(odds_ratio, pair, fill = clr)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = sprintf("%.2f%s", odds_ratio, star)),
              hjust = -0.1, size = 2.5, color = DARK) +
    scale_fill_identity() +
    geom_vline(xintercept = 1, linewidth = 0.4, color = DARK, alpha = 0.5, linetype = "dashed") +
    labs(x = "Co-occurrence OR\n(>1 = co-occur)", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    annotate("text", x = Inf, y = Inf,
             label = "All OR > 1\nSL at dependency level\n*Fisher p < 0.05",
             size = 2.5, color = DARK, hjust = 1.05, vjust = 1.3) +
    theme_sci
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== Fig3 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b(); pc <- panel_c(); pd <- panel_d()

save_panel(pa, "a"); save_panel(pb, "b"); save_panel(pc, "c"); save_panel(pd, "d")

# plot_grid (not ggdraw/draw_plot) with a wider left column: draw_plot clipped
# panel c's widest y-label ("SMARCA2") at the figure's left edge
p <- cowplot::plot_grid(pa, pb, pc, pd, ncol = 2, rel_widths = c(0.53, 0.47),
                        labels = c("a","b","c","d"),
                        label_size = 9, label_fontface = "bold",
                        label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "Fig3_Clinical.pdf"), p,
       width = 180, height = 180, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig3_Clinical.svg"), p,
       width = 180, height = 180, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "Fig3_Clinical.tiff"), p,
       width = 180, height = 180, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("Fig3_Clinical.pdf (180×180mm) ✓")
