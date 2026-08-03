# Fig3 — Clinical Stratification (R)
# Purpose: panels a+b (90×90mm) side by side → 180×90mm composite
# Note:    former panel c (TCGA survival forest) moved to Supplementary
#          Fig. S10 (R_figS10_survival.R) after manuscript restructuring.
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

save_panel <- function(p, name, w = PANEL_W, h = PANEL_H) {
  ggsave(file.path(OUT_DIR, paste0("Fig3_panel_", name, ".pdf")), p,
         width = w, height = h, units = "mm", device = cairo_pdf)
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
  df$msi_status <- factor(df$msi_status, levels = c("MSS","MSI-H"))

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

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== Fig3 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b()

save_panel(pa, "a"); save_panel(pb, "b")

# 2-panel composite: a+b side by side (rel_widths keep the two panels
# visually balanced — review feedback 2026-07-28)
p <- cowplot::plot_grid(pa, pb, ncol = 2, rel_widths = c(0.53, 0.47),
                        labels = c("a","b"),
                        label_size = 9, label_fontface = "bold",
                        label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "Fig3_Clinical.pdf"), p,
       width = 180, height = 90, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig3_Clinical.svg"), p,
       width = 180, height = 90, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "Fig3_Clinical.tiff"), p,
       width = 180, height = 90, units = "mm", device = ragg::agg_tiff, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig3_Clinical.png"), p,
       width = 180, height = 90, units = "mm", device = ragg::agg_png, dpi = 300)
REVIEW_DIR <- "figure_review"
dir.create(REVIEW_DIR, showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(OUT_DIR, "Fig3_Clinical.png"),
          file.path(REVIEW_DIR, "Fig3_Clinical.png"), overwrite = TRUE)
message("Fig3_Clinical.pdf (180×90mm) ✓")
