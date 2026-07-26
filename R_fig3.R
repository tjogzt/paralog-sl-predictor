# Fig3 — Clinical Stratification (R)
# Purpose: panels a+b (90×90mm) on top row, panel c (180×90mm) full width below
#          → 180×180mm composite
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

# ═══════════════════════════════════════════════════════════════
# PANEL C — TCGA Survival Forest Plot (multivariable v2, full width)
# ═══════════════════════════════════════════════════════════════
panel_c <- function() {
  # Single source of truth: output/tcga_survival_v2.json, written by
  # tcga_survival_v2.py (Cox PH on continuous z-scored log2 expression;
  # multivariable model adjusting for age + AJCC stage; BH FDR across the
  # 32-gene family). Never hardcoded literals.
  v2_path <- "paralog_sl_predictor/output/tcga_survival_v2.json"
  if (!file.exists(v2_path))
    stop("tcga_survival_v2.json not found — run tcga_survival_v2.py first; ",
         "simulated fallbacks are forbidden")
  v2 <- jsonlite::fromJSON(v2_path)
  # The four FDR-significant genes (BH q<0.05 in the multivariable family)
  # plus the compensating paralogs of the lead candidate pairs and ARID1A
  # for direct contrast (see manuscript text).
  genes8 <- c("PIK3CA","ARID1B","RBL1","BRCA2","PIK3CB","ARID1A","SMARCA2","CREBBP")
  pg <- v2$per_gene
  df <- tibble(
    gene = pg$gene,
    hr   = pg$multivar_age_stage$hr_multivar,
    lo   = vapply(pg$multivar_age_stage$ci_multivar, `[`, numeric(1), 1),
    hi   = vapply(pg$multivar_age_stage$ci_multivar, `[`, numeric(1), 2),
    p    = pg$multivar_age_stage$p_multivar,
    q    = pg$multivar_age_stage$q_fdr_multivar
  ) %>% filter(.data$gene %in% genes8) %>%
    arrange(desc(.data$hr)) %>%
    mutate(gene = factor(gene, levels = rev(gene)),   # highest HR on top
           fdr = .data$q < 0.05,
           clr = ifelse(fdr, RED, GRAY),
           lab = paste0(gene, ifelse(fdr, "*", "")),
           txt = sprintf("%.2f (%.2f\u2013%.2f)", hr, lo, hi),
           ypos = as.numeric(gene))

  ci_max <- max(df$hi); ci_min <- min(df$lo)
  x_left  <- max(0.4, floor(ci_min * 10) / 10 - 0.1)
  x_ann   <- ci_max * 1.12          # left edge of the HR (CI) text column
  x_right <- x_ann * 1.5            # axis right limit

  ggplot(df, aes(hr, gene)) +
    geom_vline(xintercept = 1, linewidth = 0.4, color = DARK, alpha = 0.5) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, color = clr),
                   height = 0.18, linewidth = 0.8) +
    geom_point(aes(color = clr), size = 2) +
    geom_text(aes(y = gene, label = txt), x = x_ann, hjust = 0,
              size = 2.5, family = "Arial", color = DARK) +
    annotate("text", x = x_ann, y = nrow(df) + 0.9, label = "HR (95% CI)",
             hjust = 0, size = 2.5, fontface = "bold", family = "Arial",
             color = DARK) +
    annotate("text", x = x_left, y = 0.2, hjust = 0, vjust = 0,
             label = "*FDR q < 0.05 (BH, 32-gene family)",
             size = 2.5, family = "Arial", color = GRAY) +
    scale_color_identity() +
    scale_y_discrete(labels = setNames(df$lab, df$gene),
                     expand = expansion(add = c(0.8, 1.6))) +
    scale_x_continuous(limits = c(x_left, x_right),
                       breaks = pretty(c(x_left, ci_max), n = 4)) +
    labs(x = "Multivariable hazard ratio per SD\n(Cox PH, adjusted for age + AJCC stage)",
         y = NULL) +
    theme_sci +
    theme(axis.text.y = element_text(size = TICK_FS, face = "italic"))
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== Fig3 Panel Generation (R) ===")
pa <- panel_a(); pb <- panel_b(); pc <- panel_c()

save_panel(pa, "a"); save_panel(pb, "b"); save_panel(pc, "c", w = PANEL_W * 2)

# 3-panel composite: a+b top row, c full width below
top_row <- cowplot::plot_grid(pa, pb, ncol = 2, rel_widths = c(0.53, 0.47),
                              labels = c("a","b"),
                              label_size = 9, label_fontface = "bold",
                              label_fontfamily = "Arial")
p <- cowplot::plot_grid(top_row, pc, ncol = 1,
                        labels = c("", "c"),
                        label_size = 9, label_fontface = "bold",
                        label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "Fig3_Clinical.pdf"), p,
       width = 180, height = 180, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig3_Clinical.svg"), p,
       width = 180, height = 180, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "Fig3_Clinical.tiff"), p,
       width = 180, height = 180, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("Fig3_Clinical.pdf (180×180mm) ✓")
