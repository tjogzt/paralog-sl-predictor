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
  df <- tibble(
    cancer = factor(c("Colorectal","Colorectal","Endometrial","Endometrial"),
                    levels = c("Colorectal","Endometrial")),
    msi_status = factor(c("MSI-H","MSS","MSI-H","MSS"),
                        levels = c("MSI-H","MSS")),
    auroc = c(0.631, 0.804, 0.592, NA),
    n     = c(17, 42, 20, 8))

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
  # Values recomputed from output/muttype_{cancer}_results.csv (WT − MUT,
  # manuscript Eq. 1; positive = stronger dependency in the mutant subgroup)
  df <- bind_rows(
    tibble(pair = "ARID1A->ARID1B\n(Ovarian)",    type = "Truncating", dd = 0.388),
    tibble(pair = "ARID1A->ARID1B\n(Ovarian)",    type = "Missense",   dd = 0.020),
    tibble(pair = "EP300->CREBBP\n(Colorectal)",   type = "Truncating", dd = 0.464),
    tibble(pair = "EP300->CREBBP\n(Colorectal)",   type = "Missense",   dd = 0.150),
    tibble(pair = "BRCA1->BRCA2\n(Ovarian)",       type = "Truncating", dd = 0.080),
    tibble(pair = "BRCA1->BRCA2\n(Ovarian)",       type = "Missense",   dd = 0.000),
    tibble(pair = "BRCA1->BRCA2\n(Breast)",        type = "Truncating", dd = -0.136),
    tibble(pair = "BRCA1->BRCA2\n(Breast)",        type = "Missense",   dd = 0.000))
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
  df <- tibble(
    gene = c("BRCA2","ATR","ARID1B","CRKL","SMARCA2","CREBBP","PIK3CB","HRAS"),
    hr   = c(1.116,1.112,1.084,1.084,1.045,1.040,0.981,0.971),
    se   = c(0.054,0.054,0.054,0.054,0.054,0.054,0.054,0.054),
    sig  = c("p=0.032","p=0.039","","","","","",""))
  df$gene <- factor(df$gene, levels = rev(df$gene))
  df$clr  <- ifelse(df$hr > 1 & df$sig != "", RED,
                    ifelse(df$hr > 1, BLUE, GREEN))
  df$ypos <- as.numeric(df$gene) + 0.15

  ggplot(df, aes(hr, gene)) +
    geom_vline(xintercept = 1, linewidth = 0.4, color = DARK, alpha = 0.5) +
    geom_point(aes(color = clr), size = 2) +
    geom_errorbarh(aes(xmin = hr - 1.96*se, xmax = hr + 1.96*se, color = clr),
                   height = 0.15, linewidth = 0.8) +
    geom_text(aes(label = sig, color = clr, y = ypos),
              size = 2.5, fontface = "bold", hjust = 0.5, vjust = 0) +
    scale_color_identity() +
    labs(x = "Hazard Ratio (High vs Low paralog expression)", y = NULL) +
    annotate("text", x = -Inf, y = Inf, label = "TCGA BRCA\n(n=1,082)",
             size = 2.5, color = GRAY, hjust = -0.05, vjust = 1.2, lineheight = 0.9) +
    theme_sci
}

# ═══════════════════════════════════════════════════════════════
# PANEL D — Mutational Co-occurrence
# ═══════════════════════════════════════════════════════════════
panel_d <- function() {
  df <- tibble(
    pair = c("ARID1A/ARID1B","PIK3CA/PIK3CB","BRCA1/BRCA2","EP300/CREBBP","SMARCA4/SMARCA2"),
    or   = c(6.147, 5.157, 3.422, 4.548, 4.753),
    clr  = c(RED, BLUE, ORANGE, "#0D7377", PURPLE))
  df$pair <- factor(df$pair, levels = rev(df$pair))

  ggplot(df, aes(or, pair, fill = clr)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = sprintf("%.2f", or)), hjust = -0.1, size = 2.5, color = DARK) +
    scale_fill_identity() +
    geom_vline(xintercept = 1, linewidth = 0.4, color = DARK, alpha = 0.5, linetype = "dashed") +
    labs(x = "Co-occurrence OR\n(>1 = co-occur)", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    annotate("text", x = Inf, y = -Inf, label = "All OR > 1\nSL at dependency level",
             size = 2.5, color = DARK, hjust = 1.05, vjust = -0.5) +
    theme_sci
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== Fig3 Panel Generation (R) ===")
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

ggsave(file.path(OUT_DIR, "Fig3_Clinical.pdf"), p,
       width = 180, height = 180, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig3_Clinical.svg"), p,
       width = 180, height = 180, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "Fig3_Clinical.tiff"), p,
       width = 180, height = 180, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("Fig3_Clinical.pdf (180×180mm) ✓")
