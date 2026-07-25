# Fig4 — Drug + Translational Prioritization (R)
# Purpose: Generate individual 90×90mm panels → review → patchwork → 180×180mm composite
# Usage:   Rscript R_fig4.R

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)

# ── Constants ──
BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7; ANNOT_FS <- 5.5
PANEL_W   <- 90; PANEL_H <- 90  # mm
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
    pr <- tibble(
      drug = c("AZD8330","Everolimus","Panobinostat","Triptolide","Ipatasertib","AT-9283"),
      driver = c("KRAS","PTEN","EP300","STK11","PTEN","KRAS"),
      paralog = c("HRAS","TNS2","CREBBP","SIK1","TNS2","HRAS"),
      abs_delta = c(0.352,0.315,0.339,1.221,3.24,0.331),
      drug_class = c("MEKi","mTOR/AKTi","HDACi","Other","mTOR/AKTi","Other"))
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
    tw <- tibble(
      driver = c("ARID1A","NF1","KMT2D","ATR","PPP2R1A","EP300","PIK3CA","SMARCA4",
                 "TP53","FBXW7","STK11","BRAF","KRAS","BRCA1","BRCA2","PIK3R1"),
      paralog = c("ARID1B","RASA2","KMT2C","ATM","PPP2R1B","CREBBP","PIK3CB","SMARCA2",
                  "TP63","FBXW2","SIK1","RAF1","HRAS","BRCA2","BRCA1","CRKL"),
      mean_ti = c(4.13,3.31,2.16,2.01,1.46,1.15,1.26,1.16,0.95,1.23,0.91,0.95,0.14,0.29,0.16,0.24),
      mean_selectivity = c(.237,.003,.097,.027,0,.115,-.051,.032,.022,-.0002,-.01,-.049,-.046,.045,-.039,-.148),
      classification = c("HIGH","MODERATE","MODERATE","MODERATE","LOW","MODERATE","LOW","MODERATE",
                         "MODERATE","LOW","LOW","LOW","LOW","PAN","PAN","PAN"))
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
              min.segment.length = 0, box.padding = 0.3) +
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
# PANEL C — Structural Similarity
# ═══════════════════════════════════════════════════════════════
panel_c <- function() {
  struct_path <- "paralog_sl_predictor/output/alphafold_structural_analysis.csv"
  if (file.exists(struct_path)) {
    st <- read_csv(struct_path, show_col_types = FALSE) %>%
      filter(domain_similarity > 0) %>%  # exclude zero-domain pairs (e.g. BRCA1/BRCA2)
      arrange(desc(structural_similarity)) %>% head(8)
  } else {
    st <- tibble(
      gene_a = c("EP300","SMARCA4","ARID1A","PPP2R1A","PIK3CA","KMT2D","TP53","FBXW7"),
      gene_b = c("CREBBP","SMARCA2","ARID1B","PPP2R1B","PIK3CB","KMT2C","TP63","FBXW2"),
      structural_similarity = c(.976,.966,.952,.921,.847,.829,.803,.666),
      domain_similarity     = c(1,1,.8,1,1,.6,.5,.333))
  }
  st$pair <- factor(paste0(st$gene_a, "/", st$gene_b),
                    levels = rev(paste0(st$gene_a, "/", st$gene_b)))
  df <- st %>% select(pair, structural_similarity, domain_similarity) %>%
    pivot_longer(-pair) %>%
    mutate(name = recode(name, structural_similarity = "Structural", domain_similarity = "Domain"))

  ggplot(df, aes(value, pair, fill = name)) +
    geom_col(position = position_dodge(0.7), width = 0.55) +
    scale_fill_manual(values = c(Structural = TEAL, Domain = ORANGE)) +
    labs(x = "Score", y = NULL) +
    # Legend moved out of the panel to the right: inside the panel it
    # overlapped the bottom bars (TP53/TP63 Domain, STK11/SIK1 Structural)
    coord_cartesian(clip = "off") +
    theme_sci + theme(legend.position = c(1.02, 0), legend.justification = c(0, 0),
                      legend.margin = margin(0, 0, 0, 0),
                      plot.margin = margin(4, 58, 8, 20, "pt"))
}

# (panel_c end — bottom bar: FBXW7/FBXW2 has both Structural & Domain, short values)

# ═══════════════════════════════════════════════════════════════
# PANEL D — Targetability Ranking
# ═══════════════════════════════════════════════════════════════
panel_d <- function() {
  struct_path <- "paralog_sl_predictor/output/alphafold_structural_analysis.csv"
  if (file.exists(struct_path)) {
    st <- read_csv(struct_path, show_col_types = FALSE)
    if ("clinical_targetability" %in% names(st)) {
      cand <- st %>% arrange(desc(clinical_targetability)) %>% head(10) %>%
        mutate(label = paste0(driver, "->", paralog), score = clinical_targetability)
    } else { cand <- NULL }
  } else { cand <- NULL }
  if (is.null(cand) || nrow(cand) == 0) {
    cand <- tibble(
      label = c("ARID1A>ARID1B","NF1>RASA2","KMT2D>KMT2C","PPP2R1A>PPP2R1B",
                "EP300>CREBBP","PIK3CA>PIK3CB","FBXW7>FBXW2","TP53>TP63",
                "STK11>SIK1","KRAS>HRAS"),
      score = c(.817,.615,.561,.525,.513,.495,.424,.423,.394,.390))
  }
  cand$label <- factor(cand$label, levels = rev(cand$label))
  cand$is_top <- c(TRUE, rep(FALSE, nrow(cand) - 1))
  cand$txt_col <- ifelse(cand$is_top, RED, BLUE)

  ggplot(cand, aes(score, label, fill = is_top)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = sprintf("%.3f", score), color = txt_col),
              hjust = -0.1, size = 2.5, fontface = "bold", show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = BLUE),
                      labels = c(`TRUE` = "Top candidate", `FALSE` = "Others")) +
    geom_vline(xintercept = 0.5, linewidth = 0.3, color = GRAY, linetype = "dashed", alpha = 0.3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Targetability Score", y = NULL) +
    theme_sci +
    theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))
}

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
message("=== Fig4 Panel Generation (R) ===")
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

ggsave(file.path(OUT_DIR, "Fig4_Translational.pdf"), p,
       width = 180, height = 180, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig4_Translational.svg"), p,
       width = 180, height = 180, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "Fig4_Translational.tiff"), p,
       width = 180, height = 180, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("Fig4_Translational.pdf (180×180mm) ✓")
