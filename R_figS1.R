# FigS1 — Panoramic Baseline: Cell Lines and Evaluation Summary (23 cancer types)
# Source: solid_tumor_summary.csv (real data)
library(ggplot2)
library(cowplot)
library(dplyr)
library(readr)

OUT_DIR <- "paralog_sl_predictor/output/figures"
BASE_FS <- 7; TICK_FS <- 6
RED   <- "#CB181D"; BLUE  <- "#2171B5"; GREEN <- "#238B45"
GRAY  <- "#636363"; DARK <- "#252525"

theme_sci <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.line  = element_line(linewidth = 0.4),
  axis.ticks = element_line(linewidth = 0.3),
  axis.text  = element_text(size = TICK_FS),
  axis.title = element_text(size = BASE_FS),
  plot.margin = margin(4, 6, 4, 4, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

# Load real data
d <- read_csv("paralog_sl_predictor/output/solid_tumor_summary.csv", show_col_types = FALSE)
d <- d %>% mutate(
  auroc_tier = case_when(
    is.na(dd_auroc) ~ "Not evaluable",
    dd_auroc > 0.7  ~ "AUROC > 0.7",
    dd_auroc > 0.5  ~ "AUROC 0.5-0.7",
    TRUE            ~ "AUROC < 0.5"),
  auroc_tier = factor(auroc_tier, levels = c("AUROC > 0.7","AUROC 0.5-0.7",
                                              "AUROC < 0.5","Not evaluable")),
  cancer = reorder(cancer, n_lines))

tier_colors <- c("AUROC > 0.7" = RED, "AUROC 0.5-0.7" = BLUE,
                  "AUROC < 0.5" = GRAY, "Not evaluable" = "#D9D9D9")

# Panel A: Cell lines per cancer type
pa <- ggplot(d, aes(n_lines, cancer, fill = auroc_tier)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = tier_colors, name = NULL) +
  geom_text(aes(label = n_lines), hjust = -0.15, size = 2.2, color = DARK) +
  labs(x = "Cell lines in DepMap 26Q1", y = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_sci +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.text = element_text(size = 5.5),
        legend.key.size = unit(3, "mm"))

# Panel B: Tested pairs and known positives (AUROC-tier colored + red known overlay)
pb <- ggplot(d, aes(n_pairs, cancer, fill = auroc_tier)) +
  geom_col(width = 0.7, alpha = 0.65) +
  geom_col(aes(x = n_known), width = 0.7, fill = RED, alpha = 0.9) +
  scale_fill_manual(values = tier_colors, name = NULL, guide = "none") +
  geom_text(aes(x = n_pairs, label = sprintf("%d/%d", n_known, n_pairs)),
            hjust = -0.05, size = 2, color = DARK) +
  annotate("text", x = Inf, y = -Inf,
           label = "Background: AUROC tier (see a)\nRed bars: known positives",
           hjust = 1.05, vjust = -0.3, size = 2.2, color = DARK, fontface = "italic") +
  labs(x = "Tested pairs", y = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_sci +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Save individual panels
ggsave(file.path(OUT_DIR, "FigS1_panel_a.pdf"), pa, width = 90, height = 120, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS1_panel_b.pdf"), pb, width = 90, height = 120, units = "mm", device = cairo_pdf)
message("  Individual panels saved ✓")

# Composite
p <- ggdraw() +
  draw_plot(pa, x = 0,    y = 0, width = 0.55, height = 1) +
  draw_plot(pb, x = 0.55, y = 0, width = 0.45, height = 1) +
  draw_plot_label(c("a","b"), x = c(0, 0.52), y = c(1, 1), size = 9, fontface = "bold")

ggsave(file.path(OUT_DIR, "FigS1_CellLine_Landscape.pdf"), p,
       width = 180, height = 120, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS1_CellLine_Landscape.svg"), p,
       width = 180, height = 120, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS1_CellLine_Landscape.tiff"), p,
       width = 180, height = 120, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS1_CellLine_Landscape.pdf (180×120mm, 23 cancer types, real data) ✓")
