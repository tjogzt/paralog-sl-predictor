# FigS1 — Evaluation landscape: cell lines, tested pairs, and DD AUROC by cancer
# Merged 2026-07-29: former FigS1 (panels a/b) + former FigS2 (per-lineage DD
# AUROC detail, now panel c). One shared AUROC-tier colour key (legend in a).
# Source: solid_tumor_summary.csv (real data)
library(ggplot2)
library(cowplot)
library(dplyr)
library(readr)

OUT_DIR <- "paralog_sl_predictor/output/figures"
BASE_FS <- 7; TICK_FS <- 7
RED   <- "#CB181D"; BLUE  <- "#2171B5"
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

# Panel a: Cell lines per cancer type
pa <- ggplot(d, aes(n_lines, cancer, fill = auroc_tier)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = tier_colors, name = NULL) +
  geom_text(aes(label = n_lines), hjust = -0.15, size = 2.5, color = DARK) +
  labs(x = "Cell lines in DepMap 26Q1", y = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  theme_sci +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"))

# Panel b: Tested pairs and known positives (AUROC-tier colored + red known overlay)
pb <- ggplot(d, aes(n_pairs, cancer, fill = auroc_tier)) +
  geom_col(width = 0.7, alpha = 0.65) +
  geom_col(aes(x = n_known), width = 0.7, fill = RED, alpha = 0.9) +
  scale_fill_manual(values = tier_colors, name = NULL, guide = "none") +
  geom_text(aes(x = n_pairs, label = sprintf("%d/%d", n_known, n_pairs)),
            hjust = -0.05, size = 2.5, color = DARK) +
  annotate("text", x = Inf, y = -Inf,
           label = "Background: AUROC tier (see a)\nRed bars: known positives",
           hjust = 1.05, vjust = -0.3, size = 2.5, color = DARK, fontface = "italic") +
  labs(x = "Tested pairs", y = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  theme_sci +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Panel c: DD AUROC for the evaluable lineages (absorbed from former FigS2,
# single shared 0-1 axis; tier colours as in panel a, so legend is off here)
dc <- d %>% filter(!is.na(dd_auroc)) %>% arrange(dd_auroc)
dc$cancer <- factor(as.character(dc$cancer), levels = as.character(dc$cancer))

pc <- ggplot(dc, aes(dd_auroc, cancer, fill = auroc_tier)) +
  geom_col(width = 0.62) +
  scale_fill_manual(values = tier_colors, name = NULL, guide = "none") +
  geom_text(aes(label = sprintf("%.3f", dd_auroc)),
            hjust = -0.15, size = 2.5, color = DARK) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             linewidth = 0.3, color = GRAY, alpha = 0.6) +
  # No 0.7 threshold line: it crossed the bar-end value labels of the
  # 0.65-0.68 lineages; the tier colours already encode it (see panel a).
  scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "DD AUROC", y = NULL) +
  theme_sci

# Save individual panels
ggsave(file.path(OUT_DIR, "FigS1_panel_a.pdf"), pa, width = 90, height = 120, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS1_panel_b.pdf"), pb, width = 90, height = 120, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS1_panel_c.pdf"), pc, width = 90, height = 70, units = "mm", device = cairo_pdf)
message("  Individual panels saved ✓")

# Composite: 3 panels, one row, 180x80mm (2026-07-30: 120 -> 101.25 -> 80mm;
# reviewer found even 16:9 too tall for the 23-row bar charts)
p <- ggdraw() +
  draw_plot(pa, x = 0,    y = 0, width = 0.40, height = 1) +
  draw_plot(pb, x = 0.40, y = 0, width = 0.28, height = 1) +
  draw_plot(pc, x = 0.68, y = 0, width = 0.32, height = 1) +
  draw_plot_label(c("a","b","c"), x = c(0, 0.38, 0.66), y = c(1, 1, 1),
                  size = 9, fontface = "bold", fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "FigS1_CellLine_Landscape.pdf"), p,
       width = 180, height = 80, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS1_CellLine_Landscape.svg"), p,
       width = 180, height = 80, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS1_CellLine_Landscape.tiff"), p,
       width = 180, height = 80, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS1_CellLine_Landscape.pdf (180×80mm, 3 panels: lines/pairs/AUROC, real data) ✓")
