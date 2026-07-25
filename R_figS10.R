# FigS10 — Structure-based druggability assessment
# Regenerated from output/druggability_analysis.json (alphafold_analysis.py).
# This script restores full provenance for the previously source-less legacy PDF.
library(ggplot2)
library(cowplot)
library(jsonlite)

OUT_DIR <- "paralog_sl_predictor/output/figures"
BASE_FS <- 7; TICK_FS <- 7
BLUE  <- "#2171B5"; RED   <- "#CB181D"; GRAY  <- "#636363"; DARK <- "#252525"

theme_sci <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.line  = element_line(linewidth = 0.4),
  axis.ticks = element_line(linewidth = 0.3),
  axis.text  = element_text(size = TICK_FS),
  axis.title = element_text(size = BASE_FS),
  plot.margin = margin(4, 4, 4, 4, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

d <- fromJSON("paralog_sl_predictor/output/druggability_analysis.json")
stopifnot(nrow(d) == 15)
d$is_paralog <- grepl("^Paralog", d$role)
# Pair-grouped order (paralog first within each pair), top-to-bottom as in the
# legacy figure; reverse for ggplot's bottom-up y axis.
d$gene <- factor(d$gene, levels = rev(d$gene))
d$kp <- sprintf("K=%d P=%d", d$lys_count, d$pocket_regions)

p <- ggplot(d, aes(druggability_score, gene, fill = is_paralog)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = kp), hjust = 1.08, size = 2.5, color = "white") +
  geom_text(aes(label = sprintf("%.3f", druggability_score)),
            hjust = -0.08, size = 2.5, color = DARK) +
  scale_fill_manual(name = NULL,
                    values = c(`TRUE` = RED, `FALSE` = BLUE),
                    labels = c(`TRUE` = "Paralog (target)", `FALSE` = "Driver")) +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Druggability score", y = NULL) +
  theme_sci +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        legend.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUT_DIR, "FigS10_Druggability.pdf"), p,
       width = 120, height = 120, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS10_Druggability.svg"), p,
       width = 120, height = 120, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS10_Druggability.tiff"), p,
       width = 120, height = 120, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS10_Druggability.pdf (120×120mm) ✓  15 proteins")
