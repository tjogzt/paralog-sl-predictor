# FigS11 — External check against mutation-agnostic gold standards (in4mer)
# All values computed live from reproducible artifacts:
#   output/in4mer_benchmark.csv          (pair-level DD, in4mer_benchmark.py)
#   output/in4mer_benchmark_summary.json (AUROC/CI/p, in4mer_benchmark.py)
# No simulated, random, or hardcoded data.
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

bench <- read.csv("paralog_sl_predictor/output/in4mer_benchmark.csv")
summ  <- fromJSON("paralog_sl_predictor/output/in4mer_benchmark_summary.json")

# ---- min3 framework: 10 evaluable gold pairs vs evaluable unlabeled controls
gold <- bench[bench$label == "in4mer_gold"      & !is.na(bench$dd_min3), ]
ctrl <- bench[bench$label == "unlabeled_control" & !is.na(bench$dd_min3), ]
gold$absdd <- abs(gold$dd_min3)
ctrl$absdd <- abs(ctrl$dd_min3)
p90 <- as.numeric(quantile(ctrl$absdd, 0.9))          # R type-7 == np.percentile
hi  <- gold[gold$absdd > p90, ]                        # pairs above control p90

# histogram geometry for tick placement
hh   <- hist(ctrl$absdd, breaks = 40, plot = FALSE)
ymax <- max(hh$counts)
tick_h <- ymax * 0.07
gold$lab_y <- ymax * 0.16

# stagger labels of highlighted pairs to avoid overlap
hi <- hi[order(hi$absdd, decreasing = TRUE), ]
hi$lab_y <- ymax * c(0.42, 0.66, 0.90)[seq_len(nrow(hi))]

pa <- ggplot(ctrl, aes(absdd)) +
  geom_histogram(bins = 40, fill = GRAY, alpha = 0.4, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = p90, color = DARK, linewidth = 0.5, linetype = "dashed") +
  geom_segment(data = gold,
               aes(x = absdd, xend = absdd, y = 0, yend = tick_h),
               color = RED, linewidth = 0.5) +
  geom_segment(data = hi,
               aes(x = absdd, xend = absdd, y = tick_h, yend = lab_y - ymax * 0.05),
               color = RED, linewidth = 0.3, linetype = "dotted") +
  geom_text(data = hi,
            aes(x = absdd + 0.012, y = lab_y, label = pair),
            size = 2.2, color = RED, fontface = "bold", family = "Arial", hjust = 0) +
  annotate("text", x = p90 + 0.008, y = ymax * 0.72,
           label = sprintf("control p90 = %.3f", p90),
           size = 2.5, color = DARK, hjust = 0) +
  annotate("text", x = Inf, y = Inf, label = sprintf("in4mer gold (n = %d)", nrow(gold)),
           size = 2.5, color = RED, hjust = 1.1, vjust = 1.6) +
  annotate("text", x = Inf, y = Inf, label = sprintf("unlabeled controls (n = %d)", nrow(ctrl)),
           size = 2.5, color = GRAY, hjust = 1.1, vjust = 3.4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = expression("|"*DD*"| (" >= "3-mutant framework)"), y = "Count") +
  theme_sci

# ---- Panel b: AUROC point estimates + bootstrap CI for both frameworks
fb <- data.frame(
  framework = factor(c(sprintf("\u22655 mutants (n = %d pairs)", summ$min5$n_pos),
                       sprintf("\u22653 mutants (n = %d pairs)", summ$min3$n_pos)),
                     levels = c(sprintf("\u22653 mutants (n = %d pairs)", summ$min3$n_pos),
                                sprintf("\u22655 mutants (n = %d pairs)", summ$min5$n_pos))),
  auroc = c(summ$min5$auroc, summ$min3$auroc),
  lo    = c(summ$min5$bootstrap_ci_low,  summ$min3$bootstrap_ci_low),
  hi    = c(summ$min5$bootstrap_ci_high, summ$min3$bootstrap_ci_high),
  p     = c(summ$min5$permutation_p,     summ$min3$permutation_p))

pb <- ggplot(fb, aes(y = framework)) +
  geom_vline(xintercept = 0.5, color = DARK, linewidth = 0.5, linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.15, color = BLUE, linewidth = 0.5) +
  geom_point(aes(x = auroc), color = BLUE, size = 1.8) +
  geom_text(aes(x = 1.04, label = sprintf("AUROC = %.2f\np = %.2f", auroc, p)),
            size = 2.3, color = DARK, hjust = 0, vjust = 0.5, lineheight = 0.9,
            family = "Arial") +
  annotate("text", x = 0.5, y = 2.4, label = "chance (0.50)", size = 2.3,
           color = DARK, hjust = -0.05, family = "Arial") +
  scale_x_continuous(limits = c(0, 1.45), breaks = seq(0, 1, 0.25)) +
  labs(x = "AUROC (bootstrap 95% CI)", y = NULL) +
  theme_sci + theme(axis.text.y = element_text(size = TICK_FS))

p <- ggdraw() +
  draw_plot(pa, x = 0,    y = 0, width = 0.52, height = 1) +
  draw_plot(pb, x = 0.52, y = 0, width = 0.48, height = 1) +
  draw_plot_label(c("a", "b"), x = c(0, 0.52), y = c(1, 1),
                  size = 9, fontface = "bold", fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "FigS11_in4mer_ExternalCheck.pdf"), p,
       width = 180, height = 85, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS11_in4mer_ExternalCheck.svg"), p,
       width = 180, height = 85, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS11_in4mer_ExternalCheck.tiff"), p,
       width = 180, height = 85, units = "mm", device = ragg::agg_tiff, dpi = 600)
message(sprintf("FigS11 done: %d evaluable gold, %d evaluable controls, p90 = %.3f, %d pairs above p90",
                nrow(gold), nrow(ctrl), p90, nrow(hi)))
