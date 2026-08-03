# FigS2 — Evaluation robustness: bootstrap CI, label-shuffle negative control,
#         and external check against mutation-agnostic gold standards (in4mer)
# Merged 2026-07-29: former FigS8 (panels a/b) + former FigS11 (panels c/d).
# All values computed live from reproducible artifacts:
#   output/validation_report.json          (observed AUROC, null mean, bootstrap CI)
#   output/permutation_10000.rds           (REAL 10,000 label-shuffle null draws)
#   output/bootstrap_perpair_1000.csv      (REAL 1,000 bootstrap resample draws)
#   output/in4mer_benchmark.csv            (pair-level DD, in4mer_benchmark.py)
#   output/in4mer_benchmark_summary.json   (AUROC/CI/p, in4mer_benchmark.py)
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

vr <- fromJSON("paralog_sl_predictor/output/validation_report.json")
perm_data <- readRDS("paralog_sl_predictor/output/permutation_10000.rds")
bs_real <- read.csv("paralog_sl_predictor/output/bootstrap_perpair_1000.csv")$bootstrap_auroc

obs_auroc <- vr$negative_control$observed_auroc
emp_p     <- vr$negative_control$empirical_p_value
null_vals <- perm_data$null
null_mean <- vr$negative_control$null_auroc_mean
bs_ci_lo  <- vr$bootstrap$auroc_ci_low
bs_ci_hi  <- vr$bootstrap$auroc_ci_high

# ═══════════════════════════════════════════════════════════════
# PANEL A — Bootstrap (REAL 1,000 resample draws from run_full_validation)
# ═══════════════════════════════════════════════════════════════
pa <- ggplot(data.frame(x = bs_real), aes(x)) +
  geom_histogram(bins = 30, fill = BLUE, alpha = 0.5, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = obs_auroc, color = RED, linewidth = 1.2) +
  geom_vline(xintercept = c(bs_ci_lo, bs_ci_hi), color = GRAY, linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = obs_auroc - 0.02, y = Inf,
           label = sprintf("%.3f", obs_auroc),
           size = 3, color = RED, fontface = "bold", hjust = 1, vjust = 1.5) +
  annotate("text", x = bs_ci_lo, y = Inf,
           label = sprintf("95%% CI\n[%.3f, %.3f]", bs_ci_lo, bs_ci_hi),
           size = 2.5, color = GRAY, hjust = 0, vjust = 3) +
  labs(x = "DD AUROC", y = "Frequency (1,000 iterations)") +
  theme_sci

# ═══════════════════════════════════════════════════════════════
# PANEL B — Negative control (REAL 10,000 permutations)
# ═══════════════════════════════════════════════════════════════
pb <- ggplot(data.frame(x = null_vals), aes(x)) +
  geom_histogram(bins = 40, fill = GRAY, alpha = 0.4, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = obs_auroc, color = RED, linewidth = 1.2) +
  geom_vline(xintercept = null_mean, color = DARK, linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = obs_auroc - 0.02, y = Inf,
           label = sprintf("Observed %.3f", obs_auroc),
           size = 2.5, color = RED, hjust = 1, vjust = 1.5) +
  annotate("text", x = obs_auroc - 0.02, y = Inf,
           label = sprintf("Null mean %.3f", null_mean),
           size = 2.5, color = DARK, hjust = 1, vjust = 3.3) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("empirical p = %.3f\n(10,000 permutations)", emp_p),
           size = 2.8, color = RED, fontface = "bold", hjust = 1.05, vjust = 1.5,
           lineheight = 0.9) +
  labs(x = "AUROC (shuffled labels)", y = "Frequency (10,000 permutations)") +
  theme_sci

# ═══════════════════════════════════════════════════════════════
# PANEL C — in4mer external check: |DD| of gold pairs vs unlabeled controls
# ═══════════════════════════════════════════════════════════════
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

pc <- ggplot(ctrl, aes(absdd)) +
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

# ═══════════════════════════════════════════════════════════════
# PANEL D — in4mer AUROC point estimates + bootstrap CI, both frameworks
# ═══════════════════════════════════════════════════════════════
fb <- data.frame(
  framework = factor(c(sprintf("≥5 mutants (n = %d pairs)", summ$min5$n_pos),
                       sprintf("≥3 mutants (n = %d pairs)", summ$min3$n_pos)),
                     levels = c(sprintf("≥3 mutants (n = %d pairs)", summ$min3$n_pos),
                                sprintf("≥5 mutants (n = %d pairs)", summ$min5$n_pos))),
  auroc = c(summ$min5$auroc, summ$min3$auroc),
  lo    = c(summ$min5$bootstrap_ci_low,  summ$min3$bootstrap_ci_low),
  hi    = c(summ$min5$bootstrap_ci_high, summ$min3$bootstrap_ci_high),
  p     = c(summ$min5$permutation_p,     summ$min3$permutation_p))

pd <- ggplot(fb, aes(y = framework)) +
  geom_vline(xintercept = 0.5, color = DARK, linewidth = 0.5, linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.15, color = BLUE, linewidth = 0.5) +
  geom_point(aes(x = auroc), color = BLUE, size = 1.8) +
  geom_text(aes(x = 1.03, label = sprintf("AUROC = %.2f\np = %.2f", auroc, p)),
            size = 2.3, color = DARK, hjust = 0, vjust = 0.5, lineheight = 0.9,
            family = "Arial") +
  annotate("text", x = 0.5, y = 2.4, label = "chance (0.50)", size = 2.3,
           color = DARK, hjust = -0.05, family = "Arial") +
  scale_x_continuous(limits = c(0, 1.65), breaks = seq(0, 1, 0.25)) +
  labs(x = "AUROC (bootstrap 95% CI)", y = NULL) +
  theme_sci + theme(axis.text.y = element_text(size = TICK_FS))

# ═══════════════════════════════════════════════════════════════
# COMPOSITE — 2x2 grid, 180x170mm
# ═══════════════════════════════════════════════════════════════
p <- cowplot::plot_grid(pa, pb, pc, pd, ncol = 2,
                        labels = c("a", "b", "c", "d"),
                        label_size = 9, label_fontface = "bold",
                        label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "FigS2_Evaluation_Robustness.pdf"), p,
       width = 180, height = 170, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS2_Evaluation_Robustness.svg"), p,
       width = 180, height = 170, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS2_Evaluation_Robustness.tiff"), p,
       width = 180, height = 170, units = "mm", device = ragg::agg_tiff, dpi = 600)
message(sprintf("FigS2_Evaluation_Robustness.pdf done: bootstrap+negctrl+in4mer; %d evaluable gold, %d controls, p90 = %.3f",
                nrow(gold), nrow(ctrl), p90))
