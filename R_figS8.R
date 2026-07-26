# FigS8 — Bootstrap + Negative Control (REAL 10,000 permutations)
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
bs_mean   <- vr$bootstrap$auroc_mean
bs_ci_lo  <- vr$bootstrap$auroc_ci_low
bs_ci_hi  <- vr$bootstrap$auroc_ci_high

# Panel A: Bootstrap (REAL 1,000 resample draws from run_full_validation —
# dumped to output/bootstrap_perpair_1000.csv; no simulated shape)
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

# Panel B: Negative Control (REAL 10,000 permutations)
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

p <- ggdraw() +
  draw_plot(pa, x = 0,   y = 0, width = 0.5, height = 1) +
  draw_plot(pb, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot_label(c("a","b"), x = c(0, 0.5), y = c(1, 1), size = 9, fontface = "bold", fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "FigS8_Bootstrap_NegCtrl.pdf"), p,
       width = 180, height = 90, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS8_Bootstrap_NegCtrl.svg"), p,
       width = 180, height = 90, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS8_Bootstrap_NegCtrl.tiff"), p,
       width = 180, height = 90, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS8 updated with 10,000 real permutations ✓")
