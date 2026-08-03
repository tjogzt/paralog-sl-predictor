# FigS10 — TCGA BRCA Survival Forest Plot (R)
# Purpose: single-panel forest plot (moved from main-text Fig. 3c after
#          manuscript restructuring) → 160×90mm, no panel letter
# Usage:   Rscript R_figS10_survival.R
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# ── Constants ──
BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7
FIG_W <- 160; FIG_H <- 90   # mm (single panel, 2/3-page width)
OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Colors ──
BLUE  <- "#2171B5"; RED   <- "#CB181D"
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

# ── Forest plot (identical data path to the former Fig. 3c) ──
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
x_right <- x_ann * 1.28           # axis right limit (tightened trailing blank)

p <- ggplot(df, aes(hr, gene)) +
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

ggsave(file.path(OUT_DIR, "FigS10_Survival.pdf"), p,
       width = FIG_W, height = FIG_H, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS10_Survival.svg"), p,
       width = FIG_W, height = FIG_H, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS10_Survival.tiff"), p,
       width = FIG_W, height = FIG_H, units = "mm", device = ragg::agg_tiff, dpi = 300)
ggsave(file.path(OUT_DIR, "FigS10_Survival.png"), p,
       width = FIG_W, height = FIG_H, units = "mm", device = ragg::agg_png, dpi = 300)
REVIEW_DIR <- "figure_review"
dir.create(REVIEW_DIR, showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(OUT_DIR, "FigS10_Survival.png"),
          file.path(REVIEW_DIR, "FigS10_Survival.png"), overwrite = TRUE)
message("FigS10_Survival.pdf (160×90mm) ✓")
