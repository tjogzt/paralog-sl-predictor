# FigS5 — PRISM Drug Selectivity Heatmap (R)
# Purpose: Heatmap of top drug × paralog-SL pair ΔAUC values
# Canvas: 180×135mm
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(reshape2)

OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FS <- 7; TICK_FS <- 6

theme_sci <- theme_classic(base_size = BASE_FS) + theme(
  panel.grid = element_blank(),
  axis.line    = element_blank(),
  axis.ticks   = element_line(linewidth = 0.3),
  axis.text.x  = element_text(size = TICK_FS, angle = 45, hjust = 1, vjust = 1),
  axis.text.y  = element_text(size = TICK_FS, hjust = 1),
  axis.title   = element_blank(),
  legend.position = "right",
  legend.key.height = unit(0.5, "cm"),
  legend.key.width  = unit(0.3, "cm"),
  plot.margin  = margin(4, 4, 4, 4, "pt"))

# ── Load data ──
prism_path <- "paralog_sl_predictor/output/prism_top_hits.csv"
if (file.exists(prism_path)) {
  pr <- read_csv(prism_path, show_col_types = FALSE)
} else {
  stop("prism_top_hits.csv not found")
}

# Select top drugs and top pairs by effect size
top_drugs <- pr %>%
  group_by(drug) %>%
  summarise(min_delta = min(delta_auc, na.rm = TRUE)) %>%
  arrange(min_delta) %>%
  head(20) %>%
  pull(drug)

top_pairs <- pr %>%
  mutate(pair_label = paste0(driver, "->", paralog)) %>%
  group_by(pair_label) %>%
  summarise(min_delta = min(delta_auc, na.rm = TRUE)) %>%
  arrange(min_delta) %>%
  head(15) %>%
  pull(pair_label)

# Filter to intersection
pr_sub <- pr %>%
  mutate(pair_label = paste0(driver, "->", paralog)) %>%
  filter(drug %in% top_drugs, pair_label %in% top_pairs) %>%
  group_by(drug, pair_label) %>%
  summarise(delta_auc = min(delta_auc, na.rm = TRUE), .groups = "drop")

# Shorten drug names
pr_sub$drug_short <- substr(pr_sub$drug, 1, 22)
pr_sub$drug_short <- factor(pr_sub$drug_short, levels = rev(sort(unique(pr_sub$drug_short))))
pr_sub$pair_label <- factor(pr_sub$pair_label, levels = sort(unique(pr_sub$pair_label)))

# Cap delta_auc for visualization
pr_sub$delta_capped <- pmax(pr_sub$delta_auc, -0.7)

# ── Heatmap ──
p <- ggplot(pr_sub, aes(pair_label, drug_short, fill = delta_capped)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-0.7, 0.1),
                       name = "ΔAUC\n(MUT-WT)") +
  theme_sci

ggsave(file.path(OUT_DIR, "FigS5_PRISM_Selectivity.pdf"), p,
       width = 180, height = 120, units = "mm", device = cairo_pdf)
message("FigS5_PRISM_Selectivity.pdf (180×120mm) ✓")
