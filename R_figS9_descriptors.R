# FigS9 — Sequence- and structure-derived descriptors (exploratory)
# Merged 2026-07-29: former FigS9 (panel a: k-mer Jaccard vs Needleman-Wunsch
# validation), former FigS10 (panel b: druggability), former FigS13
# (panels c/d: structural/domain conservation + composite prioritization score).
# Sources: UniProt reviewed sequences (cached), druggability_analysis.json,
#          alphafold_structural_analysis.csv — all real pipeline artifacts.
# Usage:   Rscript R_figS9_descriptors.R   (run from the project root)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)
library(httr)
library(jsonlite)

OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FS <- 7; TICK_FS <- 7; LEGEND_FS <- 7
BLUE  <- "#2171B5"; RED   <- "#CB181D"; ORANGE <- "#E6550D"
TEAL  <- "#0D7377"; GRAY  <- "#636363"; DARK   <- "#252525"

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
  legend.key.size   = unit(3, "mm"),
  plot.margin  = margin(4, 4, 4, 4, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

# ═══════════════════════════════════════════════════════════════
# PANEL A — k-mer Jaccard (k=3) vs Needleman-Wunsch identity validation
# ═══════════════════════════════════════════════════════════════
struct_path <- "paralog_sl_predictor/output/alphafold_structural_analysis.csv"

if (file.exists(struct_path)) {
  st <- read_csv(struct_path, show_col_types = FALSE)
  pair_genes <- unique(c(st$gene_a, st$gene_b))
} else {
  stop("paralog_sl_predictor/output/alphafold_structural_analysis.csv not found — ",
       "run the pipeline first; simulated fallbacks are forbidden")
}

# ── Download UniProt sequences (cached) ──
message("Loading UniProt sequences (cache preferred)...")
get_uniprot_seq <- function(gene_name) {
  url <- sprintf("https://rest.uniprot.org/uniprotkb/search?query=gene_exact:%s+AND+organism_id:9606+AND+reviewed:true&format=fasta&size=1",
                 gene_name)
  tryCatch({
    resp <- GET(url, timeout(10))
    if (status_code(resp) == 200) {
      txt <- content(resp, "text", encoding = "UTF-8")
      lines <- strsplit(txt, "\n")[[1]]
      seq_lines <- lines[!grepl("^>", lines)]
      seq <- paste(seq_lines, collapse = "")
      if (nchar(seq) > 10) return(seq)
    }
    return(NA_character_)
  }, error = function(e) NA_character_)
}

cache_file <- "paralog_sl_predictor/data/uniprot_sequences.rds"
if (file.exists(cache_file)) {
  seqs <- readRDS(cache_file)
  missing <- setdiff(pair_genes, names(seqs))
  if (length(missing) > 0) {
    for (g in missing) {
      seqs[[g]] <- get_uniprot_seq(g)
      Sys.sleep(0.3)
    }
    saveRDS(seqs, cache_file)
  }
} else {
  seqs <- list()
  for (g in pair_genes) {
    message(sprintf("  %s...", g))
    seqs[[g]] <- get_uniprot_seq(g)
    Sys.sleep(0.3)
  }
  saveRDS(seqs, cache_file)
}
message(sprintf("  %d/%d sequences obtained", sum(!is.na(seqs)), length(seqs)))

kmer_jaccard <- function(seq_a, seq_b, k = 3) {
  if (is.na(seq_a) || is.na(seq_b)) return(NA)
  get_kmers <- function(s) {
    n <- nchar(s)
    if (n < k) return(character(0))
    unique(sapply(1:(n - k + 1), function(i) substr(s, i, i + k - 1)))
  }
  ka <- get_kmers(seq_a)
  kb <- get_kmers(seq_b)
  if (length(ka) == 0 || length(kb) == 0) return(0)
  length(intersect(ka, kb)) / length(union(ka, kb))
}

# REAL Needleman-Wunsch via Bioconductor: pwalign (BioC >= 3.19) with
# Biostrings fallback
nw_identity <- function(seq_a, seq_b) {
  if (is.na(seq_a) || is.na(seq_b)) return(NA)
  aln_ns <- if (requireNamespace("pwalign", quietly = TRUE)) "pwalign" else "Biostrings"
  tryCatch({
    pa_fun  <- get("pairwiseAlignment", envir = asNamespace(aln_ns))
    pid_fun <- get("pid", envir = asNamespace(aln_ns))
    aln <- pa_fun(
      Biostrings::AAString(seq_a),
      Biostrings::AAString(seq_b),
      type = "global",
      substitutionMatrix = "BLOSUM62",
      gapOpening = 10, gapExtension = 0.5)
    pid_fun(aln, type = "PID1") / 100
  }, error = function(e) NA)
}

message("Computing k-mer Jaccard and NW identity for pairs...")
available_genes <- names(seqs)[!sapply(seqs, function(x) is.null(x) || identical(x, NA_character_))]
all_pairs <- t(combn(available_genes, 2))
known_pairs_set <- paste0(st$gene_a, "_", st$gene_b)

results <- data.frame()
for (i in 1:nrow(all_pairs)) {
  ga <- all_pairs[i, 1]; gb <- all_pairs[i, 2]
  sa <- seqs[[ga]]; sb <- seqs[[gb]]

  kj <- kmer_jaccard(sa, sb, k = 3)
  nw <- nw_identity(sa, sb)
  if (is.na(kj) || is.na(nw)) next

  is_paralog <- paste0(ga, "_", gb) %in% known_pairs_set |
                paste0(gb, "_", ga) %in% known_pairs_set

  results <- rbind(results, data.frame(
    gene_a = ga, gene_b = gb,
    kmer_jaccard = kj, nw_identity = nw,
    is_paralog = is_paralog,
    pair_label = paste0(ga, "/", gb)))
}

# Sample to ~50 pairs: keep all known paralogs + random non-paralogs
paralog_rows <- results[results$is_paralog, ]
nonparalog_rows <- results[!results$is_paralog, ]
set.seed(42)
n_sample <- max(0, 50 - nrow(paralog_rows))
if (nrow(nonparalog_rows) > n_sample) {
  nonparalog_rows <- nonparalog_rows[sample(nrow(nonparalog_rows), n_sample), ]
}
results <- rbind(paralog_rows, nonparalog_rows)

message(sprintf("  %d pairs (%d paralog, %d non-paralog)", nrow(results),
                sum(results$is_paralog), sum(!results$is_paralog)))
stopifnot(nrow(results) >= 5)

r_val <- cor(results$kmer_jaccard, results$nw_identity)
t_stat <- r_val * sqrt((nrow(results) - 2) / (1 - r_val^2))
p_val <- 2 * pt(abs(t_stat), nrow(results) - 2, lower.tail = FALSE)

# Subgroup correlations — single source of truth for the three r values
# quoted in manuscript Methods; artifact verified by audit_manuscript_numbers.py
pr <- results[results$is_paralog, ]
nr <- results[!results$is_paralog, ]
cor_p <- function(r, n) {
  if (n < 3 || is.na(r)) return(NA_real_)
  t <- r * sqrt((n - 2) / (1 - r^2))
  2 * pt(abs(t), n - 2, lower.tail = FALSE)
}
r_paralog <- if (nrow(pr) >= 3) cor(pr$kmer_jaccard, pr$nw_identity) else NA_real_
r_nonparalog <- if (nrow(nr) >= 3) cor(nr$kmer_jaccard, nr$nw_identity) else NA_real_
cor_artifact <- data.frame(
  group = c("all", "paralog", "non_paralog"),
  n = c(nrow(results), nrow(pr), nrow(nr)),
  pearson_r = c(r_val, r_paralog, r_nonparalog),
  p_value = c(p_val, cor_p(r_paralog, nrow(pr)), cor_p(r_nonparalog, nrow(nr))))
write.csv(cor_artifact,
          "paralog_sl_predictor/output/kmer_nw_correlation.csv",
          row.names = FALSE)
message(sprintf("  r: all=%.3f (n=%d), paralog=%.3f (n=%d), non-paralog=%.3f (n=%d)",
                r_val, nrow(results), r_paralog, nrow(pr), r_nonparalog, nrow(nr)))

pa <- ggplot(results, aes(nw_identity, kmer_jaccard, color = is_paralog)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(name = NULL, values = c(`TRUE` = RED, `FALSE` = BLUE),
                     labels = c(`TRUE` = "Known paralog", `FALSE` = "Non-paralog")) +
  geom_smooth(method = "lm", se = TRUE, color = GRAY, linewidth = 0.5, alpha = 0.15,
              inherit.aes = FALSE, aes(x = nw_identity, y = kmer_jaccard)) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("n = %d pairs\nr = %.3f\np = %.1e", nrow(results), r_val, p_val),
           hjust = -0.05, vjust = 1.2, size = 2.8, color = RED, fontface = "bold") +
  labs(x = "Needleman-Wunsch alignment identity",
       y = "k-mer Jaccard similarity (k=3)") +
  theme_sci +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))

# ═══════════════════════════════════════════════════════════════
# PANEL B — Structure-based druggability assessment
# ═══════════════════════════════════════════════════════════════
dg <- fromJSON("paralog_sl_predictor/output/druggability_analysis.json")
stopifnot(nrow(dg) == 15)
dg$is_paralog <- grepl("^Paralog", dg$role)
# Pair-grouped order (paralog first within each pair); reverse for ggplot's
# bottom-up y axis.
dg$gene <- factor(dg$gene, levels = rev(dg$gene))
dg$kp <- sprintf("K=%d P=%d", dg$lys_count, dg$pocket_regions)

pb <- ggplot(dg, aes(druggability_score, gene, fill = is_paralog)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = kp), hjust = 1.08, size = 2.5, color = "white") +
  geom_text(aes(label = sprintf("%.3f", druggability_score)),
            hjust = -0.08, size = 2.5, color = DARK) +
  scale_fill_manual(name = NULL,
                    values = c(`TRUE` = RED, `FALSE` = BLUE),
                    labels = c(`TRUE` = "Paralog (target)", `FALSE` = "Driver")) +
  scale_x_continuous(limits = c(0, 1.18), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Druggability score", y = NULL) +
  theme_sci +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = NA))

# ═══════════════════════════════════════════════════════════════
# PANEL C — Structural similarity and domain conservation (top 8 pairs)
# ═══════════════════════════════════════════════════════════════
stc <- st %>%
  filter(domain_similarity > 0) %>%  # exclude zero-domain pairs (e.g. BRCA1/BRCA2)
  arrange(desc(structural_similarity)) %>% head(8)
stc$pair <- factor(paste0(stc$gene_a, "/", stc$gene_b),
                   levels = rev(paste0(stc$gene_a, "/", stc$gene_b)))
dfc <- stc %>% select(pair, structural_similarity, domain_similarity) %>%
  pivot_longer(-pair) %>%
  mutate(name = recode(name, structural_similarity = "Structural",
                       domain_similarity = "Domain"))

pc <- ggplot(dfc, aes(value, pair, fill = name)) +
  geom_col(position = position_dodge(0.7), width = 0.55) +
  scale_fill_manual(values = c(Structural = TEAL, Domain = ORANGE)) +
  labs(x = "Score", y = NULL) +
  # legend below the panel: inside placement overlapped the short bottom bars;
  # reverse legend so entry order matches top-to-bottom bar order
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_sci + theme(legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.margin = margin(0, 0, 0, 0),
                    plot.margin = margin(4, 4, 4, 20, "pt"))

# ═══════════════════════════════════════════════════════════════
# PANEL D — Composite prioritization score ranking (top 10)
# ═══════════════════════════════════════════════════════════════
if (!("clinical_targetability" %in% names(st))) {
  stop("clinical_targetability column missing in ",
       "alphafold_structural_analysis.csv — run the structural/targetability ",
       "pipeline first; simulated fallbacks are forbidden")
}
cand <- st %>% arrange(desc(clinical_targetability)) %>% head(10) %>%
  mutate(label = paste0(driver, "->", paralog), score = clinical_targetability)
cand$label <- factor(cand$label, levels = rev(cand$label))
# Highlight matches the signed-DWS narrative: with the primary signed
# max(DD,0) DWS formula, ARID1A->ARID1B ranks first on the composite score
# and is the leading SELECTIVE candidate; NF1->RASA2 retains a high DWS
# through a near-zero pan-essential denominator (selectivity ~ 0) but no
# longer ranks first
cand$cat <- "Others"
cand$cat[cand$driver == "NF1"    & cand$paralog == "RASA2"]  <- "High DWS (non-selective)"
cand$cat[cand$driver == "ARID1A" & cand$paralog == "ARID1B"] <- "Rank 1 (selective candidate)"
cand$cat <- factor(cand$cat, levels = c("Rank 1 (selective candidate)",
                                        "High DWS (non-selective)", "Others"))
cat_colors <- c("Rank 1 (selective candidate)" = RED,
                "High DWS (non-selective)" = ORANGE, "Others" = BLUE)
cand$txt_col <- cat_colors[as.character(cand$cat)]

pd <- ggplot(cand, aes(score, label, fill = cat)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.3f", score), color = txt_col),
            hjust = -0.1, size = 2.5, fontface = "bold", show.legend = FALSE) +
  # direct in-bar category labels instead of a legend: a 3-entry legend
  # overlapped the bottom bars and their value labels in this 90mm panel
  geom_text(data = cand %>% filter(cat != "Others"),
            aes(x = score - 0.02, label = as.character(cat)),
            hjust = 1, size = 2.5, color = "white", fontface = "bold",
            family = "Arial", show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_manual(values = cat_colors, guide = "none") +
  geom_vline(xintercept = 0.5, linewidth = 0.3, color = GRAY,
             linetype = "dashed", alpha = 0.3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Composite prioritization score", y = NULL) +
  theme_sci

# ═══════════════════════════════════════════════════════════════
# COMPOSITE — 2x2 grid, 180x170mm
# ═══════════════════════════════════════════════════════════════
p <- cowplot::plot_grid(pa, pb, pc, pd, ncol = 2,
                        labels = c("a", "b", "c", "d"),
                        label_size = 9, label_fontface = "bold",
                        label_fontfamily = "Arial")

ggsave(file.path(OUT_DIR, "FigS9_Sequence_Structure_Descriptors.pdf"), p,
       width = 180, height = 170, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS9_Sequence_Structure_Descriptors.svg"), p,
       width = 180, height = 170, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS9_Sequence_Structure_Descriptors.tiff"), p,
       width = 180, height = 170, units = "mm", device = ragg::agg_tiff, dpi = 600)
message("FigS9_Sequence_Structure_Descriptors.pdf (180x170mm, 4 panels) ✓")
