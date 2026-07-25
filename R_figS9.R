# FigS9 — k-mer Jaccard vs Needleman-Wunsch Validation (R)
# Purpose: Validate k-mer (k=3) Jaccard as proxy for global alignment identity
# Data:    Real UniProt protein sequences for paralog genes
# Output:  90×90mm scatter plot
library(ggplot2)
library(cowplot)
library(dplyr)
library(readr)
library(httr)

OUT_DIR <- "paralog_sl_predictor/output/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_FS <- 7; TICK_FS <- 7
BLUE  <- "#2171B5"; RED   <- "#CB181D"; GRAY  <- "#636363"

theme_sci <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.line    = element_line(linewidth = 0.4),
  axis.ticks   = element_line(linewidth = 0.3),
  axis.text    = element_text(size = TICK_FS),
  axis.title   = element_text(size = BASE_FS),
  plot.margin  = margin(4, 4, 4, 4, "pt"),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA))

# ── Gene list: all unique genes from structural analysis ──
struct_path <- "paralog_sl_predictor/output/alphafold_structural_analysis.csv"
feat_path   <- "paralog_sl_predictor/data/protein_features.csv"

if (file.exists(struct_path)) {
  st <- read_csv(struct_path, show_col_types = FALSE)
  pair_genes <- unique(c(st$gene_a, st$gene_b))
} else {
  pair_genes <- c("ARID1A","ARID1B","EP300","CREBBP","PIK3CA","PIK3CB",
                  "SMARCA4","SMARCA2","PPP2R1A","PPP2R1B","KRAS","HRAS","NRAS",
                  "BRCA1","BRCA2","TP53","TP63","FBXW7","FBXW2","STK11","SIK1",
                  "KMT2D","KMT2C","NF1","RASA2","AKT1","AKT2","PTEN","TNS2",
                  "CDH1","CDH2","MAP2K1","MAP2K2","BRAF","RAF1","RB1","RBL1",
                  "ATR","ATM")
}

# ── Download UniProt sequences ──
message("Downloading UniProt sequences...")
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

# Cache sequences
cache_file <- "paralog_sl_predictor/data/uniprot_sequences.rds"
if (file.exists(cache_file)) {
  seqs <- readRDS(cache_file)
  # Download missing genes
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

# ── k-mer Jaccard (k=3) ──
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

# ── Needleman-Wunsch global alignment identity (simple implementation) ──
nw_identity <- function(seq_a, seq_b) {
  if (is.na(seq_a) || is.na(seq_b)) return(NA)
  # Use REAL Needleman-Wunsch via Bioconductor: pwalign (BioC >= 3.19)
  # with Biostrings fallback for environments where pwalign is absent.
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

# ── Build ALL pairwise combinations from genes with sequences ──
message("Computing k-mer Jaccard and NW identity for pairs...")
available_genes <- names(seqs)[!sapply(seqs, function(x) is.null(x) || identical(x, NA_character_))]
all_pairs <- t(combn(available_genes, 2))

# Mark known paralog pairs
if (file.exists(struct_path)) {
  known_pairs_set <- paste0(st$gene_a, "_", st$gene_b)
} else {
  known_pairs_set <- character(0)
}

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

# ── Plot ──
if (nrow(results) >= 5) {
  r_val <- cor(results$kmer_jaccard, results$nw_identity)
  t_stat <- r_val * sqrt((nrow(results) - 2) / (1 - r_val^2))
  p_val <- 2 * pt(abs(t_stat), nrow(results) - 2, lower.tail = FALSE)

  p <- ggplot(results, aes(nw_identity, kmer_jaccard, color = is_paralog)) +
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

  ggsave(file.path(OUT_DIR, "FigS9_kmer_validation.pdf"), p,
         width = 90, height = 90, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS9_kmer_validation.svg"), p,
       width = 90, height = 90, units = "mm", device = svglite::svglite)
ggsave(file.path(OUT_DIR, "FigS9_kmer_validation.tiff"), p,
       width = 90, height = 90, units = "mm", device = ragg::agg_tiff, dpi = 600)
  message("FigS9_kmer_validation.pdf (90×90mm) ✓")
} else {
  message("ERROR: too few pairs for validation plot")
}
