# compute_sequence_identity.R
# ============================
# Computes k-mer Jaccard sequence identity (k=3, amino acids) for every
# driver-paralog pair in TableS2 from bundled UniProt canonical sequences
# (data/uniprot_sequences.rds), and writes output/paralog_identity.csv.
#
# Motivation: data/ensembl_paralogs.csv has identity_pct = NA for all
# 66,595 pairs, so the manuscript's "DD + sequence-identity filter" AUROC
# values were previously not reproducible from any artifact. This helper
# recomputes the proxy exactly as described in Methods (k-mer Jaccard,
# k=3, validated against Needleman-Wunsch identity at r=0.88).
#
# Usage: Rscript compute_sequence_identity.R   (run from repo root)

pairs_path   <- "output/tables/TableS2_FullResults.tsv"
seqs_path    <- "data/uniprot_sequences.rds"
out_path     <- "output/paralog_identity.csv"

seqs  <- readRDS(seqs_path)
pairs <- read.delim(pairs_path, stringsAsFactors = FALSE)
pairs <- unique(pairs[, c("driver_gene", "paralog_gene")])

kmers <- function(seq, k = 3L) {
  s <- strsplit(toupper(seq), "", fixed = TRUE)[[1]]
  if (length(s) < k) return(character(0))
  starts <- seq_len(length(s) - k + 1L)
  unique(vapply(starts, function(i) paste0(s[i:(i + k - 1L)], collapse = ""), character(1)))
}

jaccard <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0L) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

kmer_cache <- lapply(seqs, kmers)

res <- pairs
res$kmer_jaccard <- mapply(function(a, b) {
  if (!a %in% names(kmer_cache) || !b %in% names(kmer_cache)) return(NA_real_)
  jaccard(kmer_cache[[a]], kmer_cache[[b]])
}, res$driver_gene, res$paralog_gene)

write.csv(res, out_path, row.names = FALSE, quote = FALSE)
message(sprintf("paralog_identity.csv written: %d pairs, %d with identity values (%.0f%% coverage)",
                nrow(res), sum(!is.na(res$kmer_jaccard)),
                100 * mean(!is.na(res$kmer_jaccard))))
missing <- unique(unlist(res[is.na(res$kmer_jaccard), c("driver_gene", "paralog_gene")]))
if (length(missing)) message("  genes without sequences: ", paste(sort(missing), collapse = ", "))
