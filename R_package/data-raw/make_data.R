library(data.table)

# Regenerate all built-in datasets from the canonical artifacts of the
# companion analysis repository (paralog-sl-predictor). Fails loudly when a
# source file is missing — a silent fallback here was the root cause of the
# v1.0.x stale-data drift.
CANON <- "../paralog_sl_predictor/output"
need <- function(...) {
  p <- file.path(...)
  if (!file.exists(p)) stop("canonical artifact not found: ", p,
                            " (run the paralog-sl-predictor pipeline first)")
  p
}

# --- gyn_drivers (example driver panels; not manuscript claims) ---
gyn_drivers <- list(
  Ovarian     = c("TP53", "BRCA1", "BRCA2", "ARID1A", "CCNE1",
                  "NF1", "RB1", "PTEN", "PIK3CA", "KRAS"),
  Endometrial = c("PTEN", "ARID1A", "PIK3CA", "CTNNB1", "TP53",
                  "PPP2R1A", "KRAS", "PIK3R1", "FBXW7", "KMT2D"),
  Cervical    = c("PIK3CA", "EP300", "FBXW7", "STK11", "ERBB2",
                  "MAPK1", "PTEN", "KRAS"),
  Lung        = c("TP53", "KRAS", "EGFR", "STK11", "KEAP1",
                  "NF1", "BRAF", "PIK3CA", "ALK", "MET"),
  Breast      = c("TP53", "PIK3CA", "PTEN", "BRCA1", "ERBB2",
                  "GATA3", "CDH1", "RB1", "NF1", "MAP3K1")
)

# --- known_sl_pairs: evidence-tiered gold standard (supplementary Table S3) ---
s3 <- as.data.frame(fread(need(CANON, "tables", "TableS3_GoldStandard.tsv")))
known_sl_pairs <- data.frame(
  gene_A    = s3$Driver,
  gene_B    = s3$Paralog,
  tier      = s3$Tier,
  direct_sl = s3$Direct_SL,
  inclusion = s3$Inclusion,
  key_ref   = s3$Key_Ref,
  stringsAsFactors = FALSE
)

# --- benchmark_methods: published CV3 references + this-study rows ---
bench <- as.data.frame(fread(need(CANON, "tables", "Table2_Benchmark.tsv")))
benchmark_methods <- data.frame(
  Method           = bench$Method,
  CV3_AUROC        = as.numeric(bench$CV3_AUROC),
  Reference        = bench$Source,
  Interpretability = bench$Interpretability,
  stringsAsFactors = FALSE
)

# --- solid_tumor_summary: primary (min >= 5 per group) pan-cancer frame ---
solid_tumor_summary <- as.data.frame(fread(need(CANON, "solid_tumor_summary.csv")))

# --- cross_cancer_summary: sensitivity (min >= 3 per group) frame ---
cross_cancer_summary <- as.data.frame(fread(need(CANON, "solid_tumor_summary_min3.csv")))

# --- therapeutic_window_summary: 21-pair DWS classification (Table S5) ---
therapeutic_window_summary <- as.data.frame(
  fread(need(CANON, "therapeutic_window_paralog_classification.csv")))

# Save all datasets
save_one <- function(obj, name) {
  assign(name, obj)
  save(list = name, file = file.path("data", paste0(name, ".rda")),
       compress = "gzip")
  cat(sprintf("Saved %s.rda (%d rows)\n", name, NROW(obj)))
}

save_one(gyn_drivers, "gyn_drivers")
save_one(known_sl_pairs, "known_sl_pairs")
save_one(benchmark_methods, "benchmark_methods")
save_one(cross_cancer_summary, "cross_cancer_summary")
save_one(solid_tumor_summary, "solid_tumor_summary")
save_one(therapeutic_window_summary, "therapeutic_window_summary")

cat("All datasets generated from", CANON, "\n")
