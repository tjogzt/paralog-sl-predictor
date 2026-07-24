library(data.table)

# --- gyn_drivers ---
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

# --- known_sl_pairs ---
known_sl_pairs <- data.frame(
  gene_A = c("SMARCA4", "ARID1A", "BRCA1", "EP300", "PIK3CA",
             "AKT1", "STK11", "FBXW7", "PPP2R1A", "CCNE1",
             "CDK4", "MAP2K1"),
  gene_B = c("SMARCA2", "ARID1B", "BRCA2", "CREBBP", "PIK3CB",
             "AKT2", "SIK1", "FBXW2", "PPP2R1B", "CCNE2",
             "CDK6", "MAP2K2"),
  evidence = rep("literature", 12),
  stringsAsFactors = FALSE
)

# --- benchmark_methods ---
benchmark_methods <- data.frame(
  Method = c("SLMGAE", "NSF4SL", "GCATSL", "GRSMF", "PiLSL", "KG4SL",
             "SLGNN", "PTGNN", "DD (this study)", "DD + ID >= 0.3"),
  CV3_AUROC = c(0.790, 0.683, 0.678, 0.656, 0.626, 0.563, 0.530,
                0.529, 0.794, 1.000),
  Reference = c(rep("Feng et al. 2024 (SD1)", 8), "This study", "This study"),
  Interpretability = c(rep("Low", 8), "High", "High"),
  stringsAsFactors = FALSE
)

# --- cross_cancer_summary ---
csv_path <- "../output/cross_cancer_summary.csv"
if (file.exists(csv_path)) {
  cross_cancer_summary <- as.data.frame(fread(csv_path))
} else {
  cross_cancer_summary <- data.frame(
    cancer     = c("Breast", "Ovarian", "Endometrial", "Cervical", "Lung"),
    n_lines    = c(50L, 55L, 28L, 16L, 95L),
    n_pairs    = c(11L, 43L, 68L, 7L, 77L),
    n_known    = c(5L, 4L, 4L, 3L, 5L),
    dd_auroc   = c(0.889, 0.846, 0.797, 0.667, 0.353),
    stringsAsFactors = FALSE
  )
}

# --- solid_tumor_summary ---
csv_path2 <- "../output/solid_tumor_summary.csv"
if (file.exists(csv_path2)) {
  solid_tumor_summary <- as.data.frame(fread(csv_path2))
} else {
  solid_tumor_summary <- data.frame(
    cancer   = c("Ovarian", "Endometrial", "Cervical", "Breast", "Lung",
                 "Colorectal", "Esophagogastric", "Biliary Tract",
                 "HNSCC", "NSCLC", "Pancreatic", "Melanoma",
                 "Hepatocellular", "Glioma", "Bladder Urothelial",
                 "Renal Cell", "SCLC", "Mesothelioma", "Neuroblastoma",
                 "Osteosarcoma", "Ewing Sarcoma", "Other Sarcoma",
                 "Rhabdomyosarcoma", "Thyroid"),
    n_lines  = c(55L, 28L, 16L, 50L, 95L, 59L, 62L, 39L,
                 26L, 112L, 43L, 46L, 46L, 37L, 63L,
                 11L, 123L, 9L, 10L, 14L, 5L, 9L,
                 3L, 3L),
    stringsAsFactors = FALSE
  )
}

# --- therapeutic_window_summary ---
csv_path3 <- "../output/therapeutic_window_summary.csv"
if (file.exists(csv_path3)) {
  therapeutic_window_summary <- as.data.frame(fread(csv_path3))
} else {
  therapeutic_window_summary <- data.frame(
    context        = c("Ovarian", "Endometrial", "Breast", "Colorectal"),
    n_lines        = c(55L, 28L, 50L, 59L),
    n_pairs        = c(18L, 20L, 13L, 20L),
    n_known        = c(9L, 9L, 6L, 10L),
    n_selective    = c(5L, 6L, 5L, 4L),
    n_good_window  = c(7L, 10L, 7L, 13L),
    stringsAsFactors = FALSE
  )
}

# Save all datasets
save_one <- function(obj, name) {
  assign(name, obj)
  save(list = name, file = file.path("data", paste0(name, ".rda")),
       compress = "gzip")
  cat(sprintf("Saved %s.rda\n", name))
}

save_one(gyn_drivers, "gyn_drivers")
save_one(known_sl_pairs, "known_sl_pairs")
save_one(benchmark_methods, "benchmark_methods")
save_one(cross_cancer_summary, "cross_cancer_summary")
save_one(solid_tumor_summary, "solid_tumor_summary")
save_one(therapeutic_window_summary, "therapeutic_window_summary")

cat("All datasets generated.\n")
