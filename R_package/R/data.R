#' Cross-Cancer Summary
#'
#' A data.frame of gynecological cancer cross-cancer analysis results.
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{cancer}{Cancer type name}
#'   \item{n_lines}{Number of cell lines}
#'   \item{n_pairs}{Number of paralog pairs tested}
#'   \item{n_known}{Number of known paralog-SL pairs}
#'   \item{dd_auroc}{DD AUROC against gold standard}
#' }
"cross_cancer_summary"

#' Solid Tumor Pan-Cancer Summary
#'
#' A data.frame of DD AUROC results across 24 solid tumor types.
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{cancer}{Cancer type name}
#'   \item{n_lines}{Number of available cell lines}
#' }
"solid_tumor_summary"

#' Therapeutic Window Summary
#'
#' A data.frame summarizing therapeutic index (TI) metrics across contexts.
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{context}{Cancer/tissue context}
#'   \item{n_lines}{Number of cell lines}
#'   \item{n_pairs}{Number of paralog pairs}
#'   \item{n_known}{Number of known SL pairs}
#'   \item{n_selective}{Number of selective candidates}
#'   \item{n_good_window}{Number with favorable TI}
#' }
"therapeutic_window_summary"

#' Gynecological and Breast/Lung Cancer Driver Genes
#'
#' A named list of driver genes for five cancer types, curated from
#' TCGA, COSMIC, and literature.
#'
#' @format A named list with 5 character vectors:
#' \describe{
#'   \item{Ovarian}{10 genes: TP53, BRCA1, BRCA2, ARID1A, CCNE1, NF1, RB1, PTEN, PIK3CA, KRAS}
#'   \item{Endometrial}{10 genes: PTEN, ARID1A, PIK3CA, CTNNB1, TP53, PPP2R1A, KRAS, PIK3R1, FBXW7, KMT2D}
#'   \item{Cervical}{8 genes: PIK3CA, EP300, FBXW7, STK11, ERBB2, MAPK1, PTEN, KRAS}
#'   \item{Lung}{10 genes: TP53, KRAS, EGFR, STK11, KEAP1, NF1, BRAF, PIK3CA, ALK, MET}
#'   \item{Breast}{10 genes: TP53, PIK3CA, PTEN, BRCA1, ERBB2, GATA3, CDH1, RB1, NF1, MAP3K1}
#' }
"gyn_drivers"

#' Known Paralog-SL Pairs (Gold Standard)
#'
#' A data.frame of 12 literature-curated paralog synthetic lethal pairs
#' used as the gold-standard positive set for validation.
#'
#' @format A data.frame with 12 rows and 2 columns:
#' \describe{
#'   \item{gene_A}{First gene in the SL pair}
#'   \item{gene_B}{Second gene in the SL pair}
#' }
"known_sl_pairs"

#' Published SL Method Benchmark (CV3)
#'
#' CV3 (gene-pair isolation) AUROC values from Feng et al. (2024)
#' Nature Communications benchmark of 8 SL prediction methods,
#' compared with our DD approach.
#'
#' @format A data.frame with 10 rows and 3 columns:
#' \describe{
#'   \item{Method}{Method name}
#'   \item{CV3_AUROC}{AUROC in the most stringent CV3 setting}
#'   \item{Interpretability}{Whether the method is interpretable}
#' }
"benchmark_methods"
