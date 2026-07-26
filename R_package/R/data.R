#' Cross-Cancer Summary (sensitivity frame)
#'
#' Per-lineage DD benchmark results on the manuscript's sensitivity frame
#' (minimum 3 mutant and 3 wild-type cell lines per driver-paralog group),
#' covering 23 solid tumor lineages. Regenerated from
#' `solid_tumor_summary_min3.csv` in the companion analysis repository
#' (paralog-sl-predictor); 12 lineages are evaluable (>= 2 known positives),
#' 9 of them with AUROC >= 0.7.
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{cancer}{Cancer type name}
#'   \item{n_lines}{Number of cell lines}
#'   \item{n_pairs}{Number of paralog pairs tested}
#'   \item{n_known}{Number of known paralog-SL pairs}
#'   \item{dd_auroc}{DD AUROC against gold standard (NA when < 2 known positives)}
#' }
"cross_cancer_summary"

#' Solid Tumor Pan-Cancer Summary (primary frame)
#'
#' Per-lineage DD benchmark results on the manuscript's primary frame
#' (minimum 5 mutant and 5 wild-type cell lines per driver-paralog group),
#' covering 19 solid tumor lineages; 8 are evaluable (>= 2 known positives)
#' and 7 reach AUROC >= 0.7. Regenerated from `solid_tumor_summary.csv` in
#' the companion analysis repository (paralog-sl-predictor).
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{cancer}{Cancer type name}
#'   \item{n_lines}{Number of cell lines}
#'   \item{n_pairs}{Number of paralog pairs tested}
#'   \item{n_known}{Number of known paralog-SL pairs}
#'   \item{dd_auroc}{DD AUROC against gold standard (NA when < 2 known positives)}
#' }
"solid_tumor_summary"

#' Therapeutic Window Summary (21-pair DWS classification)
#'
#' Per-pair dependency window score (DWS; reported as therapeutic index, TI,
#' in earlier versions) classification across evaluated lineages, from the
#' manuscript's therapeutic-window module (Table S6). Regenerated from
#' `therapeutic_window_paralog_classification.csv` in the companion
#' analysis repository (paralog-sl-predictor).
#'
#' @format A data.frame with 21 rows and columns:
#' \describe{
#'   \item{driver}{Driver gene}
#'   \item{paralog}{Paralog gene}
#'   \item{mean_ti}{Mean DWS/TI across lineages}
#'   \item{mean_dd}{Mean |DD| across lineages}
#'   \item{mean_selectivity}{Mean essential-fraction difference (MUT - WT)}
#'   \item{mean_pan_essential}{Mean pan-essential fraction of the paralog}
#'   \item{n_contexts}{Number of evaluable lineage contexts}
#'   \item{classification}{HIGH_SELECTIVITY / MODERATE / LOW_SELECTIVITY / PAN_ESSENTIAL}
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

#' Known Paralog-SL Pairs (evidence-tiered gold standard)
#'
#' The manuscript's 12 curated pairs with evidence tiers after verification
#' of the primary citations (supplementary Table S3). Tier A (3 pairs):
#' direct dual-gene perturbation synthetic-lethal evidence. Tier B (2 pairs):
#' natural-genotype conditional dependency with functional validation. The
#' Tier A + Tier B set constitutes the primary external benchmark used by
#' default in `compute_auroc()`. Tier C (5 pairs): indirect evidence only
#' (reciprocal-direction-only, other-driver, developmental redundancy, or
#' DepMap-derived). Comparators (2 pairs): mechanistic reference pairs that
#' are not sequence paralogs, used as specificity references only.
#'
#' @format A data.frame with 12 rows and columns:
#' \describe{
#'   \item{gene_A}{Driver gene}
#'   \item{gene_B}{Paralog gene}
#'   \item{tier}{Evidence tier: A, B, C, or Comparator}
#'   \item{direct_sl}{Nature of the direct evidence (Yes / Conditional / Reciprocal / ...)}
#'   \item{inclusion}{Benchmark role: Primary / Secondary / Comparator}
#'   \item{key_ref}{Key reference(s)}
#' }
"known_sl_pairs"

#' Published SL Method Benchmark (CV3)
#'
#' CV3 (gene-pair isolation) AUROC values from Feng et al. (2024)
#' Nature Communications benchmark of 8 SL prediction methods, plus
#' this-study DD values. Published values are contextual reference points
#' from a general SL gene-pair universe, NOT a head-to-head benchmark on
#' the same test set; DD was evaluated on paralog-SL pairs only (full
#' lineage-level frame, 110 entries / 8 positives, AUROC 0.676; DD + ID
#' >= 30% on a 3-pair high-identity subset, 2 positives, AUROC 1.000,
#' anecdotal).
#'
#' @format A data.frame with 10 rows and 4 columns:
#' \describe{
#'   \item{Method}{Method name}
#'   \item{CV3_AUROC}{AUROC in the most stringent CV3 setting}
#'   \item{Reference}{Source of the value}
#'   \item{Interpretability}{Whether the method is interpretable}
#' }
"benchmark_methods"
