#' paralogSL: Paralog-Based Synthetic Lethality Prediction via Delta Dependency
#'
#' A lightweight R package for de novo synthetic lethality discovery through
#' paralog compensation analysis using DepMap CRISPR dependency data.
#'
#' @keywords internal
"_PACKAGE"

# Shipped dataset referenced inside get_benchmark_table(); LazyData makes it
# visible at runtime. Remaining names are ggplot2 NSE column references in
# plot functions. Both silence R CMD check "no visible binding" notes.
utils::globalVariables(c("benchmark_methods", "CV3_AUROC", "Interpretability",
                         "Method", "cancer", "classification", "dd_auroc",
                         "group", "label", "mean_TI", "mechanism", "score",
                         "x", "y"))

#' Load DepMap CRISPR dependency data
#'
#' Reads a DepMap CRISPRGeneEffect.csv file and returns a numeric matrix
#' with cell lines as rows and genes as columns.
#'
#' @param path Path to CRISPRGeneEffect.csv
#' @return A numeric matrix (cell lines x genes)
#' @export
#' @importFrom data.table fread
load_dependency <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  dt <- data.table::fread(path, data.table = FALSE)
  rownames(dt) <- dt[, 1]
  dt <- dt[, -1, drop = FALSE]
  colnames(dt) <- gsub(" \\(\\d+\\)", "", colnames(dt))
  as.matrix(dt)
}

#' Load DepMap gene expression data
#'
#' @param path Path to OmicsExpressionProteinCodingGenesTPMLogp1.csv
#' @return A numeric matrix (cell lines x genes)
#' @export
load_expression <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  dt <- data.table::fread(path, data.table = FALSE)
  rownames(dt) <- dt[, 1]
  dt <- dt[, -1, drop = FALSE]
  colnames(dt) <- gsub(" \\(\\d+\\)", "", colnames(dt))
  as.matrix(dt)
}

#' Load DepMap cell line annotations
#'
#' @param path Path to Model.csv
#' @return A data.frame with cell line annotations
#' @export
load_models <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- data.table::fread(path, data.table = FALSE)
  colnames(df)[colnames(df) == "ModelID"] <- "DepMap_ID"
  df
}

#' Build binary mutation matrix
#'
#' @param mutations_df A data.frame with columns DepMap_ID and Gene
#' @param cell_lines Character vector of cell line IDs
#' @param genes Character vector of gene names
#' @return A binary matrix (cell lines x genes): 1 = mutated, 0 = wild-type
#' @export
build_mutation_matrix <- function(mutations_df, cell_lines, genes) {
  mat <- matrix(0, nrow = length(cell_lines), ncol = length(genes),
                dimnames = list(cell_lines, genes))
  sub <- mutations_df[mutations_df$DepMap_ID %in% cell_lines &
                      mutations_df$Gene %in% genes, , drop = FALSE]
  if (nrow(sub) == 0) return(mat)
  for (i in seq_len(nrow(sub))) {
    cl <- sub$DepMap_ID[i]; g <- sub$Gene[i]
    if (cl %in% cell_lines && g %in% genes) mat[cl, g] <- 1
  }
  mat
}

#' Compute Delta Dependency for a driver-paralog pair
#'
#' Delta Dependency DD(D, P, c) = mean(G_P | D wild-type) - mean(G_P | D mutant),
#' where G is the Chronos gene-effect score. Positive values indicate greater
#' paralog dependency in the mutant context (consistent with paralog compensation).
#' This matches Equation 1 of the companion manuscript and pcs.py in the
#' Python pipeline. The returned p_value is a Welch t-test on the dependency
#' scores themselves (distinct from the expression-based test in pcs.py,
#' which is exported there as expr_p_value and dd_p_value respectively).
#'
#' @param dependency Numeric matrix (cell lines x genes) of Chronos gene-effect scores
#' @param driver_gene Character, the driver gene name
#' @param paralog_gene Character, the paralog gene name
#' @param mut_lines Character vector of cell lines with driver mutation
#' @param wt_lines Character vector of wild-type cell lines
#' @param min_samples Minimum number of samples per group (default 3)
#' @return A list with DD, p_value, cohens_d, and group sizes
#' @export
compute_dd <- function(dependency, driver_gene, paralog_gene,
                        mut_lines, wt_lines, min_samples = 3) {
  mut_lines <- intersect(mut_lines, rownames(dependency))
  wt_lines  <- intersect(wt_lines,  rownames(dependency))

  if (length(mut_lines) < min_samples || length(wt_lines) < min_samples) {
    return(list(DD = NA, p_value = NA, cohens_d = NA,
                n_mut = length(mut_lines), n_wt = length(wt_lines)))
  }

  dep_mut <- dependency[mut_lines, paralog_gene]
  dep_wt  <- dependency[wt_lines,  paralog_gene]

  dd <- mean(dep_wt, na.rm = TRUE) - mean(dep_mut, na.rm = TRUE)
  # Welch t-test on dependency scores (t.test defaults to unequal variances)
  t_res <- tryCatch(t.test(dep_mut, dep_wt), error = function(e) NULL)
  p_val <- if (!is.null(t_res)) t_res$p.value else NA

  pooled_sd <- sqrt((var(dep_mut, na.rm = TRUE) + var(dep_wt, na.rm = TRUE)) / 2)
  cohens_d <- if (pooled_sd > 0) dd / pooled_sd else 0

  list(DD = dd, p_value = p_val, cohens_d = cohens_d,
       n_mut = length(mut_lines), n_wt = length(wt_lines))
}

#' Compute Paralog Compensation Score
#'
#' PCS = delta_expression(paralog, MUT vs WT) * necessity(paralog)
#' where necessity = -mean(gene_effect_paralog) globally.
#'
#' @param dependency Numeric matrix of Chronos gene-effect scores
#' @param expression Numeric matrix of expression values
#' @param driver_gene Character, the driver gene
#' @param paralog_gene Character, the paralog gene
#' @param mut_lines Mutant cell line IDs
#' @param wt_lines Wildtype cell line IDs
#' @return A list with PCS, delta_expression, necessity
#' @export
compute_pcs <- function(dependency, expression, driver_gene, paralog_gene,
                         mut_lines, wt_lines) {
  mut_lines <- intersect(mut_lines, intersect(rownames(dependency), rownames(expression)))
  wt_lines  <- intersect(wt_lines,  intersect(rownames(dependency), rownames(expression)))

  if (length(mut_lines) < 2 || length(wt_lines) < 2) {
    return(list(PCS = 0, delta_expression = NA, necessity = NA))
  }

  expr_mut <- mean(expression[mut_lines, paralog_gene], na.rm = TRUE)
  expr_wt  <- mean(expression[wt_lines,  paralog_gene], na.rm = TRUE)
  delta_expr <- expr_mut - expr_wt

  necessity <- -mean(dependency[, paralog_gene], na.rm = TRUE)
  pcs <- max(delta_expr, 0) * max(necessity, 0)

  list(PCS = pcs, delta_expression = delta_expr, necessity = necessity)
}

#' Run full paralog-SL analysis for a set of driver genes
#'
#' Computes DD and PCS for all paralog pairs involving specified driver genes.
#'
#' @param dependency Chronos gene-effect dependency matrix
#' @param expression Expression matrix
#' @param mutations_df Mutation annotations data.frame
#' @param paralogs_df Paralog pairs data.frame with columns gene_A, gene_B
#' @param driver_genes Character vector of driver genes to analyze
#' @param cell_lines Character vector of cell line IDs to use
#' @param min_samples Minimum samples per group (default 3)
#' @return A data.frame with DD, PCS, and statistics for each pair
#' @export
run_paralog_analysis <- function(dependency, expression, mutations_df,
                                   paralogs_df, driver_genes, cell_lines,
                                   min_samples = 3) {
  results <- list()

  for (driver in driver_genes) {
    if (!driver %in% colnames(dependency)) next

    driver_paralogs <- paralogs_df$gene_B[paralogs_df$gene_A == driver]
    driver_paralogs <- intersect(driver_paralogs, colnames(dependency))
    driver_paralogs <- intersect(driver_paralogs, colnames(expression))
    if (length(driver_paralogs) == 0) next

    mut_mat <- build_mutation_matrix(mutations_df, cell_lines, driver)
    mut_lines <- names(which(mut_mat[, driver] == 1))
    wt_lines  <- names(which(mut_mat[, driver] == 0))

    if (length(mut_lines) < min_samples || length(wt_lines) < min_samples) next

    for (para in driver_paralogs) {
      dd_res <- compute_dd(dependency, driver, para, mut_lines, wt_lines, min_samples)
      pcs_res <- compute_pcs(dependency, expression, driver, para, mut_lines, wt_lines)

      results[[length(results) + 1]] <- data.frame(
        driver_gene = driver,
        paralog_gene = para,
        DD = dd_res$DD,
        p_value = dd_res$p_value,
        cohens_d = dd_res$cohens_d,
        PCS = pcs_res$PCS,
        delta_expression = pcs_res$delta_expression,
        necessity = pcs_res$necessity,
        n_mut = dd_res$n_mut,
        n_wt = dd_res$n_wt,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) == 0) return(data.frame())
  res <- do.call(rbind, results)

  # Benjamini-Hochberg correction per driver
  res$q_value <- NA_real_
  for (drv in unique(res$driver_gene)) {
    idx <- which(res$driver_gene == drv)
    pvals <- res$p_value[idx]
    res$q_value[idx] <- p.adjust(pvals, method = "BH")
  }

  # Composite score
  res$composite_score <- with(res, {
    dd_norm <- (DD - min(DD, na.rm = TRUE)) / diff(range(DD, na.rm = TRUE))
    if (all(is.na(dd_norm))) dd_norm <- rep(0, length(DD))
    dd_norm[is.na(dd_norm)] <- 0
    0.6 * dd_norm + 0.4 * (1 - pmin(q_value, 1))
  })

  res <- res[order(-res$composite_score), , drop = FALSE]
  rownames(res) <- NULL
  res
}

#' Filter cell lines to specific cancer types
#'
#' @param models_df Cell line annotations data.frame
#' @param cancer_types Character vector of Oncotree disease names to match
#' @return A character vector of matching cell line IDs
#' @export
filter_cancer_cell_lines <- function(models_df, cancer_types) {
  disease_col <- "OncotreePrimaryDisease"
  if (!disease_col %in% colnames(models_df)) {
    warning("OncotreePrimaryDisease column not found")
    return(character(0))
  }

  pattern <- paste(cancer_types, collapse = "|")
  idx <- grepl(pattern, models_df[[disease_col]], ignore.case = TRUE)
  models_df$DepMap_ID[idx]
}

#' Compute AUROC for DD prediction against known SL pairs
#'
#' Ranks pairs by |DD| (absolute Delta Dependency), matching the validation
#' convention used throughout the companion manuscript and Python pipeline
#' ("DD alone, using only |DD|"). With the manuscript sign convention
#' (DD = mean WT − mean MUT), a large positive DD indicates compensatory
#' dependency in mutant lines; the absolute value additionally credits pairs
#' where dependency shifts strongly in either direction.
#'
#' @param results_df Results from run_paralog_analysis()
#' @param known_pairs data.frame with columns gene_A, gene_B; optionally a
#'   `tier` column (evidence tiers A/B/C/Comparator, see `known_sl_pairs`)
#' @param tiers character vector of evidence tiers treated as gold-standard
#'   positives when `known_pairs` carries a `tier` column. Defaults to the
#'   manuscript's primary external benchmark (Tier A + Tier B: direct
#'   dual-perturbation or genotype-conditional genetic evidence). Set to
#'   `NULL` to use every row of `known_pairs` (the pre-v1.1.0 behavior).
#' @return The AUROC value
#' @importFrom pROC roc
#' @export
compute_auroc <- function(results_df, known_pairs, tiers = c("A", "B")) {
  if (!is.null(tiers) && "tier" %in% names(known_pairs)) {
    known_pairs <- known_pairs[known_pairs$tier %in% tiers, ]
  }
  is_known <- (paste(results_df$driver_gene, results_df$paralog_gene) %in%
               paste(known_pairs$gene_A, known_pairs$gene_B)) |
              (paste(results_df$paralog_gene, results_df$driver_gene) %in%
               paste(known_pairs$gene_A, known_pairs$gene_B))

  if (sum(is_known) < 2 || sum(!is_known) < 2) return(NA_real_)

  roc_obj <- pROC::roc(is_known, abs(results_df$DD), quiet = TRUE)
  as.numeric(roc_obj$auc)
}

#' Benchmark DD against published method performance
#'
#' Returns the shipped `benchmark_methods` dataset: CV3 AUROC values for 8
#' published SL prediction methods from Feng et al. (2024), plus this-study
#' DD rows (full lineage-level frame, AUROC 0.676; DD + ID >= 0.3
#' high-identity subset, AUROC 1.000, anecdotal). Published values are
#' contextual reference points from a general SL gene-pair universe, not a
#' head-to-head benchmark. The dataset is regenerated from the canonical
#' Table2_Benchmark.tsv artifact of the paralog-sl-predictor pipeline.
#'
#' @return A data.frame with Method, CV3_AUROC, Reference, and
#'   Interpretability columns
#' @export
get_benchmark_table <- function() {
  benchmark_methods
}

#' Plot benchmark comparison
#'
#' @param benchmark_df Result from get_benchmark_table()
#' @return A ggplot2 object
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip theme_bw labs
#' @export
plot_benchmark <- function(benchmark_df = get_benchmark_table()) {
  benchmark_df$Method <- factor(benchmark_df$Method, levels = rev(benchmark_df$Method))

  ggplot2::ggplot(benchmark_df, ggplot2::aes(x = Method, y = CV3_AUROC,
        fill = Interpretability)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.4) +
    ggplot2::scale_fill_manual(values = c("High" = "#C0362C", "Low" = "#9B7FA6")) +
    ggplot2::labs(x = "", y = "AUROC (CV3 / paralog-SL)",
                  title = "DD vs Published SL Prediction Methods") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", CV3_AUROC)),
                       hjust = -0.2, size = 3.5)
}

#' Plot paralog co-variation from CPTAC proteomics data
#'
#' @param prot_df Protein abundance data.frame (genes x samples)
#' @param gene_a First gene name
#' @param gene_b Second gene name
#' @return A ggplot2 object
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_bw labs
#' @importFrom stats cor.test
#' @export
plot_protein_correlation <- function(prot_df, gene_a, gene_b) {
  pa <- as.numeric(prot_df[gene_a, ])
  pb <- as.numeric(prot_df[gene_b, ])
  valid <- !is.na(pa) & !is.na(pb)
  df <- data.frame(x = pa[valid], y = pb[valid])

  r_val <- cor.test(df$x, df$y)
  label <- sprintf("r = %.3f\np = %.2e", r_val$estimate, r_val$p.value)

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(alpha = 0.4, color = "#2B4C7E", size = 1.5) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#C0362C", linewidth = 1) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(x = gene_a, y = gene_b,
                  title = paste(gene_a, "vs", gene_b, "(CPTAC proteomics)")) +
ggplot2::annotate("text", x = min(df$x) + 0.6 * diff(range(df$x)),
                       y = max(df$y) - 0.1 * diff(range(df$y)),
                       label = label, hjust = 0, size = 3.5, color = "#C0362C")
}

#' Classify MSI/dMMR status from MMR gene mutations
#'
#' Classifies cell lines as MSI-H (dMMR) or MSS (pMMR) based on
#' damaging mutations in MMR genes (MLH1, MSH2, MSH6, PMS2) and POLE.
#' This is the manuscript's "MMR/POLE mutation proxy" definition: POLE
#' ultramutation is pooled with dMMR here for simplicity, whereas the
#' companion Python pipeline additionally reports POLE-mutant lines as a
#' separate stratum. Prefer curated DepMap/CCLE MSI annotations when
#' available.
#'
#' @param mutations_df Mutation data.frame with columns DepMap_ID, Gene
#' @param cell_lines Character vector of cell line IDs
#' @return A named character vector: "MSI_H" or "MSS"
#' @export
classify_msi_status <- function(mutations_df, cell_lines) {
  mmr_genes <- c("MLH1", "MSH2", "MSH6", "PMS2", "POLE")
  mmr_muts <- mutations_df$DepMap_ID[mutations_df$Gene %in% mmr_genes]
  status <- rep("MSS", length(cell_lines))
  names(status) <- cell_lines
  status[cell_lines %in% mmr_muts] <- "MSI_H"
  status
}

#' Compute Therapeutic Window metrics
#'
#' Therapeutic Index (TI) = |DD| / max(|mean_G|, pan_essential_fraction, 0.01).
#' Classifies paralogs into in vitro selectivity tiers. In the manuscript this
#' in vitro dependency selectivity score is reported as the dependency window
#' score (DWS); the function retains the `TI` field name for backward
#' compatibility (identical formula to the Python pipeline's
#' `therapeutic_index`).
#'
#' @param dependency Chronos gene-effect matrix
#' @param driver_gene Driver gene name
#' @param paralog_gene Paralog gene name
#' @param mut_lines Mutant cell line IDs
#' @param wt_lines Wild-type cell line IDs
#' @param all_cell_lines All cell line IDs (for pan-essentiality)
#' @param ceres_threshold Gene-effect threshold for essentiality (default -0.5)
#' @return A list with TI, selectivity, pan_essential_frac, classification
#' @export
compute_therapeutic_window <- function(dependency, driver_gene, paralog_gene,
                                        mut_lines, wt_lines, all_cell_lines,
                                        ceres_threshold = -0.5) {
  mut_vals <- dependency[intersect(mut_lines, rownames(dependency)), paralog_gene]
  wt_vals  <- dependency[intersect(wt_lines, rownames(dependency)), paralog_gene]
  all_vals <- dependency[intersect(all_cell_lines, rownames(dependency)), paralog_gene]

  if (length(mut_vals) < 3 || length(wt_vals) < 3) {
    return(list(TI = NA, selectivity = NA, pan_essential_frac = NA,
                classification = "INSUFFICIENT_DATA"))
  }

  dd <- mean(wt_vals, na.rm = TRUE) - mean(mut_vals, na.rm = TRUE)
  pan_mean <- abs(mean(all_vals, na.rm = TRUE))
  pan_essential_frac <- mean(all_vals < ceres_threshold, na.rm = TRUE)

  mut_essential_frac <- mean(mut_vals < ceres_threshold, na.rm = TRUE)
  wt_essential_frac  <- mean(wt_vals < ceres_threshold, na.rm = TRUE)
  selectivity <- mut_essential_frac - wt_essential_frac

  pan_essentiality <- max(pan_mean, pan_essential_frac, 0.01)
  ti <- abs(dd) / pan_essentiality

  classification <- if (pan_essential_frac > 0.5) {
    "PAN_ESSENTIAL"
  } else if (selectivity > 0.15 && ti > 1.0) {
    "HIGH_SELECTIVITY"
  } else if (selectivity > 0) {
    "MODERATE"
  } else {
    "LOW_SELECTIVITY"
  }

  list(TI = ti, selectivity = selectivity,
       pan_essential_frac = pan_essential_frac,
       DD = dd, classification = classification)
}

#' Stratify mutations by consequence type
#'
#' Classifies mutations as truncating (frameshift, nonsense, splice-site)
#' or missense. Returns a binary matrix per mutation type.
#'
#' @param mutations_df Mutation data.frame with columns DepMap_ID, Gene, VariantInfo
#' @param cell_lines Character vector of cell line IDs
#' @param driver_gene Driver gene name
#' @return A list with trunc_lines and miss_lines (character vectors of cell line IDs)
#' @export
stratify_by_mutation_type <- function(mutations_df, cell_lines, driver_gene) {
  sub <- mutations_df[mutations_df$DepMap_ID %in% cell_lines &
                      mutations_df$Gene == driver_gene, , drop = FALSE]

  if (nrow(sub) == 0) {
    return(list(trunc_lines = character(0), miss_lines = character(0)))
  }

  trunc_terms <- c("frameshift_variant", "stop_gained", "nonsense",
                   "splice_acceptor_variant", "splice_donor_variant",
                   "start_lost", "stop_lost")

  is_trunc <- sapply(sub$VariantInfo, function(v) {
    if (is.na(v)) return(FALSE)
    any(sapply(trunc_terms, function(t) grepl(t, v, ignore.case = TRUE)))
  })

  is_miss <- grepl("missense_variant|protein_altering_variant",
                   sub$VariantInfo, ignore.case = TRUE) & !is_trunc

  list(
    trunc_lines = unique(sub$DepMap_ID[is_trunc]),
    miss_lines  = unique(sub$DepMap_ID[is_miss])
  )
}
