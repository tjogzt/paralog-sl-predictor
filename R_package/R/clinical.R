
# ═══════════════════════════════════════════════════════════════
# Clinical Decision Support Functions
# ═══════════════════════════════════════════════════════════════

#' Predict trial response priority for a paralog-SL candidate
#'
#' Combines DD magnitude, therapeutic index, MSI status, and mutation type
#' into a composite clinical priority score for trial design.
#'
#' @param driver Character, driver gene name (e.g., "ARID1A")
#' @param paralog Character, paralog gene name (e.g., "ARID1B")
#' @param dd Numeric, absolute Delta Dependency value
#' @param ti Numeric, therapeutic index (default NA, computed if possible)
#' @param msi_status Character, "MSS" or "MSI-H" (default "MSS")
#' @param mutation_type Character, "truncating" or "missense" (default "truncating")
#' @param selectivity Numeric, selectivity score (default 0)
#' @return A list with priority_score, tier, and rationale
#' @export
#' @examples
#' # Manuscript Table S6 values for the leading candidate ARID1A->ARID1B
#' predict_trial_response("ARID1A", "ARID1B", dd = 0.270, ti = 2.82,
#'                        msi_status = "MSS", mutation_type = "truncating")
predict_trial_response <- function(driver, paralog, dd, ti = NA,
                                    msi_status = "MSS",
                                    mutation_type = "truncating",
                                    selectivity = 0) {
  # Base score from DD magnitude
  score <- abs(dd) * 10  # scale to 0-3 range

  # TI bonus
  if (!is.na(ti) && ti > 1.0) score <- score + log2(ti)

  # Selectivity bonus
  if (selectivity > 0.15) score <- score + 1.0

  # MSI penalty (MSI-H shows weaker signal)
  if (toupper(msi_status) == "MSI-H") score <- score * 0.7

  # Mutation type adjustment (truncating > missense)
  if (tolower(mutation_type) == "missense") score <- score * 0.8

  # Classify tier
  tier <- if (score >= 4) "HIGH_PRIORITY" else
          if (score >= 2) "MODERATE_PRIORITY" else "LOW_PRIORITY"

  rationale <- sprintf(
    "%s->%s: |DD|=%.3f, TI=%s, %s, %s mutation. Score=%.2f (%s).",
    driver, paralog, abs(dd),
    if (is.na(ti)) "NA" else sprintf("%.2f", ti),
    msi_status, mutation_type, score, tier)

  list(
    driver = driver,
    paralog = paralog,
    priority_score = round(score, 3),
    tier = tier,
    rationale = rationale,
    biomarkers = list(msi = msi_status, mutation = mutation_type,
                      dd = abs(dd), ti = ti, selectivity = selectivity)
  )
}

#' Visualize dependency shift for a driver-paralog pair
#'
#' Creates a boxplot comparing Chronos gene-effect scores of the paralog
#' between driver-mutant and wild-type cell lines.
#'
#' @param dependency Numeric matrix of Chronos gene-effect scores (cell lines x genes)
#' @param driver_gene Character, driver gene name
#' @param paralog_gene Character, paralog gene name
#' @param mut_lines Character vector of mutant cell line IDs
#' @param wt_lines Character vector of wild-type cell line IDs
#' @return A ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' p <- visualize_dependency_shift(dep, "ARID1A", "ARID1B", mut_ids, wt_ids)
#' print(p)
#' }
visualize_dependency_shift <- function(dependency, driver_gene, paralog_gene,
                                        mut_lines, wt_lines) {
  mut_lines <- intersect(mut_lines, rownames(dependency))
  wt_lines  <- intersect(wt_lines,  rownames(dependency))

  dep_mut <- dependency[mut_lines, paralog_gene]
  dep_wt  <- dependency[wt_lines,  paralog_gene]

  dd_val <- mean(dep_wt, na.rm = TRUE) - mean(dep_mut, na.rm = TRUE)
  p_val  <- tryCatch(t.test(dep_mut, dep_wt)$p.value, error = function(e) NA)

  df <- data.frame(
    score = c(dep_mut, dep_wt),
    group = factor(c(rep(paste0(driver_gene, "\nmutant"), length(dep_mut)),
                     rep(paste0(driver_gene, "\nwild-type"), length(dep_wt))),
                   levels = c(paste0(driver_gene, "\nwild-type"),
                              paste0(driver_gene, "\nmutant"))))

  p <- ggplot2::ggplot(df, ggplot2::aes(group, score, fill = group)) +
    ggplot2::geom_boxplot(alpha = 0.4, outlier.shape = NA, linewidth = 0.4) +
    ggplot2::geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
    ggplot2::scale_fill_manual(values = c("#2171B5", "#CB181D")) +
    ggplot2::annotate("text", x = 1.5, y = max(df$score, na.rm = TRUE),
                      label = sprintf("DD = %.3f\np = %.2e", dd_val,
                                      ifelse(is.na(p_val), 1, p_val)),
                      size = 3, hjust = 0.5) +
    ggplot2::labs(x = NULL,
                  y = paste0(paralog_gene, " gene-effect score"),
                  title = paste0(driver_gene, " -> ", paralog_gene, " dependency shift")) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "none",
                   plot.background = ggplot2::element_rect(fill = "white", color = NA))
  return(p)
}
