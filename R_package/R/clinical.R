
# ═══════════════════════════════════════════════════════════════
# Candidate Prioritization (deprecated clinical alias retained)
# ═══════════════════════════════════════════════════════════════

#' Predict trial response priority for a paralog-SL candidate
#'
#' \strong{Deprecated.} This name overstated the scope of the score: the
#' composite is an in vitro, association-based \emph{prioritization}
#' heuristic, not a predictor of clinical trial response. Please use
#' \code{\link{rank_dependency_candidates}} instead. This wrapper is
#' retained for backward compatibility and delegates to the new function.
#'
#' @param driver Character, driver gene name (e.g., "ARID1A")
#' @param paralog Character, paralog gene name (e.g., "ARID1B")
#' @param dd Numeric, absolute Delta Dependency value
#' @param ti Numeric, dependency window score (DWS; default NA)
#' @param msi_status Character, "MSS" or "MSI-H" (default "MSS")
#' @param mutation_type Character, "truncating" or "missense" (default "truncating")
#' @param selectivity Numeric, selectivity score (default 0)
#' @return A list with priority_score, tier, and rationale
#' @export
#' @examples
#' # Manuscript Table S5 values for the leading candidate ARID1A->ARID1B
#' rank_dependency_candidates("ARID1A", "ARID1B", dd = 0.270, ti = 2.82,
#'                            msi_status = "MSS", mutation_type = "truncating")
predict_trial_response <- function(driver, paralog, dd, ti = NA,
                                    msi_status = "MSS",
                                    mutation_type = "truncating",
                                    selectivity = 0) {
  .Deprecated("rank_dependency_candidates")
  rank_dependency_candidates(driver, paralog, dd, ti = ti,
                             msi_status = msi_status,
                             mutation_type = mutation_type,
                             selectivity = selectivity)
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
