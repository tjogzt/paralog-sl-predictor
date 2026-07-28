
# ═══════════════════════════════════════════════════════════════
# Candidate Prioritization Functions
# ═══════════════════════════════════════════════════════════════

#' Rank a paralog-SL candidate for experimental follow-up
#'
#' Combines DD magnitude, dependency window score (DWS), MSI status, and
#' mutation type into a composite \emph{prioritization} score used to rank
#' candidate paralog dependencies for experimental follow-up. The score is
#' an in vitro, association-based heuristic: it does not predict clinical
#' trial response or a therapeutic window.
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
rank_dependency_candidates <- function(driver, paralog, dd, ti = NA,
                                       msi_status = "MSS",
                                       mutation_type = "truncating",
                                       selectivity = 0) {
  # Base score from DD magnitude
  score <- abs(dd) * 10  # scale to 0-3 range

  # DWS bonus
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
    "%s->%s: |DD|=%.3f, DWS=%s, %s, %s mutation. Score=%.2f (%s).",
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
