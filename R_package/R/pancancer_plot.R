#' Plot pan-cancer DD AUROC summary
#'
#' Creates a horizontal bar chart showing Delta Dependency (DD) predictive
#' performance (AUROC) across all solid tumor types, colored by mechanistic
#' classification (TSG-driven vs oncogene-driven).
#'
#' @param summary_df A data.frame with columns: cancer, dd_auroc, mechanism
#'   (optional: "TSG", "Oncogene", "Mixed", "Insufficient").
#' @param min_pairs Minimum number of paralog pairs to include (default 5).
#' @return A ggplot2 object suitable for direct use or further customization.
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip theme_bw labs
#'   scale_fill_manual geom_hline
#' @export
plot_pancancer_summary <- function(summary_df, min_pairs = 5) {
  df <- summary_df[summary_df$n_pairs >= min_pairs, , drop = FALSE]
  if (nrow(df) == 0) stop("No cancer types with at least ", min_pairs, " pairs")

  # Add mechanism labels if not present (manuscript classification:
  # oncogene-driven = {Melanoma, NSCLC, Pancreatic}; TSG-driven = the other
  # evaluable lineages incl. SCLC; anything else = "Mixed")
  if (!"mechanism" %in% colnames(df)) {
    tsg_driven <- c("Ovarian", "Endometrial", "Cervical", "Breast",
                    "Colorectal", "Esophagogastric", "Biliary Tract",
                    "HNSCC", "Bladder Urothelial", "Hepatocellular",
                    "Glioma", "Renal Cell", "Mesothelioma", "SCLC")
    onc_driven <- c("Melanoma", "NSCLC", "Pancreatic")
    df$mechanism <- ifelse(df$cancer %in% tsg_driven, "TSG",
                    ifelse(df$cancer %in% onc_driven, "Oncogene", "Mixed"))
  }

  # Fill NAs with "Mixed"
  df$mechanism[is.na(df$mechanism)] <- "Mixed"

  # Use pre-computed AUROC if available
  if (!"dd_auroc" %in% colnames(df)) {
    stop("summary_df must contain a 'dd_auroc' column")
  }

  df <- df[order(df$dd_auroc), , drop = FALSE]
  df$cancer <- factor(df$cancer, levels = df$cancer)

  cols <- c("TSG" = "#2B4C7E", "Oncogene" = "#C0362C",
            "Mixed" = "#9B7FA6", "Insufficient" = "#808080")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = cancer, y = dd_auroc,
                                         fill = mechanism)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                        color = "grey40", alpha = 0.6) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::scale_fill_manual(values = cols, name = "Mechanism") +
    ggplot2::labs(
      x = "",
      y = "AUROC (Delta Dependency)",
      title = "Pan-cancer Paralog-SL Predictive Performance",
      subtitle = paste("Across", nrow(df), "solid tumor types"),
      caption = "Dashed line: random classifier (0.5)"
    )
  return(p)
}

#' Plot therapeutic window summary
#'
#' Visualizes the dependency window score (DWS; reported as therapeutic
#' index, TI, in earlier versions) and selectivity for paralog-SL candidates
#' ranked by preclinical prioritization score.
#'
#' @param tw_df A data.frame with columns: driver_gene (or driver),
#'   paralog_gene (or paralog), mean_TI (or mean_ti), classification -- e.g.
#'   the shipped `therapeutic_window_summary` dataset.
#' @return A ggplot2 object.
#' @export
plot_therapeutic_window <- function(tw_df) {
  if (nrow(tw_df) == 0) stop("Empty therapeutic window data.frame")

  # Accept both the canonical artifact schema (driver/paralog/mean_ti) and
  # the legacy package schema (driver_gene/paralog_gene/mean_TI)
  n <- names(tw_df)
  if (!"mean_TI" %in% n && "mean_ti" %in% n) tw_df$mean_TI <- tw_df$mean_ti
  if (!"driver_gene" %in% n && "driver" %in% n) tw_df$driver_gene <- tw_df$driver
  if (!"paralog_gene" %in% n && "paralog" %in% n) tw_df$paralog_gene <- tw_df$paralog
  if (is.null(tw_df$mean_TI)) stop("tw_df must contain a mean_TI (or mean_ti) column")

  tw_df <- tw_df[order(tw_df$mean_TI), , drop = FALSE]
  tw_df$label <- paste(tw_df$driver_gene, tw_df$paralog_gene, sep = "\u2192")
  tw_df$label <- factor(tw_df$label, levels = tw_df$label)

  tier_cols <- c("HIGH_SELECTIVITY" = "#1B7837", "MODERATE" = "#5AAE61",
                 "LOW_SELECTIVITY" = "#A6DBA0", "PAN_ESSENTIAL" = "#C0362C",
                 "INSUFFICIENT_DATA" = "#808080")

  p <- ggplot2::ggplot(tw_df, ggplot2::aes(x = label, y = mean_TI,
                                            fill = classification)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 1.0, linetype = "dashed",
                        color = "grey40", alpha = 0.6) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::scale_fill_manual(values = tier_cols, name = "Safety Tier") +
    ggplot2::labs(
      x = "",
      y = "Mean dependency window score (DWS/TI)",
      title = "Therapeutic Window for Paralog-SL Candidates",
      subtitle = "DWS > 1 indicates selective essentiality in mutant context",
      caption = "DWS = |DD| / max(|mean gene-effect|, pan-essential fraction)"
    )
  return(p)
}
