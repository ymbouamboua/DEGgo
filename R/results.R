#' Process differential expression results
#'
#' Classifies genes according to statistical significance and
#' log2 fold-change thresholds, computes \code{-log10(padj)}
#' values for visualization, and extracts significant DEGs.
#'
#' Genes are assigned to one of three categories:
#' \itemize{
#'   \item \code{"Upregulated"}
#'   \item \code{"Downregulated"}
#'   \item \code{"Not Significant"}
#' }
#'
#' A gene is considered significant when:
#' \itemize{
#'   \item adjusted p-value < \code{padj_cutoff}
#'   \item absolute log2 fold-change > \code{logfc_cutoff}
#' }
#'
#' @param res_df Differential expression result table containing at least:
#'   \itemize{
#'     \item \code{padj}
#'     \item \code{log2FoldChange}
#'   }
#' @param padj_cutoff Adjusted p-value significance threshold.
#' @param logfc_cutoff Absolute log2 fold-change threshold.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{res_df}: Updated differential expression table with:
#'   \itemize{
#'     \item \code{negLog10Padj}
#'     \item \code{significance}
#'   }
#'   \item \code{sig_deg}: Subset of significant differentially
#'   expressed genes.
#' }
#'
#' @details
#' This helper is used internally by \code{run_deggo()} after
#' differential expression analysis to standardize downstream
#' visualization, reporting, volcano plots, and enrichment analyses.
#'
#' @keywords internal
#' @noRd
process_deg_results <- function(
    res_df,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25
) {

  res_df$negLog10Padj <- -log10(res_df$padj + 1e-300)

  res_df$significance <- "Ns"

  res_df$significance[
    !is.na(res_df$padj) &
      res_df$padj < padj_cutoff &
      res_df$log2FoldChange > logfc_cutoff
  ] <- "Up"

  res_df$significance[
    !is.na(res_df$padj) &
      res_df$padj < padj_cutoff &
      res_df$log2FoldChange < -logfc_cutoff
  ] <- "Down"

  sig_deg <- res_df[
    res_df$significance != "Ns",
    ,
    drop = FALSE
  ]

  list(
    res_df = res_df,
    sig_deg = sig_deg
  )
}


