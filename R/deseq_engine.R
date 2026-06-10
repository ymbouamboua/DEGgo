
# =========================================================
# RUN DESEQ2
# =========================================================
#' Run DESeq2 differential expression analysis
#'
#' @param counts Raw count matrix with genes as rows and samples as columns.
#' @param metadata Sample metadata data frame. Row names must match sample names.
#' @param design_formula Formula specifying the experimental design.
#'
#' @return A list with the DESeq2 object and result table.
#'
#' @keywords internal
#' @noRd
run_deseq <- function(
    counts,
    metadata,
    design_formula = ~ condition
) {

  log_msg("Running DESeq2 analysis...")

  metadata$condition <- as.factor(
    metadata$condition
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(

    countData = counts,

    colData = metadata,

    design = design_formula
  )

  dds <- DESeq2::DESeq(dds)

  res <- DESeq2::results(dds)

  res_df <- as.data.frame(res)

  res_df$ENSEMBL <- rownames(res_df)

  return(
    list(
      dds = dds,
      res_df = res_df
    )
  )
}


# =========================================================
# PROCESS DEG RESULTS
# =========================================================
#' Process differential expression results
#'
#' @param res_df Differential expression results data frame.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param logfc_cutoff Absolute log2 fold-change cutoff.
#'
#' @return A list with processed results and significant DEGs.
#'
#' @keywords internal
#' @noRd
process_deg_results <- function(
    res_df,
    padj_cutoff = 0.05,
    logfc_cutoff = 1
) {

  log_msg("Processing DEG results...")

  res_df <- res_df |>

    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(padj)
    )

  res_df$padj[
    res_df$padj == 0
  ] <- 1e-300


  # -------------------------------------------------------
  # Significance Classification
  # -------------------------------------------------------

  res_df$significance <- dplyr::case_when(

    res_df$padj < padj_cutoff &
      res_df$log2FoldChange > logfc_cutoff
    ~ "Upregulated",

    res_df$padj < padj_cutoff &
      res_df$log2FoldChange < -logfc_cutoff
    ~ "Downregulated",

    TRUE ~ "Not Significant"
  )


  # -------------------------------------------------------
  # Volcano Coordinates
  # -------------------------------------------------------

  res_df$negLog10Padj <- -log10(
    res_df$padj
  )


  # -------------------------------------------------------
  # Significant DEGs
  # -------------------------------------------------------

  sig_deg <- subset(

    res_df,

    padj < padj_cutoff &
      abs(log2FoldChange) > logfc_cutoff
  )

  return(
    list(
      res_df = res_df,
      sig_deg = sig_deg
    )
  )
}
