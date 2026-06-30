# ============================================================ #
# Single mode
# ============================================================ #

#' Run DEGgo single-comparison differential expression analysis
#'
#' Internal helper that runs one differential expression analysis and returns
#' annotated results without generating plots, GO enrichment or exports.
#'
#' @param counts Integer count matrix.
#' @param metadata Sample metadata.
#' @param method Differential expression method.
#' @param design_formula Design formula.
#' @param contrast Contrast vector.
#' @param orgdb Annotation database.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param logfc_cutoff Absolute log2 fold-change cutoff.
#' @param log DEGgo logger.
#'
#' @return DEGgo result object for one comparison.
#'
#' @keywords internal
#' @noRd
.run_deggo_single <- function(
    counts,
    metadata,
    method,
    design_formula,
    contrast,
    orgdb,
    padj_cutoff,
    logfc_cutoff,
    log
) {

  log("[7/11] Running differential expression analysis", type = "step")

  clean_deg_tables = TRUE

  de_results <- run_de(
    counts = counts,
    metadata = metadata,
    method = method,
    design_formula = design_formula,
    contrast = contrast
  )

  dds <- de_results$dds %||% de_results$object %||% NULL

  if (method == "DESeq2" && !is.null(dds)) {
    dds <- .annotate_dds(
      dds = dds,
      orgdb = orgdb
    )
  }

  log("[8/11] Annotating DE results", type = "step")

  processed <- .deggo_annotate_result_table(
    res_df = de_results$res_df,
    orgdb = orgdb,
    keytype = "ENSEMBL",
    gene_col = "ENSEMBL",
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff
  )

  res_df <- processed$res_df
  sig_deg <- processed$sig_deg

  de_results$dds <- dds
  de_results$res_df <- res_df
  de_results$metadata <- metadata

  de_results$results <- list(
    single_comparison = res_df
  )

  de_results$sig_deg <- list(
    single_comparison = sig_deg
  )

  de_results$sig_deg_clean <- list(
    single_comparison = sig_deg
  )

  de_results$summary <- .deggo_build_summary(
    results = de_results$results,
    sig_deg = de_results$sig_deg,
    logfc_cutoff = logfc_cutoff
  )

  if (method == "DESeq2" && !is.null(dds)) {

    de_results$gene_annotation <- data.frame(
      ENSEMBL = rownames(dds),
      SYMBOL = SummarizedExperiment::rowData(dds)$SYMBOL,
      ENTREZID = SummarizedExperiment::rowData(dds)$ENTREZID,
      GENENAME = SummarizedExperiment::rowData(dds)$GENENAME,
      stringsAsFactors = FALSE
    )
  }

  de_results
}


# ============================================================ #
# Pairwise mode
# ============================================================ #

#' Run DEGgo pairwise differential expression analysis
#'
#' Internal helper that runs DESeq2 pairwise contrasts and returns annotated
#' result tables without generating plots, GO enrichment or exports.
#'
#' @param counts Integer count matrix.
#' @param metadata Sample metadata.
#' @param orgdb Annotation database.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param logfc_cutoff Absolute log2 fold-change cutoff.
#' @param pairwise_group_cols Metadata columns defining pairwise groups.
#' @param pairwise_contrast_col Column storing generated contrast labels.
#' @param pairwise_contrasts Optional named pairwise contrast list.
#' @param pairwise_mode Pairwise comparison mode.
#' @param log DEGgo logger.
#'
#' @return DEGgo result object containing all pairwise comparisons.
#'
#' @keywords internal
#' @noRd
.run_deggo_pairwise <- function(
    counts,
    metadata,
    orgdb,
    padj_cutoff,
    logfc_cutoff,
    pairwise_group_cols,
    pairwise_contrast_col,
    pairwise_contrasts,
    pairwise_mode,
    log
) {

  log("[7/11] Running pairwise DESeq2 contrasts", type = "step")

  clean_deg_tables = TRUE

  de_results <- run_deseq2_pairwise(
    counts = counts,
    metadata = metadata,
    group_cols = pairwise_group_cols,
    contrast_col = pairwise_contrast_col,
    contrast_list = pairwise_contrasts,
    pairwise_mode = pairwise_mode
  )

  de_results$dds <- .annotate_dds(
    dds = de_results$dds,
    orgdb = orgdb
  )

  de_results$metadata <- metadata

  log("[8/11] Annotating pairwise DE tables", type = "step")

  results_annotated <- list()
  sig_deg <- list()

  for (nm in names(de_results$results)) {

    processed <- .deggo_annotate_result_table(
      res_df = de_results$results[[nm]],
      orgdb = orgdb,
      keytype = "ENSEMBL",
      gene_col = "ENSEMBL",
      padj_cutoff = padj_cutoff,
      logfc_cutoff = logfc_cutoff
    )

    results_annotated[[nm]] <- processed$res_df
    sig_deg[[nm]] <- processed$sig_deg
  }

  de_results$results <- results_annotated
  de_results$sig_deg <- sig_deg
  de_results$sig_deg_clean <- sig_deg

  log("[9/11] Summarizing pairwise results", type = "step")

  de_results$summary <- .deggo_build_summary(
    results = de_results$results,
    sig_deg = de_results$sig_deg,
    logfc_cutoff = logfc_cutoff
  )

  if (!is.null(de_results$dds)) {

    dds <- de_results$dds

    de_results$gene_annotation <- data.frame(
      ENSEMBL = rownames(dds),
      SYMBOL = SummarizedExperiment::rowData(dds)$SYMBOL,
      ENTREZID = SummarizedExperiment::rowData(dds)$ENTREZID,
      GENENAME = SummarizedExperiment::rowData(dds)$GENENAME,
      stringsAsFactors = FALSE
    )
  }

  de_results
}
