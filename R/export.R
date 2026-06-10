# =========================================================
# EXPORT RESULTS
# =========================================================

export_deg_results <- function(
    res_df,
    sig_deg,
    output_dir = "DEGgo_out"
) {
  
  log_msg("Exporting results...")
  
  cols <- c(
    "ENSEMBL",
    "SYMBOL",
    "ENTREZID",
    "significance",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "negLog10Padj"
  )
  
  res_df <- res_df[, cols, drop = FALSE]
  
  sig_deg <- sig_deg[
    order(sig_deg$padj),
    cols,
    drop = FALSE
  ]
  
  utils::write.table(
    res_df,
    file.path(output_dir, "All_DEGs.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  utils::write.table(
    sig_deg,
    file.path(output_dir, "Significant_DEGs.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}