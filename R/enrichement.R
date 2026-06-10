#' Run Gene Ontology enrichment analysis
#'
#' Performs GO enrichment analysis on
#' significant differentially expressed genes.
#'
#' @param sig_deg Significant DEG dataframe.
#'
#' @param ontology GO ontology category.
#' One of "BP", "MF", or "CC".
#' @param orgdb OrgDb annotation database used for gene ID conversion.
#' @param output_dir Directory for exported results.
#'
#' @return A list containing:
#' \itemize{
#'   \item GO enrichment results
#'   \item GO barplot
#'   \item GO dotplot
#' }
#'
#' @export

# =========================================================
# GO ENRICHMENT
# =========================================================

run_go_enrichment <- function(
    sig_deg,
    ontology = "BP",
    orgdb,
    output_dir = "DEGgo_out"
) {
  
  log_msg("Running GO enrichment...")
  
  if (is.null(sig_deg) || nrow(sig_deg) == 0) {
    log_msg("No significant DEGs found. Skipping GO enrichment.", type = "warn")
    return(NULL)
  }
  
  if (!"ENTREZID" %in% colnames(sig_deg)) {
    log_msg("ENTREZID column not found. Skipping GO enrichment.", type = "warn")
    return(NULL)
  }
  
  entrez <- unique(stats::na.omit(as.character(sig_deg$ENTREZID)))
  
  if (!length(entrez)) {
    log_msg("No ENTREZ IDs available for enrichment. Skipping GO enrichment.", type = "warn")
    return(NULL)
  }
  
  ego <- clusterProfiler::enrichGO(
    gene = entrez,
    OrgDb = orgdb,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  ego_df <- as.data.frame(ego)
  
  if (!nrow(ego_df)) {
    log_msg("No enriched GO terms found. Skipping GO plot.", type = "warn")
    return(list(go_results = ego_df, go_plot = NULL))
  }
  
  utils::write.csv(
    ego_df,
    file.path(output_dir, paste0("GO_", ontology, "_enrichment.csv")),
    row.names = FALSE
  )
  
  go_plot <- barplot(
    ego,
    showCategory = 20,
    title = paste0("GO enrichment: ", ontology)
  )

  
  ggplot2::ggsave(
    filename = file.path(output_dir, "GO_Barplot.png"),
    plot = go_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  go_plot <- enrichplot::dotplot(
    ego,
    showCategory = 20,
    title = paste0("GO enrichment: ", ontology)
  )
  
  ggplot2::ggsave(
    filename = file.path(output_dir, "GO_Dotplot.png"),
    plot = go_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  list(
    go_results = ego_df,
    go_object = ego,
    go_plot = go_plot
  )
}