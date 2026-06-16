#' Run Gene Ontology enrichment analysis
#'
#' Performs Gene Ontology enrichment analysis on significant differentially
#' expressed genes using Entrez identifiers.
#'
#' @param sig_deg Significant DEG data frame. Must contain an \code{ENTREZID}
#'   column.
#' @param ontology Character. GO ontology to test. One of \code{"BP"},
#'   \code{"MF"}, or \code{"CC"}. Default is \code{"BP"}.
#' @param orgdb OrgDb annotation database used by
#'   \code{clusterProfiler::enrichGO()}.
#' @param output_dir Character. Directory where enrichment tables and plots
#'   are exported. Default is \code{"DEGgo_out"}.
#' @param show_category Integer. Number of GO terms to display in plots.
#'   Default is \code{20}.
#' @param pvalue_cutoff Numeric. P-value cutoff for enrichment.
#'   Default is \code{0.05}.
#' @param qvalue_cutoff Numeric. Q-value cutoff for enrichment.
#'   Default is \code{0.2}.
#' @param p_adjust_method Character. Multiple-testing correction method.
#'   Default is \code{"BH"}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{go_results}}{GO enrichment result table as a data frame.}
#'   \item{\code{go_object}}{Original \code{enrichResult} object returned by
#'   \code{clusterProfiler::enrichGO()}.}
#'   \item{\code{bar_plot}}{GO enrichment bar plot.}
#'   \item{\code{dot_plot}}{GO enrichment dot plot.}
#'   \item{\code{ontology}}{Ontology tested.}
#'   \item{\code{n_input_genes}}{Number of unique Entrez IDs used.}
#' }
#'
#' @details
#' The function expects that Entrez identifiers have already been added to
#' the DEG table, for example using \code{map_entrez_ids()}.
#'
#' If no significant genes, no Entrez IDs, or no enriched terms are found,
#' the function returns \code{NULL} or an empty result object without failing.
#'
#' @examples
#' \dontrun{
#' orgdb <- org.Mm.eg.db::org.Mm.eg.db
#'
#' go <- run_go_enrichment(
#'   sig_deg = sig_deg,
#'   ontology = "BP",
#'   orgdb = orgdb,
#'   output_dir = "DEGgo_out"
#' )
#'
#' go$go_results
#' go$bar_plot
#' go$dot_plot
#' }
#'
#' @export
run_go_enrichment <- function(
    sig_deg,
    ontology = c("BP", "MF", "CC"),
    orgdb,
    output_dir = "DEGgo_out",
    show_category = 20,
    pvalue_cutoff = 0.05,
    qvalue_cutoff = 0.2,
    p_adjust_method = "BH"
) {

  ontology <- match.arg(ontology)

  log <- .msg(verbose = TRUE)
  log("Running GO enrichment...")

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required.", call. = FALSE)
  }

  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Package 'enrichplot' is required.", call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (is.null(sig_deg) || nrow(sig_deg) == 0) {
    log("No significant DEGs found. Skipping GO enrichment.", type = "warn")
    return(NULL)
  }

  if (!"ENTREZID" %in% colnames(sig_deg)) {
    log("ENTREZID column not found. Skipping GO enrichment.", type = "warn")
    return(NULL)
  }

  entrez <- unique(stats::na.omit(as.character(sig_deg$ENTREZID)))

  if (!length(entrez)) {
    log("No ENTREZ IDs available for enrichment. Skipping GO enrichment.",
        type = "warn")
    return(NULL)
  }

  ego <- clusterProfiler::enrichGO(
    gene = entrez,
    OrgDb = orgdb,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = p_adjust_method,
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    readable = TRUE
  )

  ego_df <- as.data.frame(ego)

  utils::write.csv(
    ego_df,
    file = file.path(output_dir, paste0("GO_", ontology, "_enrichment.csv")),
    row.names = FALSE
  )

  if (!nrow(ego_df)) {
    log("No enriched GO terms found. Skipping GO plots.", type = "warn")

    return(list(
      go_results = ego_df,
      go_object = ego,
      bar_plot = NULL,
      dot_plot = NULL,
      ontology = ontology,
      n_input_genes = length(entrez)
    ))
  }

  dot_plot <- enrichplot::dotplot(
    ego,
    showCategory = show_category,
    title = paste0("GO enrichment: ", ontology)
  )

  ggplot2::ggsave(
    filename = file.path(output_dir, paste0("GO_", ontology, "_dotplot.png")),
    plot = dot_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  list(
    go_results = ego_df,
    go_object = ego,
    dot_plot = dot_plot,
    ontology = ontology,
    n_input_genes = length(entrez)
  )
}
