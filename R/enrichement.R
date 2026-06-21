#' Run Gene Ontology enrichment analysis
#'
#' Performs Gene Ontology enrichment analysis separately on
#' upregulated and downregulated significant differentially
#' expressed genes using Entrez identifiers.
#'
#' @param sig_deg Significant DEG data frame containing Entrez IDs
#'   and log2 fold-change values.
#' @param comparison Character. Comparison name used to label GO
#'   enrichment results and exported files.
#' @param ontology Character. GO ontology to test. One of
#'   \code{"BP"}, \code{"MF"} or \code{"CC"}.
#' @param orgdb OrgDb annotation database used by
#'   \code{clusterProfiler::enrichGO()}.
#' @param output_dir Character. Output directory.
#' @param pvalue_cutoff Numeric. Raw p-value cutoff.
#' @param qvalue_cutoff Numeric. q-value cutoff.
#' @param p_adjust_method Character. Multiple testing correction method.
#' @param logfc_col Character. Column containing log2 fold-change values.
#' @param entrez_col Character. Column containing Entrez gene identifiers.
#' @param min_genes Integer. Minimum number of genes required to
#'   perform enrichment.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{go_results}}{GO enrichment results table.}
#'   \item{\code{ontology}}{Ontology tested.}
#'   \item{\code{comparison}}{Comparison identifier.}
#'   \item{\code{n_up_genes}}{Number of upregulated genes tested.}
#'   \item{\code{n_down_genes}}{Number of downregulated genes tested.}
#' }
#'
#' @export
#'
run_go_enrichment <- function(
    sig_deg,
    comparison = NULL,
    ontology = c("BP", "MF", "CC"),
    orgdb,
    output_dir = "DEGgo_out",
    pvalue_cutoff = 0.05,
    qvalue_cutoff = 0.2,
    p_adjust_method = "BH",
    logfc_col = "log2FoldChange",
    entrez_col = "ENTREZID",
    min_genes = 10
) {

  ontology <- match.arg(ontology)
  log <- .deggo_msg(verbose = TRUE, prefix = "DEGgo")

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  comparison <- comparison %||% "comparison"

  if (is.null(sig_deg) || !is.data.frame(sig_deg) || nrow(sig_deg) == 0) {
    log(paste0("No significant DEGs for ", comparison, ". Skipping GO."), type = "warn")
    return(NULL)
  }

  required <- c(entrez_col, logfc_col)
  missing <- setdiff(required, colnames(sig_deg))

  if (length(missing) > 0) {
    log(
      paste0(
        "Missing column(s): ",
        paste(missing, collapse = ", "),
        ". Skipping GO for ",
        comparison,
        "."
      ),
      type = "warn"
    )
    return(NULL)
  }

  sig_deg[[logfc_col]] <- suppressWarnings(as.numeric(sig_deg[[logfc_col]]))

  deg_list <- list(
    Up = sig_deg[
      !is.na(sig_deg[[logfc_col]]) &
        sig_deg[[logfc_col]] > 0,
      ,
      drop = FALSE
    ],
    Down = sig_deg[
      !is.na(sig_deg[[logfc_col]]) &
        sig_deg[[logfc_col]] < 0,
      ,
      drop = FALSE
    ]
  )

  run_one_go <- function(df, Regulation) {

    entrez <- unique(stats::na.omit(as.character(df[[entrez_col]])))

    if (length(entrez) < min_genes) {
      log(
        paste0(
          comparison,
          " - ",
          Regulation,
          ": fewer than ",
          min_genes,
          " genes."
        ),
        type = "warn"
      )
      return(data.frame())
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

    if (nrow(ego_df) == 0) {
      log(
        paste0("No enriched GO terms for ", comparison, " - ", Regulation),
        type = "warn"
      )
      return(data.frame())
    }

    ego_df$Description <- toupper(ego_df$Description)
    ego_df$comparison <- comparison
    ego_df$Regulation <- Regulation
    ego_df$ontology <- ontology
    ego_df$n_input_genes <- length(entrez)

    ego_df
  }

  go_table <- do.call(
    rbind,
    lapply(names(deg_list), function(reg) {
      run_one_go(deg_list[[reg]], reg)
    })
  )

  if (is.null(go_table) || nrow(go_table) == 0) {
    return(NULL)
  }

  rownames(go_table) <- NULL

  utils::write.table(
    go_table,
    file = file.path(
      output_dir,
      paste0(comparison, "_GO_", ontology, "_enrichment.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  list(
    go_results = go_table,
    ontology = ontology,
    comparison = comparison,
    n_up_genes = nrow(deg_list$Up),
    n_down_genes = nrow(deg_list$Down)
  )
}
