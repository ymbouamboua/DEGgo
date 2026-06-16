#' Run Gene Ontology enrichment analysis
#'
#' Performs Gene Ontology enrichment analysis separately on upregulated and
#' downregulated significant differentially expressed genes using Entrez IDs.
#'
#' @param sig_deg Significant DEG data frame. Must contain Entrez IDs and
#'   log2 fold-change values.
#' @param ontology Character. GO ontology to test. One of \code{"BP"},
#'   \code{"MF"}, or \code{"CC"}.
#' @param orgdb OrgDb annotation database used by
#'   \code{clusterProfiler::enrichGO()}.
#' @param output_dir Character. Directory where enrichment tables and plots
#'   are exported.
#' @param show_category Integer. Number of GO terms to display per regulation.
#' @param pvalue_cutoff Numeric. P-value cutoff for enrichment.
#' @param qvalue_cutoff Numeric. Q-value cutoff for enrichment.
#' @param p_adjust_method Character. Multiple-testing correction method.
#' @param logfc_col Character. Column containing log2 fold-change values.
#' @param entrez_col Character. Column containing Entrez IDs.
#' @param min_genes Integer. Minimum number of genes required for enrichment.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{go_results}}{Merged GO enrichment table with Regulation column.}
#'   \item{\code{go_by_regulation}}{Separate Up and Down enrichGO results.}
#'   \item{\code{go_plot}}{Combined Up/Down GO plot.}
#'   \item{\code{ontology}}{Ontology tested.}
#'   \item{\code{n_up_genes}}{Number of upregulated genes tested.}
#'   \item{\code{n_down_genes}}{Number of downregulated genes tested.}
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
    p_adjust_method = "BH",
    logfc_col = "log2FoldChange",
    entrez_col = "ENTREZID",
    min_genes = 10
) {

  ontology <- match.arg(ontology)

  log <- .msg(verbose = TRUE)
  log("Running GO enrichment for Up and Down DEGs...")

  for (pkg in c("clusterProfiler", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required.", call. = FALSE)
    }
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(sig_deg) || !is.data.frame(sig_deg) || nrow(sig_deg) == 0) {
    log("No significant DEGs found. Skipping GO enrichment.", type = "warn")
    return(NULL)
  }

  required <- c(entrez_col, logfc_col)
  missing <- setdiff(required, colnames(sig_deg))

  if (length(missing) > 0) {
    log(
      paste0("Missing column(s): ", paste(missing, collapse = ", "),
             ". Skipping GO enrichment."),
      type = "warn"
    )
    return(NULL)
  }

  sig_deg[[logfc_col]] <- as.numeric(sig_deg[[logfc_col]])

  deg_list <- list(
    Up = sig_deg[!is.na(sig_deg[[logfc_col]]) & sig_deg[[logfc_col]] > 0, , drop = FALSE],
    Down = sig_deg[!is.na(sig_deg[[logfc_col]]) & sig_deg[[logfc_col]] < 0, , drop = FALSE]
  )

  run_one_go <- function(df, regulation) {

    entrez <- unique(stats::na.omit(as.character(df[[entrez_col]])))

    if (length(entrez) < min_genes) {
      log(paste0("Skipping ", regulation, ": fewer than ", min_genes, " genes."), type = "warn")

      return(list(
        go_results = data.frame(),
        go_object = NULL,
        regulation = regulation,
        n_input_genes = length(entrez)
      ))
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

    if (!nrow(ego_df)) {
      log(paste0("No enriched GO terms for ", regulation, "."), type = "warn")

      return(list(
        go_results = data.frame(),
        go_object = ego,
        regulation = regulation,
        n_input_genes = length(entrez)
      ))
    }

    ego_df$Description <- toupper(ego_df$Description)
    ego_df$Regulation <- regulation
    ego_df$Ontology <- ontology
    ego_df$n_input_genes <- length(entrez)

    out_prefix <- paste0("GO_", ontology, "_", regulation)

    utils::write.csv(
      ego_df,
      file = file.path(output_dir, paste0(out_prefix, "_enrichment.csv")),
      row.names = FALSE
    )

    utils::write.table(
      ego_df,
      file = file.path(output_dir, paste0(out_prefix, "_enrichment.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    saveRDS(
      ego,
      file = file.path(output_dir, paste0(out_prefix, "_enrichGO.rds"))
    )

    list(
      go_results = ego_df,
      go_object = ego,
      regulation = regulation,
      n_input_genes = length(entrez)
    )
  }

  go_by_regulation <- lapply(names(deg_list), function(reg) {
    run_one_go(deg_list[[reg]], reg)
  })

  names(go_by_regulation) <- names(deg_list)

  go_table <- do.call(
    rbind,
    lapply(go_by_regulation, function(x) {
      if (is.null(x$go_results) || nrow(x$go_results) == 0) return(NULL)
      x$go_results
    })
  )

  go_plot <- NULL

  if (is.null(go_table) || nrow(go_table) == 0) {
    log("No significant GO terms found for Up or Down genes.", type = "warn")
    go_table <- data.frame()
  } else {
    rownames(go_table) <- NULL

    utils::write.csv(
      go_table,
      file = file.path(output_dir, paste0("GO_", ontology, "_UpDown_enrichment.csv")),
      row.names = FALSE
    )

    utils::write.table(
      go_table,
      file = file.path(output_dir, paste0("GO_", ontology, "_UpDown_enrichment.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    go_plot <- plot_go_terms(
      go_df = go_table,
      comparison = paste0("GO ", ontology, " enrichment"),
      top_n = show_category
    )

    ggplot2::ggsave(
      filename = file.path(output_dir, paste0("GO_", ontology, "_UpDown_plot.png")),
      plot = go_plot,
      width = 8,
      height = 6,
      dpi = 300
    )

    ggplot2::ggsave(
      filename = file.path(output_dir, paste0("GO_", ontology, "_UpDown_plot.pdf")),
      plot = go_plot,
      width = 8,
      height = 6
    )
  }

  list(
    go_results = go_table,
    go_by_regulation = go_by_regulation,
    go_plot = go_plot,
    ontology = ontology,
    n_up_genes = nrow(deg_list$Up),
    n_down_genes = nrow(deg_list$Down)
  )
}
