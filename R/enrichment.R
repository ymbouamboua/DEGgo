
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
    entrez_col = NULL,
    min_genes = 10
) {


  ontology <- match.arg(ontology)
  log <- .deggo_msg(verbose = TRUE, prefix = "DEGgo")

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  comparison <- comparison %||% "comparison"

  if (is.null(entrez_col)) {
    entrez_col <- if ("entrez_id" %in% colnames(sig_deg)) {
      "entrez_id"
    } else if ("ENTREZID" %in% colnames(sig_deg)) {
      "ENTREZID"
    } else if ("entrez" %in% colnames(sig_deg)) {
      "entrez"
    } else {
      NA_character_
    }
  }
  
  if (is.na(entrez_col)) {
    log(
      paste0("No Entrez ID column found. Skipping GO for ", comparison, "."),
      type = "warn"
    )
    return(NULL)
  }
  
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



# ========================================================= #
# PLOT GO TERMS
# ========================================================= #
#' Plot GO terms by regulation status
#'
#' @param go_df GO enrichment data frame from \code{run_go_enrichment()}.
#' @param go_terms_of_interest Optional GO terms to display.
#' @param comparison Optional comparison name used as plot title.
#' @param top_n Number of GO terms to show per regulation.
#' @param wrap_width Width for wrapping GO term labels.
#' @param color_values Named colors for Up and Down.
#' @param size_range Point size range.
#' @param style Theme style passed to \code{.deggo_theme()}.
#' @param x_angle X-axis text angle.
#' @param txtsize Base font size.
#'
#' @return A ggplot object.
#' @export
plot_go_terms <- function(
    go_df,
    go_terms_of_interest = NULL,
    comparison = NULL,
    top_n = 10,
    wrap_width = 50,
    color_values = c(
      Up = "#740001",
      Down = "#6497b1"
    ),
    size_range = c(2, 6),
    style = "bw",
    x_angle = 0,
    txtsize = 12
) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.", call. = FALSE)

  if (!requireNamespace("stringr", quietly = TRUE))
    stop("Package 'stringr' is required.", call. = FALSE)

  if (!requireNamespace("forcats", quietly = TRUE))
    stop("Package 'forcats' is required.", call. = FALSE)

  if (is.null(go_df) || !is.data.frame(go_df) || nrow(go_df) == 0)
    stop("go_df is empty.", call. = FALSE)

  required <- c(
    "Description",
    "FoldEnrichment",
    "p.adjust",
    "Count",
    "Regulation"
  )

  miss <- setdiff(required, colnames(go_df))

  if (length(miss))
    stop(
      "Missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )

  ## ----------------------------- #
  ## Clean
  ## ----------------------------- #

  df <- go_df

  df$Description <- toupper(trimws(df$Description))
  df$Regulation <- factor(df$Regulation,
                          levels = c("Up", "Down"))

  df$FoldEnrichment <- as.numeric(df$FoldEnrichment)
  df$p.adjust <- as.numeric(df$p.adjust)
  df$Count <- as.integer(df$Count)

  keep <-
    !is.na(df$Regulation) &
    is.finite(df$FoldEnrichment) &
    is.finite(df$p.adjust) &
    df$p.adjust > 0

  df <- df[keep, , drop = FALSE]

  if (!nrow(df))
    stop("No valid GO terms after cleaning.", call. = FALSE)

  ## ----------------------------- #
  ## Keyword filtering
  ## ----------------------------- #

  if (!is.null(go_terms_of_interest)) {

    pattern <- paste(
      toupper(go_terms_of_interest),
      collapse = "|"
    )

    keep <- grepl(
      pattern,
      df$Description,
      ignore.case = TRUE
    )

    df <- df[keep, , drop = FALSE]

    if (!nrow(df))
      stop(
        "No GO terms matched the requested keywords.",
        call. = FALSE
      )
  }

  ## ----------------------------- #
  ## Keep best term per regulation
  ## ----------------------------- #

  ord <- order(
    df$Regulation,
    df$Description,
    df$p.adjust,
    -df$FoldEnrichment
  )

  df <- df[ord, ]

  dup <- duplicated(
    paste(df$Regulation, df$Description)
  )

  df <- df[!dup, , drop = FALSE]

  ## ----------------------------- #
  ## Top N
  ## ----------------------------- #

  if (is.null(go_terms_of_interest)) {

    split_df <- split(df, df$Regulation)

    split_df <- lapply(split_df, function(x) {

      x <- x[order(
        x$p.adjust,
        -x$FoldEnrichment
      ), ]

      head(x, top_n)

    })

    df <- do.call(rbind, split_df)
    rownames(df) <- NULL
  }

  ## ----------------------------- #
  ## Plot variables
  ## ----------------------------- #

  df$log10FDR <- -log10(df$p.adjust)

  df$Description_wrapped <-
    stringr::str_wrap(
      df$Description,
      width = wrap_width
    )

  df$Description_wrapped <-
    forcats::fct_reorder(
      df$Description_wrapped,
      df$log10FDR
    )

  regs <- unique(as.character(df$Regulation))

  color_values <- color_values[
    names(color_values) %in% regs
  ]

  ## ----------------------------- #
  ## Plot
  ## ----------------------------- #

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = log10FDR,
      y = Description_wrapped
    )
  ) +

    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = log10FDR,
        yend = Description_wrapped,
        colour = Regulation
      ),
      linewidth = 0.8
    ) +

    ggplot2::geom_point(
      ggplot2::aes(
        size = FoldEnrichment,
        colour = Regulation
      )
    ) +

    ggplot2::scale_colour_manual(
      values = color_values
    ) +

    ggplot2::scale_size(
      range = size_range
    ) +

    ggplot2::labs(
      x = expression(-log[10](FDR)),
      y = NULL,
      colour = "Regulation",
      size = "Fold enrichment",
      title = ifelse(
        is.null(comparison),
        "GO enrichment",
        comparison
      )
    ) +

    .deggo_theme(
      style = style,
      x.ang = x_angle,
      txtsize = txtsize
    )
}


# ========================================================= #
# PLOT ALL GO TERMS
# ========================================================= #
#' Plot GO enrichment terms across all DEGgo comparisons
#'
#' Generate Gene Ontology enrichment plots for all comparisons stored in a
#' DEGgo results object. This function iterates over `results$go_results`,
#' creates one GO plot per comparison using `plot_go_terms()`, and optionally
#' saves each plot as PNG and PDF.
#'
#' @param results A DEGgo results object returned by `run_deggo()`. Must contain
#'   a `go_results` list.
#' @param top_n Integer. Number of top enriched GO terms to display per plot.
#'   Default is `10`.
#' @param txtsize Numeric. Base text size used in the GO plots. Default is `8`.
#' @param style Character. Plot theme style passed to `plot_go_terms()`.
#'   Default is `"bw"`.
#' @param output_dir Optional character. Directory where GO plots are saved.
#'   If `NULL`, plots are returned but not written to disk.
#' @param width Numeric. Plot width in inches for saved files. Default is `8`.
#' @param height Numeric. Plot height in inches for saved files. Default is `6`.
#' @param dpi Integer. Resolution in dots per inch for saved PNG files.
#'   Default is `300`.
#'
#' @return A named list of `ggplot` objects, one per comparison. Comparisons
#'   without enriched GO terms return `NULL`.
#'
#' @examples
#' \dontrun{
#' go_plots <- plot_all_go_terms(
#'   results = results,
#'   top_n = 10,
#'   output_dir = "GO_plots"
#' )
#' }
#'
#' @export
#'
plot_all_go_terms <- function(
    results,
    top_n = 10,
    txtsize = 8,
    style = "bw",
    output_dir = NULL,
    width = 8,
    height = 6,
    dpi = 300
) {

  plots <- lapply(names(results$go_results), function(comp) {

    go_df <- results$go_results[[comp]]$go_results

    if (is.null(go_df) || nrow(go_df) == 0) {
      return(NULL)
    }

    p <- plot_go_terms(
      go_df = go_df,
      comparison = comp,
      top_n = top_n,
      txtsize = txtsize,
      style = style
    )

    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

      ggplot2::ggsave(
        filename = file.path(output_dir, paste0(comp, "_GO_terms.png")),
        plot = p,
        width = width,
        height = height,
        dpi = dpi,
        bg = "white"
      )

      ggplot2::ggsave(
        filename = file.path(output_dir, paste0(comp, "_GO_terms.pdf")),
        plot = p,
        width = width,
        height = height,
        bg = "white"
      )
    }

    p
  })

  names(plots) <- names(results$go_results)

  plots
}




.run_go_pipeline <- function(
    sig_deg,
    comparison,
    ontology,
    orgdb,
    dirs,
    txtsize
){

  go <- run_go_enrichment(
    sig_deg,
    comparison,
    ontology,
    dirs$go,
    orgdb
  )

  plot <- NULL

  if(!is.null(go$go_results) &&
     nrow(go$go_results)>0){

    plot <- plot_go_terms(
      go$go_results,
      comparison=paste0(comparison," GO"),
      txtsize=txtsize
    )

  }

  list(

    go_results = go,

    go_plot = plot

  )

}
