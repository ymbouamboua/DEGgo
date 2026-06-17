# ========================================================= #
# Null helper
# ========================================================= #
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# ========================================================= #
# Create deggo directories
# ========================================================= #
#' Create deggo directories
#' @keywords internal
#' @noRd
.deggo_dirs <- function(output_dir, analysis_mode = "pairwise") {

  if (analysis_mode == "pairwise") {
    dirs <- list(
      results     = file.path(output_dir, "pairwise_results"),
      significant = file.path(output_dir, "pairwise_significant"),
      volcano     = file.path(output_dir, "pairwise_volcano"),
      heatmaps    = file.path(output_dir, "pairwise_heatmaps"),
      go          = file.path(output_dir, "pairwise_GO"),
      pca         = file.path(output_dir, "pairwise_PCA")
    )
  } else {
    dirs <- list(
      results     = file.path(output_dir, "results"),
      significant = file.path(output_dir, "significant"),
      volcano     = file.path(output_dir, "volcano_plots"),
      heatmaps    = file.path(output_dir, "heatmaps"),
      go          = file.path(output_dir, "GO"),
      pca         = file.path(output_dir, "PCA")
    )
  }

  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}



# ========================================================= #
# Merge pairwise GO tables
# ========================================================= #
#' Merge pairwise GO
#' @keywords internal
#' @noRd
.merge_pairwise_go <- function(go_results, go_dir, ontology) {

  go_merged <- do.call(
    rbind,
    lapply(go_results, function(x) {
      if (is.null(x) || is.null(x$go_results) || nrow(x$go_results) == 0) {
        return(NULL)
      }
      x$go_results
    })
  )

  if (is.null(go_merged) || nrow(go_merged) == 0) {
    return(NULL)
  }

  rownames(go_merged) <- NULL

  utils::write.table(
    go_merged,
    file.path(go_dir, paste0("ALL_pairwise_GO_", ontology, "_merged.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  go_merged
}

# ========================================================= #
# FORMAT RUNTIME
# ========================================================= #

#' Format runtime
#' @keywords internal
#' @noRd
.format_duration <- function(seconds) {

  if (is.null(seconds) || is.na(seconds)) {
    return(NA_character_)
  }

  seconds <- as.numeric(seconds)

  hrs  <- floor(seconds / 3600)
  mins <- floor((seconds %% 3600) / 60)
  secs <- seconds %% 60

  if (hrs > 0) {
    sprintf("%dh %dm %.1fs", hrs, mins, secs)
  } else if (mins > 0) {
    sprintf("%dm %.1fs", mins, secs)
  } else {
    sprintf("%.1fs", secs)
  }
}


# ========================================================= #
# LOGGER
# ========================================================= #

#' Create package logger
#'
#' Lightweight logger used throughout GnRHcell.
#'
#' @param verbose Logical.
#' @param prefix Package prefix.
#'
#' @return Logging function.
#'
#' @keywords internal
#' @noRd
.msg <- function(
    verbose = TRUE,
    prefix = "INFO"
) {

  t0 <- Sys.time()

  function(
    ...,
    type = c(
      "info",
      "step",
      "done",
      "warn",
      "error",
      "header"
    ),
    duration = NULL
  ) {

    if (!isTRUE(verbose)) {
      return(invisible(NULL))
    }

    type <- match.arg(type)

    txt <- paste(..., collapse = " ")

    if (type == "done") {

      if (is.null(duration)) {
        duration <- as.numeric(
          difftime(
            Sys.time(),
            t0,
            units = "secs"
          )
        )
      }

      txt <- paste0(
        txt,
        " Duration: ",
        .format_duration(duration)
      )
    }

    tag <- switch(
      type,
      info   = "[INFO]",
      step   = "[STEP]",
      done   = "[DONE]",
      warn   = "[WARN]",
      error  = "[ERROR]",
      header = paste0("[", prefix, "]")
    )

    cat(
      tag,
      " ",
      txt,
      "\n",
      sep = ""
    )

    invisible(NULL)
  }
}


# ========================================================= #
# SAVE SESSION INFO
# ========================================================= #

#' Save session info
#'
#' @param output_dir outdir directory
#'
#' @keywords internal
#' @noRd
.safe_write_session_info <- function(output_dir) {
  x <- tryCatch(
    utils::capture.output(utils::sessionInfo()),
    error = function(e) {
      c(
        "sessionInfo() failed.",
        paste("Reason:", conditionMessage(e)),
        "",
        "Fallback:",
        utils::capture.output(Sys.info()),
        utils::capture.output(R.version)
      )
    }
  )

  writeLines(
    x,
    file.path(output_dir, "sessionInfo.txt")
  )
}


# ========================================================= #
# .flag_low_mad
# ========================================================= #

#' Flag low outliers using MAD
#'
#' @param x Numeric vector.
#' @param nmad Number of MADs below the median.
#'
#' @return Logical vector.
#' @keywords internal
#' @noRd
.flag_low_mad <- function(x, nmad = 3) {
  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, na.rm = TRUE)
  x < med - nmad * madv
}

# ========================================================= #
# .flag_high_mad
# ========================================================= #

#' Flag high outliers using MAD
#'
#' @param x Numeric vector.
#' @param nmad Number of MADs above the median.
#'
#' @return Logical vector.
#' @keywords internal
#' @noRd
.flag_high_mad <- function(x, nmad = 3) {
  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, na.rm = TRUE)
  x > med + nmad * madv
}



# ========================================================= #
# EXTRACT NORMALIZED EXPRESSION
# ========================================================= #

#' Extract normalized gene expression
#'
#' Extract normalized or variance-stabilized expression values from a
#' DESeq2 object and return a long-format table.
#'
#' @param dds DESeqDataSet object.
#' @param metadata Sample metadata.
#' @param genes Character vector of gene symbols or identifiers.
#' @param assay One of \code{"vst"}, \code{"normalized"}, \code{"log2_normalized"} or \code{"raw"}.
#' @param gene_col Character. Column containing gene names or identifiers in rowData(dds).
#' @param collapse_fun Function used to collapse duplicated gene rows.
#'
#' @return Long-format expression table.
#' @export
#'
extract_expression <- function(
    dds,
    metadata,
    genes,
    assay = c("vst", "normalized", "log2_normalized", "raw"),
    gene_col = "SYMBOL",
    collapse_fun = sum
) {

  assay <- match.arg(assay)

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!"sample" %in% colnames(metadata)) {
    stop("metadata must contain a 'sample' column.", call. = FALSE)
  }

  rd <- as.data.frame(SummarizedExperiment::rowData(dds))

  mat <- switch(
    assay,
    raw = DESeq2::counts(dds, normalized = FALSE),
    normalized = DESeq2::counts(dds, normalized = TRUE),
    log2_normalized = log2(DESeq2::counts(dds, normalized = TRUE) + 1),
    vst = SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
  )

  gene_ids <- if (gene_col %in% colnames(rd)) {
    rd[[gene_col]]
  } else {
    rownames(mat)
  }

  gene_ids <- as.character(gene_ids)

  keep <- gene_ids %in% genes & !is.na(gene_ids) & gene_ids != ""

  if (!any(keep)) {
    stop("None of the requested genes were found.", call. = FALSE)
  }

  mat <- mat[keep, , drop = FALSE]
  gene_ids <- gene_ids[keep]

  mat <- stats::aggregate(
    mat,
    by = list(gene = gene_ids),
    FUN = collapse_fun
  )

  rownames(mat) <- mat$gene
  mat$gene <- NULL

  found_genes <- intersect(genes, rownames(mat))
  mat <- mat[found_genes, , drop = FALSE]

  expr_wide <- as.data.frame(t(mat), check.names = FALSE)
  expr_wide$sample <- rownames(expr_wide)
  rownames(expr_wide) <- NULL

  metadata <- metadata[!duplicated(metadata$sample), , drop = FALSE]

  expr_wide <- merge(
    expr_wide,
    metadata,
    by = "sample",
    all.x = TRUE,
    sort = FALSE
  )

  expr_long <- stats::reshape(
    expr_wide,
    varying = found_genes,
    v.names = "expression",
    timevar = "gene",
    times = found_genes,
    direction = "long"
  )

  rownames(expr_long) <- NULL

  expr_long$assay <- assay

  expr_long <- expr_long[
    ,
    c(
      "sample",
      "gene",
      "expression",
      setdiff(colnames(expr_long), c("sample", "gene", "expression", "id", "assay")),
      "assay"
    ),
    drop = FALSE
  ]

  expr_long$id <- NULL

  expr_long
}
