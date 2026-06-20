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
.deggo_dirs <- function(output_dir, analysis_mode) {
  dirs <- list(
    qc = file.path(output_dir, paste0(analysis_mode, "_QC")),
    pca = file.path(output_dir, paste0(analysis_mode, "_PCA")),
    volcano = file.path(output_dir, paste0(analysis_mode, "_volcano")),
    heatmaps = file.path(output_dir, paste0(analysis_mode, "_heatmaps")),
    results = file.path(output_dir, paste0(analysis_mode, "_results")),
    significant = file.path(output_dir, paste0(analysis_mode, "_significant")),
    go = file.path(output_dir, paste0(analysis_mode, "_GO"))
  )

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



# ======================================================= #
# .write_deggo_manifest
# ======================================================= #
#' Write DEGgo output manifest
#'
#' Internal helper used by DEGgo reports to summarize generated output files.
#'
#' @param output_dir DEGgo output directory.
#' @param analysis_mode Analysis mode: "single" or "pairwise".
#' @param method Differential expression method.
#'
#' @return Invisibly returns a data frame.
#'
#' @keywords internal
#' @noRd
.write_deggo_manifest <- function(
    output_dir,
    dirs,
    analysis_mode = c("single", "pairwise")
) {

  analysis_mode <- match.arg(analysis_mode)

  manifest <- data.frame(
    Folder = c(
      basename(dirs$results),
      basename(dirs$significant),
      basename(dirs$go),
      basename(dirs$heatmaps),
      basename(dirs$pca),
      basename(dirs$qc),
      paste0(analysis_mode, "_summary.tsv"),
      "sessionInfo.txt",
      "reproducibility"
    ),
    Content = c(
      "Full differential expression results",
      "Significant DEGs only",
      "GO enrichment tables and GO plots",
      "DEG heatmaps",
      "Sample-level PCA plots",
      "Sample correlation and hierarchical clustering",
      "Summary of differential expression results",
      "Reproducibility and software versions",
      "Run parameters and cleaned input files"
    ),
    stringsAsFactors = FALSE
  )

  utils::write.table(
    manifest,
    file.path(output_dir, "DEGgo_output_manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(manifest)
}


# ======================================================= #
# .safe_report
# ======================================================= #
#' Safely generate a DEGgo report
#'
#' Internal wrapper around [generate_deggo_report()] that prevents report
#' generation errors from stopping the main DEGgo workflow.
#'
#' @param res DEGgo result object.
#'
#' @return Report output invisibly, or `NULL` if report generation is disabled
#'   or fails.
#'
#' @keywords internal
#' @noRd
.safe_report <- function(
    res,
    output_dir,
    generate_report = TRUE,
    report_formats = "html",
    report_template = NULL
) {
  if (!isTRUE(generate_report)) return(NULL)

  tryCatch(
    generate_deggo_report(
      results = res,
      output_dir = output_dir,
      formats = report_formats,
      report_template = report_template
    ),
    error = function(e) {
      warning("Report generation failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}


# ======================================================= #
# .make_run_params
# ======================================================= #
#' Create DEGgo run parameter metadata
#'
#' Builds a named list describing the main parameters used in the current
#' DEGgo run. This is used for reproducibility tracking and reporting.
#'
#' @return A named list of run parameters.
#'
#' @keywords internal
#' @noRd
.make_run_params <- function(
    deggo_version,
    organism,
    method,
    analysis_mode,
    ontology,
    padj_cutoff,
    logfc_cutoff,
    filter_method,
    min_count,
    min_samples,
    min_total,
    top_n_heatmap,
    top_n_labels,
    prepare_input,
    gene_col,
    feature_col,
    sample_col,
    design_formula,
    pairwise_group_cols,
    pairwise_contrast_col,
    pairwise_mode,
    output_dir,
    repro_dir
) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  list(
    deggo_version = deggo_version,
    date = as.character(Sys.Date()),
    organism = organism,
    method = method,
    analysis_mode = analysis_mode,
    ontology = ontology,
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff,
    filter_method = filter_method,
    min_count = min_count,
    min_samples = min_samples,
    min_total = min_total,
    top_n_heatmap = top_n_heatmap,
    top_n_labels = top_n_labels,
    prepare_input = prepare_input,
    gene_col = paste(gene_col, collapse = ", "),
    feature_col = paste(feature_col, collapse = ", "),
    sample_col = paste(sample_col, collapse = ", "),
    design_formula = paste(deparse(design_formula), collapse = ""),
    pairwise_group_cols = paste(pairwise_group_cols %||% NA, collapse = ", "),
    pairwise_contrast_col = pairwise_contrast_col,
    pairwise_mode = pairwise_mode,
    output_dir = output_dir,
    reproducibility_dir = repro_dir
  )
}

# ======================================================= #
# .save_repro
# ======================================================= #
#' Save DEGgo reproducibility files
#'
#' Writes DEGgo run parameters to the reproducibility directory.
#'
#' @param res DEGgo result object containing a `run_params` element.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
.save_repro <- function(
    res,
    repro_dir,
    save_reproducibility = TRUE
) {
  if (!isTRUE(save_reproducibility)) return(invisible(NULL))

  params_df <- data.frame(
    parameter = names(res$run_params),
    value = vapply(
      res$run_params,
      function(x) paste(x, collapse = ", "),
      character(1)
    ),
    stringsAsFactors = FALSE
  )

  utils::write.table(
    params_df,
    file.path(repro_dir, "run_parameters.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(NULL)
}

# ======================================================= #
# .save_clean_input_files
# ======================================================= #
#' Save cleaned input files
#'
#' Saves the filtered count matrix, matched metadata, and a summary of the
#' cleaned input data to the reproducibility directory.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
.save_clean_input_files <- function(
    counts,
    metadata,
    repro_dir,
    save_clean_inputs = TRUE
) {
  if (!isTRUE(save_clean_inputs)) return(invisible(NULL))

  utils::write.table(
    counts,
    file.path(repro_dir, "clean_filtered_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  utils::write.table(
    metadata,
    file.path(repro_dir, "clean_matched_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  lib_sizes <- colSums(counts, na.rm = TRUE)

  input_summary <- data.frame(
    metric = c(
      "samples",
      "genes_after_filtering",
      "total_counts",
      "median_library_size",
      "min_library_size",
      "max_library_size"
    ),
    value = c(
      ncol(counts),
      nrow(counts),
      sum(counts, na.rm = TRUE),
      stats::median(lib_sizes, na.rm = TRUE),
      min(lib_sizes, na.rm = TRUE),
      max(lib_sizes, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  utils::write.table(
    input_summary,
    file.path(repro_dir, "clean_input_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(NULL)
}


# ======================================================= #
# .make_pca_list
# ======================================================= #
#' Generate PCA plots for DEGgo results
#'
#' Creates PCA plots colored by sample and available metadata variables such as
#' tissue, treatment, sex, and condition.
#'
#' @param dds A DESeq2 dataset object.
#' @param md A metadata data frame matched to the count matrix.
#'
#' @return A named list of PCA plot objects.
#'
#' @keywords internal
#' @noRd
.make_pca_list <- function(dds, md, pca_dir) {
  pca_list <- list()

  pca_list$sample <- plot_pca(
    dds = dds,
    metadata = md,
    output_dir = pca_dir,
    filename = "PCA_by_sample",
    color_by = "sample",
    title = "PCA by sample"
  )

  for (col in intersect(c("tissue", "treatment", "sex", "condition"), colnames(md))) {
    pca_list[[col]] <- plot_pca(
      dds = dds,
      metadata = md,
      output_dir = pca_dir,
      filename = paste0("PCA_by_", col),
      color_by = col,
      title = paste("PCA by", col)
    )
  }

  if (all(c("tissue", "treatment") %in% colnames(md))) {
    pca_list$tissue_treatment <- plot_pca(
      dds = dds,
      metadata = md,
      output_dir = pca_dir,
      filename = "PCA_tissue_treatment",
      color_by = "tissue",
      shape_by = "treatment",
      title = "PCA tissue + treatment"
    )
  }

  pca_list
}


