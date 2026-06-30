# ============================================================ #
# DEGgo main workflow helpers
# ============================================================ #

`%||%` <- function(x, y) if (is.null(x)) y else x

# ============================================================ #
# Main run_deggo
# ============================================================ #
#' Run DEGgo bulk RNA-seq downstream analysis
#'
#' Run a complete and automated bulk RNA-seq downstream analysis workflow.
#'
#' DEGgo can perform optional raw sample quality control, optional removal of
#' flagged samples, marker-based biological validation, input preparation,
#' sample matching, gene identifier cleaning, low-expression filtering,
#' differential expression analysis, result annotation, PCA, heatmap and volcano
#' visualization, Gene Ontology enrichment, reproducibility export, and automated
#' HTML/PDF/PPTX report generation.
#'
#' The function supports either a single differential expression analysis or
#' multiple pairwise contrasts. Pairwise mode currently uses DESeq2.
#'
#' @param counts Raw count table, data frame or matrix. Rows should represent
#'   genes/features and columns should represent samples. If `prepare_input =
#'   TRUE`, DEGgo attempts to detect gene and sample columns automatically using
#'   `gene_col`, `feature_col`, and `sample_col`.
#' @param metadata Sample metadata data frame. Must contain a sample identifier
#'   column matching the count matrix sample names. For single analysis, it must
#'   contain `condition` unless this is encoded in the design/contrast workflow.
#'   For pairwise analysis, `condition` can be automatically generated from
#'   `pairwise_group_cols`.
#' @param gene_col Character vector of possible gene identifier column names.
#'   Used during input preparation.
#' @param feature_col Character vector of possible feature/gene symbol column
#'   names. Used during input preparation and marker-based checks.
#' @param sample_col Character vector of possible sample identifier column names
#'   in `metadata`.
#' @param prepare_input Logical. If `TRUE`, automatically prepares and matches
#'   the count table and metadata using `prepare_counts_metadata()`. If `FALSE`,
#'   `counts` must already be a numeric matrix with sample names matching
#'   metadata row names or a sample column.
#' @param raw_qc Logical. If `TRUE`, run exploratory raw sample QC before input
#'   preparation and differential expression analysis using
#'   `explore_bulk_rnaseq()`.
#' @param remove_flagged Logical. If `TRUE`, remove samples flagged by the raw QC
#'   table before continuing the workflow. Default is `FALSE` to keep sample
#'   exclusion under user control.
#' @param qc_markers Optional character vector of marker genes used during raw
#'   and clean QC.
#' @param marker_sets Optional named list of marker gene sets used for
#'   marker-based biological validation with `marker_score_check()`.
#' @param qc_sample_col Optional sample column used by `remove_flagged_samples()`.
#'   If `NULL`, the first value of `sample_col` is used.
#' @param qc_output_prefix Character prefix used for QC output directories.
#' @param output_dir Output directory. If `NULL`, results are written to
#'   `"DEGgo_results"`. A dated DEGgo subdirectory is created automatically.
#' @param padj_cutoff Adjusted p-value cutoff used to define significant DEGs.
#' @param logfc_cutoff Absolute log2 fold-change cutoff used to define
#'   significant DEGs.
#' @param top_n_heatmap Number of top genes to display in DEG heatmaps.
#' @param top_n_labels Number of top genes to label in volcano plots.
#' @param ontology Gene Ontology namespace. One of `"BP"`, `"MF"`, or `"CC"`.
#' @param organism Organism used for gene annotation. One of `"human"`,
#'   `"mouse"`, `"rat"`, or `"custom"`.
#' @param orgdb Optional AnnotationDbi OrgDb object. Required when
#'   `organism = "custom"`.
#' @param method Differential expression method. One of `"DESeq2"`, `"edgeR"`,
#'   or `"limma"`. Pairwise mode currently supports DESeq2 only.
#' @param analysis_mode Analysis mode. `"single"` runs one differential
#'   expression analysis. `"pairwise"` runs multiple pairwise DESeq2 contrasts.
#' @param contrast Optional contrast passed to the differential expression
#'   engine in single-analysis mode.
#' @param design_formula Design formula used by the differential expression
#'   engine, for example `~ condition` or `~ batch + condition`.
#' @param pairwise_group_cols Character vector of metadata columns used to define
#'   pairwise groups. Required when `analysis_mode = "pairwise"`.
#' @param pairwise_contrast_col Name of the metadata column used to store or
#'   generate pairwise contrast groups.
#' @param pairwise_contrasts Optional named list of pairwise contrasts. If
#'   provided, these contrasts are used instead of automatically generated
#'   contrasts.
#' @param filter_method Low-expression filtering method. One of `"count"`,
#'   `"cpm"`, or `"none"`.
#' @param pairwise_mode Pairwise comparison mode. One of `"all"`,
#'   `"within_first"`, or `"within_second"`.
#' @param min_count Minimum count threshold used by count-based filtering.
#' @param min_samples Minimum number of samples required to pass the count/CPM
#'   threshold.
#' @param min_total Minimum total count required for a gene to be retained.
#' @param generate_report Logical. If `TRUE`, generate DEGgo report files.
#' @param report_formats Character vector of report formats, for example
#'   `"html"`, `"pdf"`, or `c("html", "pdf")`.
#' @param report_template Optional path to a custom R Markdown report template.
#'   If `NULL`, the package DEGgo template is used.
#' @param generate_pptx Logical. If `TRUE`, generate a PowerPoint report.
#' @param pptx_file Optional path to the output PowerPoint file. If `NULL`,
#'   `"DEGgo_Report.pptx"` is written inside `output_dir`.
#' @param save_reproducibility Logical. If `TRUE`, save reproducibility files,
#'   including run objects and session information.
#' @param save_clean_inputs Logical. If `TRUE`, save cleaned count and metadata
#'   tables after sample matching and filtering.
#' @param txtsize Base text size used in DEGgo visualizations.
#' @param seed Random seed for reproducibility.
#'
#' @return A named list containing differential expression results, significant
#'   DEGs, summary tables, plots, GO enrichment results, QC outputs, cleaned
#'   counts and metadata, report paths, PowerPoint path, output directories,
#'   run parameters, manifest, and DEGgo version.
#'
#' @details
#' In single mode, DEGgo runs one differential expression analysis using the
#' selected method and design formula. In pairwise mode, DEGgo builds or uses
#' predefined pairwise contrasts from metadata columns and runs DESeq2 for each
#' contrast.
#'
#' The raw QC step is performed before input preparation and filtering. The
#' sample QC step is performed after gene filtering and uses the final matched
#' count matrix and metadata.
#'
#' Automatic sample removal is disabled by default because exclusion of samples
#' should usually be reviewed by the analyst.
#'
#' @examples
#' \dontrun{
#' results <- run_deggo(
#'   counts = counts,
#'   metadata = metadata,
#'   organism = "mouse",
#'   method = "DESeq2",
#'   analysis_mode = "single",
#'   design_formula = ~ condition,
#'   contrast = c("condition", "treated", "control"),
#'   sample_col = "sample"
#' )
#'
#' pairwise_results <- run_deggo(
#'   counts = counts,
#'   metadata = metadata,
#'   organism = "mouse",
#'   method = "DESeq2",
#'   analysis_mode = "pairwise",
#'   pairwise_group_cols = c("treatment", "sex", "tissue"),
#'   pairwise_contrasts = pairwise_contrasts,
#'   sample_col = "sample",
#'   raw_qc = TRUE,
#'   remove_flagged = FALSE
#' )
#' }
#'
#' @seealso
#' `prepare_counts_metadata`,
#' `explore_bulk_rnaseq`,
#' `remove_flagged_samples`,
#' `marker_score_check`,
#' `run_de`,
#' `run_deseq2_pairwise`,
#' `run_go_enrichment`,
#' `generate_deggo_report`,
#' `generate_deggo_pptx`
#'
#' @export
#'
run_deggo <- function(
    counts,
    metadata,
    gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_id"),
    feature_col = c("gene_name", "SYMBOL", "symbol", "gene_symbol", "external_gene_name"),
    sample_col = c("sample", "Sample", "SAMPLE"),
    prepare_input = TRUE,
    raw_qc = TRUE,
    remove_flagged = FALSE,
    qc_markers = NULL,
    marker_sets = NULL,
    qc_sample_col = NULL,
    qc_output_prefix = "DEGgo_QC",
    output_dir = NULL,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25,
    top_n_heatmap = 50,
    top_n_labels = 10,
    ontology = c("BP", "MF", "CC"),
    organism = c("human", "mouse", "rat", "custom"),
    orgdb = NULL,
    method = c("DESeq2", "edgeR", "limma"),
    analysis_mode = c("single", "pairwise"),
    contrast = NULL,
    design_formula = ~ condition,
    pairwise_group_cols = NULL,
    pairwise_contrast_col = "comparison_group",
    pairwise_contrasts = NULL,
    filter_method = c("count", "cpm", "none"),
    pairwise_mode = c("all", "within_first", "within_second"),
    min_count = 5,
    min_samples = 2,
    min_total = 10,
    generate_report = TRUE,
    report_formats = "html",
    report_template = NULL,
    generate_pptx = FALSE,
    pptx_file = NULL,
    save_reproducibility = TRUE,
    save_clean_inputs = TRUE,
    txtsize = 12,
    seed = 123
) {

  set.seed(seed)

  method <- match.arg(method)
  filter_method <- match.arg(filter_method)
  organism <- match.arg(organism)
  analysis_mode <- match.arg(analysis_mode)
  pairwise_mode <- match.arg(pairwise_mode)
  ontology <- match.arg(ontology)

  log <- .deggo_msg(verbose = TRUE, prefix = "DEGgo")
  t_start <- Sys.time()

  log("==== STARTING DEGgo ANALYSIS ====", type = "header")

  if (is.null(output_dir)) {
    output_dir <- "DEGgo_results"
  }

  output_dir <- file.path(
    output_dir,
    paste0("DEGgo_", format(Sys.Date(), "%Y_%m_%d"), "_", method, "_", analysis_mode)
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  dirs <- .deggo_dirs(output_dir, analysis_mode)

  repro_dir <- file.path(output_dir, "reproducibility")
  dir.create(repro_dir, recursive = TRUE, showWarnings = FALSE)

  deggo_version <- tryCatch(
    as.character(utils::packageVersion("DEGgo")),
    error = function(e) "development"
  )

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (is.null(report_template)) {

    candidates <- c(
      system.file(
        "rmarkdown/templates/deggo_report/skeleton/skeleton.Rmd",
        package = "DEGgo"
      ),
      file.path(
        getwd(),
        "inst/rmarkdown/templates/deggo_report/skeleton/skeleton.Rmd"
      ),
      file.path(
        dirname(getwd()),
        "inst/rmarkdown/templates/deggo_report/skeleton/skeleton.Rmd"
      ),
      file.path(
        dirname(dirname(getwd())),
        "inst/rmarkdown/templates/deggo_report/skeleton/skeleton.Rmd"
      )
    )

    candidates <- candidates[nzchar(candidates)]
    candidates <- candidates[file.exists(candidates)]

    if (!length(candidates)) {
      warning(
        "DEGgo report template not found. Report generation will be skipped."
      )
      report_template <- NULL
    } else {
      report_template <- candidates[1]
    }
  }

  qc_raw <- NULL
  qc_clean <- NULL
  marker_check <- NULL

  # ---------------------------------------------------------- #
  # Raw QC
  # ---------------------------------------------------------- #

  if (isTRUE(raw_qc)) {

    log("[QC] Raw sample quality control", type = "step")

    qc_raw <- explore_bulk_rnaseq(
      counts = counts,
      metadata = metadata,
      gene_col = gene_col,
      feature_col = feature_col,
      sample_col = sample_col,
      markers = qc_markers,
      output_dir = file.path(output_dir, paste0(qc_output_prefix, "_raw"))
    )

    if (isTRUE(remove_flagged)) {

      log("[QC] Removing flagged samples", type = "step")

      scol <- qc_sample_col %||% sample_col[1]

      counts_for_removal <- .deggo_counts_with_gene_col(
        counts = counts,
        gene_col = gene_col
      )

      cleaned <- remove_flagged_samples(
        counts = counts_for_removal,
        metadata = metadata,
        qc_table = qc_raw$qc,
        sample_col = scol,
        remove_col = "recommend_remove",
        gene_cols = unique(c(gene_col, feature_col)),
        verbose = TRUE
      )

      counts <- cleaned$counts
      metadata <- cleaned$metadata

      log("[QC] Re-running QC after sample removal", type = "step")

      qc_clean <- explore_bulk_rnaseq(
        counts = counts,
        metadata = metadata,
        gene_col = gene_col,
        feature_col = feature_col,
        sample_col = sample_col,
        markers = qc_markers,
        output_dir = file.path(output_dir, paste0(qc_output_prefix, "_clean"))
      )
    }
  }

  # ---------------------------------------------------------- #
  # Input preparation
  # ---------------------------------------------------------- #

  log("[1/11] Matching counts and metadata", type = "step")

  if (isTRUE(prepare_input)) {

    counts_input <- .deggo_counts_with_gene_col(
      counts = counts,
      gene_col = gene_col
    )

    prep <- prepare_counts_metadata(
      counts = counts_input,
      metadata = metadata,
      gene_col = gene_col,
      feature_col = feature_col,
      sample_col = sample_col,
      verbose = FALSE
    )

    counts <- prep$counts
    metadata <- prep$metadata

  } else {

    counts <- as.matrix(counts)
    suppressWarnings(storage.mode(counts) <- "numeric")
    counts <- round(counts)
    storage.mode(counts) <- "integer"

    scol <- sample_col[sample_col %in% colnames(metadata)][1]

    if (!is.na(scol)) {
      rownames(metadata) <- metadata[[scol]]
    } else if ("sample" %in% colnames(metadata)) {
      rownames(metadata) <- metadata$sample
    } else {
      stop("No sample column found in metadata.", call. = FALSE)
    }
  }

  if (!is.matrix(counts)) {
    stop("counts must be a matrix after input preparation.", call. = FALSE)
  }

  # ---------------------------------------------------------- #
  # Marker validation
  # ---------------------------------------------------------- #

  if (!is.null(marker_sets)) {

    log("[QC] Marker-based biological validation", type = "step")

    marker_genes <- unique(unlist(marker_sets))

    plot_gene_heatmap(
      counts = counts,
      metadata = metadata,
      genes = marker_genes,
      gene_col = gene_col[1],
      feature_col = feature_col[1],
      sample_col = sample_col[1],
      annotation_cols = intersect(
        c("tissue", "sex", "treatment", "condition"),
        colnames(metadata)
      ),
      order_by = intersect(
        c("tissue", "sex", "treatment", "condition"),
        colnames(metadata)
      ),
      output_dir = file.path(output_dir, paste0(qc_output_prefix, "_marker_check")),
      filename = "Marker_Check_Heatmap",
      main = "Marker-based sample validation"
    )

    marker_check <- marker_score_check(
      counts = counts,
      metadata = metadata,
      marker_sets = marker_sets,
      sample_col = sample_col[1],
      group_col = "tissue",
      feature_col = feature_col[1]
    )
  }

  # ---------------------------------------------------------- #
  # Annotation and validation
  # ---------------------------------------------------------- #

  log("[2/11] Loading annotation database", type = "step")
  orgdb <- .get_orgdb(organism = organism, orgdb = orgdb)

  log("[3/11] Validating inputs", type = "step")

  sample_col_use <- sample_col[sample_col %in% colnames(metadata)][1]

  if (is.na(sample_col_use)) {
    stop("No valid sample column found in metadata.", call. = FALSE)
  }

  if (is.null(rownames(metadata)) || any(rownames(metadata) == "")) {
    rownames(metadata) <- metadata[[sample_col_use]]
  }

  if (analysis_mode == "single") {

    if (is.null(design_formula)) {
      stop("'design_formula' is required in single mode.", call. = FALSE)
    }

    if (is.null(contrast)) {
      stop("'contrast' is required in single mode.", call. = FALSE)
    }

    required_metadata_cols <- unique(c(
      sample_col_use,
      all.vars(design_formula),
      contrast[1]
    ))

    missing_cols <- setdiff(required_metadata_cols, colnames(metadata))

    if (length(missing_cols) > 0) {
      stop(
        "Metadata is missing required column(s) for single mode: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    contrast_var <- contrast[1]
    contrast_levels <- contrast[2:3]

    if (!all(contrast_levels %in% metadata[[contrast_var]])) {
      stop(
        "Contrast levels not found in metadata column '", contrast_var, "': ",
        paste(setdiff(contrast_levels, metadata[[contrast_var]]), collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (analysis_mode == "pairwise") {

    if (method != "DESeq2") {
      stop("Pairwise mode is currently available only with DESeq2.", call. = FALSE)
    }

    if (is.null(pairwise_group_cols)) {
      stop("'pairwise_group_cols' is required for pairwise mode.", call. = FALSE)
    }

    required_metadata_cols <- unique(c(
      sample_col_use,
      pairwise_group_cols
    ))

    missing_cols <- setdiff(required_metadata_cols, colnames(metadata))

    if (length(missing_cols) > 0) {
      stop(
        "Metadata is missing required column(s) for pairwise mode: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    metadata[[pairwise_contrast_col]] <- apply(
      metadata[, pairwise_group_cols, drop = FALSE],
      1,
      paste,
      collapse = "_"
    )

    metadata[[pairwise_contrast_col]] <- make.names(
      metadata[[pairwise_contrast_col]]
    )
  }

  sample_ids <- colnames(counts)

  metadata <- metadata[match(sample_ids, rownames(metadata)), , drop = FALSE]
  rownames(metadata) <- sample_ids
  metadata[[sample_col_use]] <- sample_ids

  # ---------------------------------------------------------- #
  # Cleaning and filtering
  # ---------------------------------------------------------- #

  log("[4/11] Cleaning gene identifiers", type = "step")

  counts <- clean_ensembl_ids(counts)

  log("[5/11] Filtering low-expression genes", type = "step")

  counts <- preprocess_counts(
    counts = counts,
    metadata = metadata,
    filter_method = filter_method,
    min_count = min_count,
    min_samples = min_samples,
    min_total = min_total
  )

  colnames(counts) <- sample_ids
  rownames(metadata) <- sample_ids
  metadata[[sample_col_use]] <- sample_ids

  .save_clean_input_files(
    counts = counts,
    metadata = metadata,
    repro_dir = repro_dir,
    save_clean_inputs = save_clean_inputs
  )

  # ---------------------------------------------------------- #
  # Sample QC after filtering
  # ---------------------------------------------------------- #

  log("[6/11] Sample quality control", type = "step")

  sample_qc <- run_sample_qc(
    counts = counts,
    metadata = metadata,
    txtsize = txtsize,
    output_dir = dirs$qc,
    annotation_cols = intersect(
      c("condition", "treatment", "sex", "tissue"),
      colnames(metadata)
    )
  )

  # ---------------------------------------------------------- #
  # Differential expression
  # ---------------------------------------------------------- #

  if (analysis_mode == "pairwise") {

    de_results <- .run_deggo_pairwise(
      counts = counts,
      metadata = metadata,
      orgdb = orgdb,
      dirs = dirs,
      ontology = ontology,
      padj_cutoff = padj_cutoff,
      logfc_cutoff = logfc_cutoff,
      top_n_heatmap = top_n_heatmap,
      top_n_labels = top_n_labels,
      txtsize = txtsize,
      pairwise_group_cols = pairwise_group_cols,
      pairwise_contrast_col = pairwise_contrast_col,
      pairwise_contrasts = pairwise_contrasts,
      pairwise_mode = pairwise_mode,
      log = log
    )

  } else {

    de_results <- .run_deggo_single(
      counts = counts,
      metadata = metadata,
      method = method,
      design_formula = design_formula,
      contrast = contrast,
      orgdb = orgdb,
      dirs = dirs,
      ontology = ontology,
      padj_cutoff = padj_cutoff,
      logfc_cutoff = logfc_cutoff,
      top_n_heatmap = top_n_heatmap,
      top_n_labels = top_n_labels,
      txtsize = txtsize,
      log = log
    )
  }

  # ---------------------------------------------------------- #
  # Final metadata and report/PPTX compatibility
  # ---------------------------------------------------------- #

  de_results$counts <- counts
  de_results$metadata <- de_results$metadata %||% metadata
  de_results$output_dir <- output_dir
  de_results$output_dirs <- dirs
  de_results$version <- deggo_version
  de_results$sample_qc <- sample_qc
  de_results$qc_raw <- qc_raw
  de_results$qc_clean <- qc_clean
  de_results$marker_check <- marker_check

  if (is.null(de_results$summary)) {
    stop("DEGgo internal error: summary table was not created.", call. = FALSE)
  }

  .write_session_info(output_dir)

  de_results$run_params <- .make_run_params(
    deggo_version = deggo_version,
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
    txtsize = txtsize,
    prepare_input = prepare_input,
    gene_col = gene_col,
    feature_col = feature_col,
    sample_col = sample_col,
    design_formula = design_formula,
    pairwise_group_cols = pairwise_group_cols,
    pairwise_contrast_col = pairwise_contrast_col,
    pairwise_mode = pairwise_mode,
    output_dir = output_dir,
    repro_dir = repro_dir
  )

  .save_repro(
    res = de_results,
    repro_dir = repro_dir,
    save_reproducibility = save_reproducibility
  )

  # ---------------------------------------------------------- #
  # Write run-level summary
  # ---------------------------------------------------------- #

  utils::write.table(
    de_results$summary,
    file.path(output_dir, paste0(analysis_mode, "_summary.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # ---------------------------------------------------------- #
  # Write output manifest
  # ---------------------------------------------------------- #

  de_results$output_manifest <- .write_deggo_manifest(
    output_dir = output_dir,
    dirs = dirs,
    analysis_mode = analysis_mode
  )

  # ---------------------------------------------------------- #
  # Organize run-level files BEFORE report/PPTX
  # ---------------------------------------------------------- #

  run_files <- .deggo_organize_run_files(
    output_dir = output_dir,
    analysis_mode = analysis_mode
  )

  de_results$summary_file <- run_files$summary_file
  de_results$session_file <- run_files$session_file
  de_results$manifest_file <- run_files$manifest_file

  # ---------------------------------------------------------- #
  # Report
  # ---------------------------------------------------------- #

  de_results$report_files <- NULL

  if (isTRUE(generate_report)) {

    log("[REPORT] Generating DEGgo report", type = "step")

    de_results$report_files <- .deggo_report(
      res = de_results,
      output_dir = output_dir,
      generate_report = generate_report,
      report_formats = report_formats,
      report_template = report_template
    )

  }

  # ---------------------------------------------------------- #
  # PowerPoint
  # ---------------------------------------------------------- #

  de_results$pptx_file <- NULL

  if (isTRUE(generate_pptx)) {

    report_dir <- file.path(output_dir, "DEGgo_Report")
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

    if (is.null(pptx_file)) {
      pptx_file <- file.path(report_dir, "DEGgo_Report.pptx")
    }

    de_results$pptx_file <- generate_deggo_pptx(
      results = de_results,
      output_file = pptx_file
    )

  }


  log(
    paste0("==== DEGgo ", toupper(analysis_mode), " ANALYSIS COMPLETE ===="),
    type = "done",
    duration = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  )

  de_results
}
