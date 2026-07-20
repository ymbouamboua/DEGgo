
# ============================================================ #
# Main run_deggo
# ============================================================ #
#' Run DEGgo bulk RNA-seq downstream analysis
#'
#' Automated bulk RNA-seq differential expression workflow including QC,
#' preprocessing, differential expression, annotation, visualization,
#' GO enrichment, reporting, and reproducibility exports.
#' @param counts Raw count matrix or data frame.
#' @param metadata Sample metadata data frame.
#' @param project_name Optional project name shown in HTML, PDF and PowerPoint reports.
#' @param gene_col Candidate gene identifier column names.
#' @param feature_col Candidate gene symbol/name column names.
#' @param sample_col Candidate sample identifier column names.
#' @param prepare_input Logical. Prepare and match input tables.
#' @param raw_qc Logical. Run raw count QC.
#' @param remove_flagged Logical. Remove QC-flagged samples.
#' @param qc_markers Optional marker genes for QC.
#' @param marker_sets Optional named marker gene sets.
#' @param qc_sample_col Optional sample column for QC.
#' @param qc_output_prefix QC output prefix.
#' @param output_dir Output directory.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param logfc_cutoff Absolute log2 fold-change cutoff.
#' @param top_n_heatmap Number of genes shown in heatmaps.
#' @param top_n_labels Number of genes labelled in volcano plots.
#' @param min_expr_count Minimum expression count for clean DEG filtering.
#' @param min_expr_samples Minimum samples passing expression threshold.
#' @param min_prevalence Minimum prevalence for clean DEG filtering.
#' @param max_sample_fraction Maximum single-sample fraction.
#' @param max_group_sample_fraction Maximum fraction of group expression contributed by a single sample.
#' @param min_group_mean Minimum group mean expression.
#' @param expr_filter_groups Grouping variables for expression filtering.
#' @param clean_deg_tables Logical. Apply post-DEG expression cleaning.
#' @param ontology GO ontology.
#' @param organism Organism name.
#' @param orgdb Optional custom OrgDb object.
#' @param method Differential expression method.
#' @param analysis_mode Single or pairwise analysis mode.
#' @param contrast Contrast vector for single analysis.
#' @param design_formula Design formula.
#' @param pairwise_group_cols Metadata columns used to build pairwise groups.
#' @param pairwise_contrast_col Name of pairwise contrast column.
#' @param pairwise_contrasts Optional named list of pairwise contrasts.
#' @param filter_method Gene filtering method.
#' @param pairwise_mode Pairwise contrast generation mode.
#' @param min_count Minimum count threshold.
#' @param min_samples Minimum number of samples passing `min_count`.
#' @param min_total Minimum total count.
#' @param generate_report Logical. Generate report.
#' @param report_formats Report formats.
#' @param report_template Optional report template path.
#' @param generate_pptx Logical. Generate PowerPoint.
#' @param pptx_file Optional PowerPoint output file.
#' @param save_reproducibility Logical. Save reproducibility bundle.
#' @param save_clean_inputs Logical. Save cleaned input tables.
#' @param txtsize Base text size.
#' @param seed Random seed.
#'
#' @return A DEGgo results object.
#'
#' @export
#'
run_deggo <- function(
    counts,
    metadata,
    project_name = NULL,
    gene_col = c(
      "gene_id", "GeneID", "gene", "Gene",
      "ENSEMBL", "ensembl", "ensembl_id"
    ),
    feature_col = c(
      "gene_name", "SYMBOL", "symbol",
      "gene_symbol", "external_gene_name"
    ),
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
    min_expr_count = 20,
    min_expr_samples = 3,
    min_prevalence = 0.6,
    max_sample_fraction = 0.45,
    max_group_sample_fraction = 0.45,
    min_group_mean = 10,
    min_group_median = 20,
    max_group_cv = NULL,
    expr_filter_groups = "auto",
    clean_deg_tables = TRUE,
    ontology = c("BP", "MF", "CC"),
    organism = c("human", "mouse", "rat", "custom"),
    orgdb = NULL,
    method = c("DESeq2", "edgeR", "limma"),
    analysis_mode = c("single", "pairwise"),
    contrast = NULL,
    design_formula = ~ condition,
    pairwise_group_cols = NULL,
    pairwise_contrast_col = "group",
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
    heatmap_annotation_cols = "auto",
    palette = "default",
    seed = 123
) {
  set.seed(seed)

  method <- match.arg(method)
  filter_method <- match.arg(filter_method)
  organism <- match.arg(organism)
  analysis_mode <- match.arg(analysis_mode)
  pairwise_mode <- match.arg(pairwise_mode)
  ontology <- match.arg(ontology)
  palette <- match.arg(
    palette,
    choices = c(
      "default",
      "nature",
      "jama",
      "nejm",
      "lancet",
      "viridis",
      "okabe"
    )
  )

  log <- .deggo_msg(verbose = TRUE, prefix = "DEGgo")
  t_start <- Sys.time()

  log("==== STARTING DEGgo ANALYSIS ====", type = "header")

  # ---------------------------------------------------------- #
  # 1. Initialize run
  # ---------------------------------------------------------- #

  cfg <- .deggo_initialize(
    output_dir = output_dir,
    method = method,
    analysis_mode = analysis_mode,
    organism = organism,
    ontology = ontology,
    report_template = report_template,
    log = log
  )

  project_name <- project_name %||% "DEGgo RNA-seq analysis"
  project_name <- as.character(project_name)[1]

  output_dir <- cfg$output_dir
  dirs <- cfg$dirs
  repro_dir <- cfg$repro_dir
  deggo_version <- cfg$deggo_version
  report_template <- cfg$report_template

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  # ---------------------------------------------------------- #
  # 2. Raw QC
  # ---------------------------------------------------------- #

  qc <- .deggo_raw_qc(
    counts = counts,
    metadata = metadata,
    raw_qc = raw_qc,
    remove_flagged = remove_flagged,
    gene_col = gene_col,
    feature_col = feature_col,
    sample_col = sample_col,
    qc_markers = qc_markers,
    qc_sample_col = qc_sample_col,
    qc_output_prefix = qc_output_prefix,
    output_dir = output_dir
  )

  counts <- qc$counts
  metadata <- qc$metadata

  # ---------------------------------------------------------- #
  # 3. Prepare inputs
  # ---------------------------------------------------------- #

  prep <- .deggo_prepare_inputs(
    counts = counts,
    metadata = metadata,
    prepare_input = prepare_input,
    gene_col = gene_col,
    feature_col = feature_col,
    sample_col = sample_col,
    filter_method = filter_method,
    min_count = min_count,
    min_samples = min_samples,
    min_total = min_total,
    repro_dir = repro_dir,
    save_clean_inputs = save_clean_inputs,
    txtsize = txtsize,
    dirs = dirs
  )

  counts <- prep$counts
  metadata <- prep$metadata
  sample_qc <- prep$sample_qc

  # ---------------------------------------------------------- #
  # 4. Marker validation
  # ---------------------------------------------------------- #

  marker_check <- NULL

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
  # 5. Annotation database
  # ---------------------------------------------------------- #

  log("[2/11] Loading annotation database", type = "step")

  orgdb <- .get_orgdb(
    organism = organism,
    orgdb = orgdb
  )

  # ---------------------------------------------------------- #
  # 6. Validate analysis inputs
  # ---------------------------------------------------------- #

  log("[3/11] Validating inputs", type = "step")

  sample_col_use <- sample_col[sample_col %in% colnames(metadata)][1]

  if (is.na(sample_col_use)) {
    stop("No valid sample column found in metadata.", call. = FALSE)
  }

  if (is.null(rownames(metadata)) || any(rownames(metadata) == "")) {
    rownames(metadata) <- as.character(metadata[[sample_col_use]])
  }

  sample_ids <- colnames(counts)

  metadata <- metadata[match(sample_ids, rownames(metadata)), , drop = FALSE]
  rownames(metadata) <- sample_ids
  metadata[[sample_col_use]] <- sample_ids

  if (analysis_mode == "single") {

    if (is.null(design_formula)) {
      stop("'design_formula' is required in single mode.", call. = FALSE)
    }

    if (is.null(contrast)) {
      stop("'contrast' is required in single mode.", call. = FALSE)
    }

    required_cols <- unique(c(
      sample_col_use,
      all.vars(design_formula),
      contrast[1]
    ))

    missing_cols <- setdiff(required_cols, colnames(metadata))

    if (length(missing_cols)) {
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
        "Contrast levels not found in metadata column '",
        contrast_var,
        "': ",
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

    required_cols <- unique(c(
      sample_col_use,
      pairwise_group_cols
    ))

    missing_cols <- setdiff(required_cols, colnames(metadata))

    if (length(missing_cols)) {
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

  # ---------------------------------------------------------- #
  # 7. Differential expression only
  # ---------------------------------------------------------- #

  if (analysis_mode == "single") {

    de_results <- .run_deggo_single(
      counts = counts,
      metadata = metadata,
      method = method,
      design_formula = design_formula,
      contrast = contrast,
      orgdb = orgdb,
      padj_cutoff = padj_cutoff,
      logfc_cutoff = logfc_cutoff,
      log = log
    )

  } else {

    de_results <- .run_deggo_pairwise(
      counts = counts,
      metadata = metadata,
      orgdb = orgdb,
      pairwise_group_cols = pairwise_group_cols,
      pairwise_contrast_col = pairwise_contrast_col,
      pairwise_contrasts = pairwise_contrasts,
      pairwise_mode = pairwise_mode,
      padj_cutoff = padj_cutoff,
      logfc_cutoff = logfc_cutoff,
      log = log
    )
  }

  # ---------------------------------------------------------- #
  # 8. Attach run-level objects
  # ---------------------------------------------------------- #

  de_results$counts <- counts
  de_results$metadata <- de_results$metadata %||% metadata
  de_results$project_name <- project_name
  de_results$output_dir <- output_dir
  de_results$output_dirs <- dirs
  de_results$version <- deggo_version
  de_results$sample_qc <- sample_qc
  de_results$qc_raw <- qc$qc_raw
  de_results$qc_clean <- qc$qc_clean
  de_results$marker_check <- marker_check

  # ---------------------------------------------------------- #
  # 9. Optional DEG cleaning
  # ---------------------------------------------------------- #

  if (isTRUE(clean_deg_tables)) {

    log("[POST] Cleaning significant DEG tables", type = "step")

    de_results <- .deggo_clean_deg_tables(
      de_results = de_results,
      counts = counts,
      metadata = de_results$metadata,
      min_expr_count = min_expr_count,
      min_expr_samples = min_expr_samples,
      min_prevalence = min_prevalence,
      max_sample_fraction = max_sample_fraction,
      max_group_sample_fraction = max_group_sample_fraction,
      min_group_mean = min_group_mean,
      min_group_median = min_group_median,
      max_group_cv = max_group_cv,
      expr_filter_groups = expr_filter_groups,
      log = log
    )

  } else {

    de_results$sig_deg_clean <- de_results$sig_deg
  }

  de_results$summary <- .deggo_build_summary(
    results = de_results$results,
    sig_deg = de_results$sig_deg_clean,
    go_results = NULL,
    method = method,
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff
  )


  # ---------------------------------------------------------- #
  # 10. PCA variables
  # ---------------------------------------------------------- #

  if (analysis_mode == "single") {
    pca_vars <- unique(c(
      all.vars(design_formula),
      contrast[1]
    ))
  } else {
    pca_vars <- unique(pairwise_group_cols)
  }

  pca_vars <- setdiff(pca_vars, sample_col_use)
  pca_vars <- intersect(pca_vars, colnames(metadata))

  # ---------------------------------------------------------- #
  # 11. Visualizations
  # ---------------------------------------------------------- #

  if (analysis_mode == "single") {
    plot_order_vars <- unique(c(
      setdiff(all.vars(design_formula), contrast[1]),
      contrast[1]
    ))
  } else {
    plot_order_vars <- unique(pairwise_group_cols)
  }

  plot_order_vars <- intersect(plot_order_vars, colnames(metadata))

  de_results <- .deggo_make_plots(
    de_results = de_results,
    counts = counts,
    metadata = de_results$metadata,
    dirs = dirs,
    analysis_mode = analysis_mode,
    method = method,
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff,
    top_n_heatmap = top_n_heatmap,
    top_n_labels = top_n_labels,
    txtsize = txtsize,
    log = log,
    contrast = contrast,
    pairwise_contrast_col = pairwise_contrast_col,
    pairwise_contrasts = pairwise_contrasts,
    sample_col = sample_col_use,
    pca_vars = pca_vars,
    plot_order_vars = plot_order_vars,
    design_formula = design_formula,
    pairwise_group_cols = pairwise_group_cols,
    heatmap_annotation_cols = heatmap_annotation_cols,
    palette = palette
  )

  # ---------------------------------------------------------- #
  # 12. GO enrichment
  # ---------------------------------------------------------- #

  de_results <- .deggo_make_go(
    de_results = de_results,
    dirs = dirs,
    ontology = ontology,
    orgdb = orgdb,
    txtsize = txtsize,
    log = log
  )

  summary_go <- .deggo_build_summary(
    results = de_results$results,
    sig_deg = de_results$sig_deg_clean,
    go_results = de_results$go_results,
    method = method,
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff
  )

  de_results$summary$go_terms <- summary_go$go_terms[
    match(de_results$summary$comparison, summary_go$comparison)
  ]

  # ---------------------------------------------------------- #
  # 13. Run parameters
  # ---------------------------------------------------------- #

  de_results$run_params <- .make_run_params(
    deggo_version = deggo_version,
    project_name = project_name,
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

  # ---------------------------------------------------------- #
  # 14. Export
  # ---------------------------------------------------------- #

  de_results <- .deggo_export_all(
    de_results = de_results,
    dirs = dirs,
    output_dir = output_dir,
    analysis_mode = analysis_mode,
    repro_dir = repro_dir,
    save_reproducibility = save_reproducibility
  )

  de_results$output_manifest <- .write_deggo_manifest(
    output_dir = output_dir,
    dirs = dirs,
    analysis_mode = analysis_mode
  )

  run_files <- .deggo_organize_run_files(
    output_dir = output_dir,
    analysis_mode = analysis_mode
  )

  de_results$summary_file <- run_files$summary_file
  de_results$session_file <- run_files$session_file
  de_results$manifest_file <- run_files$manifest_file

  # ---------------------------------------------------------- #
  # 15. Report
  # ---------------------------------------------------------- #

  de_results$report_files <- NULL

  if (isTRUE(generate_report)) {

    log("[REPORT] Generating DEGgo report", type = "step")

    de_results$report_files <- .deggo_report(
      res = de_results,
      output_dir = output_dir,
      generate_report = generate_report,
      report_formats = report_formats,
      report_template = report_template,
      project_name = project_name
    )
  }

  # ---------------------------------------------------------- #
  # 16. PowerPoint
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
      output_file = pptx_file,
      project_name = project_name
    )
  }

  log(
    paste0("==== DEGgo ", toupper(analysis_mode), " ANALYSIS COMPLETE ===="),
    type = "done",
    duration = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  )

  de_results
}
