# ======================================================= #
# DEGgo.R
# Automated bulk RNA-seq downstream analysis
# ======================================================= #

#' Run DEGgo bulk RNA-seq downstream analysis
#'
#' Run an automated bulk RNA-seq downstream analysis workflow including
#' input preparation, sample matching, differential expression analysis,
#' result annotation, visualization, heatmap generation,
#' PCA analysis, Gene Ontology enrichment, and optional HTML/PDF report
#' generation.
#'
#' @param counts Raw count table or matrix.
#' @param metadata Sample metadata data frame.
#' @param gene_col Character vector of possible gene identifier columns.
#' @param feature_col Character vector of possible feature annotation columns.
#' @param sample_col Character vector of possible sample identifier columns.
#' @param prepare_input Logical. If TRUE, automatically matches counts and metadata.
#' @param output_dir Output directory. If NULL, a dated directory is created.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param logfc_cutoff Absolute log2 fold-change cutoff.
#' @param top_n_heatmap Number of top genes for heatmap.
#' @param top_n_labels Number of top genes labeled in volcano plot.
#' @param ontology GO ontology: BP, MF, or CC.
#' @param organism Organism: human, mouse, or rat.
#' @param orgdb Optional OrgDb object.
#' @param method DE method: DESeq2, edgeR, or limma.
#' @param analysis_mode Analysis mode: single or pairwise.
#' @param contrast Optional contrast for single analysis.
#' @param design_formula Design formula.
#' @param pairwise_group_cols Metadata columns for pairwise groups.
#' @param pairwise_contrast_col Name of pairwise contrast column.
#' @param pairwise_contrasts Optional named list of contrasts.
#' @param filter_method Filtering method: count, cpm, or none.
#' @param pairwise_auto Automatically generate pairwise contrasts.
#' @param pairwise_mode Pairwise mode: all, within_first, or within_second.
#' @param min_count Minimum count for filtering.
#' @param min_samples Minimum samples for filtering.
#' @param min_total Minimum total count for filtering.
#' @param generate_report Logical. Generate DEGgo report.
#' @param report_formats Character vector. Report formats: html and/or pdf.
#' @param report_template Optional path to report R Markdown template.
#' @param seed Random seed.
#'
#' @return A list containing DEG results, plots, GO results, reports, and output paths.
#' @export
run_deggo <- function(
    counts,
    metadata,
    gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_id"),
    feature_col = c("gene_name", "SYMBOL", "symbol", "gene_symbol", "external_gene_name"),
    sample_col = c("sample", "Sample", "SAMPLE"),
    prepare_input = TRUE,
    output_dir = NULL,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25,
    top_n_heatmap = 50,
    top_n_labels = 10,
    ontology = c("BP", "MF", "CC"),
    organism = c("human", "mouse", "rat"),
    orgdb = NULL,
    method = c("DESeq2", "edgeR", "limma"),
    analysis_mode = c("single", "pairwise"),
    contrast = NULL,
    design_formula = ~ condition,
    pairwise_group_cols = NULL,
    pairwise_contrast_col = "comparison_group",
    pairwise_contrasts = NULL,
    filter_method = c("count", "cpm", "none"),
    pairwise_auto = TRUE,
    pairwise_mode = c("all", "within_first", "within_second"),
    min_count = 5,
    min_samples = 2,
    min_total = 10,
    generate_report = TRUE,
    report_formats = "html",
    report_template = NULL,
    seed = 123
) {

  set.seed(seed)

  method <- match.arg(method)
  filter_method <- match.arg(filter_method)
  organism <- match.arg(organism)
  analysis_mode <- match.arg(analysis_mode)
  pairwise_mode <- match.arg(pairwise_mode)
  ontology <- match.arg(ontology)

  log <- .msg(verbose = TRUE, prefix = "DEGgo")
  t_start <- Sys.time()

  log("==== STARTING DEGgo ANALYSIS ====", type = "header")

  if (is.null(output_dir)) {
    output_dir <- "DEGgo_results"
  }

  output_dir <- file.path(
    output_dir,
    paste0(
      "DEGgo_",
      format(Sys.Date(), "%Y_%m_%d"),
      "_",
      method,
      "_",
      analysis_mode
    )
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  dirs <- .deggo_dirs(output_dir, analysis_mode)

  deggo_version <- tryCatch(
    as.character(utils::packageVersion("DEGgo")),
    error = function(e) "development"
  )

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  log("[1/11] Matching counts and metadata", type = "step")

  if (isTRUE(prepare_input)) {

    prep <- prepare_counts_metadata(
      counts = counts,
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

    if (!is.numeric(counts)) {
      stop("counts must be numeric when prepare_input = FALSE.", call. = FALSE)
    }

    counts <- round(counts)
    storage.mode(counts) <- "integer"

    if ("sample" %in% colnames(metadata)) {
      rownames(metadata) <- metadata$sample
    }
  }

  if (!is.matrix(counts)) {
    stop("counts must be a matrix after input preparation.", call. = FALSE)
  }

  if (analysis_mode == "pairwise") {

    if (method != "DESeq2") {
      stop("Pairwise mode is currently available only with DESeq2.", call. = FALSE)
    }

    if (is.null(pairwise_group_cols)) {
      stop("'pairwise_group_cols' is required for pairwise mode.", call. = FALSE)
    }

    missing_pairwise_cols <- setdiff(pairwise_group_cols, colnames(metadata))

    if (length(missing_pairwise_cols) > 0) {
      stop(
        "Missing pairwise grouping columns: ",
        paste(missing_pairwise_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  log("[2/11] Loading annotation database", type = "step")

  orgdb <- .get_orgdb(
    organism = organism,
    orgdb = orgdb
  )

  log("[3/11] Validating inputs", type = "step")

  if (!"condition" %in% colnames(metadata)) {
    if (analysis_mode == "pairwise") {
      metadata$condition <- apply(
        metadata[, pairwise_group_cols, drop = FALSE],
        1,
        paste,
        collapse = "_"
      )
    } else {
      stop("Metadata must contain 'condition'.", call. = FALSE)
    }
  }

  validate_inputs(
    counts = counts,
    metadata = metadata,
    condition_col = "condition"
  )

  sample_ids <- colnames(counts)

  metadata <- metadata[
    match(sample_ids, rownames(metadata)),
    ,
    drop = FALSE
  ]

  rownames(metadata) <- sample_ids

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

  # ======================================================= #
  # Pairwise mode
  # ======================================================= #

  if (analysis_mode == "pairwise") {

    log("[6/11] Running pairwise DESeq2 contrasts", type = "step")

    de_results <- run_deseq2_pairwise(
      counts = counts,
      metadata = metadata,
      group_cols = pairwise_group_cols,
      contrast_col = pairwise_contrast_col,
      contrast_list = pairwise_contrasts,
      pairwise_mode = pairwise_mode
    )

    de_results$dds <- .annotate_dds(
      dds = de_results$dds,
      orgdb = orgdb
    )

    log("[7/11] PCA analysis", type = "step")

    de_results$pca <- list(
      sample = plot_pca(
        de_results$dds,
        de_results$metadata,
        dirs$pca,
        "PCA_by_sample",
        color_by = "sample",
        title = "PCA by sample"
      ),
      tissue = plot_pca(
        de_results$dds,
        de_results$metadata,
        dirs$pca,
        "PCA_by_tissue",
        color_by = "tissue",
        title = "PCA by tissue"
      ),
      treatment = plot_pca(
        de_results$dds,
        de_results$metadata,
        dirs$pca,
        "PCA_by_treatment",
        color_by = "treatment",
        title = "PCA by treatment"
      ),
      sex = plot_pca(
        de_results$dds,
        de_results$metadata,
        dirs$pca,
        "PCA_by_sex",
        color_by = "sex",
        title = "PCA by sex"
      ),
      tissue_treatment = plot_pca(
        de_results$dds,
        de_results$metadata,
        dirs$pca,
        "PCA_tissue_treatment",
        color_by = "tissue",
        shape_by = "treatment",
        title = "PCA tissue + treatment"
      )
    )

    log("[8/11] Annotating DE tables", type = "step")

    de_results$results <- lapply(
      de_results$results,
      annotate_de_results,
      orgdb = orgdb,
      padj_cutoff = padj_cutoff,
      logfc_cutoff = logfc_cutoff
    )

    de_results$sig_deg <- list()
    de_results$volcano_plots <- list()
    de_results$heatmaps <- list()
    de_results$go_results <- list()
    de_results$go_plots <- list()

    log("[9/11] Exporting plots, heatmaps and GO", type = "step")

    vsd_pairwise <- DESeq2::vst(de_results$dds, blind = FALSE)

    for (nm in names(de_results$results)) {

      log("Processing comparison:", nm, type = "info")

      res_df_i <- de_results$results[[nm]]

      sig_i <- res_df_i[
        !is.na(res_df_i$padj) &
          res_df_i$padj < padj_cutoff &
          abs(res_df_i$log2FoldChange) > logfc_cutoff,
        ,
        drop = FALSE
      ]

      de_results$sig_deg[[nm]] <- sig_i

      utils::write.table(
        res_df_i,
        file.path(dirs$results, paste0(nm, ".tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      utils::write.table(
        sig_i,
        file.path(dirs$significant, paste0(nm, "_significant.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      de_results$volcano_plots[[nm]] <- plot_volcano(
        res_df = res_df_i,
        top_n_labels = top_n_labels,
        output_dir = dirs$volcano,
        filename = paste0(nm, "_Volcano_Plot"),
        title = nm,
        logfc_cutoff = logfc_cutoff,
        padj_cutoff = padj_cutoff
      )

      de_results$heatmaps[[nm]] <- plot_heatmap(
        vsd = vsd_pairwise,
        res_df = res_df_i,
        metadata = de_results$metadata,
        sample_subset = de_results$samples[[nm]],
        top_n_heatmap = top_n_heatmap,
        padj_cutoff = padj_cutoff,
        output_dir = dirs$heatmaps,
        main = nm,
        filename = paste0(nm, "_Heatmap"),
        annotation_cols = intersect(
          c("treatment", "sex", "tissue"),
          colnames(de_results$metadata)
        ),
        order_by = intersect(
          c("tissue", "sex", "treatment"),
          colnames(de_results$metadata)
        )
      )

      de_results$go_results[[nm]] <- run_go_enrichment(
        sig_deg = sig_i,
        comparison = nm,
        ontology = ontology,
        output_dir = dirs$go,
        orgdb = orgdb
      )

      go_df_i <- de_results$go_results[[nm]]$go_results

      if (!is.null(go_df_i) && nrow(go_df_i) > 0) {

        de_results$go_plots[[nm]] <- plot_go_terms(
          go_df = go_df_i,
          comparison = paste0(nm, " GO enrichment"),
          top_n = 10,
          style = "bw"
        )

        ggplot2::ggsave(
          filename = file.path(dirs$go, paste0(nm, "_GO_terms.png")),
          plot = de_results$go_plots[[nm]],
          width = 8,
          height = 6,
          dpi = 300
        )

        ggplot2::ggsave(
          filename = file.path(dirs$go, paste0(nm, "_GO_terms.pdf")),
          plot = de_results$go_plots[[nm]],
          width = 8,
          height = 6
        )
      }
    }

    de_results$go_merged <- .merge_pairwise_go(
      go_results = de_results$go_results,
      go_dir = dirs$go,
      ontology = ontology
    )

    log("[10/11] Summarizing pairwise results", type = "step")

    de_results$summary <- do.call(
      rbind,
      lapply(names(de_results$results), function(nm) {

        x <- de_results$results[[nm]]

        data.frame(
          comparison = nm,
          total_genes = nrow(x),
          significant = sum(
            !is.na(x$padj) &
              x$padj < padj_cutoff &
              abs(x$log2FoldChange) > logfc_cutoff,
            na.rm = TRUE
          ),
          up = sum(
            !is.na(x$padj) &
              x$padj < padj_cutoff &
              x$log2FoldChange > logfc_cutoff,
            na.rm = TRUE
          ),
          down = sum(
            !is.na(x$padj) &
              x$padj < padj_cutoff &
              x$log2FoldChange < -logfc_cutoff,
            na.rm = TRUE
          ),
          stringsAsFactors = FALSE
        )
      })
    )

    utils::write.table(
      de_results$summary,
      file.path(output_dir, "pairwise_summary.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    de_results$output_dir <- output_dir
    de_results$output_dirs <- dirs
    de_results$version <- deggo_version

    .safe_write_session_info(output_dir)

    log("[11/11] Generating DEGgo report", type = "step")

    de_results$run_params <- list(
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
      output_dir = output_dir
    )

    de_results$report_files <- NULL

    if (isTRUE(generate_report)) {
      de_results$report_files <- generate_deggo_report(
        results = de_results,
        output_dir = output_dir,
        formats = report_formats,
        report_template = report_template
      )
    }

    log(
      "==== DEGgo PAIRWISE ANALYSIS COMPLETE ====",
      type = "done",
      duration = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    )

    return(de_results)
  }

  # ======================================================= #
  # Single mode
  # ======================================================= #

  log("[6/11] Running differential expression analysis", type = "step")

  de_results <- run_de(
    counts = counts,
    metadata = metadata,
    method = method,
    design_formula = design_formula,
    contrast = contrast
  )

  dds <- de_results$dds %||% de_results$object

  dds <- .annotate_dds(
    dds = dds,
    orgdb = orgdb
  )

  res_df <- de_results$res_df

  log("[7/11] Annotating DE results", type = "step")

  res_df$SYMBOL <- map_ensembl_to_feature(
    ensembl_vec = res_df$ENSEMBL,
    orgdb = orgdb
  )

  res_df <- map_entrez_ids(
    res_df = res_df,
    orgdb = orgdb
  )

  processed <- process_deg_results(
    res_df = res_df,
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff
  )

  res_df <- processed$res_df
  sig_deg <- processed$sig_deg

  log("[8/11] Generating volcano, PCA and heatmap", type = "step")

  volcano_plot <- plot_volcano(
    res_df = res_df,
    top_n_labels = top_n_labels,
    output_dir = dirs$volcano,
    filename = "Volcano_Plot",
    logfc_cutoff = logfc_cutoff,
    padj_cutoff = padj_cutoff
  )

  pca_plot <- NULL
  heatmap_matrix <- NULL
  vsd <- NULL

  if (method == "DESeq2") {

    pca_results <- plot_pca(
      dds = dds,
      metadata = metadata,
      output_dir = dirs$pca
    )

    vsd <- pca_results$vsd
    pca_plot <- pca_results$pca_plot

    heatmap_matrix <- plot_heatmap(
      vsd = vsd,
      res_df = res_df,
      metadata = metadata,
      top_n_heatmap = top_n_heatmap,
      padj_cutoff = padj_cutoff,
      output_dir = dirs$heatmaps,
      main = "Top DEG heatmap",
      filename = "Top_DEG_Heatmap",
      annotation_cols = intersect(
        c("condition", "treatment", "sex", "tissue"),
        colnames(metadata)
      ),
      order_by = intersect(
        c("condition", "tissue", "sex", "treatment"),
        colnames(metadata)
      )
    )
  }

  log("[9/11] Running GO enrichment", type = "step")

  go_results <- list()

  go_results <- run_go_enrichment(
    sig_deg = sig_deg,
    comparison = "single_comparison",
    ontology = ontology,
    output_dir = dirs$go,
    orgdb = orgdb
  )

  go_df_i <- go_results$go_results

  go_plot <- NULL

  if (!is.null(go_df_i) && nrow(go_df_i) > 0) {

    go_plot <- plot_go_terms(
      go_df = go_df_i,
      comparison = "single_comparison GO enrichment",
      top_n = 10,
      style = "bw"
    )

    ggplot2::ggsave(
      filename = file.path(dirs$go, "single_comparison_GO_terms.png"),
      plot = go_plot,
      width = 8,
      height = 6,
      dpi = 300
    )

    ggplot2::ggsave(
      filename = file.path(dirs$go, "single_comparison_GO_terms.pdf"),
      plot = go_plot,
      width = 8,
      height = 6
    )
  }

  de_results$go_plot <- go_plot

  log("[10/11] Exporting final results", type = "step")

  export_deg_results(
    res_df,
    sig_deg,
    dirs$results
  )

  de_results$dds <- dds
  de_results$res_df <- res_df
  de_results$sig_deg <- sig_deg
  de_results$volcano_plot <- volcano_plot
  de_results$pca_plot <- pca_plot
  de_results$heatmap_matrix <- heatmap_matrix
  de_results$go_results <- go_results
  de_results$output_dir <- output_dir
  de_results$output_dirs <- dirs
  de_results$version <- deggo_version

  de_results$gene_annotation <- data.frame(
    ENSEMBL = rownames(dds),
    SYMBOL = SummarizedExperiment::rowData(dds)$SYMBOL,
    ENTREZID = SummarizedExperiment::rowData(dds)$ENTREZID,
    GENENAME = SummarizedExperiment::rowData(dds)$GENENAME,
    stringsAsFactors = FALSE
  )

  .safe_write_session_info(output_dir)

  log("[11/11] Generating DEGgo report", type = "step")

  de_results$run_params <- list(
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
    output_dir = output_dir
  )

  de_results$report_files <- NULL

  if (isTRUE(generate_report)) {
    de_results$report_files <- generate_deggo_report(
      results = de_results,
      output_dir = output_dir,
      formats = report_formats,
      report_template = report_template
    )
  }

  log(
    "==== DEGgo SINGLE ANALYSIS COMPLETE ====",
    type = "done",
    duration = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  )

  return(de_results)
}
