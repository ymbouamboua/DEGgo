#' Run DEGgo bulk RNA-seq downstream analysis
#'
#' Run an automated bulk RNA-seq downstream analysis workflow including
#' input preparation, sample matching, differential expression analysis,
#' result annotation, visualization, heatmap generation,
#' PCA analysis, and Gene Ontology enrichment.
#'
#' DEGgo automatically detects common gene identifier columns
#' (e.g. gene_id, ENSEMBL, gene) and optional gene annotation
#' columns (e.g. gene_name, SYMBOL), allowing flexible input
#' count tables from multiple RNA-seq pipelines.
#'
#' @param counts Raw count table or matrix.
#' @param metadata Sample metadata data frame.
#' @param gene_col Character vector of possible gene identifier columns.
#'   Used to identify the column containing unique gene IDs
#'   (e.g. ENSEMBL IDs) in the counts table.
#'   The first matching column is used automatically.
#' @param feature_col Character vector of possible feature annotation
#'   columns containing gene symbols or gene names.
#'   If present, this information is retained and used for
#'   annotation, plotting, and reporting.
#' @param sample_col Character vector of possible sample identifier
#'   columns in the metadata table.
#'   The first matching column is used to match metadata
#'   with count matrix columns.
#' @param prepare_input Logical. If TRUE, automatically matches counts and
#'   metadata using \code{prepare_counts_metadata()}.
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
#' @param seed Random seed.
#'
#' @return A list containing DEG results, plots, GO results, and output path.
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
  log("[1/10] Resolving analysis parameters", type = "step")

  if (is.null(output_dir)) {
    output_dir <- file.path(
      "DEGgo_results",
      paste0(format(Sys.Date(), "%Y-%m-%d"), "_", method, "_", analysis_mode)
    )
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  deggo_version <- tryCatch(
    as.character(utils::packageVersion("DEGgo")),
    error = function(e) "development"
  )

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  # ------------------------------------------------------- #
  # Step 2: input preparation
  # ------------------------------------------------------- #

  if (isTRUE(prepare_input)) {

    log("[2/10] Matching counts and metadata", type = "step")

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

  if (!is.numeric(counts) && !is.integer(counts)) {
    stop("counts must be numeric/integer after input preparation.", call. = FALSE)
  }

  # ------------------------------------------------------- #
  # Pairwise validation
  # ------------------------------------------------------- #

  if (analysis_mode == "pairwise" && method != "DESeq2") {
    stop(
      "Pairwise mode is currently available only with method = 'DESeq2'.",
      call. = FALSE
    )
  }

  if (analysis_mode == "pairwise") {

    if (is.null(pairwise_group_cols)) {
      stop(
        "'pairwise_group_cols' must be provided when analysis_mode = 'pairwise'.",
        call. = FALSE
      )
    }

    missing_pairwise_cols <- setdiff(pairwise_group_cols, colnames(metadata))

    if (length(missing_pairwise_cols) > 0) {
      stop(
        "Missing pairwise grouping columns in metadata: ",
        paste(missing_pairwise_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  log("DE method:", method, type = "info")
  log("Analysis mode:", analysis_mode, type = "info")
  log("Design:", paste(deparse(design_formula), collapse = ""), type = "info")
  log("Organism:", organism, type = "info")
  log("Adjusted p-value cutoff:", padj_cutoff, type = "info")
  log("Log2 fold-change cutoff:", logfc_cutoff, type = "info")
  log("GO ontology:", ontology, type = "info")
  log("Output directory:", output_dir, type = "info")

  log("[3/10] Loading annotation database", type = "step")

  orgdb <- .get_orgdb(
    organism = organism,
    orgdb = orgdb
  )

  log("[4/10] Validating thresholds and inputs", type = "step")

  if (padj_cutoff <= 0 || padj_cutoff >= 1) {
    stop("padj_cutoff must be between 0 and 1.", call. = FALSE)
  }

  if (logfc_cutoff < 0) {
    stop("logfc_cutoff must be positive.", call. = FALSE)
  }

  if (top_n_heatmap <= 0) {
    stop("top_n_heatmap must be positive.", call. = FALSE)
  }

  if (top_n_labels <= 0) {
    stop("top_n_labels must be positive.", call. = FALSE)
  }

  if (!"condition" %in% colnames(metadata)) {
    if (analysis_mode == "pairwise" && !is.null(pairwise_group_cols)) {
      metadata$condition <- apply(
        metadata[, pairwise_group_cols, drop = FALSE],
        1,
        paste,
        collapse = "_"
      )
    } else {
      stop(
        "Metadata must contain 'condition' or provide pairwise_group_cols.",
        call. = FALSE
      )
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

  log("[5/10] Cleaning gene identifiers", type = "step")

  counts <- clean_ensembl_ids(counts)

  log("[6/10] Filtering low-expression genes", type = "step")

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

  # ------------------------------------------------------- #
  # Pairwise mode
  # ------------------------------------------------------- #

  if (analysis_mode == "pairwise") {

    log("[7/10] Running pairwise DESeq2 contrasts", type = "step")

    de_results <- run_deseq2_pairwise(
      counts = counts,
      metadata = metadata,
      group_cols = pairwise_group_cols,
      contrast_col = pairwise_contrast_col,
      contrast_list = pairwise_contrasts,
      pairwise_mode = pairwise_mode
    )

    pca_dir <- file.path(output_dir, "pairwise_PCA")
    dir.create(pca_dir, showWarnings = FALSE, recursive = TRUE)

    de_results$pca <- list(
      sample = plot_pca(
        de_results$dds,
        de_results$metadata,
        pca_dir,
        "PCA_by_sample",
        color_by = "sample",
        title = "PCA by sample"
      ),
      tissue = plot_pca(
        de_results$dds,
        de_results$metadata,
        pca_dir,
        "PCA_by_tissue",
        color_by = "tissue",
        title = "PCA by tissue"
      ),
      treatment = plot_pca(
        de_results$dds,
        de_results$metadata,
        pca_dir,
        "PCA_by_treatment",
        color_by = "treatment",
        title = "PCA by treatment"
      ),
      sex = plot_pca(
        de_results$dds,
        de_results$metadata,
        pca_dir,
        "PCA_by_sex",
        color_by = "sex",
        title = "PCA by sex"
      ),
      tissue_treatment = plot_pca(
        de_results$dds,
        de_results$metadata,
        pca_dir,
        "PCA_tissue_treatment",
        color_by = "tissue",
        shape_by = "treatment",
        title = "PCA tissue + treatment"
      )
    )

    log("[8/10] Annotating pairwise DEG tables", type = "step")

    de_results$results <- lapply(
      de_results$results,
      annotate_de_results,
      orgdb = orgdb,
      padj_cutoff = padj_cutoff,
      logfc_cutoff = logfc_cutoff
    )

    result_dir <- file.path(output_dir, "pairwise_results")
    sig_dir <- file.path(output_dir, "pairwise_significant")
    plot_dir <- file.path(output_dir, "pairwise_plots")
    go_dir <- file.path(output_dir, "pairwise_GO")
    heatmap_dir <- file.path(output_dir, "pairwise_heatmaps")

    dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(sig_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

    de_results$sig_deg <- list()
    de_results$volcano_plots <- list()
    de_results$heatmaps <- list()
    de_results$go_results <- list()

    log("[9/10] Exporting pairwise plots and GO enrichment", type = "step")

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
        file.path(result_dir, paste0(nm, ".tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      utils::write.table(
        sig_i,
        file.path(sig_dir, paste0(nm, "_significant.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      de_results$volcano_plots[[nm]] <- plot_volcano(
        res_df = res_df_i,
        top_n_labels = top_n_labels,
        output_dir = plot_dir,
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
        output_dir = heatmap_dir,
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

      go_i <- run_go_enrichment(
        sig_deg = sig_i,
        ontology = ontology,
        output_dir = go_dir,
        orgdb = orgdb
      )

      de_results$go_results[[nm]] <- go_i

      if (!is.null(go_i) && !is.null(go_i$go_results)) {
        utils::write.table(
          go_i$go_results,
          file.path(go_dir, paste0(nm, "_GO_", ontology, "_enrichment.tsv")),
          sep = "\t",
          row.names = FALSE,
          quote = FALSE
        )
      }
    }

    log("[10/10] Summarizing pairwise results", type = "step")

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

    .safe_write_session_info(output_dir)

    de_results$output_dir <- output_dir
    de_results$version <- deggo_version

    log(
      "==== DEGgo PAIRWISE ANALYSIS COMPLETE ====",
      type = "done",
      duration = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    )

    return(de_results)
  }

  # ------------------------------------------------------- #
  # Single mode
  # ------------------------------------------------------- #

  log("[7/10] Running differential expression analysis", type = "step")

  de_results <- run_de(
    counts = counts,
    metadata = metadata,
    method = method,
    design_formula = design_formula,
    contrast = contrast
  )

  dds <- de_results$dds %||% de_results$object
  res_df <- de_results$res_df

  log("[8/10] Annotating DEG results", type = "step")

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

  up_n <- sum(
    !is.na(res_df$padj) &
      res_df$padj < padj_cutoff &
      res_df$log2FoldChange > logfc_cutoff,
    na.rm = TRUE
  )

  down_n <- sum(
    !is.na(res_df$padj) &
      res_df$padj < padj_cutoff &
      res_df$log2FoldChange < -logfc_cutoff,
    na.rm = TRUE
  )

  log("[9/10] Generating plots and GO enrichment", type = "step")

  volcano_plot <- plot_volcano(
    res_df = res_df,
    top_n_labels = top_n_labels,
    output_dir = output_dir,
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
      output_dir = output_dir
    )

    vsd <- pca_results$vsd
    pca_plot <- pca_results$pca_plot

    heatmap_matrix <- plot_heatmap(
      vsd = vsd,
      res_df = res_df,
      metadata = metadata,
      top_n_heatmap = top_n_heatmap,
      padj_cutoff = padj_cutoff,
      output_dir = output_dir,
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

  } else {
    log("PCA and heatmap currently available only for DESeq2.", type = "warn")
  }

  go_results <- run_go_enrichment(
    sig_deg = sig_deg,
    ontology = ontology,
    output_dir = output_dir,
    orgdb = orgdb
  )

  log("[10/10] Exporting final results", type = "step")

  export_deg_results(
    res_df,
    sig_deg,
    output_dir
  )

  summary_text <- paste0(
    "DEGgo Analysis Summary\n",
    "======================\n\n",
    "DEGgo Version: ", deggo_version, "\n",
    "DE method: ", method, "\n",
    "Analysis mode: ", analysis_mode, "\n",
    "Organism: ", organism, "\n",
    "GO Ontology: ", ontology, "\n\n",
    "Adjusted p-value cutoff: ", padj_cutoff, "\n",
    "Log2 fold-change cutoff: ", logfc_cutoff, "\n\n",
    "Total genes analyzed: ", nrow(res_df), "\n",
    "Significant DEGs: ", nrow(sig_deg), "\n",
    "Upregulated genes: ", up_n, "\n",
    "Downregulated genes: ", down_n, "\n"
  )

  writeLines(summary_text, file.path(output_dir, "Analysis_Summary.txt"))

  .safe_write_session_info(output_dir)

  generate_report(
    sig_deg = sig_deg,
    go_results = go_results,
    output_dir = output_dir,
    report_template = file.path("inst", "report_template.Rmd")
  )

  log("Significant DEGs:", nrow(sig_deg), type = "info")
  log("Upregulated genes:", up_n, type = "info")
  log("Downregulated genes:", down_n, type = "info")

  log(
    "==== DEGgo ANALYSIS COMPLETE ====",
    type = "done",
    duration = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  )

  list(
    dds = dds,
    res_df = res_df,
    sig_deg = sig_deg,
    volcano_plot = volcano_plot,
    pca_plot = pca_plot,
    heatmap_matrix = heatmap_matrix,
    go_results = go_results,
    output_dir = output_dir,
    version = deggo_version
  )
}
