#' Run DEGgo bulk RNA-seq downstream analysis
#'
#' Run an automated bulk RNA-seq downstream analysis workflow including
#' differential expression analysis, result annotation, visualization,
#' heatmap generation, PCA analysis, and Gene Ontology enrichment.
#'
#' The function supports two analysis modes:
#' \itemize{
#'   \item \code{"single"}: run one differential expression analysis using
#'   \code{design_formula} and optional \code{contrast}.
#'   \item \code{"pairwise"}: run multiple pairwise contrasts generated
#'   automatically or supplied by the user.
#' }
#'
#' @param counts Raw count matrix with genes as rows and samples as columns.
#'   Column names must match sample identifiers in \code{metadata}.
#' @param metadata Sample metadata data frame. Row names should correspond to
#'   sample names in \code{counts}. If a \code{sample} column is present, it can
#'   be used upstream by \code{prepare_counts_metadata()} to define row names.
#' @param output_dir Character. Directory where results will be exported.
#'   If \code{NULL}, a timestamped output directory is created automatically.
#' @param padj_cutoff Numeric. Adjusted p-value cutoff used to define
#'   significant differentially expressed genes.
#' @param logfc_cutoff Numeric. Absolute log2 fold-change cutoff used to define
#'   significant differentially expressed genes.
#' @param top_n_heatmap Integer. Number of top genes displayed in the heatmap.
#' @param top_n_labels Integer. Number of top genes labeled in the volcano plot.
#' @param ontology Character. Gene Ontology category used for enrichment.
#'   One of \code{"BP"}, \code{"MF"}, or \code{"CC"}.
#' @param organism Character. Organism used for gene annotation.
#'   One of \code{"human"}, \code{"mouse"}, or \code{"rat"}.
#' @param orgdb Optional OrgDb annotation object. If \code{NULL}, the OrgDb is
#'   selected automatically from \code{organism}.
#' @param method Character. Differential expression engine.
#'   One of \code{"DESeq2"}, \code{"edgeR"}, or \code{"limma"}.
#' @param analysis_mode Character. Analysis mode. Use \code{"single"} for one
#'   differential expression analysis or \code{"pairwise"} for multiple
#'   pairwise comparisons.
#' @param contrast Optional contrast passed to the differential expression
#'   method in \code{"single"} mode. For DESeq2, this is typically a character
#'   vector such as \code{c("condition", "treated", "control")}.
#' @param design_formula Formula. Design formula used for differential
#'   expression analysis in \code{"single"} mode.
#' @param pairwise_group_cols Character vector. Metadata columns used to build
#'   combined groups for pairwise analysis. For example,
#'   \code{c("condition", "sex")} creates groups such as
#'   \code{"WT_F"} and \code{"KO_M"}.
#' @param pairwise_contrast_col Character. Name of the combined metadata column
#'   created for pairwise analysis.
#' @param pairwise_contrasts Optional named list of pairwise contrasts. Each
#'   contrast must be a character vector of length three:
#'   \code{c(pairwise_contrast_col, level_1, level_2)}. If \code{NULL} and
#'   \code{pairwise_auto = TRUE}, contrasts are generated automatically.
#' @param filter_method Character. Gene filtering method.
#'   One of \code{"count"}, \code{"cpm"}, or \code{"none"}.
#' @param pairwise_auto Logical. If \code{TRUE}, automatically generate
#'   pairwise contrasts when \code{pairwise_contrasts = NULL}.
#' @param pairwise_mode Character. Pairwise contrast generation mode:
#'   \itemize{
#'     \item \code{"all"}: all pairwise comparisons among combined groups.
#'     \item \code{"within_first"}: comparisons within levels of the first
#'     column in \code{pairwise_group_cols}.
#'     \item \code{"within_second"}: comparisons within levels of the second
#'     column in \code{pairwise_group_cols}.
#'   }
#' @param min_count Integer. Minimum count required for count-based filtering.
#' @param min_samples Integer. Minimum number of samples required to pass
#'   the count threshold.
#' @param min_total Integer. Minimum total count required across all samples.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return
#' If \code{analysis_mode = "single"}, returns a list containing:
#' \itemize{
#'   \item \code{dds}: DESeq2 dataset object or equivalent model object.
#'   \item \code{res_df}: annotated differential expression results.
#'   \item \code{sig_deg}: significant differentially expressed genes.
#'   \item \code{volcano_plot}: volcano plot object.
#'   \item \code{pca_plot}: PCA plot object, when available.
#'   \item \code{heatmap_matrix}: expression matrix used for heatmap,
#'   when available.
#'   \item \code{go_results}: Gene Ontology enrichment results.
#' }
#'
#' If \code{analysis_mode = "pairwise"}, returns a list containing:
#' \itemize{
#'   \item \code{dds}: DESeq2 dataset object.
#'   \item \code{results}: named list of pairwise differential expression
#'   result tables.
#'   \item \code{summary}: summary table for all pairwise comparisons.
#'   \item \code{metadata}: metadata used in the pairwise analysis.
#'   \item \code{contrast_col}: combined contrast column name.
#'   \item \code{contrast_levels}: available contrast levels.
#' }
#'
#' @details
#' In \code{"pairwise"} mode, the function first creates a combined group
#' variable from \code{pairwise_group_cols}. For example, if
#' \code{pairwise_group_cols = c("condition", "sex")}, the combined groups may
#' be \code{"WT_F"}, \code{"WT_M"}, \code{"KO_F"}, and \code{"KO_M"}.
#'
#' When \code{pairwise_contrasts = NULL} and \code{pairwise_auto = TRUE},
#' pairwise contrasts are generated automatically according to
#' \code{pairwise_mode}. Users can also provide a custom named list of
#' contrasts for full control.
#'
#' @examples
#' \dontrun{
#' # Single comparison
#' results <- run_deggo(
#'   counts = counts,
#'   metadata = metadata,
#'   organism = "mouse",
#'   method = "DESeq2",
#'   design_formula = ~ condition,
#'   contrast = c("condition", "treated", "control")
#' )
#'
#' # Automatic pairwise comparisons among condition-sex groups
#' results_pairwise <- run_deggo(
#'   counts = counts,
#'   metadata = metadata,
#'   organism = "mouse",
#'   method = "DESeq2",
#'   analysis_mode = "pairwise",
#'   pairwise_group_cols = c("condition", "sex"),
#'   pairwise_mode = "all"
#' )
#'
#' # User-defined pairwise contrasts
#' my_contrasts <- list(
#'   treated_male_vs_control_male =
#'     c("comparison_group", "treated_M", "control_M"),
#'   treated_female_vs_control_female =
#'     c("comparison_group", "treated_F", "control_F")
#' )
#'
#' results_custom <- run_deggo(
#'   counts = counts,
#'   metadata = metadata,
#'   organism = "mouse",
#'   method = "DESeq2",
#'   analysis_mode = "pairwise",
#'   pairwise_group_cols = c("condition", "sex"),
#'   pairwise_contrasts = my_contrasts
#' )
#' }
#'
#' @export
#' 
run_deggo <- function(
    counts,
    metadata,
    output_dir = NULL,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25,
    top_n_heatmap = 50,
    top_n_labels = 10,
    ontology = "BP",
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
  
  log <- .msg(verbose = TRUE)
  t_start <- Sys.time()
  
  log("==== STARTING DEGgo ANALYSIS ====", type = "header")
  
  log("Resolving analysis parameters", type = "step")
  
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  
  if ("sample" %in% colnames(metadata)) {
    rownames(metadata) <- metadata$sample
  }
  
  if (analysis_mode == "pairwise" && method != "DESeq2") {
    stop("Pairwise mode is currently available only with method = 'DESeq2'.",
         call. = FALSE)
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
  
  if (is.null(output_dir)) {
    output_dir <- paste0(
      "DEGgo_",
      method,
      "_",
      analysis_mode,
      "_",
      format(Sys.time(), "%Y%m%d_%H%M%S")
    )
  }
  
  log("DE method:", method, type = "info")
  log("Analysis mode:", analysis_mode, type = "info")
  log("Design:", paste(deparse(design_formula), collapse = ""), type = "info")
  log("Organism:", organism, type = "info")
  log("Adjusted p-value cutoff:", padj_cutoff, type = "info")
  log("Log2 fold-change cutoff:", logfc_cutoff, type = "info")
  log("GO ontology:", ontology, type = "info")
  log("Output directory:", output_dir, type = "info")
  
  log("Loading annotation database", type = "step")
  orgdb <- .get_orgdb(organism = organism, orgdb = orgdb)
  
  cellDEverse_version <- "0.1.0"
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  
  log("Validating thresholds and inputs", type = "step")
  
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
  
  log("Cleaning gene identifiers", type = "step")
  counts <- clean_ensembl_ids(counts)
  
  log("Filtering low-expression genes", type = "step")
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
  
  if (analysis_mode == "pairwise") {
    
    log("Running pairwise DESeq2 contrasts", type = "step")
    
    de_results <- run_deseq2_pairwise(
      counts = counts,
      metadata = metadata,
      group_cols = pairwise_group_cols,
      contrast_col = pairwise_contrast_col,
      contrast_list = pairwise_contrasts,
      pairwise_mode = pairwise_mode
    )
    
    log("Annotating pairwise DEG tables", type = "step")
    
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
    
    dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(sig_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)
    
    de_results$sig_deg <- list()
    de_results$volcano_plots <- list()
    de_results$go_results <- list()
    
    log("Exporting pairwise results, plots and GO enrichment", type = "step")
    
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
        output_dir = plot_dir
      )
      
      if (file.exists(file.path(plot_dir, "Volcano_Plot.png"))) {
        file.rename(
          file.path(plot_dir, "Volcano_Plot.png"),
          file.path(plot_dir, paste0(nm, "_Volcano_Plot.png"))
        )
      }
      
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
      
      if (file.exists(file.path(go_dir, "GO_Dotplot.png"))) {
        file.rename(
          file.path(go_dir, "GO_Dotplot.png"),
          file.path(go_dir, paste0(nm, "_GO_Dotplot.png"))
        )
      }
      
      if (file.exists(file.path(go_dir, "GO_Barplot.png"))) {
        file.rename(
          file.path(go_dir, "GO_Barplot.png"),
          file.path(go_dir, paste0(nm, "_GO_Barplot.png"))
        )
      }
    }
    
    log("Summarizing pairwise results", type = "step")
    
    de_results$summary <- do.call(
      rbind,
      lapply(names(de_results$results), function(nm) {
        x <- de_results$results[[nm]]
        data.frame(
          comparison = nm,
          total_genes = nrow(x),
          significant = sum(x$significance != "Not Significant", na.rm = TRUE),
          up = sum(x$significance == "Upregulated", na.rm = TRUE),
          down = sum(x$significance == "Downregulated", na.rm = TRUE),
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
    
    writeLines(
      utils::capture.output(utils::sessionInfo()),
      file.path(output_dir, "sessionInfo.txt")
    )
    
    log(
      "Pairwise analysis completed successfully.",
      type = "done",
      duration = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    )
    
    return(de_results)
  }
  
  log("Running differential expression analysis", type = "step")
  
  de_results <- run_de(
    counts = counts,
    metadata = metadata,
    method = method,
    design_formula = design_formula,
    contrast = contrast
  )
  
  dds <- de_results$dds %||% de_results$object
  res_df <- de_results$res_df
  
  log("Annotating DEG results", type = "step")
  
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
  
  up_n <- if (nrow(sig_deg)) {
    sum(sig_deg$log2FoldChange > 0, na.rm = TRUE)
  } else {
    0
  }
  
  down_n <- if (nrow(sig_deg)) {
    sum(sig_deg$log2FoldChange < 0, na.rm = TRUE)
  } else {
    0
  }
  
  log("Generating volcano plot", type = "step")
  
  volcano_plot <- plot_volcano(
    res_df = res_df,
    top_n_labels = top_n_labels,
    output_dir = output_dir
  )
  
  pca_plot <- NULL
  heatmap_matrix <- NULL
  vsd <- NULL
  
  if (method == "DESeq2") {
    
    log("Generating PCA plot", type = "step")
    
    pca_results <- plot_pca(
      dds = dds,
      metadata = metadata,
      output_dir = output_dir
    )
    
    vsd <- pca_results$vsd
    pca_plot <- pca_results$pca_plot
    
    log("Generating heatmap", type = "step")
    
    heatmap_matrix <- plot_heatmap(
      vsd = vsd,
      res_df = res_df,
      metadata = metadata,
      top_n_heatmap = top_n_heatmap,
      padj_cutoff = padj_cutoff,
      output_dir = output_dir
    )
    
  } else {
    
    log(
      "PCA and heatmap currently available only for DESeq2.",
      type = "warn"
    )
  }
  
  log("Running GO enrichment", type = "step")
  
  go_results <- run_go_enrichment(
    sig_deg = sig_deg,
    ontology = ontology,
    output_dir = output_dir,
    orgdb = orgdb
  )
  
  log("Exporting final results", type = "step")
  
  export_deg_results(
    res_df,
    sig_deg,
    output_dir
  )
  
  summary_text <- paste0(
    "DEGgo Analysis Summary\n",
    "======================\n\n",
    "DEGgo Version: ", cellDEverse_version, "\n",
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
  writeLines(
    utils::capture.output(utils::sessionInfo()),
    file.path(output_dir, "sessionInfo.txt")
  )
  
  generate_report(
    sig_deg = sig_deg,
    go_results = go_results,
    output_dir = output_dir,
    report_template = file.path("inst", "report_template.Rmd")
  )
  
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
    go_results = go_results
  )
}