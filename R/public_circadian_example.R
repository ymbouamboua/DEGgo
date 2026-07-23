# ============================================================
# Public GSE98965 circadian RNA-seq example
# ============================================================


#' Download the GSE98965 baboon transcriptomic atlas
#'
#' Download the processed gene-level FPKM expression matrix from the
#' Gene Expression Omnibus accession GSE98965.
#'
#' The dataset corresponds to the baboon diurnal transcriptome atlas
#' described by Mure et al. and contains gene-level expression values
#' across multiple tissues and Zeitgeber Time points.
#'
#' If the expected supplementary file is already present in
#' `data_dir`, it is reused and is not downloaded again.
#'
#' @param data_dir Character string. Directory in which the GSE98965
#'   supplementary files are stored. The directory is created when it
#'   does not exist.
#'
#' @return A `data.frame` containing the complete GSE98965 processed
#'   expression table.
#'
#' @details
#' This function requires the Bioconductor package \pkg{GEOquery}.
#' Supplementary files are retrieved with
#' [GEOquery::getGEOSuppFiles()].
#'
#' @references
#' Mure LS et al. Diurnal transcriptome atlas of a primate across major
#' neural and peripheral tissues. Science. 2018.
#'
#' GEO accession: GSE98965.
#'
#' @examples
#' \dontrun{
#' baboon_expression <- download_gse98965(
#'   data_dir = "GSE98965"
#' )
#'
#' dim(baboon_expression)
#' }
#'
#' @export
download_gse98965 <- function(
    data_dir = "GSE98965"
) {
  
  if (
    length(data_dir) != 1L ||
    is.na(data_dir) ||
    !nzchar(data_dir)
  ) {
    stop(
      "`data_dir` must be a non-empty character string.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop(
      "The Bioconductor package 'GEOquery' is required. ",
      "Install it with:\n",
      "BiocManager::install('GEOquery')",
      call. = FALSE
    )
  }
  
  dir.create(
    data_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  expected_file <- file.path(
    data_dir,
    "GSE98965_baboon_tissue_expression_FPKM.csv.gz"
  )
  
  if (!file.exists(expected_file)) {
    
    GEOquery::getGEOSuppFiles(
      GEO = "GSE98965",
      makeDirectory = FALSE,
      baseDir = data_dir
    )
  }
  
  candidates <- list.files(
    path = data_dir,
    pattern = "FPKM.*\\.csv\\.gz$|\\.csv\\.gz$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  expression_file <- if (file.exists(expected_file)) {
    
    expected_file
    
  } else if (length(candidates) > 0L) {
    
    candidates[1L]
    
  } else {
    
    stop(
      "The GSE98965 expression matrix was not found in: ",
      normalizePath(
        data_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      call. = FALSE
    )
  }
  
  expression_data <- utils::read.csv(
    expression_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (nrow(expression_data) == 0L) {
    stop(
      "The downloaded GSE98965 expression table is empty.",
      call. = FALSE
    )
  }
  
  if (ncol(expression_data) < 2L) {
    stop(
      "The downloaded GSE98965 expression table does not contain ",
      "sufficient columns.",
      call. = FALSE
    )
  }
  
  attr(
    expression_data,
    "source_file"
  ) <- normalizePath(
    expression_file,
    winslash = "/",
    mustWork = FALSE
  )
  
  expression_data
}


#' Prepare one GSE98965 tissue for circadian analysis
#'
#' Extract one tissue from the GSE98965 baboon transcriptomic atlas,
#' construct a numeric genes-by-samples expression matrix, filter lowly
#' expressed genes, generate numeric Zeitgeber Time metadata, and prepare
#' gene-symbol annotations.
#'
#' @param expression_data A `data.frame` containing the processed
#'   GSE98965 FPKM expression table, usually returned by
#'   [DEGgo::download_gse98965()].
#'
#' @param tissue_code Character string specifying the GSE98965 tissue
#'   code to extract. The default is `"LIV"` for liver.
#'
#' @param min_expression Numeric value. Minimum FPKM expression required
#'   for a gene to be considered expressed in a sample.
#'
#' @param min_samples Positive integer. Minimum number of samples in
#'   which a gene must reach `min_expression`.
#'
#' @param log_transform Logical. If `TRUE`, expression values are
#'   transformed as `log2(FPKM + 1)`.
#'
#' @return A named list containing:
#'
#' \describe{
#'   \item{expr}{Numeric genes-by-samples expression matrix.}
#'   \item{metadata}{Sample metadata with `sample`, `tissue`, `ZT`,
#'     and numeric `time` columns.}
#'   \item{gene_annotation}{Gene identifier and symbol annotation table.}
#'   \item{filter_summary}{Summary of expression filtering.}
#'   \item{tissue_code}{Selected tissue code.}
#'   \item{assay}{Expression representation returned by the function.}
#' }
#'
#' @details
#' Sample columns are detected using the pattern
#' `TISSUE.ZTXX`, for example `LIV.ZT02`.
#'
#' Duplicate gene identifiers are preserved using
#' [base::make.unique()]. The resulting unique identifiers are stored
#' in the `matrix_gene_id` annotation column and correspond exactly to
#' the expression matrix row names.
#'
#' @examples
#' \dontrun{
#' baboon_expression <- download_gse98965()
#'
#' liver <- prepare_circadian_tissue(
#'   expression_data = baboon_expression,
#'   tissue_code = "LIV",
#'   min_expression = 1,
#'   min_samples = 3
#' )
#'
#' dim(liver$expr)
#' head(liver$metadata)
#' head(liver$gene_annotation)
#' }
#'
#' @export
prepare_circadian_tissue <- function(
    expression_data,
    tissue_code = "LIV",
    min_expression = 1,
    min_samples = 3,
    log_transform = TRUE
) {
  
  if (!is.data.frame(expression_data)) {
    stop(
      "`expression_data` must be a data.frame.",
      call. = FALSE
    )
  }
  
  if (nrow(expression_data) == 0L) {
    stop(
      "`expression_data` contains no rows.",
      call. = FALSE
    )
  }
  
  if (
    length(tissue_code) != 1L ||
    is.na(tissue_code) ||
    !nzchar(tissue_code)
  ) {
    stop(
      "`tissue_code` must be a non-empty character string.",
      call. = FALSE
    )
  }
  
  if (
    length(min_expression) != 1L ||
    !is.numeric(min_expression) ||
    is.na(min_expression) ||
    !is.finite(min_expression) ||
    min_expression < 0
  ) {
    stop(
      "`min_expression` must be a finite non-negative number.",
      call. = FALSE
    )
  }
  
  if (
    length(min_samples) != 1L ||
    !is.numeric(min_samples) ||
    is.na(min_samples) ||
    !is.finite(min_samples) ||
    min_samples < 1 ||
    min_samples != as.integer(min_samples)
  ) {
    stop(
      "`min_samples` must be a positive integer.",
      call. = FALSE
    )
  }
  
  if (
    length(log_transform) != 1L ||
    is.na(log_transform) ||
    !is.logical(log_transform)
  ) {
    stop(
      "`log_transform` must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  
  tissue_code <- as.character(tissue_code)
  min_samples <- as.integer(min_samples)
  
  sample_pattern <- paste0(
    "^",
    tissue_code,
    "\\.ZT[0-9]{2}$"
  )
  
  sample_columns <- grep(
    sample_pattern,
    colnames(expression_data),
    value = TRUE
  )
  
  if (length(sample_columns) == 0L) {
    
    available_sample_columns <- grep(
      "^[A-Za-z0-9]+\\.ZT[0-9]{2}$",
      colnames(expression_data),
      value = TRUE
    )
    
    available_tissues <- sort(
      unique(
        sub(
          "\\.ZT[0-9]{2}$",
          "",
          available_sample_columns
        )
      )
    )
    
    available_text <- if (length(available_tissues) > 0L) {
      paste(
        utils::head(
          available_tissues,
          20L
        ),
        collapse = ", "
      )
    } else {
      "none detected"
    }
    
    stop(
      "No samples were found for tissue code '",
      tissue_code,
      "'. Available tissue codes include: ",
      available_text,
      ".",
      call. = FALSE
    )
  }
  
  if (min_samples > length(sample_columns)) {
    stop(
      "`min_samples` cannot exceed the number of selected samples (",
      length(sample_columns),
      ").",
      call. = FALSE
    )
  }
  
  gene_id_candidates <- c(
    "EnsemblID",
    "gene_id",
    "GeneID",
    "gene",
    "Gene",
    "ENSEMBL",
    "ensembl",
    "ensembl_id"
  )
  
  symbol_candidates <- c(
    "Symbol",
    "symbol",
    "SYMBOL",
    "gene_symbol",
    "gene_name",
    "GeneName"
  )
  
  gene_id_column <- intersect(
    gene_id_candidates,
    colnames(expression_data)
  )
  
  symbol_column <- intersect(
    symbol_candidates,
    colnames(expression_data)
  )
  
  if (length(gene_id_column) == 0L) {
    stop(
      "No supported gene identifier column was detected. ",
      "Supported names include: ",
      paste(
        gene_id_candidates,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  
  gene_id_column <- gene_id_column[1L]
  
  symbol_column <- if (length(symbol_column) > 0L) {
    symbol_column[1L]
  } else {
    NA_character_
  }
  
  gene_ids <- as.character(
    expression_data[[gene_id_column]]
  )
  
  valid_gene_ids <- !is.na(gene_ids) &
    nzchar(
      trimws(gene_ids)
    )
  
  if (!any(valid_gene_ids)) {
    stop(
      "No valid gene identifiers were found in column '",
      gene_id_column,
      "'.",
      call. = FALSE
    )
  }
  
  expr_data <- expression_data[
    valid_gene_ids,
    sample_columns,
    drop = FALSE
  ]
  
  expr_data[] <- lapply(
    expr_data,
    function(x) {
      suppressWarnings(
        as.numeric(x)
      )
    }
  )
  
  expr <- as.matrix(expr_data)
  
  storage.mode(expr) <- "numeric"
  
  gene_ids <- gene_ids[
    valid_gene_ids
  ]
  
  matrix_gene_ids <- make.unique(
    gene_ids
  )
  
  rownames(expr) <- matrix_gene_ids
  
  symbols <- if (!is.na(symbol_column)) {
    
    as.character(
      expression_data[[symbol_column]][valid_gene_ids]
    )
    
  } else {
    
    rep(
      NA_character_,
      length(gene_ids)
    )
  }
  
  symbols[
    !is.na(symbols) &
      !nzchar(trimws(symbols))
  ] <- NA_character_
  
  annotation <- data.frame(
    gene_id = gene_ids,
    matrix_gene_id = matrix_gene_ids,
    symbol = symbols,
    stringsAsFactors = FALSE
  )
  
  finite_genes <- rowSums(
    is.finite(expr)
  ) == ncol(expr)
  
  n_non_finite <- sum(
    !finite_genes
  )
  
  expr <- expr[
    finite_genes,
    ,
    drop = FALSE
  ]
  
  annotation <- annotation[
    finite_genes,
    ,
    drop = FALSE
  ]
  
  if (nrow(expr) == 0L) {
    stop(
      "No genes remained after removing non-finite expression rows.",
      call. = FALSE
    )
  }
  
  keep_expressed <- rowSums(
    expr >= min_expression
  ) >= min_samples
  
  n_before_filtering <- nrow(expr)
  n_retained <- sum(keep_expressed)
  n_removed <- sum(!keep_expressed)
  
  if (n_retained == 0L) {
    stop(
      "No genes passed the expression filter. Consider reducing ",
      "`min_expression` or `min_samples`.",
      call. = FALSE
    )
  }
  
  expr <- expr[
    keep_expressed,
    ,
    drop = FALSE
  ]
  
  annotation <- annotation[
    keep_expressed,
    ,
    drop = FALSE
  ]
  
  if (isTRUE(log_transform)) {
    expr <- log2(
      expr + 1
    )
  }
  
  metadata <- data.frame(
    sample = colnames(expr),
    tissue = tissue_code,
    ZT = sub(
      paste0(
        "^",
        tissue_code,
        "\\."
      ),
      "",
      colnames(expr)
    ),
    stringsAsFactors = FALSE
  )
  
  metadata$time <- suppressWarnings(
    as.numeric(
      sub(
        "^ZT",
        "",
        metadata$ZT
      )
    )
  )
  
  if (anyNA(metadata$time)) {
    stop(
      "One or more Zeitgeber Time values could not be converted ",
      "to numeric values.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(metadata$sample)) {
    stop(
      "Duplicate sample identifiers were detected.",
      call. = FALSE
    )
  }
  
  metadata <- metadata[
    order(
      metadata$time,
      metadata$sample
    ),
    ,
    drop = FALSE
  ]
  
  expr <- expr[
    ,
    metadata$sample,
    drop = FALSE
  ]
  
  rownames(metadata) <- metadata$sample
  
  if (!identical(colnames(expr), metadata$sample)) {
    stop(
      "Expression and metadata sample order could not be aligned.",
      call. = FALSE
    )
  }
  
  filter_summary <- data.frame(
    metric = c(
      "Genes before expression filtering",
      "Genes with non-finite expression removed",
      "Genes retained",
      "Genes removed by expression filter",
      "Minimum expression",
      "Minimum samples",
      "Expression transformation"
    ),
    value = c(
      as.character(n_before_filtering),
      as.character(n_non_finite),
      as.character(n_retained),
      as.character(n_removed),
      as.character(min_expression),
      as.character(min_samples),
      if (isTRUE(log_transform)) {
        "log2(FPKM + 1)"
      } else {
        "FPKM"
      }
    ),
    stringsAsFactors = FALSE
  )
  
  list(
    expr = expr,
    metadata = metadata,
    gene_annotation = annotation,
    filter_summary = filter_summary,
    tissue_code = tissue_code,
    assay = if (isTRUE(log_transform)) {
      "log2_normalized"
    } else {
      "normalized"
    }
  )
}


#' Run the public GSE98965 circadian RNA-seq example
#'
#' Run a complete, reproducible circadian rhythmicity analysis using one
#' tissue from the public GSE98965 baboon transcriptomic atlas.
#'
#' The workflow downloads the public processed expression matrix,
#' selects one tissue, filters lowly expressed genes, converts Zeitgeber
#' Time labels into numeric values, validates expression-metadata
#' correspondence, runs MetaCycle and cosinor analysis, adds gene-symbol
#' annotations, generates diagnostic plots, and exports result tables.
#'
#' @param tissue_code Character string specifying the GSE98965 tissue
#'   code. The default `"LIV"` corresponds to liver.
#'
#' @param data_dir Character string. Directory used to store the
#'   downloaded GSE98965 supplementary files.
#'
#' @param output_dir Character string. Parent output directory for the
#'   analysis. A tissue-specific subdirectory is created automatically.
#'
#' @param assay Character string indicating the expression representation
#'   supplied to [DEGgo::run_deggo_rhythmicity()]. One of `"log2_normalized"`,
#'   `"vst"`, `"raw"`, or `"normalized"`.
#'
#'   The default is `"log2_normalized"` because
#'   [DEGgo::prepare_circadian_tissue()] returns `log2(FPKM + 1)` values.
#'
#' @param min_expression Numeric value. Minimum FPKM expression used to
#'   define an expressed gene.
#'
#' @param min_samples Positive integer. Minimum number of tissue samples
#'   in which a gene must reach `min_expression`.
#'
#' @param padj_cutoff Numeric value between zero and one. Adjusted
#'   p-value threshold used to classify rhythmic genes.
#'
#' @param n_top_plots Positive integer. Number of top rhythmic genes for
#'   which individual expression and fitted-rhythm plots are generated.
#'
#' @param show_gene_id Logical. If `TRUE`, rhythmicity plot labels include
#'   both the gene symbol and matrix gene identifier.
#'
#' @param seed Integer random seed used for reproducibility.
#'
#' @return A named list containing:
#'
#' \describe{
#'   \item{summary}{Combined rhythmicity result table.}
#'   \item{rhythmic_genes}{Genes classified as rhythmic by at least one
#'     selected method.}
#'   \item{rhythmicity_summary}{Counts of genes by rhythmicity class.}
#'   \item{analysis_summary}{Summary of the public dataset analysis.}
#'   \item{prepared}{Prepared expression matrix, metadata, and gene
#'     annotation.}
#'   \item{deggo}{Complete object returned by
#'     [run_deggo_rhythmicity()].}
#'   \item{output_dir}{Tissue-specific output directory.}
#' }
#'
#' @details
#' The default analysis uses:
#'
#' \itemize{
#'   \item MetaCycle algorithms ARS, JTK, and Lomb-Scargle;
#'   \item a tested period range of 20 to 28 hours;
#'   \item a 24-hour single-component cosinor model;
#'   \item Benjamini-Hochberg adjusted p-values;
#'   \item gene symbols for output tables and figure labels.
#' }
#'
#' The workflow is intended as a reproducible demonstration and as a
#' starting point for custom circadian RNA-seq analyses.
#'
#' @references
#' Mure LS et al. Diurnal transcriptome atlas of a primate across major
#' neural and peripheral tissues. Science. 2018.
#'
#' GEO accession: GSE98965.
#'
#' @examples
#' \dontrun{
#' results <- run_public_circadian_example(
#'   tissue_code = "LIV",
#'   data_dir = "GSE98965",
#'   output_dir = "DEGgo_rhythmicity_results",
#'   min_expression = 1,
#'   min_samples = 3,
#'   padj_cutoff = 0.05,
#'   n_top_plots = 25,
#'   show_gene_id = FALSE,
#'   seed = 4173
#' )
#'
#' head(results$summary)
#' head(results$rhythmic_genes)
#' results$rhythmicity_summary
#' }
#'
#' @seealso
#' [DEGgo::download_gse98965()],
#' [DEGgo::prepare_circadian_tissue()],
#' [DEGgo::run_deggo_rhythmicity()]
#'
#' @export
run_public_circadian_example <- function(
    tissue_code = "LIV",
    data_dir = "GSE98965",
    output_dir = "DEGgo_rhythmicity_results",
    assay = c(
      "log2_normalized",
      "vst",
      "raw",
      "normalized"
    ),
    min_expression = 1,
    min_samples = 3,
    padj_cutoff = 0.05,
    n_top_plots = 25,
    show_gene_id = TRUE,
    seed = 4173
) {
  
  assay <- match.arg(assay)
  
  if (
    length(output_dir) != 1L ||
    is.na(output_dir) ||
    !nzchar(output_dir)
  ) {
    stop(
      "`output_dir` must be a non-empty character string.",
      call. = FALSE
    )
  }
  
  if (
    length(padj_cutoff) != 1L ||
    !is.numeric(padj_cutoff) ||
    is.na(padj_cutoff) ||
    !is.finite(padj_cutoff) ||
    padj_cutoff <= 0 ||
    padj_cutoff >= 1
  ) {
    stop(
      "`padj_cutoff` must be a numeric value between 0 and 1.",
      call. = FALSE
    )
  }
  
  if (
    length(n_top_plots) != 1L ||
    !is.numeric(n_top_plots) ||
    is.na(n_top_plots) ||
    !is.finite(n_top_plots) ||
    n_top_plots < 1 ||
    n_top_plots != as.integer(n_top_plots)
  ) {
    stop(
      "`n_top_plots` must be a positive integer.",
      call. = FALSE
    )
  }
  
  if (
    length(show_gene_id) != 1L ||
    is.na(show_gene_id) ||
    !is.logical(show_gene_id)
  ) {
    stop(
      "`show_gene_id` must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  
  if (
    length(seed) != 1L ||
    !is.numeric(seed) ||
    is.na(seed) ||
    !is.finite(seed)
  ) {
    stop(
      "`seed` must be a finite numeric value.",
      call. = FALSE
    )
  }
  
  set.seed(
    as.integer(seed)
  )
  
  cat("[1/4] Loading GSE98965\n")
  
  baboon_fpkm <- download_gse98965(
    data_dir = data_dir
  )
  
  cat(
    "[2/4] Preparing tissue: ",
    tissue_code,
    "\n",
    sep = ""
  )
  
  prepared <- prepare_circadian_tissue(
    expression_data = baboon_fpkm,
    tissue_code = tissue_code,
    min_expression = min_expression,
    min_samples = min_samples,
    log_transform = identical(
      assay,
      "log2_normalized"
    )
  )
  
  if (
    identical(assay, "vst") ||
    identical(assay, "raw")
  ) {
    warning(
      "The GSE98965 public matrix contains processed FPKM values. ",
      "`assay = \"log2_normalized\"` is recommended. ",
      "The requested assay was: ",
      assay,
      ".",
      call. = FALSE
    )
  }
  
  analysis_output_dir <- file.path(
    output_dir,
    tissue_code
  )
  
  dir.create(
    analysis_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  required_annotation_columns <- c(
    "matrix_gene_id",
    "symbol"
  )
  
  missing_annotation_columns <- setdiff(
    required_annotation_columns,
    colnames(prepared$gene_annotation)
  )
  
  if (length(missing_annotation_columns) > 0L) {
    stop(
      "`prepared$gene_annotation` is missing required columns: ",
      paste(
        missing_annotation_columns,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  
  prepared$gene_annotation$matrix_gene_id <- as.character(
    prepared$gene_annotation$matrix_gene_id
  )
  
  prepared$gene_annotation$symbol <- as.character(
    prepared$gene_annotation$symbol
  )
  
  valid_symbol <- !is.na(
    prepared$gene_annotation$symbol
  ) &
    nzchar(
      trimws(
        prepared$gene_annotation$symbol
      )
    ) &
    prepared$gene_annotation$symbol != "NA"
  
  annotation_order <- order(
    !valid_symbol
  )
  
  prepared$gene_annotation <- prepared$gene_annotation[
    annotation_order,
    ,
    drop = FALSE
  ]
  
  prepared$gene_annotation <- prepared$gene_annotation[
    !duplicated(
      prepared$gene_annotation$matrix_gene_id
    ),
    ,
    drop = FALSE
  ]
  
  cat("[3/4] Running DEGgo rhythmicity analysis\n")
  
  rhythm_results <- run_deggo_rhythmicity(
    expr = prepared$expr,
    metadata = prepared$metadata,
    sample_col = "sample",
    time_col = "time",
    group_col = NULL,
    assay = assay,
    methods = c(
      "meta2d",
      "cosinor"
    ),
    period_range = c(
      20,
      28
    ),
    cycle_length = 24,
    cycMethod = c(
      "ARS",
      "JTK",
      "LS"
    ),
    padj_cutoff = padj_cutoff,
    cosinor_engine = "auto",
    output_dir = analysis_output_dir,
    project_name = paste(
      "GSE98965",
      tissue_code,
      "circadian RNA-seq analysis"
    ),
    generate_plots = TRUE,
    n_top_plots = as.integer(
      n_top_plots
    ),
    gene_annotation = prepared$gene_annotation,
    gene_id_col = "matrix_gene_id",
    gene_symbol_col = "symbol",
    show_gene_id = show_gene_id,
    seed = as.integer(seed)
  )
  
  cat("[4/4] Preparing and exporting results\n")
  
  combined_results <- as.data.frame(
    rhythm_results$summary,
    stringsAsFactors = FALSE
  )
  
  if (nrow(combined_results) == 0L) {
    stop(
      "The DEGgo rhythmicity summary is empty.",
      call. = FALSE
    )
  }
  
  result_gene_column <- intersect(
    c(
      "gene",
      "gene_id",
      "matrix_gene_id",
      "feature",
      "Gene",
      "EnsemblID",
      "ENSEMBL"
    ),
    colnames(combined_results)
  )
  
  if (length(result_gene_column) > 0L) {
    
    result_gene_column <- result_gene_column[1L]
    
    combined_results$matrix_gene_id <- as.character(
      combined_results[[result_gene_column]]
    )
    
  } else if (
    !is.null(rownames(combined_results)) &&
    !identical(
      rownames(combined_results),
      as.character(
        seq_len(
          nrow(combined_results)
        )
      )
    )
  ) {
    
    combined_results$matrix_gene_id <- rownames(
      combined_results
    )
    
  } else {
    
    stop(
      "Unable to identify gene identifiers in the DEGgo ",
      "rhythmicity summary.",
      call. = FALSE
    )
  }
  
  annotation_match <- match(
    combined_results$matrix_gene_id,
    prepared$gene_annotation$matrix_gene_id
  )
  
  annotated_symbol <- prepared$gene_annotation$symbol[
    annotation_match
  ]
  
  if ("symbol" %in% colnames(combined_results)) {
    
    existing_symbol <- as.character(
      combined_results$symbol
    )
    
    replace_symbol <- is.na(existing_symbol) |
      !nzchar(trimws(existing_symbol)) |
      existing_symbol == "NA"
    
    existing_symbol[
      replace_symbol
    ] <- annotated_symbol[
      replace_symbol
    ]
    
    combined_results$symbol <- existing_symbol
    
  } else {
    
    combined_results$symbol <- annotated_symbol
  }
  
  invalid_symbol <- is.na(
    combined_results$symbol
  ) |
    !nzchar(
      trimws(
        combined_results$symbol
      )
    ) |
    combined_results$symbol == "NA"
  
  combined_results$gene_label <- combined_results$symbol
  
  combined_results$gene_label[
    invalid_symbol
  ] <- combined_results$matrix_gene_id[
    invalid_symbol
  ]
  
  if (isTRUE(show_gene_id)) {
    
    symbol_differs_from_id <- !invalid_symbol &
      combined_results$symbol !=
      combined_results$matrix_gene_id
    
    combined_results$gene_label[
      symbol_differs_from_id
    ] <- paste0(
      combined_results$symbol[
        symbol_differs_from_id
      ],
      " (",
      combined_results$matrix_gene_id[
        symbol_differs_from_id
      ],
      ")"
    )
  }
  
  if (!"rhythmic_by" %in% colnames(combined_results)) {
    stop(
      "The DEGgo summary table does not contain a ",
      "'rhythmic_by' column.",
      call. = FALSE
    )
  }
  
  rhythmic_genes <- combined_results[
    !is.na(combined_results$rhythmic_by) &
      combined_results$rhythmic_by != "none",
    ,
    drop = FALSE
  ]
  
  padj_columns <- intersect(
    c(
      "metacycle_padj",
      "cosinor_padj",
      "meta2d_BH.Q",
      "meta2d_pvalue"
    ),
    colnames(rhythmic_genes)
  )
  
  if (
    nrow(rhythmic_genes) > 0L &&
    length(padj_columns) > 0L
  ) {
    
    rhythmic_genes$best_adjusted_pvalue <- apply(
      rhythmic_genes[
        ,
        padj_columns,
        drop = FALSE
      ],
      1L,
      function(x) {
        
        x <- suppressWarnings(
          as.numeric(x)
        )
        
        if (all(is.na(x))) {
          return(NA_real_)
        }
        
        min(
          x,
          na.rm = TRUE
        )
      }
    )
    
    rhythmic_genes <- rhythmic_genes[
      order(
        rhythmic_genes$best_adjusted_pvalue,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
  }
  
  rhythmicity_summary <- as.data.frame(
    table(
      combined_results$rhythmic_by,
      useNA = "ifany"
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(rhythmicity_summary) <- c(
    "rhythmic_by",
    "n_genes"
  )
  
  expression_export <- data.frame(
    gene_id = rownames(
      prepared$expr
    ),
    symbol = prepared$gene_annotation$symbol[
      match(
        rownames(prepared$expr),
        prepared$gene_annotation$matrix_gene_id
      )
    ],
    prepared$expr,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  utils::write.csv(
    expression_export,
    file = file.path(
      analysis_output_dir,
      paste0(
        "GSE98965_",
        tissue_code,
        "_log2_FPKM_expression.csv"
      )
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    prepared$metadata,
    file = file.path(
      analysis_output_dir,
      paste0(
        "GSE98965_",
        tissue_code,
        "_metadata.csv"
      )
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    prepared$gene_annotation,
    file = file.path(
      analysis_output_dir,
      paste0(
        "GSE98965_",
        tissue_code,
        "_gene_annotation.csv"
      )
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    prepared$filter_summary,
    file = file.path(
      analysis_output_dir,
      paste0(
        "GSE98965_",
        tissue_code,
        "_filter_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    combined_results,
    file = file.path(
      analysis_output_dir,
      paste0(
        "GSE98965_",
        tissue_code,
        "_combined_rhythmicity_results.csv"
      )
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    rhythmic_genes,
    file = file.path(
      analysis_output_dir,
      paste0(
        "GSE98965_",
        tissue_code,
        "_rhythmic_genes.csv"
      )
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    rhythmicity_summary,
    file = file.path(
      analysis_output_dir,
      paste0(
        "GSE98965_",
        tissue_code,
        "_rhythmicity_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  analysis_summary <- data.frame(
    metric = c(
      "Dataset",
      "Tissue",
      "Number of samples",
      "Number of unique time points",
      "Time points",
      "First time point",
      "Last time point",
      "Genes analysed",
      "Rhythmic genes",
      "Adjusted p-value cutoff",
      "Output directory"
    ),
    value = c(
      "GSE98965",
      tissue_code,
      nrow(prepared$metadata),
      length(
        unique(
          prepared$metadata$time
        )
      ),
      paste(
        sort(
          unique(
            prepared$metadata$time
          )
        ),
        collapse = ", "
      ),
      min(
        prepared$metadata$time,
        na.rm = TRUE
      ),
      max(
        prepared$metadata$time,
        na.rm = TRUE
      ),
      nrow(prepared$expr),
      nrow(rhythmic_genes),
      padj_cutoff,
      normalizePath(
        analysis_output_dir,
        winslash = "/",
        mustWork = FALSE
      )
    ),
    stringsAsFactors = FALSE
  )
  
  utils::write.csv(
    analysis_summary,
    file = file.path(
      analysis_output_dir,
      paste0(
        "GSE98965_",
        tissue_code,
        "_analysis_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  structure(
    list(
      summary = combined_results,
      rhythmic_genes = rhythmic_genes,
      rhythmicity_summary = rhythmicity_summary,
      analysis_summary = analysis_summary,
      prepared = prepared,
      deggo = rhythm_results,
      output_dir = analysis_output_dir
    ),
    class = c(
      "deggo_public_circadian",
      "list"
    )
  )
}