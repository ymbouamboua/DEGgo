
#' Prepare count matrix and metadata with matched samples
#'
#' @param counts A data.frame or matrix of raw counts.
#' @param metadata A data.frame containing sample information.
#' @param gene_col Name of the gene ID column in counts.
#' @param sample_col Name of the sample column in metadata.
#' @param keep_common_only Logical. If TRUE, keep only samples present in both.
#' @param verbose Logical. Print summary messages.
#'
#' @return A list with matched count matrix and metadata.
#' @export
prepare_counts_metadata <- function(
    counts,
    metadata,
    gene_col = c("GeneID", "gene", "Gene", "ENSEMBL"),
    sample_col = c("sample", "Sample", "SAMPLE"),
    keep_common_only = TRUE,
    verbose = TRUE
) {
  
  if (!is.data.frame(counts) && !is.matrix(counts)) {
    stop("counts must be a data.frame or matrix.", call. = FALSE)
  }
  
  if (!is.data.frame(metadata)) {
    metadata <- as.data.frame(metadata)
  }
  
  counts <- as.data.frame(counts, check.names = FALSE)
  
  gene_col <- gene_col[gene_col %in% colnames(counts)][1]
  sample_col <- sample_col[sample_col %in% colnames(metadata)][1]
  
  if (is.na(gene_col)) {
    stop("No gene column found in counts.", call. = FALSE)
  }
  
  if (is.na(sample_col)) {
    stop("No sample column found in metadata.", call. = FALSE)
  }
  
  gene_ids <- counts[[gene_col]]
  
  counts_mat <- counts[, setdiff(colnames(counts), gene_col), drop = FALSE]
  rownames(counts_mat) <- gene_ids
  
  counts_mat <- as.matrix(counts_mat)
  storage.mode(counts_mat) <- "integer"
  
  metadata[[sample_col]] <- as.character(metadata[[sample_col]])
  
  common_samples <- intersect(colnames(counts_mat), metadata[[sample_col]])
  
  if (!length(common_samples)) {
    stop("No common samples between counts and metadata.", call. = FALSE)
  }
  
  missing_in_meta <- setdiff(colnames(counts_mat), metadata[[sample_col]])
  missing_in_counts <- setdiff(metadata[[sample_col]], colnames(counts_mat))
  
  if (!keep_common_only) {
    if (length(missing_in_meta) || length(missing_in_counts)) {
      stop(
        "Sample mismatch.\n",
        "Missing in metadata: ", paste(missing_in_meta, collapse = ", "), "\n",
        "Missing in counts: ", paste(missing_in_counts, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  counts_mat <- counts_mat[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, metadata[[sample_col]]), , drop = FALSE]
  
  rownames(metadata) <- metadata[[sample_col]]
  
  if (!all(colnames(counts_mat) == rownames(metadata))) {
    stop("Final sample order mismatch.", call. = FALSE)
  }
  
  if (verbose) {
    message("Counts genes   : ", nrow(counts_mat))
    message("Counts samples : ", ncol(counts_mat))
    message("Metadata rows  : ", nrow(metadata))
    
    if (length(missing_in_meta)) {
      message("Removed from counts because missing in metadata: ",
              paste(missing_in_meta, collapse = ", "))
    }
    
    if (length(missing_in_counts)) {
      message("Removed from metadata because missing in counts: ",
              paste(missing_in_counts, collapse = ", "))
    }
    
    message("Sample matching: OK")
  }
  
  list(
    counts = counts_mat,
    metadata = metadata,
    gene_col = gene_col,
    sample_col = sample_col,
    common_samples = common_samples,
    missing_in_meta = missing_in_meta,
    missing_in_counts = missing_in_counts
  )
}


#' Validate count matrix and metadata
#'
#' @param counts Count matrix with genes as rows and samples as columns.
#' @param metadata Sample metadata with samples as rows.
#' @param condition_col Column containing the experimental condition.
#'
#' @return Invisibly returns TRUE.
#' @keywords internal
#' @noRd
validate_inputs <- function(
    counts,
    metadata,
    condition_col = "condition"
) {
  
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("'counts' must be a matrix or data.frame.", call. = FALSE)
  }
  
  counts <- as.matrix(counts)
  
  if (!is.numeric(counts) && !is.integer(counts)) {
    stop("'counts' must contain numeric or integer values.", call. = FALSE)
  }
  
  if (anyNA(counts)) {
    stop("'counts' contains missing values.", call. = FALSE)
  }
  
  if (any(counts < 0)) {
    stop("'counts' contains negative values.", call. = FALSE)
  }
  
  if (is.null(rownames(counts))) {
    stop("'counts' must have gene IDs as rownames.", call. = FALSE)
  }
  
  if (is.null(colnames(counts))) {
    stop("'counts' must have sample names as colnames.", call. = FALSE)
  }
  
  if (!is.data.frame(metadata)) {
    metadata <- as.data.frame(metadata)
  }
  
  if (is.null(rownames(metadata))) {
    stop("'metadata' must have sample names as rownames.", call. = FALSE)
  }
  
  if (!condition_col %in% colnames(metadata)) {
    stop(
      "Metadata must contain column: ",
      condition_col,
      call. = FALSE
    )
  }
  
  if (!all(colnames(counts) %in% rownames(metadata))) {
    missing <- setdiff(colnames(counts), rownames(metadata))
    stop(
      "Samples missing in metadata: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!all(colnames(counts) == rownames(metadata))) {
    stop(
      "Sample order mismatch between counts columns and metadata rownames. ",
      "Run prepare_counts_metadata() first.",
      call. = FALSE
    )
  }
  
  if (length(unique(metadata[[condition_col]])) < 2) {
    stop(
      "Metadata column '",
      condition_col,
      "' must contain at least two groups.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


# =========================================================
# CLEAN ENSEMBL IDS
# =========================================================

clean_ensembl_ids <- function(counts) {
  
  log_msg("Cleaning ENSEMBL IDs...")
  
  counts <- as.matrix(counts)
  sample_names <- colnames(counts)
  
  clean_ids <- sub("\\..*$", "", rownames(counts))
  rownames(counts) <- clean_ids
  
  counts <- rowsum(
    x = counts,
    group = rownames(counts),
    reorder = FALSE
  )
  
  colnames(counts) <- sample_names
  storage.mode(counts) <- "integer"
  
  counts
}

# =========================================================
# PREPROCESS COUNTS
# =========================================================

preprocess_counts <- function(
    counts,
    metadata = NULL,
    filter_method = c("count", "cpm", "none"),
    min_count = 5,
    min_samples = 2,
    min_total = 10
) {
  
  log_msg("Preprocessing counts...")
  
  filter_method <- match.arg(filter_method)
  
  counts <- round(as.matrix(counts))
  mode(counts) <- "integer"
  
  if (filter_method == "none") {
    log_msg("No gene filtering applied.", type = "warn")
    return(counts)
  }
  
  if (filter_method == "cpm") {
    
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("Package 'edgeR' is required for CPM filtering.", call. = FALSE)
    }
    
    if (is.null(metadata) || !"condition" %in% colnames(metadata)) {
      stop("metadata with 'condition' column is required for CPM filtering.", call. = FALSE)
    }
    
    keep <- edgeR::filterByExpr(
      counts,
      group = metadata$condition
    )
    
  } else {
    
    keep <- (
      rowSums(counts >= min_count) >= min_samples
    ) &
      (
        rowSums(counts) >= min_total
      )
  }
  
  log_msg(
    paste0(
      "Genes before filtering: ", nrow(counts),
      " | Genes retained: ", sum(keep),
      " | Filter method: ", filter_method
    )
  )
  
  counts[keep, , drop = FALSE]
}