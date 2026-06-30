#=========================================================#
# CHECK RAW COUNTS
#=========================================================#

#' Validate raw count table
#'
#' Performs sanity checks on bulk RNA-seq count tables.
#'
#' @param counts Raw count table.
#' @param gene_col Possible gene ID columns.
#' @param feature_col Possible gene symbol columns.
#' @param min_library_size Minimum acceptable library size.
#' @param min_detected_genes Minimum acceptable detected genes.
#'
#' @return List of QC results.
#' @export

check_raw_counts <- function(
    counts,
    gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl"),
    feature_col = c("gene_name", "SYMBOL", "symbol"),
    min_library_size = 1e6,
    min_detected_genes = 5000
) {
  counts <- as.data.frame(counts, check.names = FALSE)

  gcol <- gene_col[gene_col %in% colnames(counts)][1]
  fcol <- feature_col[feature_col %in% colnames(counts)][1]

  if (is.na(gcol)) stop("No gene ID column found.", call. = FALSE)

  sample_cols <- colnames(counts)[!colnames(counts) %in% c(gcol, fcol)]

  mat <- as.matrix(counts[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"

  # -------------------------------------------------- #
  # Sample name checks
  # -------------------------------------------------- #

  dup_idx <- duplicated(colnames(counts)) | duplicated(colnames(counts), fromLast = TRUE)

  dup_cols <- data.frame(
    position = which(dup_idx),
    column = colnames(counts)[dup_idx],
    stringsAsFactors = FALSE
  )

  dup_samples <- unique(sample_cols[duplicated(sample_cols)])

  # -------------------------------------------------- #
  # Duplicated expression profiles
  # -------------------------------------------------- #

  dup_genes <- unique(counts[[gcol]][duplicated(counts[[gcol]])])

  # -------------------------------------------------- #
  # Count checks
  # -------------------------------------------------- #

  missing_values <- sum(is.na(mat))
  negative_values <- sum(mat < 0, na.rm = TRUE)
  non_integer <- sum(abs(mat - round(mat)) > 1e-6, na.rm = TRUE)

  # -------------------------------------------------- #
  # Library QC
  # -------------------------------------------------- #

  qc <- data.frame(
    sample = colnames(mat),
    library_size = colSums(mat, na.rm = TRUE),
    detected_genes = colSums(mat > 0, na.rm = TRUE),
    low_library = colSums(mat, na.rm = TRUE) < min_library_size,
    low_detected = colSums(mat > 0, na.rm = TRUE) < min_detected_genes,
    stringsAsFactors = FALSE
  )

  # -------------------------------------------------- #
  # Report
  # -------------------------------------------------- #

  cat("Genes: ", nrow(mat), "\n")
  cat("Samples: ", ncol(mat), "\n")
  cat("Duplicated genes: ", length(dup_genes), "\n")
  cat("Duplicated sample names: ", length(dup_samples), "\n")
  cat("Missing values: ", missing_values, "\n")
  cat("Negative counts: ", negative_values, "\n")
  cat("Non-integer counts: ", non_integer, "\n")

  invisible(list(
    qc = qc,
    gene_col = gcol,
    feature_col = fcol,
    duplicated_columns = dup_cols,
    duplicated_samples = dup_samples,
    duplicated_gene_ids = dup_genes,
    missing_values = missing_values,
    negative_counts = negative_values,
    non_integer_counts = non_integer
  ))
}



# ========================================================= #
# PREPARE COUNTS MATRIX
# ========================================================= #
#' Prepare count matrix for DEGgo analysis
#'
#' Converts a count table containing gene annotation columns and sample
#' count columns into a clean numeric matrix suitable for downstream
#' differential expression analysis.
#'
#' The function:
#' \enumerate{
#'   \item Extracts sample columns defined in the metadata.
#'   \item Verifies that all metadata samples are present in the count table.
#'   \item Converts count columns to numeric values.
#'   \item Uses a gene identifier column as row names when available.
#'   \item Returns a numeric count matrix compatible with DESeq2,
#'   edgeR, and limma workflows.
#' }
#'
#' @param counts A count table containing gene annotation columns and
#'   sample count columns.
#' @param metadata A metadata data frame containing sample information.
#' @param sample_col Column name in \code{metadata} containing sample IDs.
#'   Default is \code{"sample"}.
#' @param gene_cols Character vector of possible gene identifier columns.
#'   The first matching column found in \code{counts} is used as row names.
#'   Default includes common gene ID and gene symbol columns.
#'
#' @return A numeric matrix with:
#' \describe{
#'   \item{rows}{Genes.}
#'   \item{columns}{Samples.}
#' }
#'
#' @details
#' This helper is used internally by \code{\link{run_deggo}} to ensure
#' that only numeric sample count columns are retained before differential
#' expression analysis.
#'
#' If one of the specified \code{gene_cols} is found, it is used as the
#' row names of the returned matrix. Duplicate gene identifiers are made
#' unique using \code{\link[base]{make.unique}}.
#'
#' The function stops if:
#' \itemize{
#'   \item Metadata samples are missing from the count table.
#'   \item Non-numeric values remain after coercion.
#' }
#'
#' @examples
#' counts_mat <- prepare_deggo_counts(
#'   counts = counts,
#'   metadata = metadata
#' )
#'
#' dim(counts_mat)
#' head(colnames(counts_mat))
#' head(rownames(counts_mat))
#'
#' @seealso
#' \code{\link{prepare_counts_metadata}},
#' \code{\link{run_deggo}}
#'
#' @keywords internal
#' @noRd
prepare_deggo_counts <- function(
    counts,
    metadata,
    sample_col = "sample",
    gene_cols = c(
      "gene_id",
      "gene_name",
      "GeneID",
      "SYMBOL",
      "symbol"
    )
) {

  counts <- as.data.frame(counts, check.names = FALSE)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  sample_ids <- metadata[[sample_col]]

  missing <- setdiff(sample_ids, colnames(counts))

  if (length(missing) > 0) {
    stop(
      "Samples missing from counts: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  count_mat <- counts[, sample_ids, drop = FALSE]

  count_mat[] <- lapply(count_mat, function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })

  if (anyNA(count_mat)) {
    stop(
      "Non-numeric values detected in count columns after coercion.",
      call. = FALSE
    )
  }

  present_gene_cols <- intersect(
    gene_cols,
    colnames(counts)
  )

  if (length(present_gene_cols) > 0) {
    row_ids <- counts[[present_gene_cols[1]]]
  } else {
    row_ids <- rownames(counts)
  }

  rownames(count_mat) <- make.unique(
    as.character(row_ids)
  )

  count_mat <- as.matrix(count_mat)
  storage.mode(count_mat) <- "numeric"

  count_mat
}




# ========================================================= #
# PREPARE COUNTS METADATA
# ========================================================= #
#' Prepare count matrix and metadata with matched samples
#'
#' @param counts Raw count table or matrix.
#' @param metadata Sample metadata data frame.
#' @param gene_col Possible gene identifier column names in counts.
#' @param feature_col Possible gene symbol/annotation column names in counts.
#' @param sample_col Possible sample identifier column names in metadata.
#' @param keep_common_only Logical. If TRUE, keeps only samples present
#'   in both counts and metadata.
#' @param verbose Logical. Print preparation summary.
#'
#' @return A list containing prepared counts, metadata, feature map,
#'   selected column names, common samples, and missing sample names.
#'
#' @seealso
#' See also the internal input validation used by DEGgo.
#'
#' @keywords internal
#' @noRd
prepare_counts_metadata <- function(
    counts,
    metadata,
    gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_id"),
    feature_col = c("gene_name", "SYMBOL", "symbol", "gene_symbol", "external_gene_name"),
    sample_col = c("sample", "Sample", "SAMPLE"),
    keep_common_only = TRUE,
    verbose = TRUE
) {

  if (!is.data.frame(counts) && !is.matrix(counts)) {
    stop("counts must be a data.frame or matrix.", call. = FALSE)
  }

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  counts <- as.data.frame(counts, check.names = FALSE)

  gene_col <- gene_col[gene_col %in% colnames(counts)][1]
  feature_col <- feature_col[feature_col %in% colnames(counts)][1]
  sample_col <- sample_col[sample_col %in% colnames(metadata)][1]

  if (is.na(gene_col)) {
    stop("No gene identifier column found in counts.", call. = FALSE)
  }

  if (is.na(sample_col)) {
    stop("No sample column found in metadata.", call. = FALSE)
  }

  gene_ids <- as.character(counts[[gene_col]])
  gene_ids <- make.unique(gene_ids)

  annotation_cols <- unique(c(gene_col, feature_col))
  annotation_cols <- annotation_cols[!is.na(annotation_cols)]

  count_df <- counts[
    ,
    setdiff(colnames(counts), annotation_cols),
    drop = FALSE
  ]

  metadata[[sample_col]] <- as.character(metadata[[sample_col]])

  common_samples <- intersect(
    colnames(count_df),
    metadata[[sample_col]]
  )

  if (!length(common_samples)) {
    stop("No common samples between counts and metadata.", call. = FALSE)
  }

  missing_in_meta <- setdiff(colnames(count_df), metadata[[sample_col]])
  missing_in_counts <- setdiff(metadata[[sample_col]], colnames(count_df))

  if (!keep_common_only && (length(missing_in_meta) > 0 || length(missing_in_counts) > 0)) {
    stop(
      "Sample mismatch.\n",
      "Missing in metadata: ", paste(missing_in_meta, collapse = ", "), "\n",
      "Missing in counts: ", paste(missing_in_counts, collapse = ", "),
      call. = FALSE
    )
  }

  count_df <- count_df[, common_samples, drop = FALSE]

  count_df[] <- lapply(count_df, function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })

  if (anyNA(count_df)) {
    bad_cols <- names(which(vapply(count_df, function(x) anyNA(x), logical(1))))

    stop(
      "Non-numeric values detected in count columns after coercion: ",
      paste(bad_cols, collapse = ", "),
      call. = FALSE
    )
  }

  counts_mat <- as.matrix(count_df)
  rownames(counts_mat) <- gene_ids
  counts_mat <- round(counts_mat)
  storage.mode(counts_mat) <- "integer"

  metadata <- metadata[
    match(common_samples, metadata[[sample_col]]),
    ,
    drop = FALSE
  ]

  rownames(metadata) <- metadata[[sample_col]]

  if (!"sample" %in% colnames(metadata)) {
    metadata$sample <- metadata[[sample_col]]
  }

  feature_map <- NULL

  if (!is.na(feature_col)) {
    feature_map <- data.frame(
      ENSEMBL = gene_ids,
      SYMBOL = as.character(counts[[feature_col]]),
      stringsAsFactors = FALSE
    )
  }

  if (verbose) {
    cat("Gene ID column : ", gene_col, "\n")
    if (!is.na(feature_col)) {
      cat("Feature column : ", feature_col, "\n")
    }
    cat("Sample column  : ", sample_col, "\n")
    cat("Counts genes   : ", nrow(counts_mat), "\n")
    cat("Counts samples : ", ncol(counts_mat), "\n")
    cat("Metadata rows  : ", nrow(metadata), "\n")
    cat("Sample matching: OK", "\n")
  }

  list(
    counts = counts_mat,
    metadata = metadata,
    gene_col = gene_col,
    feature_col = feature_col,
    feature_map = feature_map,
    sample_col = sample_col,
    common_samples = common_samples,
    missing_in_meta = missing_in_meta,
    missing_in_counts = missing_in_counts
  )
}


# ========================================================= #
# VALIDATE INPUT
# ========================================================= #

#' Validate count matrix and metadata
#'
#' @param counts Count matrix with genes as rows and samples as columns.
#' @param metadata Sample metadata with samples as rows.
#' @param condition_col Column containing the experimental condition.
#'
#' @return Invisibly returns TRUE.
#' @keywords internal
#' @noRd
# validate_inputs <- function(
#     counts,
#     metadata,
#     condition_col = "condition"
# ) {
#
#   if (!is.matrix(counts) && !is.data.frame(counts)) {
#     stop("'counts' must be a matrix or data.frame.", call. = FALSE)
#   }
#
#   counts <- as.matrix(counts)
#
#   if (!is.numeric(counts) && !is.integer(counts)) {
#     stop("'counts' must contain numeric or integer values.", call. = FALSE)
#   }
#
#   if (anyNA(counts)) {
#     stop("'counts' contains missing values.", call. = FALSE)
#   }
#
#   if (any(counts < 0)) {
#     stop("'counts' contains negative values.", call. = FALSE)
#   }
#
#   if (is.null(rownames(counts))) {
#     stop("'counts' must have gene IDs as rownames.", call. = FALSE)
#   }
#
#   if (is.null(colnames(counts))) {
#     stop("'counts' must have sample names as colnames.", call. = FALSE)
#   }
#
#   if (!is.data.frame(metadata)) {
#     metadata <- as.data.frame(metadata)
#   }
#
#   if (is.null(rownames(metadata))) {
#     stop("'metadata' must have sample names as rownames.", call. = FALSE)
#   }
#
#   if (!condition_col %in% colnames(metadata)) {
#     stop(
#       "Metadata must contain column: ",
#       condition_col,
#       call. = FALSE
#     )
#   }
#
#   if (!all(colnames(counts) %in% rownames(metadata))) {
#     missing <- setdiff(colnames(counts), rownames(metadata))
#     stop(
#       "Samples missing in metadata: ",
#       paste(missing, collapse = ", "),
#       call. = FALSE
#     )
#   }
#
#   if (!all(colnames(counts) == rownames(metadata))) {
#     stop(
#       "Sample order mismatch between counts columns and metadata rownames. ",
#       "Run prepare_counts_metadata() first.",
#       call. = FALSE
#     )
#   }
#
#   if (length(unique(metadata[[condition_col]])) < 2) {
#     stop(
#       "Metadata column '",
#       condition_col,
#       "' must contain at least two groups.",
#       call. = FALSE
#     )
#   }
#
#   invisible(TRUE)
# }

validate_inputs <- function(counts, metadata, design_formula = NULL, contrast = NULL) {

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame.")
  }

  design_vars <- character(0)

  if (!is.null(design_formula)) {
    design_vars <- all.vars(design_formula)
  }

  contrast_var <- NULL
  if (!is.null(contrast)) {
    contrast_var <- contrast[1]
  }

  required_cols <- unique(c(design_vars, contrast_var))

  missing_cols <- setdiff(required_cols, colnames(metadata))

  if (length(missing_cols) > 0) {
    stop(
      "Metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  TRUE
}

# ========================================================= #
# CLEAN ENSEMBL IDS
# ========================================================= #


#' Clean and collapse Ensembl gene identifiers
#'
#' Removes Ensembl version suffixes and collapses duplicated
#' gene identifiers by summing counts across rows.
#'
#' For example:
#'
#' \preformatted{
#' ENSG000001234.5 -> ENSG000001234
#' }
#'
#' After cleaning, duplicated identifiers are merged using
#' \code{\link[base]{rowsum}}.
#'
#' @param counts Count matrix with genes as rows and samples as columns.
#'
#' @return Integer count matrix with cleaned Ensembl identifiers
#'   as row names.
#'
#' @details
#' Ensembl version numbers are removed using:
#'
#' \code{sub("\\\\..*$", "", gene_id)}
#'
#' Duplicate Ensembl identifiers generated after version removal
#' are collapsed by summing counts across rows.
#'
#' @seealso
#' \code{\link{prepare_counts_metadata}},
#' \code{\link{preprocess_counts}}
#'
#' @keywords internal
#' @noRd
clean_ensembl_ids <- function(counts) {

  log <- .deggo_msg(verbose = TRUE)

  log("Cleaning ENSEMBL IDs...")

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


# ========================================================= #
# PROCESS COUNTS RNA-SEQ MATRIX
# ========================================================= #

#' Preprocess RNA-seq count matrix
#'
#' Filters lowly expressed genes prior to differential
#' expression analysis.
#'
#' Three filtering strategies are supported:
#'
#' \itemize{
#'   \item \code{"count"}: retain genes expressed above a minimum
#'   count threshold in a minimum number of samples.
#'   \item \code{"cpm"}: use
#'   \code{edgeR::filterByExpr()}.
#'   \item \code{"none"}: no filtering.
#' }
#'
#' @param counts Raw count matrix.
#' @param metadata Sample metadata. Required when
#'   \code{filter_method = "cpm"}.
#' @param filter_method Filtering method.
#'   One of \code{"count"}, \code{"cpm"}, or \code{"none"}.
#' @param min_count Minimum count threshold for count-based filtering.
#' @param min_samples Minimum number of samples passing
#'   \code{min_count}.
#' @param min_total Minimum total counts across all samples.
#'
#' @return Filtered integer count matrix.
#'
#' @details
#' Count-based filtering retains genes satisfying:
#'
#' \deqn{
#' \sum(count_i \ge min\_count) \ge min\_samples
#' }
#'
#' and
#'
#' \deqn{
#' \sum count_i \ge min\_total
#' }
#'
#' CPM-based filtering uses
#' \code{edgeR::filterByExpr()}, which accounts for
#' library sizes and experimental groups.
#'
#' @seealso
#' \code{\link{prepare_counts_metadata}},
#' \code{\link{validate_inputs}}
#'
#' @keywords internal
#' @noRd
preprocess_counts <- function(
    counts,
    metadata = NULL,
    filter_method = c("count", "cpm", "none"),
    min_count = 5,
    min_samples = 2,
    min_total = 10
) {

  log <- .deggo_msg(verbose = TRUE)

  log("Preprocessing counts...")

  filter_method <- match.arg(filter_method)

  counts <- round(as.matrix(counts))
  mode(counts) <- "integer"

  if (filter_method == "none") {
    log("No gene filtering applied.", type = "warn")
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

  log(
    paste0(
      "Genes before filtering: ", nrow(counts),
      " | Genes retained: ", sum(keep),
      " | Filter method: ", filter_method
    )
  )

  counts[keep, , drop = FALSE]
}


# ========================================================= #
# REMOVE FLAGGED SAMPLES
# ========================================================= #
#' Remove flagged samples from counts and metadata
#'
#' Removes samples flagged for exclusion in a QC table and returns
#' synchronized counts and metadata objects.
#'
#' @param counts Count matrix or data frame. Rows are genes and columns
#'   are samples. Gene annotation columns are preserved.
#' @param metadata Sample metadata.
#' @param qc_table QC table containing sample identifiers and removal flags.
#' @param sample_col Column containing sample names. Default is
#'   \code{"sample_original"}.
#' @param remove_col Logical column indicating samples to remove.
#'   Default is \code{"recommend_remove"}.
#' @param gene_cols Character vector of annotation columns to preserve.
#'   Default is \code{c("gene_id", "gene_name")}.
#' @param verbose Logical; print summary messages. Default is \code{TRUE}.
#'
#' @return A list with counts, metadata, and removed sample names.
#'
#' @examples
#' \dontrun{
#' qc <- explore_bulk_rnaseq(counts, metadata)
#'
#' cleaned <- remove_flagged_samples(
#'   counts = counts,
#'   metadata = metadata,
#'   qc_table = qc$qc
#' )
#'
#' counts <- cleaned$counts
#' metadata <- cleaned$metadata
#' }
#'
#' @export
remove_flagged_samples <- function(
    counts,
    metadata,
    qc_table,
    sample_col = "sample_original",
    remove_col = "recommend_remove",
    gene_cols = c("gene_id", "gene_name"),
    verbose = TRUE
) {

  counts <- as.data.frame(counts, check.names = FALSE)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  qc_table <- as.data.frame(qc_table, stringsAsFactors = FALSE)

  if (!sample_col %in% colnames(qc_table)) {
    stop("Column not found in qc_table: ", sample_col, call. = FALSE)
  }

  if (!remove_col %in% colnames(qc_table)) {
    stop("Column not found in qc_table: ", remove_col, call. = FALSE)
  }

  if (!sample_col %in% colnames(metadata)) {
    stop("Column not found in metadata: ", sample_col, call. = FALSE)
  }

  remove_samples <- qc_table[[sample_col]][
    qc_table[[remove_col]] %in% TRUE
  ]

  remove_samples <- unique(as.character(remove_samples))

  metadata_clean <- metadata[
    !metadata[[sample_col]] %in% remove_samples,
    ,
    drop = FALSE
  ]

  keep_cols <- !colnames(counts) %in% remove_samples
  counts_clean <- counts[, keep_cols, drop = FALSE]

  keep_gene_cols <- intersect(gene_cols, colnames(counts_clean))
  counts_clean <- counts_clean[
    ,
    unique(c(keep_gene_cols, setdiff(colnames(counts_clean), keep_gene_cols))),
    drop = FALSE
  ]

  if (verbose) {
    cat("Removed samples:", length(remove_samples), "\n")

    if (length(remove_samples) > 0) {
      cat(paste(remove_samples, collapse = ", "), "\n")
    }

    cat("Remaining samples:", nrow(metadata_clean), "\n")
  }

  list(
    counts = counts_clean,
    metadata = metadata_clean,
    removed_samples = remove_samples
  )
}



# ============================================================ #
# .deggo_qc_prepare
# ============================================================ #
#' Prepare counts and metadata for DEGgo quality control
#'
#' Internal helper that standardizes count matrices and metadata before
#' exploratory quality control analyses.
#'
#' The function:
#' \itemize{
#'   \item Ensures an explicit gene identifier column is present.
#'   \item Detects gene, feature and sample identifier columns.
#'   \item Matches samples shared between the count table and metadata.
#'   \item Extracts a numeric count matrix suitable for QC functions.
#'   \item Returns harmonized metadata aligned with the count matrix.
#' }
#'
#' This helper is used by QC-related functions such as
#' `explore_bulk_rnaseq()`, `run_sample_qc()`, and other visualization
#' utilities that require synchronized counts and metadata.
#'
#' @param counts Count table, data frame or matrix.
#' @param metadata Sample metadata data frame.
#' @param gene_col Character vector of possible gene identifier column names.
#' @param feature_col Character vector of possible feature annotation column
#'   names (gene symbols, external gene names, etc.).
#' @param sample_col Character vector of possible sample identifier column
#'   names in `metadata`.
#'
#' @return
#' A named list containing:
#' \describe{
#'   \item{mat}{Numeric count matrix with genes as rows and samples as columns.}
#'   \item{metadata}{Metadata reordered to match the count matrix columns.}
#'   \item{counts_df}{Original count table converted to a data frame.}
#'   \item{gene_col}{Detected gene identifier column.}
#'   \item{feature_col}{Detected feature annotation column, if available.}
#'   \item{sample_col}{Detected sample identifier column.}
#' }
#'
#' @details
#' The returned count matrix contains only samples present in both the count
#' table and metadata. Gene identifiers are preserved exactly as provided and
#' are assigned as row names of the returned matrix.
#'
#' This function performs no expression filtering or normalization; it simply
#' harmonizes the inputs for downstream quality control.
#'
#' @keywords internal
#' @noRd
.deggo_qc_prepare <- function(
    counts,
    metadata,
    gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_id"),
    feature_col = c("gene_name", "SYMBOL", "symbol", "gene_symbol", "external_gene_name"),
    sample_col = c("sample", "Sample", "SAMPLE")
) {

  counts <- .deggo_counts_with_gene_col(counts, gene_col = gene_col)
  counts <- as.data.frame(counts, check.names = FALSE)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  gene_col <- gene_col[gene_col %in% colnames(counts)][1]
  feature_col <- feature_col[feature_col %in% colnames(counts)][1]
  sample_col <- sample_col[sample_col %in% colnames(metadata)][1]

  if (is.na(gene_col)) stop("No gene ID column found.", call. = FALSE)
  if (is.na(sample_col)) stop("No sample column found.", call. = FALSE)

  rownames(metadata) <- as.character(metadata[[sample_col]])

  drop_cols <- c(gene_col, feature_col)
  drop_cols <- drop_cols[!is.na(drop_cols)]

  sample_cols <- intersect(setdiff(colnames(counts), drop_cols), rownames(metadata))

  if (!length(sample_cols)) {
    stop("No matching samples between counts and metadata.", call. = FALSE)
  }

  mat <- as.matrix(counts[, sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- as.character(counts[[gene_col]])

  metadata <- metadata[sample_cols, , drop = FALSE]

  list(
    mat = mat,
    metadata = metadata,
    counts_df = counts,
    gene_col = gene_col,
    feature_col = feature_col,
    sample_col = sample_col
  )
}


# ============================================================ #
# .deggo_counts_with_gene_cols
# ============================================================ #
#' Ensure a gene identifier column is present
#'
#' Internal helper that guarantees a count table contains an explicit gene
#' identifier column. If no recognized gene identifier column is found, the
#' function converts the row names into a new `"gene_id"` column while
#' preserving the original count data.
#'
#' This helper is primarily used before automatic input preparation and raw
#' quality control to standardize user inputs originating from either matrices
#' (gene IDs stored as row names) or data frames (gene IDs stored as columns).
#'
#' @param counts Count matrix or data frame.
#' @param gene_col Character vector of acceptable gene identifier column names.
#'
#' @return
#' A data frame containing an explicit gene identifier column. If one already
#' exists, the input is returned unchanged (except for coercion to a data
#' frame when necessary).
#'
#' @details
#' If none of the columns specified in `gene_col` are found, the row names are
#' used to create a new `"gene_id"` column. An error is thrown if neither a
#' valid gene identifier column nor row names are available.
#'
#' @keywords internal
#' @noRd
.deggo_counts_with_gene_col <- function(counts, gene_col = "gene_id") {

  if (!is.matrix(counts) && !is.data.frame(counts)) {
    return(counts)
  }

  counts <- as.data.frame(counts, check.names = FALSE)

  if (length(intersect(gene_col, colnames(counts))) > 0) {
    return(counts)
  }

  gene_id <- rownames(counts)

  if (is.null(gene_id) || all(is.na(gene_id)) || all(gene_id == "")) {
    stop(
      "No gene ID column found and row names are empty. ",
      "Please provide a gene identifier column using 'gene_col'.",
      call. = FALSE
    )
  }

  rownames(counts) <- NULL

  data.frame(
    gene_id = gene_id,
    counts,
    check.names = FALSE
  )
}




# ======================================================= #
# Sample QC: correlation and hierarchical clustering
# ======================================================= #
#' Run sample-level RNA-seq quality control
#'
#' Generate sample-level quality control outputs from a raw or filtered count
#' matrix, including a Spearman sample correlation heatmap, hierarchical
#' clustering dendrogram, and sample correlation matrix.
#'
#' This function is designed to be used inside the DEGgo workflow after
#' count filtering and metadata matching. It saves all QC files directly into
#' the provided output directory.
#'
#' @param counts Numeric count matrix with genes in rows and samples in columns.
#' @param metadata Sample metadata data frame. Row names must match sample names
#'   in `colnames(counts)`.
#' @param output_dir Output directory where QC files will be saved.
#' @param annotation_cols Character vector of metadata columns used as heatmap
#'   annotations.
#' @param dpi Resolution used for the hierarchical clustering PNG.
#' @param style Plot theme style.
#' @param txtsize Base text size.
#'
#' @return Invisibly returns a list containing the sample correlation matrix,
#'   output directory, and paths to generated files.
#'
#' @details
#' The function uses `log2(count + 1)` transformed counts to compute Spearman
#' sample correlations and hierarchical clustering based on Euclidean distance.
#'
#' @examples
#' \dontrun{
#' qc <- run_sample_qc(
#'   counts = counts,
#'   metadata = metadata,
#'   output_dir = "pairwise_QC"
#' )
#' }
#'
#' @export
#'
run_sample_qc <- function(
    counts,
    metadata,
    output_dir,
    annotation_cols = c("condition", "treatment", "sex", "tissue"),
    dpi = 300,
    style = "bw",
    txtsize = 12
) {
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required.", call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  mat <- as.matrix(counts)
  suppressWarnings(storage.mode(mat) <- "numeric")

  if (!is.numeric(mat)) {
    stop("counts must be numeric.", call. = FALSE)
  }

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (is.null(rownames(metadata))) {
    stop("metadata must have row names matching sample names.", call. = FALSE)
  }

  common_samples <- intersect(colnames(mat), rownames(metadata))

  if (length(common_samples) < 2) {
    stop("At least two matched samples are required for sample QC.", call. = FALSE)
  }

  mat <- mat[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  meta_cols <- intersect(annotation_cols, colnames(metadata))

  sample_metrics <- data.frame(
    sample = colnames(mat),
    library_size = colSums(mat, na.rm = TRUE),
    detected_genes = colSums(mat > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (length(meta_cols) > 0) {
    sample_metrics <- cbind(
      sample_metrics,
      metadata[sample_metrics$sample, meta_cols, drop = FALSE]
    )
  }

  color_priority <- c("condition", "treatment", "sex", "tissue")
  color_col <- intersect(color_priority, colnames(sample_metrics))[1] %||% NULL

  if (!is.null(color_col)) {
    sample_metrics[[color_col]] <- factor(sample_metrics[[color_col]])

    color_df <- sample_metrics[, color_col, drop = FALSE]
    color_list <- .deggo_annotation_colors(color_df)
    color_values <- color_list[[color_col]]

    p_qc <- ggplot2::ggplot(
      sample_metrics,
      ggplot2::aes(
        x = library_size,
        y = detected_genes,
        color = .data[[color_col]]
      )
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      ggplot2::scale_color_manual(values = color_values, drop = FALSE)

  } else {
    p_qc <- ggplot2::ggplot(
      sample_metrics,
      ggplot2::aes(
        x = library_size,
        y = detected_genes
      )
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.9)
  }

  p_qc <- p_qc +
    ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      color = "black",
      linewidth = 0.4
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    ggplot2::labs(
      title = "Detected genes vs library size",
      x = "Library size",
      y = "Detected genes",
      color = color_col
    ) +
    .deggo_theme(style = style, txtsize = txtsize)

  ggplot2::ggsave(
    filename = file.path(output_dir, "Detected_Genes_vs_Library_Size.png"),
    plot = p_qc,
    width = 7,
    height = 5,
    dpi = dpi
  )

  utils::write.table(
    sample_metrics,
    file.path(output_dir, "Sample_QC_Metrics.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  logmat <- log2(mat + 1)

  annotation_col <- metadata[, meta_cols, drop = FALSE]

  if (ncol(annotation_col) > 0) {
    annotation_col[] <- lapply(annotation_col, factor)
  } else {
    annotation_col <- NULL
  }

  ann_colors <- .deggo_annotation_colors(annotation_col)

  cor_mat <- stats::cor(
    logmat,
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  pheatmap::pheatmap(
    cor_mat,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    filename = file.path(output_dir, "Sample_Correlation_Heatmap.png"),
    main = "Sample correlation",
    color = grDevices::colorRampPalette(
      c("#6497b1", "#F7F7F7", "#740001")
    )(100),
    fontsize = 8,
    width = 12,
    height = 10
  )

  grDevices::png(
    filename = file.path(output_dir, "Hierarchical_Clustering.png"),
    width = 12,
    height = 6,
    units = "in",
    res = dpi
  )

  plot(
    stats::hclust(stats::dist(t(logmat))),
    main = "Hierarchical clustering",
    xlab = "",
    sub = ""
  )

  grDevices::dev.off()

  utils::write.table(
    cor_mat,
    file.path(output_dir, "Sample_Correlation_Matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  invisible(list(
    metrics = sample_metrics,
    correlation = cor_mat,
    qc_dir = output_dir,
    files = list(
      detected_genes_library_size = file.path(output_dir, "Detected_Genes_vs_Library_Size.png"),
      correlation_heatmap = file.path(output_dir, "Sample_Correlation_Heatmap.png"),
      clustering = file.path(output_dir, "Hierarchical_Clustering.png"),
      qc_metrics = file.path(output_dir, "Sample_QC_Metrics.tsv"),
      correlation_matrix = file.path(output_dir, "Sample_Correlation_Matrix.tsv")
    )
  ))
}



# ========================================================= #
# MARKER SCORE CHECK
# ========================================================= #
#' Validate sample identity using marker gene signatures
#'
#' Computes signature scores from predefined marker gene sets,
#' predicts the most likely sample identity, and flags potential
#' sample swaps or annotation mismatches.
#'
#' For each marker set, the function calculates the mean expression
#' of all detected marker genes across samples. The predicted group
#' corresponds to the marker set with the highest score.
#'
#' @param counts A data frame or matrix containing gene expression counts.
#'   Rows correspond to genes and columns correspond to samples.
#' @param metadata A data frame containing sample metadata.
#' @param marker_sets Named list of marker gene vectors. Each list element
#'   defines a biological group or tissue signature.
#' @param sample_col Column name in \code{metadata} containing sample IDs.
#'   Default is \code{"sample"}.
#' @param group_col Column name in \code{metadata} containing the expected
#'   sample annotation (e.g. tissue, condition). Default is \code{"tissue"}.
#' @param feature_col Column name in \code{counts} containing gene symbols.
#'   Default is \code{"gene_name"}.
#' @param log_transform Logical; if \code{TRUE}, computes scores from
#'   \code{log2(count + 1)} transformed values. Default is \code{TRUE}.
#' @param plot Logical; if \code{TRUE}, returns a barplot of marker scores.
#'   Default is \code{TRUE}.
#' @param style Theme style passed to \code{.deggo_theme()}.
#' @param txtsize Texte size passed to \code{.deggo_theme()}
#'
#' @return A list containing:
#' \describe{
#'   \item{scores}{Data frame with marker scores, predicted groups,
#'   and swap flags.}
#'   \item{long}{Long-format version of the score table.}
#'   \item{swaps}{Subset of samples where predicted and expected
#'   groups differ.}
#'   \item{plot}{A ggplot object showing marker scores per sample
#'   if \code{plot = TRUE}; otherwise \code{NULL}.}
#' }
#'
#' @details
#' The function performs a simple signature-based classification:
#' \enumerate{
#'   \item Matches marker genes against the feature column.
#'   \item Computes the mean expression of matched genes for each sample.
#'   \item Predicts the group with the highest signature score.
#'   \item Flags samples where predicted and expected annotations differ.
#' }
#'
#' Gene matching is case-insensitive.
#'
#' @examples
#' \dontrun{
#' marker_sets <- list(
#'   BAT = c("Ucp1", "Cidea", "Ppargc1a"),
#'   WAT = c("Adipoq", "Lep", "Fabp4"),
#'   TESTIS = c("Amh", "Sox9", "Ddx4"),
#'   OVARY = c("Foxl2", "Fshr", "Bmp15")
#' )
#'
#' res <- marker_score_check(
#'   counts = counts,
#'   metadata = metadata,
#'   marker_sets = marker_sets,
#'   sample_col = "sample",
#'   group_col = "tissue"
#' )
#'
#' head(res$scores)
#' res$swaps
#' print(res$plot)
#'}
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs
#'
#' @export
marker_score_check <- function(
    counts,
    metadata,
    marker_sets,
    sample_col = "sample",
    group_col = "tissue",
    feature_col = "gene_name",
    log_transform = TRUE,
    plot = TRUE,
    style = "classic",
    txtsize = 12
) {

  counts <- as.data.frame(counts, check.names = FALSE)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  sample_cols <- metadata[[sample_col]]
  sample_cols <- sample_cols[sample_cols %in% colnames(counts)]

  if (length(sample_cols) == 0) {
    stop("No metadata samples found in counts columns.", call. = FALSE)
  }

  gene_names <- toupper(counts[[feature_col]])

  compute_one_score <- function(genes) {

    keep <- gene_names %in% toupper(genes)

    if (!any(keep)) {
      return(rep(NA_real_, length(sample_cols)))
    }

    mat <- as.matrix(counts[keep, sample_cols, drop = FALSE])
    storage.mode(mat) <- "numeric"

    if (log_transform) {
      mat <- log2(mat + 1)
    }

    colMeans(mat, na.rm = TRUE)
  }

  scores <- data.frame(
    sample = sample_cols,
    stringsAsFactors = FALSE
  )

  for (nm in names(marker_sets)) {
    scores[[paste0(nm, "_score")]] <- compute_one_score(marker_sets[[nm]])
  }

  md <- metadata[, c(sample_col, group_col), drop = FALSE]
  colnames(md) <- c("sample", "true_group")

  scores <- merge(scores, md, by = "sample", all.x = TRUE)

  score_cols <- paste0(names(marker_sets), "_score")
  score_mat <- as.matrix(scores[, score_cols, drop = FALSE])

  scores$predicted_group <- names(marker_sets)[
    max.col(score_mat, ties.method = "first")
  ]

  scores$possible_swap <- (
    scores$true_group != scores$predicted_group
  )

  scores_long <- reshape2::melt(
    scores,
    id.vars = c(
      "sample",
      "true_group",
      "predicted_group",
      "possible_swap"
    ),
    measure.vars = score_cols,
    variable.name = "marker_set",
    value.name = "score"
  )

  scores_long$marker_set <- sub(
    "_score$",
    "",
    scores_long$marker_set
  )

  p <- NULL

  if (plot) {

    p <- ggplot2::ggplot(
      scores_long,
      ggplot2::aes(
        x = sample,
        y = score,
        fill = marker_set
      )
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = NULL,
        y = "Mean marker score",
        fill = "Marker set"
      ) +
      .deggo_theme(style = style, txtsize = txtsize)
  }

  list(
    scores = scores,
    long = scores_long,
    swaps = scores[scores$possible_swap, , drop = FALSE],
    plot = p
  )
}




.deggo_initialize <- function(
    output_dir,
    method,
    analysis_mode,
    organism,
    ontology,
    report_template,
    log
) {

  if (is.null(output_dir))
    output_dir <- "DEGgo_results"

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

  dir.create(output_dir,
             recursive = TRUE,
             showWarnings = FALSE)

  output_dir <- normalizePath(
    output_dir,
    winslash="/",
    mustWork=FALSE
  )

  dirs <- .deggo_dirs(
    output_dir,
    analysis_mode
  )

  repro_dir <- file.path(
    output_dir,
    "reproducibility"
  )

  dir.create(
    repro_dir,
    recursive=TRUE,
    showWarnings=FALSE
  )

  deggo_version <- tryCatch(

    as.character(
      utils::packageVersion("DEGgo")
    ),

    error=function(e)
      "development"

  )

  if (is.null(report_template)) {

    report_template <- .deggo_existing_file(

      system.file(
        "rmarkdown/templates/deggo_report/skeleton/skeleton.Rmd",
        package="DEGgo"
      ),

      file.path(
        getwd(),
        "inst/rmarkdown/templates/deggo_report/skeleton/skeleton.Rmd"
      )

    )

  }

  list(

    output_dir=output_dir,

    dirs=dirs,

    repro_dir=repro_dir,

    deggo_version=deggo_version,

    report_template=report_template

  )

}




.deggo_raw_qc <- function(

  counts,
  metadata,

  raw_qc,
  remove_flagged,

  gene_col,
  feature_col,
  sample_col,

  qc_markers,

  qc_sample_col,

  qc_output_prefix,

  output_dir

){

  qc_raw <- NULL
  qc_clean <- NULL

  if(!raw_qc){

    return(list(

      counts=counts,

      metadata=metadata,

      qc_raw=NULL,

      qc_clean=NULL

    ))

  }

  qc_raw <- explore_bulk_rnaseq(

    counts=counts,

    metadata=metadata,

    gene_col=gene_col,

    feature_col=feature_col,

    sample_col=sample_col,

    markers=qc_markers,

    output_dir=file.path(
      output_dir,
      paste0(qc_output_prefix,"_raw")
    )

  )

  if(remove_flagged){

    cleaned <- remove_flagged_samples(

      counts=.deggo_counts_with_gene_col(
        counts,
        gene_col
      ),

      metadata=metadata,

      qc_table=qc_raw$qc,

      sample_col=qc_sample_col %||% sample_col[1]

    )

    counts <- cleaned$counts
    metadata <- cleaned$metadata

    qc_clean <- explore_bulk_rnaseq(

      counts,

      metadata,

      gene_col,

      feature_col,

      sample_col,

      markers=qc_markers,

      output_dir=file.path(
        output_dir,
        paste0(qc_output_prefix,"_clean")
      )

    )

  }

  list(

    counts=counts,

    metadata=metadata,

    qc_raw=qc_raw,

    qc_clean=qc_clean

  )

}




.deggo_prepare_inputs <- function(

  counts,
  metadata,

  prepare_input,

  gene_col,
  feature_col,
  sample_col,

  filter_method,
  min_count,
  min_samples,
  min_total,

  repro_dir,
  save_clean_inputs,

  txtsize,
  dirs

){

  if(prepare_input){
    prep <- prepare_counts_metadata(
      counts=.deggo_counts_with_gene_col(
        counts,
        gene_col
      ),
      metadata=metadata,
      gene_col=gene_col,
      feature_col=feature_col,
      sample_col=sample_col
    )
    counts <- prep$counts
    metadata <- prep$metadata
  }

  counts <- clean_ensembl_ids(counts)

  counts <- preprocess_counts(
    counts,
    metadata,
    filter_method,
    min_count,
    min_samples,
    min_total
  )

  sample_qc <- run_sample_qc(
    counts,
    metadata,
    txtsize=txtsize,
    output_dir=dirs$qc
  )

  .save_clean_input_files(
    counts,
    metadata,
    repro_dir,
    save_clean_inputs
  )

  list(
    counts=counts,
    metadata=metadata,
    sample_qc=sample_qc
  )
}




.deggo_run_analysis <- function(
  analysis_mode,
  counts,
  metadata,
  ...
){

  if(analysis_mode=="single"){
    return(
      .run_deggo_single(
        counts=counts,
        metadata=metadata,
        ...
      )
    )
  }

  .run_deggo_pairwise(
    counts=counts,
    metadata=metadata,
    ...
  )
}




