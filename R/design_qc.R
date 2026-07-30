
# ----------------------------------------------------------------------- #
# Normalize count input for design_qc()
# ----------------------------------------------------------------------- #

#' Prepare count and metadata inputs for design quality control
#'
#' Internal helper used to validate, clean, reorder, and standardize a
#' gene-by-sample count table together with its associated sample metadata.
#'
#' The function detects or extracts gene identifiers, optionally removes
#' Ensembl version suffixes, validates sample correspondence, converts the
#' expression values to a numeric matrix, and optionally aggregates duplicated
#' gene identifiers by summing their counts.
#'
#' @param counts A matrix, data frame, or matrix-like object containing gene
#'   expression counts. Rows should represent genes and columns should represent
#'   samples. Gene identifiers may be stored in the row names or in a dedicated
#'   column specified by `gene_col`.
#'
#' @param metadata A data frame or data-frame-like object containing sample
#'   metadata. It must contain the sample identifier column specified by
#'   `sample_col`.
#'
#' @param sample_col A character string giving the name of the metadata column
#'   containing sample identifiers. These identifiers must match column names
#'   in `counts`.
#'
#' @param gene_col An optional character string giving the name of the column
#'   in `counts` that contains gene identifiers. When `NULL`, the function first
#'   checks whether valid gene identifiers are stored in the row names. If not,
#'   it attempts to detect a common gene identifier column automatically.
#'
#' @param strip_ensembl_version Logical. If `TRUE`, trailing numeric Ensembl
#'   version suffixes are removed from gene identifiers. For example,
#'   `"ENSG00000141510.18"` becomes `"ENSG00000141510"`.
#'
#' @param aggregate_duplicates Logical. If `TRUE`, duplicated gene identifiers
#'   are aggregated by summing counts across duplicated rows. If `FALSE`, the
#'   function stops when duplicated gene identifiers are detected.
#'
#' @return A list with the following elements:
#'
#' \describe{
#'   \item{counts}{
#'     An integer matrix containing genes in rows and samples in columns.
#'     Sample columns are ordered according to `metadata[[sample_col]]`.
#'   }
#'   \item{metadata}{
#'     The validated metadata data frame.
#'   }
#'   \item{gene_col}{
#'     The gene identifier column used, or `NULL` when gene identifiers were
#'     obtained from row names.
#'   }
#'   \item{duplicated_genes_aggregated}{
#'     The number of duplicated gene rows detected before aggregation.
#'   }
#' }
#'
#' @details
#' Sample identifiers must be unique, non-missing, and present as columns in the
#' count table. Count values must be numeric, finite, and non-negative.
#'
#' A warning is emitted when non-integer values are detected because downstream
#' methods such as DESeq2 variance-stabilizing transformation and regularized
#' logarithm transformation generally expect raw integer counts.
#'
#' @keywords internal
#'
#' @noRd
.dq_prepare_counts <- function(
    counts,
    metadata,
    sample_col,
    gene_col = NULL,
    strip_ensembl_version = TRUE,
    aggregate_duplicates = TRUE
) {
  counts <- as.data.frame(
    counts,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  metadata <- as.data.frame(
    metadata,
    stringsAsFactors = FALSE
  )

  if (!sample_col %in% colnames(metadata)) {
    stop(
      "Sample column '",
      sample_col,
      "' was not found in metadata."
    )
  }

  sample_ids <- trimws(
    as.character(metadata[[sample_col]])
  )

  if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
    stop(
      "Metadata sample identifiers contain missing or empty values."
    )
  }

  if (anyDuplicated(sample_ids)) {
    duplicated_samples <- unique(
      sample_ids[duplicated(sample_ids)]
    )

    stop(
      "Metadata contains duplicated sample identifiers: ",
      paste(duplicated_samples, collapse = ", ")
    )
  }

  # ----------------------------------------------------------------------- #
  # Detect whether gene IDs are already stored in row names
  # ----------------------------------------------------------------------- #

  has_valid_rownames <- (
    !is.null(rownames(counts)) &&
      !identical(
        rownames(counts),
        as.character(seq_len(nrow(counts)))
      ) &&
      all(!is.na(rownames(counts))) &&
      all(nzchar(rownames(counts)))
  )

  # ----------------------------------------------------------------------- #
  # Detect gene identifier column
  # ----------------------------------------------------------------------- #

  if (is.null(gene_col) && !has_valid_rownames) {
    gene_candidates <- c(
      "gene_id",
      "EnsemblID",
      "ENSEMBL",
      "ensembl_id",
      "GeneID",
      "gene",
      "Gene",
      "feature_id",
      "FeatureID"
    )

    detected <- gene_candidates[
      gene_candidates %in% colnames(counts)
    ]

    if (length(detected) > 0L) {
      gene_col <- detected[1]
    }
  }

  if (!is.null(gene_col)) {
    if (!gene_col %in% colnames(counts)) {
      stop(
        "Gene identifier column '",
        gene_col,
        "' was not found in the count table. Available columns: ",
        paste(colnames(counts), collapse = ", ")
      )
    }

    gene_ids <- trimws(
      as.character(counts[[gene_col]])
    )
  } else if (has_valid_rownames) {
    gene_ids <- trimws(
      as.character(rownames(counts))
    )
  } else {
    stop(
      "Counts must have unique gene row names or contain a gene ",
      "identifier column. Supply gene_col explicitly, for example ",
      "gene_col = 'gene_id'."
    )
  }

  if (isTRUE(strip_ensembl_version)) {
    gene_ids <- sub(
      "\\.[0-9]+$",
      "",
      gene_ids
    )
  }

  keep_gene <- (
    !is.na(gene_ids) &
      nzchar(gene_ids)
  )

  if (!all(keep_gene)) {
    counts <- counts[
      keep_gene,
      ,
      drop = FALSE
    ]

    gene_ids <- gene_ids[keep_gene]
  }

  if (length(gene_ids) == 0L) {
    stop(
      "No valid gene identifiers remain after cleaning."
    )
  }

  # ----------------------------------------------------------------------- #
  # Extract sample columns in metadata order
  # ----------------------------------------------------------------------- #

  missing_samples <- setdiff(
    sample_ids,
    colnames(counts)
  )

  if (length(missing_samples) > 0L) {
    stop(
      "Samples missing from the count table: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  count_matrix <- as.matrix(
    counts[
      ,
      sample_ids,
      drop = FALSE
    ]
  )

  suppressWarnings(
    storage.mode(count_matrix) <- "numeric"
  )

  if (anyNA(count_matrix)) {
    stop(
      "The count matrix contains missing or non-numeric values."
    )
  }

  if (any(!is.finite(count_matrix))) {
    stop(
      "The count matrix contains non-finite values."
    )
  }

  if (any(count_matrix < 0)) {
    stop(
      "The count matrix contains negative values."
    )
  }

  # Raw counts should normally be integer-like
  if (any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
    warning(
      "The count matrix contains non-integer values. ",
      "VST and rlog generally expect raw integer counts.",
      call. = FALSE
    )
  }

  # ----------------------------------------------------------------------- #
  # Aggregate duplicated gene IDs
  # ----------------------------------------------------------------------- #

  n_duplicates <- sum(duplicated(gene_ids))

  if (n_duplicates > 0L) {
    if (!isTRUE(aggregate_duplicates)) {
      duplicated_ids <- unique(
        gene_ids[duplicated(gene_ids)]
      )

      stop(
        "Counts contain duplicated gene identifiers: ",
        paste(head(duplicated_ids, 20L), collapse = ", "),
        if (length(duplicated_ids) > 20L) " ..." else ""
      )
    }

    count_matrix <- rowsum(
      count_matrix,
      group = gene_ids,
      reorder = FALSE
    )
  } else {
    rownames(count_matrix) <- gene_ids
  }

  storage.mode(count_matrix) <- "integer"

  if (anyDuplicated(rownames(count_matrix))) {
    stop(
      "Gene identifiers remain duplicated after count preparation."
    )
  }

  if (anyDuplicated(colnames(count_matrix))) {
    stop(
      "The count matrix contains duplicated sample columns."
    )
  }

  list(
    counts = count_matrix,
    metadata = metadata,
    gene_col = gene_col,
    duplicated_genes_aggregated = n_duplicates
  )
}



#' Experimental design quality control for bulk RNA-seq
#'
#' Evaluates metadata structure, repeated experimental units, PCA-associated
#' effects, optional PERMANOVA, and recommends a differential-expression model.
#'
#' @param counts Integer-like matrix/data.frame with genes in rows and samples in columns.
#' @param metadata Data frame with one row per sample.
#' @param sample_col Metadata column containing sample identifiers.
#' @param group_col Main biological condition to test.
#' @param gene_col Optional character string specifying the column in
#'   `counts` that contains gene identifiers. When `NULL`, gene identifiers
#'   are expected to be stored in the row names.
#'
#' @param strip_ensembl_version Logical. If `TRUE`, remove Ensembl version
#'   suffixes from gene identifiers, for example converting
#'   `"ENSG00000141510.18"` to `"ENSG00000141510"`.
#'
#' @param aggregate_duplicates Logical. If `TRUE`, rows with duplicated gene
#'   identifiers after preprocessing are aggregated by summing their counts.
#'   If `FALSE`, duplicated identifiers should trigger an error or be retained
#'   according to the function implementation.
#' @param unit_col Independent experimental unit (e.g. cage, donor, patient).
#'   If NULL, a conservative name-based detector is used.
#' @param batch_cols Optional technical/biological covariates to assess.
#' @param transform PCA transformation: "vst", "rlog", or "logcpm".
#' @param n_top Number of most variable genes used for PCA.
#' @param n_pcs Number of PCs used for structure assessment.
#' @param run_permanova Run vegan::adonis2 when vegan is installed.
#' @param permutations Number of PERMANOVA permutations.
#' @param dominance_ratio Minimum weighted PCA R2 ratio for declaring one
#'   factor stronger than the main group.
#' @param min_effect_r2 Minimum weighted PCA R2 for an effect to be considered relevant.
#' @param verbose Print a concise summary.
#'
#' @return An object of class "deggo_design_qc".
#' @export
#'
#' @examples
#' \dontrun{
#' qc <- design_qc(
#'   counts = counts,
#'   metadata = metadata,
#'   sample_col = "sample",
#'   group_col = "Group",
#'   unit_col = "Cage"
#' )
#' print(qc)
#' plot(qc)
#' }
design_qc <- function(
    counts,
    metadata,
    sample_col,
    group_col,
    unit_col = NULL,
    batch_cols = NULL,
    gene_col = NULL,
    strip_ensembl_version = TRUE,
    aggregate_duplicates = TRUE,
    transform = c("vst", "rlog", "logcpm"),
    n_top = 500L,
    n_pcs = 5L,
    run_permanova = TRUE,
    permutations = 999L,
    dominance_ratio = 2,
    min_effect_r2 = 0.10,
    verbose = TRUE
) {
  transform <- match.arg(transform)

  # ----------------------------------------------------------------------- #
  # Prepare counts
  # ----------------------------------------------------------------------- #

  prepared <- .dq_prepare_counts(
    counts = counts,
    metadata = metadata,
    sample_col = sample_col,
    gene_col = gene_col,
    strip_ensembl_version = strip_ensembl_version,
    aggregate_duplicates = aggregate_duplicates
  )

  counts <- prepared$counts
  metadata <- prepared$metadata

  # ----------------------------------------------------------------------- #
  # Validate aligned inputs
  # ----------------------------------------------------------------------- #

  x <- .dq_validate_inputs(
    counts = counts,
    metadata = metadata,
    sample_col = sample_col,
    group_col = group_col,
    unit_col = unit_col,
    batch_cols = batch_cols
  )

  counts <- x$counts
  metadata <- x$metadata
  unit_col <- x$unit_col
  batch_cols <- x$batch_cols

  # ----------------------------------------------------------------------- #
  # Detect repeated structure
  # ----------------------------------------------------------------------- #

  repeated <- .dq_repeated_structure(
    metadata = metadata,
    group_col = group_col,
    unit_col = unit_col
  )

  # ----------------------------------------------------------------------- #
  # Transform and run PCA
  # ----------------------------------------------------------------------- #

  transformed <- .dq_transform_counts(
    counts = counts,
    metadata = metadata,
    transform = transform
  )

  pca <- .dq_run_pca(
    transformed = transformed,
    metadata = metadata,
    n_top = n_top,
    n_pcs = n_pcs
  )

  # ----------------------------------------------------------------------- #
  # Factors assessed by PCA and PERMANOVA
  # ----------------------------------------------------------------------- #

  factors <- unique(
    c(
      group_col,
      unit_col,
      batch_cols
    )
  )

  factors <- factors[
    !is.na(factors) &
      nzchar(factors)
  ]

  pca_effects <- .dq_pca_effects(
    scores = pca$scores,
    metadata = metadata,
    factors = factors,
    variance = pca$variance
  )

  permanova <- .dq_permanova(
    transformed = transformed[
      pca$genes,
      ,
      drop = FALSE
    ],
    metadata = metadata,
    factors = factors,
    permutations = permutations,
    enabled = run_permanova
  )

  # ----------------------------------------------------------------------- #
  # Determine dominant factor
  # ----------------------------------------------------------------------- #

  dominance <- .dq_detect_dominance(
    pca_effects = pca_effects,
    group_col = group_col,
    unit_col = unit_col,
    batch_cols = batch_cols,
    dominance_ratio = dominance_ratio,
    min_effect_r2 = min_effect_r2
  )

  # ----------------------------------------------------------------------- #
  # Recommend model
  # ----------------------------------------------------------------------- #

  recommendation <- .dq_recommend_model(
    metadata = metadata,
    group_col = group_col,
    unit_col = unit_col,
    repeated = repeated,
    dominance = dominance
  )

  preparation_warnings <- character()

  if (prepared$duplicated_genes_aggregated > 0L) {
    preparation_warnings <- paste0(
      prepared$duplicated_genes_aggregated,
      " duplicated gene rows were aggregated by summing their counts."
    )
  }

  # ----------------------------------------------------------------------- #
  # Build result
  # ----------------------------------------------------------------------- #

  out <- list(
    call = match.call(),

    parameters = list(
      sample_col = sample_col,
      group_col = group_col,
      unit_col = unit_col,
      batch_cols = batch_cols,
      gene_col = prepared$gene_col,
      strip_ensembl_version = strip_ensembl_version,
      aggregate_duplicates = aggregate_duplicates,
      transform = transform,
      n_top = n_top,
      n_pcs = n_pcs,
      dominance_ratio = dominance_ratio,
      min_effect_r2 = min_effect_r2
    ),

    input = list(
      genes = nrow(counts),
      samples = ncol(counts),
      gene_col = prepared$gene_col,
      duplicated_genes_aggregated =
        prepared$duplicated_genes_aggregated
    ),

    summary = data.frame(
      samples = ncol(counts),
      groups = nlevels(metadata[[group_col]]),

      independent_units = if (is.null(unit_col)) {
        NA_integer_
      } else {
        nlevels(metadata[[unit_col]])
      },

      repeated_units = repeated$has_repeated_units,
      crossed_group_unit = repeated$unit_crosses_group,
      dominant_factor = dominance$dominant_factor,
      recommended_method = recommendation$method,
      stringsAsFactors = FALSE
    ),

    counts = counts,
    metadata = metadata,
    structure = repeated,
    transformed = transformed,
    pca = pca,
    pca_effects = pca_effects,
    permanova = permanova,
    dominance = dominance,
    recommendation = recommendation,

    warnings = unique(
      c(
        preparation_warnings,
        x$warnings,
        repeated$warnings,
        recommendation$warnings
      )
    )
  )

  class(out) <- "deggo_design_qc"

  if (isTRUE(verbose)) {
    print(out)
  }

  out
}


#' Print a DEGgo design quality control summary
#'
#' Displays a concise summary of the study design assessment returned by
#' [design_qc()]. The summary includes the number of samples and groups,
#' whether repeated measurements are present, the dominant experimental factor,
#' the recommended differential expression framework, and any detected design
#' warnings.
#'
#' @param x An object of class `"deggo_design_qc"` returned by
#'   [design_qc()].
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return The input object, invisibly.
#'
#' @seealso [design_qc()]
#'
#' @examples
#' \dontrun{
#' res <- design_qc(
#'   counts = counts,
#'   metadata = metadata,
#'   sample_col = "sample",
#'   group_col = "group"
#' )
#'
#' print(res)
#' }
#'
#' @method print deggo_design_qc
#' @export
print.deggo_design_qc <- function(x, ...) {
  cat("\nDEGgo DesignQC\n")
  cat(strrep("-", 50), "\n", sep = "")
  cat("Samples:             ", x$summary$samples, "\n", sep = "")
  cat("Groups:              ", x$summary$groups, "\n", sep = "")

  if (!is.na(x$summary$independent_units)) {
    cat("Independent units:   ", x$summary$independent_units,
        " (", x$parameters$unit_col, ")\n", sep = "")
  } else {
    cat("Independent units:   not specified/detected\n")
  }

  cat("Repeated units:      ", ifelse(x$summary$repeated_units, "yes", "no"), "\n", sep = "")
  cat("Dominant factor:     ", x$summary$dominant_factor, "\n", sep = "")
  cat("Recommended method:  ", x$summary$recommended_method, "\n", sep = "")
  cat("Reason:              ", x$recommendation$reason, "\n", sep = "")

  if (length(x$warnings)) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat("  - ", w, "\n", sep = "")
  }

  invisible(x)
}


#' Plot DesignQC PCA
#'
#' @param x A deggo_design_qc object.
#' @param color_by Metadata variable used for colour.
#' @param shape_by Optional metadata variable used for shape.
#' @param label Show sample labels.
#' @param ... Unused.
#' @export
plot.deggo_design_qc <- function(x, color_by = NULL, shape_by = NULL,
                                 label = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  color_by <- color_by %||% x$parameters$group_col
  df <- x$pca$scores

  if (!color_by %in% names(df)) {
    stop("Unknown color_by column: ", color_by, call. = FALSE)
  }
  if (!is.null(shape_by) && !shape_by %in% names(df)) {
    stop("Unknown shape_by column: ", shape_by, call. = FALSE)
  }

  aes_args <- list(x = quote(PC1), y = quote(PC2), colour = as.name(color_by))
  if (!is.null(shape_by)) aes_args$shape <- as.name(shape_by)

  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", 100 * x$pca$variance[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * x$pca$variance[2]),
      colour = color_by,
      shape = shape_by
    ) +
    ggplot2::theme_classic(base_size = 9)

  if (isTRUE(label)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data[[x$parameters$sample_col]]),
      check_overlap = TRUE,
      vjust = -0.7,
      show.legend = FALSE,
      size = 2.5
    )
  }

  p
}


#' Aggregate raw counts by an experimental unit
#'
#' @param counts Gene-by-sample count matrix.
#' @param metadata Sample metadata.
#' @param sample_col Sample identifier column.
#' @param unit_col Experimental unit column, e.g. Cage.
#' @param group_col Group column.
#' @param extra_cols Optional columns required to be constant inside each unit.
#' @return List containing aggregated counts and metadata.
#' @export
aggregate_counts_by_unit <- function(
    counts,
    metadata,
    sample_col,
    unit_col,
    group_col,
    extra_cols = NULL
) {
  x <- .dq_validate_inputs(
    counts, metadata, sample_col, group_col,
    unit_col = unit_col,
    batch_cols = extra_cols
  )
  counts <- x$counts
  metadata <- x$metadata

  cols <- unique(c(unit_col, group_col, extra_cols))
  for (cl in setdiff(cols, unit_col)) {
    n_by_unit <- tapply(metadata[[cl]], metadata[[unit_col]], function(z) {
      length(unique(as.character(z[!is.na(z)])))
    })
    if (any(n_by_unit > 1L)) {
      bad <- names(n_by_unit)[n_by_unit > 1L]
      stop(
        "Column '", cl, "' is not constant within ", unit_col,
        ". Conflicting units: ", paste(bad, collapse = ", "),
        call. = FALSE
      )
    }
  }

  units <- levels(droplevels(metadata[[unit_col]]))
  agg <- vapply(units, function(u) {
    rowSums(counts[, metadata[[unit_col]] == u, drop = FALSE])
  }, FUN.VALUE = numeric(nrow(counts)))

  rownames(agg) <- rownames(counts)
  colnames(agg) <- units

  first <- match(units, metadata[[unit_col]])
  agg_meta <- metadata[first, cols, drop = FALSE]
  rownames(agg_meta) <- units
  agg_meta[[sample_col]] <- units
  agg_meta <- agg_meta[, unique(c(sample_col, cols)), drop = FALSE]

  list(counts = agg, metadata = agg_meta)
}


# ====================================================================== #
# Internal helpers
# ====================================================================== #

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Validate DesignQC inputs
#'
#' Internal helper that validates count and metadata inputs, checks sample
#' correspondence, standardizes grouping variables, and resolves experimental
#' unit and batch columns.
#'
#' @param counts A gene-by-sample count matrix or data frame.
#' @param metadata A data frame containing sample-level metadata.
#' @param sample_col Character string naming the sample identifier column.
#' @param group_col Character string naming the biological group column.
#' @param unit_col Optional character string naming the experimental unit column.
#' @param batch_cols Optional character vector naming batch-related columns.
#'
#' @return A list containing validated counts, reordered metadata, the resolved
#'   experimental unit column, batch columns, and validation warnings.
#'
#' @keywords internal
#' @noRd
.dq_validate_inputs <- function(counts, metadata, sample_col, group_col,
                                unit_col = NULL, batch_cols = NULL) {
  warnings <- character()

  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("'counts' must be a matrix or data.frame.", call. = FALSE)
  }
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"

  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts))) {
    stop("Counts must have unique gene row names.", call. = FALSE)
  }
  if (is.null(colnames(counts)) || anyDuplicated(colnames(counts))) {
    stop("Counts must have unique sample column names.", call. = FALSE)
  }
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("Counts must contain finite non-negative values.", call. = FALSE)
  }

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  required <- c(sample_col, group_col)
  missing_required <- setdiff(required, names(metadata))
  if (length(missing_required)) {
    stop("Missing metadata columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }

  if (anyNA(metadata[[sample_col]]) || anyDuplicated(metadata[[sample_col]])) {
    stop("Metadata sample identifiers must be non-missing and unique.", call. = FALSE)
  }

  ids <- as.character(metadata[[sample_col]])
  missing_meta <- setdiff(colnames(counts), ids)
  missing_counts <- setdiff(ids, colnames(counts))
  if (length(missing_meta) || length(missing_counts)) {
    stop(
      "Counts/metadata sample mismatch. Missing in metadata: ",
      paste(missing_meta, collapse = ", "),
      "; missing in counts: ", paste(missing_counts, collapse = ", "),
      call. = FALSE
    )
  }

  metadata <- metadata[match(colnames(counts), ids), , drop = FALSE]
  rownames(metadata) <- colnames(counts)

  metadata[[group_col]] <- droplevels(factor(metadata[[group_col]]))
  if (nlevels(metadata[[group_col]]) < 2L) {
    stop("'group_col' must contain at least two groups.", call. = FALSE)
  }

  if (is.null(unit_col)) {
    detected <- .dq_detect_unit_column(metadata, sample_col, group_col)
    unit_col <- detected$unit_col
    warnings <- c(warnings, detected$warning)
  }

  if (!is.null(unit_col)) {
    if (!unit_col %in% names(metadata)) {
      stop("Unknown unit_col: ", unit_col, call. = FALSE)
    }
    metadata[[unit_col]] <- droplevels(factor(metadata[[unit_col]]))
  }

  batch_cols <- unique(batch_cols %||% character())
  bad_batch <- setdiff(batch_cols, names(metadata))
  if (length(bad_batch)) {
    stop("Unknown batch_cols: ", paste(bad_batch, collapse = ", "), call. = FALSE)
  }
  for (cl in batch_cols) metadata[[cl]] <- droplevels(factor(metadata[[cl]]))

  list(
    counts = counts,
    metadata = metadata,
    unit_col = unit_col,
    batch_cols = batch_cols,
    warnings = warnings[nzchar(warnings)]
  )
}



#' Detect an experimental unit column
#'
#' Internal helper that conservatively searches metadata for a repeated
#' experimental unit such as a donor, patient, subject, cage, or animal.
#'
#' @param metadata A sample metadata data frame.
#' @param sample_col Character string naming the sample identifier column.
#' @param group_col Character string naming the biological group column.
#'
#' @return A list with the detected `unit_col`, or `NULL` when no suitable
#'   repeated-unit column is found, together with an explanatory warning.
#'
#' @keywords internal
#' @noRd
.dq_detect_unit_column <- function(metadata, sample_col, group_col) {
  aliases <- c(
    "cage", "donor", "patient", "subject", "individual", "animal",
    "mouse", "rat", "litter", "participant", "person", "id"
  )
  nms <- names(metadata)
  canonical <- tolower(gsub("[^a-z0-9]", "", nms))
  candidate <- which(canonical %in% aliases & !nms %in% c(sample_col, group_col))

  if (!length(candidate)) {
    return(list(
      unit_col = NULL,
      warning = paste0(
        "No experimental-unit column was supplied or conservatively detected. ",
        "Set unit_col explicitly for cage/donor/patient designs."
      )
    ))
  }

  repeated <- candidate[vapply(candidate, function(i) {
    anyDuplicated(metadata[[i]]) > 0L && length(unique(metadata[[i]])) > 1L
  }, logical(1))]

  if (!length(repeated)) {
    return(list(unit_col = NULL, warning = "Candidate unit columns were unique and therefore not treated as repeated units."))
  }

  chosen <- nms[repeated[1]]
  list(
    unit_col = chosen,
    warning = paste0(
      "Experimental unit auto-detected as '", chosen,
      "'. Verify this choice; explicit unit_col is safer."
    )
  )
}




#' Characterize repeated-measures structure
#'
#' Internal helper that summarizes the number of samples and biological groups
#' represented by each experimental unit and determines whether units are
#' repeated or cross biological groups.
#'
#' @param metadata A sample metadata data frame.
#' @param group_col Character string naming the biological group column.
#' @param unit_col Optional character string naming the experimental unit column.
#'
#' @return A list describing repeated units, cross-group units, samples per
#'   unit, groups per unit, and any design warnings.
#'
#' @keywords internal
#' @noRd
.dq_repeated_structure <- function(metadata, group_col, unit_col) {
  if (is.null(unit_col)) {
    return(list(
      has_repeated_units = FALSE,
      unit_crosses_group = FALSE,
      samples_per_unit = NULL,
      groups_per_unit = NULL,
      warnings = character()
    ))
  }

  samples_per_unit <- table(metadata[[unit_col]])
  groups_per_unit <- tapply(metadata[[group_col]], metadata[[unit_col]], function(z) {
    length(unique(as.character(z)))
  })

  crossed <- any(groups_per_unit > 1L)
  warnings <- character()
  if (crossed) {
    warnings <- c(
      warnings,
      paste0(
        "Some ", unit_col, " levels occur in multiple ", group_col,
        " levels. This is a paired/repeated-measures structure, not a simple nested design."
      )
    )
  }

  list(
    has_repeated_units = any(samples_per_unit > 1L),
    unit_crosses_group = crossed,
    samples_per_unit = samples_per_unit,
    groups_per_unit = groups_per_unit,
    warnings = warnings
  )
}



#' Transform counts for DesignQC
#'
#' Internal helper that transforms raw counts using DESeq2 variance-stabilizing
#' or regularized-logarithm transformation when available, with log2-CPM as a
#' fallback.
#'
#' @param counts A numeric gene-by-sample count matrix.
#' @param metadata A sample metadata data frame ordered to match `counts`.
#' @param transform Character string specifying `"vst"`, `"rlog"`, or
#'   `"logcpm"`.
#'
#' @return A numeric transformed expression matrix with genes in rows and
#'   samples in columns.
#'
#' @keywords internal
#' @noRd
.dq_transform_counts <- function(counts, metadata, transform) {
  rounded <- round(counts)
  de_available <- requireNamespace("DESeq2", quietly = TRUE)
  transform_failed <- FALSE

  if (transform %in% c("vst", "rlog") && de_available) {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = rounded,
      colData = metadata,
      design = ~1
    )
    tx <- tryCatch(
      {
        if (transform == "vst") {
          # DESeq2::vst() estimates its trend from `nsub = 1000` genes and
          # therefore errors for small matrices commonly used in unit tests.
          # The full transformation is appropriate and stable in that case.
          if (nrow(dds) < 1000L) {
            DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
          } else {
            DESeq2::vst(dds, blind = TRUE)
          }
        } else {
          DESeq2::rlog(dds, blind = TRUE)
        }
      },
      error = function(e) {
        transform_failed <<- TRUE
        warning(
          paste0(
            "DESeq2 ", transform, " transformation failed (", conditionMessage(e),
            "); falling back to log2-CPM for DesignQC."
          ),
          call. = FALSE
        )
        NULL
      }
    )

    if (!is.null(tx)) {
      return(SummarizedExperiment::assay(tx))
    }
  }

  if (transform != "logcpm" && !de_available && !transform_failed) {
    warning(
      "DESeq2 is unavailable; falling back to log2-CPM for DesignQC.",
      call. = FALSE
    )
  }

  lib <- colSums(counts)
  if (any(lib <= 0)) stop("All samples must have a positive library size.", call. = FALSE)
  log2(t(t(counts + 0.5) / (lib + 1) * 1e6) + 1)
}


#' Run principal component analysis for DesignQC
#'
#' Internal helper that selects the most variable genes and performs PCA on
#' transformed expression values.
#'
#' @param transformed A transformed gene-by-sample expression matrix.
#' @param metadata A sample metadata data frame ordered to match the expression
#'   matrix.
#' @param n_top Integer specifying the maximum number of variable genes to use.
#' @param n_pcs Integer specifying the maximum number of principal components
#'   to retain.
#'
#' @return A list containing the PCA fit, sample scores with metadata, explained
#'   variance, and selected genes.
#'
#' @keywords internal
#' @noRd
.dq_run_pca <- function(transformed, metadata, n_top, n_pcs) {
  vars <- apply(transformed, 1L, stats::var, na.rm = TRUE)
  vars[!is.finite(vars)] <- 0
  n_top <- min(as.integer(n_top), nrow(transformed))
  genes <- names(sort(vars, decreasing = TRUE))[seq_len(n_top)]

  fit <- stats::prcomp(t(transformed[genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
  variance <- fit$sdev^2 / sum(fit$sdev^2)
  n_pcs <- min(as.integer(n_pcs), ncol(fit$x))

  scores <- as.data.frame(fit$x[, seq_len(n_pcs), drop = FALSE])
  scores <- cbind(scores, metadata[rownames(scores), , drop = FALSE])

  list(
    fit = fit,
    scores = scores,
    variance = variance[seq_len(n_pcs)],
    genes = genes
  )
}




#' Estimate metadata effects on principal components
#'
#' Internal helper that fits linear models for metadata factors across retained
#' principal components and summarizes their adjusted R-squared values.
#'
#' @param scores A data frame containing PCA scores.
#' @param metadata A sample metadata data frame.
#' @param factors Character vector naming metadata factors to evaluate.
#' @param variance Numeric vector containing the explained variance of the
#'   retained principal components.
#'
#' @return A data frame containing weighted and maximum adjusted R-squared
#'   values for each evaluated factor.
#'
#' @keywords internal
#' @noRd
.dq_pca_effects <- function(scores, metadata, factors, variance) {
  if (!length(factors)) return(data.frame())

  rows <- lapply(factors, function(factor_name) {
    z <- metadata[[factor_name]]
    if (length(unique(z[!is.na(z)])) < 2L) return(NULL)

    pc_names <- paste0("PC", seq_along(variance))
    r2 <- vapply(pc_names, function(pc) {
      df <- data.frame(y = scores[[pc]], x = z)
      df <- df[stats::complete.cases(df), , drop = FALSE]
      if (nrow(df) < 3L || length(unique(df$x)) < 2L) return(NA_real_)
      summary(stats::lm(y ~ x, data = df))$adj.r.squared
    }, numeric(1))
    r2 <- pmax(r2, 0, na.rm = FALSE)

    data.frame(
      factor = factor_name,
      weighted_r2 = sum(r2 * variance, na.rm = TRUE) / sum(variance[is.finite(r2)]),
      max_pc_r2 = max(r2, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out[order(out$weighted_r2, decreasing = TRUE), , drop = FALSE]
}



#' Run PERMANOVA for DesignQC
#'
#' Internal helper that evaluates associations between sample-level expression
#' distances and metadata factors using marginal PERMANOVA.
#'
#' @param transformed A transformed gene-by-sample expression matrix.
#' @param metadata A sample metadata data frame.
#' @param factors Character vector naming metadata factors to test.
#' @param permutations Integer specifying the number of permutations.
#' @param enabled Logical indicating whether PERMANOVA should be performed.
#'
#' @return A PERMANOVA result data frame, `NULL` when skipped, or `NULL` with an
#'   attached note when the required package is unavailable.
#'
#' @keywords internal
#' @noRd
.dq_permanova <- function(transformed, metadata, factors, permutations, enabled) {
  if (!isTRUE(enabled)) return(NULL)
  if (!requireNamespace("vegan", quietly = TRUE)) {
    return(structure(NULL, note = "Package 'vegan' is not installed; PERMANOVA skipped."))
  }
  if (!length(factors)) return(NULL)

  keep <- vapply(factors, function(cl) length(unique(metadata[[cl]])) > 1L, logical(1))
  factors <- factors[keep]
  if (!length(factors)) return(NULL)

  dat <- metadata[, factors, drop = FALSE]
  form <- stats::as.formula(paste("dist_mat ~", paste(factors, collapse = " + ")))
  dist_mat <- stats::dist(t(transformed))

  fit <- vegan::adonis2(
    formula = form,
    data = dat,
    permutations = as.integer(permutations),
    by = "margin"
  )

  tab <- as.data.frame(fit)
  tab$factor <- rownames(tab)
  rownames(tab) <- NULL
  tab <- tab[tab$factor %in% factors, , drop = FALSE]
  tab[, c("factor", setdiff(names(tab), "factor")), drop = FALSE]
}



#' Detect dominant experimental factors
#'
#' Internal helper that compares PCA-associated effect sizes for the biological
#' group, experimental unit, and batch variables.
#'
#' @param pca_effects A data frame returned by `.dq_pca_effects()`.
#' @param group_col Character string naming the biological group column.
#' @param unit_col Optional character string naming the experimental unit column.
#' @param batch_cols Character vector naming batch-related columns.
#' @param dominance_ratio Numeric threshold defining how much larger a
#'   non-group effect must be relative to the group effect.
#' @param min_effect_r2 Numeric minimum weighted adjusted R-squared required for
#'   a factor to be considered dominant.
#'
#' @return A list describing the dominant factor, group effect, dominant effect,
#'   ratio relative to the group, and whether a non-group factor dominates.
#'
#' @keywords internal
#' @noRd
.dq_detect_dominance <- function(pca_effects, group_col, unit_col, batch_cols,
                                 dominance_ratio, min_effect_r2) {
  if (!nrow(pca_effects)) {
    return(list(
      dominant_factor = "undetermined",
      group_r2 = NA_real_,
      dominant_r2 = NA_real_,
      ratio_to_group = NA_real_,
      is_non_group_dominant = FALSE
    ))
  }

  group_r2 <- pca_effects$weighted_r2[match(group_col, pca_effects$factor)]
  group_r2 <- if (length(group_r2) && is.finite(group_r2)) group_r2 else 0

  candidates <- unique(c(unit_col, batch_cols))
  candidate_rows <- pca_effects[pca_effects$factor %in% candidates, , drop = FALSE]

  if (!nrow(candidate_rows)) {
    return(list(
      dominant_factor = group_col,
      group_r2 = group_r2,
      dominant_r2 = group_r2,
      ratio_to_group = 1,
      is_non_group_dominant = FALSE
    ))
  }

  i <- which.max(candidate_rows$weighted_r2)
  best <- candidate_rows[i, , drop = FALSE]
  ratio <- best$weighted_r2 / max(group_r2, 1e-6)
  dominates <- is.finite(best$weighted_r2) &&
    best$weighted_r2 >= min_effect_r2 &&
    ratio >= dominance_ratio

  list(
    dominant_factor = if (dominates) best$factor else group_col,
    group_r2 = group_r2,
    dominant_r2 = if (dominates) best$weighted_r2 else group_r2,
    ratio_to_group = ratio,
    is_non_group_dominant = dominates
  )
}


#' Recommend a differential expression model
#'
#' Internal helper that recommends DESeq2, dream, or unit-level aggregation
#' based on the experimental-unit structure, repeated measurements, and PCA
#' dominance assessment.
#'
#' @param metadata A sample metadata data frame.
#' @param group_col Character string naming the biological group column.
#' @param unit_col Optional character string naming the experimental unit column.
#' @param repeated A repeated-structure summary returned by
#'   `.dq_repeated_structure()`.
#' @param dominance A dominance summary returned by `.dq_detect_dominance()`.
#'
#' @return A list containing the recommended method, model formula, rationale,
#'   alternative strategy, and design warnings.
#'
#' @keywords internal
#' @noRd
.dq_recommend_model <- function(metadata, group_col, unit_col, repeated, dominance) {
  warnings <- character()
  group_n <- table(metadata[[group_col]])

  if (is.null(unit_col)) {
    return(list(
      method = "DESeq2",
      formula = paste0("~ ", group_col),
      reason = "No repeated experimental unit was specified or detected.",
      alternative = NA_character_,
      warnings = warnings
    ))
  }

  unit_by_group <- tapply(metadata[[unit_col]], metadata[[group_col]], function(z) {
    length(unique(as.character(z)))
  })

  if (any(unit_by_group < 2L)) {
    warnings <- c(
      warnings,
      "At least one group has fewer than two independent units; treatment inference is not reliably estimable."
    )
  }

  if (!repeated$has_repeated_units) {
    return(list(
      method = "DESeq2",
      formula = paste0("~ ", group_col),
      reason = paste0("Each ", unit_col, " occurs once; no within-unit correlation is present."),
      alternative = NA_character_,
      warnings = warnings
    ))
  }

  if (repeated$unit_crosses_group) {
    return(list(
      method = "dream",
      formula = paste0("~ ", group_col, " + (1|", unit_col, ")"),
      reason = paste0(
        "The same ", unit_col,
        " contributes observations to multiple groups, requiring a repeated-measures model."
      ),
      alternative = paste0("Paired DESeq2 with design ~ ", unit_col, " + ", group_col,
                           " may be used when the fixed-effect design matrix has full rank."),
      warnings = warnings
    ))
  }

  reason <- paste0(
    "Multiple samples share each ", unit_col,
    "; treating all samples as independent would create pseudoreplication"
  )
  if (dominance$is_non_group_dominant && identical(dominance$dominant_factor, unit_col)) {
    reason <- paste0(reason, ", and ", unit_col, " dominates the assessed PCA structure.")
  } else {
    reason <- paste0(reason, ".")
  }

  list(
    method = "dream",
    formula = paste0("~ ", group_col, " + (1|", unit_col, ")"),
    reason = reason,
    alternative = paste0(
      "Aggregate raw counts by ", unit_col,
      " and run DESeq2 (~ ", group_col,
      ") when the scientific estimand is the unit-level effect."
    ),
    warnings = warnings
  )
}
