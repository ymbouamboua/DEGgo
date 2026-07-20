# ============================================================ #
# DEGgo cleaning pipeline
# ============================================================ #

`%||%` <- function(x, y) if (is.null(x)) y else x


#' Clean one DEG table using contrast-specific expression support
#'
#' Removes significant genes with weak, sparse, or outlier-driven expression.
#' Filtering is performed on raw counts using only samples belonging to the
#' current comparison.
#'
#' @param res_df Differential-expression result table.
#' @param counts Raw count matrix or count data frame.
#' @param metadata Sample metadata.
#' @param sample_subset Samples belonging to the current comparison.
#' @param comparison Comparison name.
#' @param min_expr_count Minimum raw count defining an expressed sample.
#' @param min_expr_samples Minimum expressed samples required within one group.
#' @param min_prevalence Minimum expressed-sample proportion within one group.
#' @param max_sample_fraction Maximum contribution of one sample to expression
#'   across all contrast samples.
#' @param max_group_sample_fraction Maximum contribution of one sample within
#'   its biological group.
#' @param min_group_mean Minimum mean raw count within one group.
#' @param min_group_median Minimum median raw count within one group. When NULL,
#'   defaults to `min_expr_count`.
#' @param max_group_cv Optional maximum coefficient of variation within a group.
#'   Use NULL to disable this filter.
#' @param expr_filter_groups Grouping metadata column or `"auto"`.
#' @param sample_col Sample identifier column.
#' @param gene_col Candidate gene identifier columns.
#'
#' @return Filtered DEG table with diagnostic expression-support columns.
#'
#' @keywords internal
#' @noRd
.deggo_clean_deg_table <- function(
    res_df,
    counts,
    metadata,
    sample_subset = NULL,
    comparison = NULL,
    min_expr_count = 20,
    min_expr_samples = 3,
    min_prevalence = 0.6,
    max_sample_fraction = 0.45,
    max_group_sample_fraction = 0.45,
    min_group_mean = 10,
    min_group_median = NULL,
    max_group_cv = NULL,
    expr_filter_groups = "auto",
    sample_col = "sample",
    gene_col = c(
      "gene_id",
      "ENSEMBL",
      "ensembl",
      "ensembl_id",
      "GeneID"
    )
) {

  # ---------------------------------------------------------- #
  # Validate DEG table
  # ---------------------------------------------------------- #

  if (
    is.null(res_df) ||
    !is.data.frame(res_df) ||
    nrow(res_df) == 0L
  ) {
    return(res_df)
  }

  res_df <- as.data.frame(
    res_df,
    stringsAsFactors = FALSE
  )

  if (!"ENSEMBL" %in% colnames(res_df)) {
    if ("gene_id" %in% colnames(res_df)) {
      res_df$ENSEMBL <- res_df$gene_id
    } else {
      stop(
        "`res_df` must contain an `ENSEMBL` or `gene_id` column.",
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------- #
  # Validate filtering parameters
  # ---------------------------------------------------------- #

  .check_scalar <- function(
    x,
    name,
    minimum = -Inf,
    maximum = Inf,
    integer = FALSE,
    allow_null = FALSE
  ) {

    if (is.null(x)) {
      if (allow_null) {
        return(invisible(TRUE))
      }

      stop(
        "`", name, "` cannot be NULL.",
        call. = FALSE
      )
    }

    if (
      length(x) != 1L ||
      is.na(x) ||
      !is.numeric(x) ||
      !is.finite(x) ||
      x < minimum ||
      x > maximum ||
      (integer && x != as.integer(x))
    ) {
      stop(
        "`", name, "` contains an invalid value.",
        call. = FALSE
      )
    }

    invisible(TRUE)
  }

  .check_scalar(
    min_expr_count,
    "min_expr_count",
    minimum = 0
  )

  .check_scalar(
    min_expr_samples,
    "min_expr_samples",
    minimum = 1,
    integer = TRUE
  )

  if (!is.null(min_prevalence)) {
    .check_scalar(
      min_prevalence,
      "min_prevalence",
      minimum = 0,
      maximum = 1
    )
  }

  .check_scalar(
    max_sample_fraction,
    "max_sample_fraction",
    minimum = 0,
    maximum = 1
  )

  .check_scalar(
    max_group_sample_fraction,
    "max_group_sample_fraction",
    minimum = 0,
    maximum = 1
  )

  .check_scalar(
    min_group_mean,
    "min_group_mean",
    minimum = 0
  )

  if (is.null(min_group_median)) {
    min_group_median <- min_expr_count
  }

  .check_scalar(
    min_group_median,
    "min_group_median",
    minimum = 0
  )

  if (!is.null(max_group_cv)) {
    .check_scalar(
      max_group_cv,
      "max_group_cv",
      minimum = 0
    )
  }

  min_expr_samples <- as.integer(
    min_expr_samples
  )

  # ---------------------------------------------------------- #
  # Prepare metadata
  # ---------------------------------------------------------- #

  metadata <- as.data.frame(
    metadata,
    stringsAsFactors = FALSE
  )

  if (!sample_col %in% colnames(metadata)) {

    if (
      is.null(rownames(metadata)) ||
      any(!nzchar(rownames(metadata)))
    ) {
      stop(
        "Metadata must contain sample identifiers.",
        call. = FALSE
      )
    }

    metadata[[sample_col]] <- rownames(metadata)
  }

  metadata[[sample_col]] <- as.character(
    metadata[[sample_col]]
  )

  if (
    anyNA(metadata[[sample_col]]) ||
    any(!nzchar(metadata[[sample_col]]))
  ) {
    stop(
      "Metadata contains missing or empty sample identifiers.",
      call. = FALSE
    )
  }

  if (anyDuplicated(metadata[[sample_col]])) {
    stop(
      "Metadata sample identifiers must be unique.",
      call. = FALSE
    )
  }

  rownames(metadata) <- metadata[[sample_col]]

  # ---------------------------------------------------------- #
  # Prepare count matrix
  # ---------------------------------------------------------- #

  counts_df <- as.data.frame(
    counts,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  gene_col_use <- gene_col[
    gene_col %in% colnames(counts_df)
  ][1]

  if (!is.na(gene_col_use)) {

    gene_ids <- as.character(
      counts_df[[gene_col_use]]
    )

    sample_columns <- intersect(
      rownames(metadata),
      colnames(counts_df)
    )

    if (length(sample_columns) < 2L) {
      stop(
        "Fewer than two metadata samples were found in counts.",
        call. = FALSE
      )
    }

    count_mat <- as.matrix(
      counts_df[, sample_columns, drop = FALSE]
    )

    rownames(count_mat) <- gene_ids

  } else {

    count_mat <- as.matrix(counts)

    if (
      is.null(rownames(count_mat)) ||
      any(!nzchar(rownames(count_mat)))
    ) {
      stop(
        paste0(
          "Counts must contain a gene ID column or gene IDs ",
          "as row names."
        ),
        call. = FALSE
      )
    }
  }

  suppressWarnings(
    storage.mode(count_mat) <- "numeric"
  )

  if (anyNA(count_mat)) {
    stop(
      paste0(
        "Count matrix contains NA values after numeric conversion. ",
        "Only numeric sample columns should be supplied."
      ),
      call. = FALSE
    )
  }

  if (any(count_mat < 0)) {
    stop(
      "Count matrix contains negative values.",
      call. = FALSE
    )
  }

  rownames(count_mat) <- sub(
    "\\..*$",
    "",
    as.character(rownames(count_mat))
  )

  res_df$ENSEMBL <- sub(
    "\\..*$",
    "",
    as.character(res_df$ENSEMBL)
  )

  if (anyDuplicated(rownames(count_mat))) {
    count_mat <- rowsum(
      count_mat,
      group = rownames(count_mat),
      reorder = FALSE
    )
  }

  # ---------------------------------------------------------- #
  # Select contrast-specific samples
  # ---------------------------------------------------------- #

  samples <- intersect(
    colnames(count_mat),
    rownames(metadata)
  )

  if (!is.null(sample_subset)) {

    sample_subset <- unique(
      as.character(sample_subset)
    )

    sample_subset <- sample_subset[
      !is.na(sample_subset) &
        nzchar(sample_subset)
    ]

    samples <- intersect(
      samples,
      sample_subset
    )
  }

  if (length(samples) < 2L) {
    stop(
      "Fewer than two valid samples are available for comparison '",
      comparison %||% "unknown",
      "'.",
      call. = FALSE
    )
  }

  count_mat <- count_mat[
    ,
    samples,
    drop = FALSE
  ]

  metadata <- metadata[
    samples,
    ,
    drop = FALSE
  ]

  # ---------------------------------------------------------- #
  # Determine grouping variable
  # ---------------------------------------------------------- #

  group_col <- expr_filter_groups

  if (identical(group_col, "auto")) {

    group_col <- NULL

    if (!is.null(comparison)) {

      if (grepl(
        "Male_vs_Female|Female_vs_Male",
        comparison,
        ignore.case = TRUE
      )) {
        group_col <- "sex"

      } else if (grepl(
        "PAMH_vs_PBS|PBS_vs_PAMH",
        comparison,
        ignore.case = TRUE
      )) {
        group_col <- "treatment"
      }
    }

    if (
      is.null(group_col) ||
      !group_col %in% colnames(metadata)
    ) {

      candidates <- c(
        "group",
        "condition",
        "treatment",
        "sex"
      )

      group_col <- candidates[
        candidates %in% colnames(metadata)
      ][1]
    }
  }

  if (
    is.null(group_col) ||
    length(group_col) != 1L ||
    is.na(group_col) ||
    !group_col %in% colnames(metadata)
  ) {
    stop(
      paste0(
        "Cannot determine `expr_filter_groups`. ",
        "Provide a valid metadata column."
      ),
      call. = FALSE
    )
  }

  groups <- as.character(
    metadata[colnames(count_mat), group_col]
  )

  valid_group <-
    !is.na(groups) &
    nzchar(groups)

  count_mat <- count_mat[
    ,
    valid_group,
    drop = FALSE
  ]

  metadata <- metadata[
    valid_group,
    ,
    drop = FALSE
  ]

  groups <- groups[valid_group]

  if (ncol(count_mat) < 2L) {
    stop(
      "Fewer than two samples remain after group filtering.",
      call. = FALSE
    )
  }

  group_indices <- split(
    seq_along(groups),
    groups
  )

  group_indices <- group_indices[
    lengths(group_indices) > 0L
  ]

  if (!length(group_indices)) {
    stop(
      "No valid expression-filter groups were found.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------- #
  # Match DEG genes to count matrix
  # ---------------------------------------------------------- #

  gene_index <- match(
    res_df$ENSEMBL,
    rownames(count_mat)
  )

  matched <- !is.na(gene_index)
  matched_rows <- which(matched)

  keep <- rep(
    FALSE,
    nrow(res_df)
  )

  # Diagnostic columns
  res_df$expr_filter_n_samples <- NA_integer_
  res_df$expr_filter_prevalence <- NA_real_
  res_df$expr_filter_max_fraction <- NA_real_
  res_df$expr_filter_max_group_fraction <- NA_real_
  res_df$expr_filter_max_group_mean <- NA_real_
  res_df$expr_filter_max_group_median <- NA_real_
  res_df$expr_filter_min_supporting_cv <- NA_real_
  res_df$expr_filter_support_groups <- NA_integer_

  if (!any(matched)) {
    return(
      res_df[0, , drop = FALSE]
    )
  }

  expr <- count_mat[
    gene_index[matched],
    ,
    drop = FALSE
  ]

  # ---------------------------------------------------------- #
  # Overall expression diagnostics
  # ---------------------------------------------------------- #

  total_expression <- rowSums(
    expr,
    na.rm = TRUE
  )

  n_expr_overall <- rowSums(
    expr >= min_expr_count,
    na.rm = TRUE
  )

  prevalence_overall <- rowMeans(
    expr >= min_expr_count,
    na.rm = TRUE
  )

  max_fraction <- apply(
    expr,
    1,
    function(x) {

      total <- sum(x, na.rm = TRUE)

      if (!is.finite(total) || total <= 0) {
        return(Inf)
      }

      max(x, na.rm = TRUE) / total
    }
  )

  # ---------------------------------------------------------- #
  # Group-specific diagnostics
  # ---------------------------------------------------------- #

  n_genes <- nrow(expr)
  n_groups <- length(group_indices)

  support_matrix <- matrix(
    FALSE,
    nrow = n_genes,
    ncol = n_groups,
    dimnames = list(
      NULL,
      names(group_indices)
    )
  )

  group_mean_matrix <- matrix(
    NA_real_,
    nrow = n_genes,
    ncol = n_groups,
    dimnames = list(
      NULL,
      names(group_indices)
    )
  )

  group_median_matrix <- group_mean_matrix
  group_fraction_matrix <- group_mean_matrix
  group_cv_matrix <- group_mean_matrix
  group_n_supported_matrix <- group_mean_matrix
  group_prevalence_matrix <- group_mean_matrix

  for (j in seq_along(group_indices)) {

    ii <- group_indices[[j]]
    x <- expr[, ii, drop = FALSE]

    group_n <- ncol(x)

    required_samples <- min_expr_samples

    if (!is.null(min_prevalence)) {
      required_samples <- max(
        required_samples,
        ceiling(group_n * min_prevalence)
      )
    }

    n_supported <- rowSums(
      x >= min_expr_count,
      na.rm = TRUE
    )

    prevalence <- rowMeans(
      x >= min_expr_count,
      na.rm = TRUE
    )

    group_mean <- rowMeans(
      x,
      na.rm = TRUE
    )

    group_median <- apply(
      x,
      1,
      stats::median,
      na.rm = TRUE
    )

    group_total <- rowSums(
      x,
      na.rm = TRUE
    )

    group_max_fraction <- vapply(
      seq_len(nrow(x)),
      function(i) {

        total <- group_total[i]

        if (!is.finite(total) || total <= 0) {
          return(Inf)
        }

        max(
          x[i, ],
          na.rm = TRUE
        ) / total
      },
      numeric(1)
    )

    group_sd <- apply(
      x,
      1,
      stats::sd,
      na.rm = TRUE
    )

    group_cv <- ifelse(
      is.finite(group_mean) &
        group_mean > 0,
      group_sd / group_mean,
      Inf
    )

    sample_support_pass <-
      n_supported >= required_samples

    prevalence_pass <- if (is.null(min_prevalence)) {
      rep(TRUE, n_genes)
    } else {
      prevalence >= min_prevalence
    }

    mean_pass <-
      is.finite(group_mean) &
      group_mean >= min_group_mean

    median_pass <-
      is.finite(group_median) &
      group_median >= min_group_median

    group_dominance_pass <-
      is.finite(group_max_fraction) &
      group_max_fraction <= max_group_sample_fraction

    cv_pass <- if (is.null(max_group_cv)) {
      rep(TRUE, n_genes)
    } else {
      is.finite(group_cv) &
        group_cv <= max_group_cv
    }

    support_matrix[, j] <-
      sample_support_pass &
      prevalence_pass &
      mean_pass &
      median_pass &
      group_dominance_pass &
      cv_pass

    group_mean_matrix[, j] <- group_mean
    group_median_matrix[, j] <- group_median
    group_fraction_matrix[, j] <- group_max_fraction
    group_cv_matrix[, j] <- group_cv
    group_n_supported_matrix[, j] <- n_supported
    group_prevalence_matrix[, j] <- prevalence
  }

  # ---------------------------------------------------------- #
  # Final filtering decision
  # ---------------------------------------------------------- #

  support_groups <- rowSums(
    support_matrix,
    na.rm = TRUE
  )

  overall_dominance_pass <-
    is.finite(max_fraction) &
    max_fraction <= max_sample_fraction

  nonzero_pass <-
    is.finite(total_expression) &
    total_expression > 0

  keep_matched <-
    support_groups >= 1L &
    overall_dominance_pass &
    nonzero_pass

  keep[matched_rows] <- keep_matched

  # ---------------------------------------------------------- #
  # Diagnostic summaries
  # ---------------------------------------------------------- #

  max_group_mean <- apply(
    group_mean_matrix,
    1,
    function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) NA_real_ else max(x)
    }
  )

  max_group_median <- apply(
    group_median_matrix,
    1,
    function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) NA_real_ else max(x)
    }
  )

  max_group_fraction <- apply(
    group_fraction_matrix,
    1,
    function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) Inf else max(x)
    }
  )

  min_supporting_cv <- vapply(
    seq_len(n_genes),
    function(i) {

      supporting <- which(
        support_matrix[i, ]
      )

      if (!length(supporting)) {
        return(NA_real_)
      }

      min(
        group_cv_matrix[i, supporting],
        na.rm = TRUE
      )
    },
    numeric(1)
  )

  # ---------------------------------------------------------- #
  # Attach diagnostics
  # ---------------------------------------------------------- #

  res_df$expr_filter_n_samples[matched_rows] <-
    n_expr_overall

  res_df$expr_filter_prevalence[matched_rows] <-
    prevalence_overall

  res_df$expr_filter_max_fraction[matched_rows] <-
    max_fraction

  res_df$expr_filter_max_group_fraction[matched_rows] <-
    max_group_fraction

  res_df$expr_filter_max_group_mean[matched_rows] <-
    max_group_mean

  res_df$expr_filter_max_group_median[matched_rows] <-
    max_group_median

  res_df$expr_filter_min_supporting_cv[matched_rows] <-
    min_supporting_cv

  res_df$expr_filter_support_groups[matched_rows] <-
    support_groups

  res_df[
    keep,
    ,
    drop = FALSE
  ]
}


#' Export cleaned DEG tables
#'
#' @keywords internal
#' @noRd
.deggo_clean_deg_tables <- function(
    de_results,
    counts,
    metadata,
    min_expr_count,
    min_expr_samples,
    min_prevalence,
    max_sample_fraction,
    max_group_sample_fraction,
    min_group_mean,
    min_group_median,
    max_group_cv,
    expr_filter_groups,
    log = NULL
) {

  if (
    is.null(de_results$sig_deg) ||
    !is.list(de_results$sig_deg)
  ) {
    return(de_results)
  }

  de_results$sig_deg_clean <- vector(
    mode = "list",
    length = length(de_results$sig_deg)
  )

  names(de_results$sig_deg_clean) <-
    names(de_results$sig_deg)

  is_pairwise <- !is.null(de_results$pairwise_contrasts) ||
    !is.null(de_results$contrast_list)

  for (nm in names(de_results$sig_deg)) {

    sample_subset <- NULL

    if (
      !is.null(de_results$samples) &&
      nm %in% names(de_results$samples)
    ) {
      sample_subset <- as.character(
        de_results$samples[[nm]]
      )
    }

    if (
      isTRUE(is_pairwise) &&
      (
        is.null(sample_subset) ||
        length(sample_subset) < 2L
      )
    ) {
      stop(
        "Missing contrast-specific samples for pairwise comparison '",
        nm,
        "'.",
        call. = FALSE
      )
    }

    clean <- .deggo_clean_deg_table(
      res_df = de_results$sig_deg[[nm]],
      counts = counts,
      metadata = metadata,
      sample_subset = sample_subset,
      comparison = nm,
      min_expr_count = min_expr_count,
      min_expr_samples = min_expr_samples,
      min_prevalence = min_prevalence,
      max_sample_fraction = max_sample_fraction,
      max_group_sample_fraction = max_group_sample_fraction,
      min_group_mean = min_group_mean,
      min_group_median = min_group_median,
      max_group_cv = max_group_cv,
      expr_filter_groups = expr_filter_groups
    )

    de_results$sig_deg_clean[[nm]] <- clean

    if (!is.null(log)) {
      log(
        paste0(
          nm,
          ": ",
          nrow(clean),
          "/",
          nrow(de_results$sig_deg[[nm]]),
          " retained using ",
          if (is.null(sample_subset)) {
            "all matched samples"
          } else {
            paste0(
              length(sample_subset),
              " contrast-specific samples"
            )
          }
        ),
        type = "info"
      )
    }
  }

  de_results
}




#' Export cleaned DEG tables
#'
#' @keywords internal
#' @noRd
.deggo_export_clean_deg_tables <- function(
    de_results,
    dirs
) {

  if (is.null(de_results$sig_deg_clean))
    return(invisible())

  for (nm in names(de_results$sig_deg_clean)) {

    utils::write.table(

      de_results$sig_deg_clean[[nm]],

      file.path(
        dirs$significant,
        paste0(
          nm,
          "_significant_clean.tsv"
        )
      ),

      sep = "\t",
      quote = FALSE,
      row.names = FALSE

    )

  }

  invisible()

}
