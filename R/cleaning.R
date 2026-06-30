# ============================================================ #
# DEGgo cleaning pipeline
# ============================================================ #

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Clean one DEG table using expression support
#'
#' Remove significant DEGs with weak biological support.
#'
#' @keywords internal
#' @noRd
#'
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
    min_group_mean = 10,
    expr_filter_groups = "auto"
) {

  if (is.null(res_df) || !nrow(res_df))
    return(res_df)

  if (!"ENSEMBL" %in% names(res_df)) {

    if ("gene_id" %in% names(res_df)) {

      res_df$ENSEMBL <- res_df$gene_id

    } else {

      stop("res_df must contain ENSEMBL or gene_id.")

    }
  }

  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"

  rownames(counts) <- sub("\\..*$", "", rownames(counts))
  res_df$ENSEMBL <- sub("\\..*$", "", res_df$ENSEMBL)

  metadata <- as.data.frame(metadata)

  if (!"sample" %in% names(metadata))
    metadata$sample <- rownames(metadata)

  rownames(metadata) <- metadata$sample

  samples <- intersect(colnames(counts), rownames(metadata))

  if (!is.null(sample_subset))
    samples <- intersect(samples, sample_subset)

  counts <- counts[, samples, drop = FALSE]
  metadata <- metadata[samples, , drop = FALSE]

  ## ----------------------------- #
  ## determine grouping variable
  ## ----------------------------- #

  group_col <- expr_filter_groups

  if (identical(group_col, "auto")) {

    if (!is.null(comparison)) {

      if (grepl("Male_vs_Female|Female_vs_Male", comparison)) {

        group_col <- "sex"

      } else if (grepl("PAMH_vs_PBS|PBS_vs_PAMH", comparison)) {

        group_col <- "treatment"

      }
    }

    if (is.null(group_col) || !group_col %in% names(metadata)) {

      for (x in c("condition",
                  "treatment",
                  "sex")) {

        if (x %in% names(metadata)) {

          group_col <- x
          break

        }

      }

    }

  }

  if (!group_col %in% names(metadata)) {

    stop(
      "Cannot determine expr_filter_groups."
    )

  }

  groups <- metadata[colnames(counts), group_col]

  keep <- logical(nrow(res_df))

  ids <- res_df$ENSEMBL

  for (i in seq_along(ids)) {

    id <- ids[i]

    if (!id %in% rownames(counts))
      next

    x <- counts[id, ]

    if (sum(x) == 0)
      next

    ## one-sample domination

    frac <- max(x) / sum(x)

    if (frac > max_sample_fraction)
      next

    ok <- FALSE

    for (g in unique(groups)) {

      z <- x[groups == g]

      required <- min_expr_samples

      if (!is.null(min_prevalence)) {

        required <- max(
          required,
          ceiling(length(z) * min_prevalence)
        )

      }

      if (

        sum(z >= min_expr_count) >= required &&
        mean(z) >= min_group_mean

      ) {

        ok <- TRUE

      }

    }

    keep[i] <- ok

  }

  res_df[keep, , drop = FALSE]

}


#' Clean all DEG tables
#'
#' @keywords internal
#' @noRd
#'
.deggo_clean_deg_tables <- function(
    de_results,
    counts,
    metadata,
    min_expr_count,
    min_expr_samples,
    min_prevalence,
    max_sample_fraction,
    min_group_mean,
    expr_filter_groups,
    log = NULL
) {

  if (is.null(de_results$sig_deg))
    return(de_results)

  de_results$sig_deg_clean <- list()

  for (nm in names(de_results$sig_deg)) {

    subset <- NULL

    if (!is.null(de_results$samples))
      subset <- de_results$samples[[nm]]

    clean <- .deggo_clean_deg_table(

      res_df = de_results$sig_deg[[nm]],
      counts = counts,
      metadata = metadata,
      sample_subset = subset,
      comparison = nm,
      min_expr_count = min_expr_count,
      min_expr_samples = min_expr_samples,
      min_prevalence = min_prevalence,
      max_sample_fraction = max_sample_fraction,
      min_group_mean = min_group_mean,
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
          " retained"
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
