# ======================================================= #===== #
# DEGgo rhythmicity analysis (MetaCycle + cosinor)
# ======================================================= #===== #
#
# Standalone module for detecting rhythmic (e.g. circadian) genes in
# time-course expression data using MetaCycle (meta2d) and single-component
# cosinor regression (cosinor2), plus an optional cosinor-based test for
# differential rhythmicity between two groups.
#
# This module is intentionally independent from run_deggo(): it can be run
# on any expression matrix / DESeq2 dds object with a numeric time column in
# metadata.


# ======================================================= #== #
# CREATE RHYTHMICITY DIRECTORIES
# ======================================================= #== #
#' Create DEGgo rhythmicity output directories
#' @keywords internal
#' @noRd
.deggo_rhythm_dirs <- function(output_dir) {
  dirs <- list(
    results = file.path(output_dir, "rhythmicity_results"),
    plots = file.path(output_dir, "rhythmicity_plots")
  )

  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}


# ======================================================= #== #
# PREPARE RHYTHMICITY INPUT
# ======================================================= #== #
#' Prepare expression matrix and metadata for rhythmicity analysis
#'
#' Accepts either a plain genes x samples numeric matrix/data.frame, or a
#' DESeq2 `DESeqDataSet`, and returns a matched matrix and metadata table
#' ordered by time.
#'
#' @keywords internal
#' @noRd
.prepare_rhythm_input <- function(
    expr,
    metadata,
    sample_col = "sample",
    time_col = "time",
    group_col = NULL,
    assay = c(
      "vst",
      "raw",
      "normalized",
      "log2_normalized"
    )
) {

  assay <- match.arg(assay)

  metadata <- as.data.frame(
    metadata,
    stringsAsFactors = FALSE
  )

  # ======================================================= #
  # 1. Validate metadata
  # ======================================================= #

  if (!sample_col %in% colnames(metadata)) {
    stop(
      "metadata must contain a '",
      sample_col,
      "' column.",
      call. = FALSE
    )
  }

  if (!time_col %in% colnames(metadata)) {
    stop(
      "metadata must contain a '",
      time_col,
      "' column.",
      call. = FALSE
    )
  }

  if (!is.numeric(metadata[[time_col]])) {
    stop(
      "metadata[['",
      time_col,
      "']] must be numeric, for example time in hours.",
      call. = FALSE
    )
  }

  if (
    !is.null(group_col) &&
    !group_col %in% colnames(metadata)
  ) {
    stop(
      "metadata must contain a '",
      group_col,
      "' column.",
      call. = FALSE
    )
  }

  metadata[[sample_col]] <- as.character(
    metadata[[sample_col]]
  )

  if (anyNA(metadata[[sample_col]])) {
    stop(
      "metadata[['",
      sample_col,
      "']] contains missing sample identifiers.",
      call. = FALSE
    )
  }

  if (anyDuplicated(metadata[[sample_col]])) {
    duplicated_samples <- unique(
      metadata[[sample_col]][
        duplicated(metadata[[sample_col]])
      ]
    )

    stop(
      "metadata contains duplicated sample identifiers: ",
      paste(
        duplicated_samples,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  if (anyNA(metadata[[time_col]])) {
    stop(
      "metadata[['",
      time_col,
      "']] contains missing time values.",
      call. = FALSE
    )
  }

  # ======================================================= #
  # 2. Extract or validate expression matrix
  # ======================================================= #

  is_dds <- inherits(
    expr,
    "DESeqDataSet"
  )

  if (is_dds) {

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop(
        "Package 'DESeq2' is required to extract expression ",
        "from a DESeqDataSet object.",
        call. = FALSE
      )
    }

    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop(
        "Package 'SummarizedExperiment' is required to extract ",
        "the VST assay.",
        call. = FALSE
      )
    }

    mat <- switch(
      assay,

      raw = DESeq2::counts(
        expr,
        normalized = FALSE
      ),

      normalized = DESeq2::counts(
        expr,
        normalized = TRUE
      ),

      log2_normalized = log2(
        DESeq2::counts(
          expr,
          normalized = TRUE
        ) + 1
      ),

      vst = SummarizedExperiment::assay(
        DESeq2::vst(
          expr,
          blind = FALSE
        )
      )
    )

  } else {

    mat <- as.matrix(expr)

    if (!is.numeric(mat)) {
      stop(
        "For matrix input, expr must contain numeric expression values.",
        call. = FALSE
      )
    }

    if (assay == "raw") {

      if (any(mat < 0, na.rm = TRUE)) {
        stop(
          "assay = 'raw' requires non-negative expression values.",
          call. = FALSE
        )
      }

      mat <- log2(
        mat + 1
      )

    } else if (assay == "normalized") {

      if (any(mat < 0, na.rm = TRUE)) {
        stop(
          "assay = 'normalized' requires non-negative expression values.",
          call. = FALSE
        )
      }

      mat <- log2(
        mat + 1
      )

    } else if (assay == "log2_normalized") {

      # Matrix is already log2 transformed.
      mat <- mat

    } else if (assay == "vst") {

      # A plain matrix cannot be VST-transformed without a DESeqDataSet.
      # It is therefore assumed to already contain transformed values.
      warning(
        "assay = 'vst' was used with a matrix input. ",
        "The matrix is assumed to already contain VST-like or ",
        "otherwise transformed expression values.",
        call. = FALSE
      )
    }
  }

  # ======================================================= #
  # 3. Validate expression matrix
  # ======================================================= #

  if (is.null(rownames(mat))) {
    stop(
      "Expression input must have gene identifiers as row names.",
      call. = FALSE
    )
  }

  if (is.null(colnames(mat))) {
    stop(
      "Expression input must have column names matching sample IDs.",
      call. = FALSE
    )
  }

  if (anyDuplicated(colnames(mat))) {
    duplicated_samples <- unique(
      colnames(mat)[
        duplicated(colnames(mat))
      ]
    )

    stop(
      "Expression input contains duplicated sample names: ",
      paste(
        duplicated_samples,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  storage.mode(mat) <- "double"

  if (any(!is.finite(mat), na.rm = TRUE)) {
    stop(
      "Expression matrix contains non-finite values.",
      call. = FALSE
    )
  }

  # ======================================================= #
  # 4. Match expression and metadata samples
  # ======================================================= #

  common <- intersect(
    colnames(mat),
    metadata[[sample_col]]
  )

  if (length(common) < 4L) {
    stop(
      "At least 4 matched samples with time points are required ",
      "for rhythmicity analysis; found ",
      length(common),
      ".",
      call. = FALSE
    )
  }

  missing_in_metadata <- setdiff(
    colnames(mat),
    metadata[[sample_col]]
  )

  missing_in_expression <- setdiff(
    metadata[[sample_col]],
    colnames(mat)
  )

  if (length(missing_in_metadata) > 0L) {
    warning(
      length(missing_in_metadata),
      " expression sample(s) were not present in metadata and ",
      "will be excluded.",
      call. = FALSE
    )
  }

  if (length(missing_in_expression) > 0L) {
    warning(
      length(missing_in_expression),
      " metadata sample(s) were not present in the expression ",
      "matrix and will be excluded.",
      call. = FALSE
    )
  }

  metadata <- metadata[
    metadata[[sample_col]] %in% common,
    ,
    drop = FALSE
  ]

  metadata <- metadata[
    order(
      metadata[[time_col]],
      metadata[[sample_col]]
    ),
    ,
    drop = FALSE
  ]

  rownames(metadata) <- NULL

  mat <- mat[
    ,
    metadata[[sample_col]],
    drop = FALSE
  ]

  # ======================================================= #
  # 5. Validate final correspondence
  # ======================================================= #

  if (!identical(
    colnames(mat),
    metadata[[sample_col]]
  )) {
    stop(
      "Internal error: expression columns and metadata samples ",
      "are not in the same order.",
      call. = FALSE
    )
  }

  list(
    mat = mat,
    metadata = metadata,
    assay = assay,
    n_genes = nrow(mat),
    n_samples = ncol(mat)
  )
}



# ======================================================= #== #
# RUN METACYCLE
# ======================================================= #== #
#' Run MetaCycle (meta2d) rhythmicity detection
#' @keywords internal
#' @noRd
.run_metacycle <- function(
    mat,
    metadata,
    time_col,
    minper,
    maxper,
    cycMethod,
    log
) {

  time_vec <- metadata[[time_col]]

  in_df <- data.frame(CycID = rownames(mat), mat, check.names = FALSE)

  res <- tryCatch(
    MetaCycle::meta2d(
      infile = "deggo_rhythm_input.csv",
      inDF = in_df,
      filestyle = "csv",
      timepoints = time_vec,
      minper = minper,
      maxper = maxper,
      cycMethod = cycMethod,
      outputFile = FALSE,
      parallelize = FALSE
    ),
    error = function(e) {
      log(paste("MetaCycle failed:", conditionMessage(e)), type = "error")
      NULL
    }
  )

  if (is.null(res) || is.null(res$meta) || !nrow(res$meta)) {
    return(NULL)
  }

  out <- res$meta
  colnames(out)[colnames(out) == "CycID"] <- "gene"

  keep <- c("gene", grep("^meta2d_", colnames(out), value = TRUE))
  out <- out[, intersect(keep, colnames(out)), drop = FALSE]

  colnames(out) <- sub("^meta2d_", "metacycle_", colnames(out))
  colnames(out) <- sub("^metacycle_pvalue$", "metacycle_pvalue", colnames(out))
  colnames(out) <- sub("^metacycle_BH\\.Q$", "metacycle_padj", colnames(out))

  out$gene <- as.character(out$gene)
  out
}


# ======================================================= #== #
# MANUAL LM()-BASED COSINOR FALLBACK
# ======================================================= #== #
#' Fit a single-component cosinor model via lm(), without cosinor/cosinor2
#'
#' Fits `y ~ cos(2*pi*t/period) + sin(2*pi*t/period)` and back-transforms the
#' coefficients to mesor/amplitude/acrophase using the same parameterization
#' (`y = M + A*cos(wt - phi)`) as `cosinor::cosinor.lm()`'s transformed
#' table. Mesor, amplitude magnitude, and the zero-amplitude p-value match
#' the package engine exactly (both ultimately fit the same `lm()` model).
#'
#' Acrophase can occasionally differ by pi between engines: this
#' implementation uses `atan2(b2, b1)`, which is quadrant-aware, whereas
#' `cosinor:::summary.cosinor.lm()` derives its "acr" via plain
#' `atan(b2 / b1)`, which cannot distinguish a coefficient pair `(b1, b2)`
#' from `(-b1, -b2)` (e.g. when both raw coefficients are negative). The
#' `atan2()`-based acrophase here is the statistically correct one.
#'
#' @keywords internal
#' @noRd
.cosinor_fit_manual <- function(y, time, period) {

  w <- 2 * pi / period
  df <- data.frame(y = y, cc = cos(w * time), ss = sin(w * time))

  fit <- tryCatch(stats::lm(y ~ cc + ss, data = df), error = function(e) NULL)

  if (is.null(fit)) {
    return(
      list(
        mesor = NA_real_,
        amplitude = NA_real_,
        acrophase = NA_real_,
        pvalue = NA_real_
      )
    )
  }

  co <- stats::coef(fit)
  b1 <- unname(co["cc"])
  b2 <- unname(co["ss"])

  pvalue <- tryCatch({
    null_fit <- stats::lm(y ~ 1, data = df)
    stats::anova(null_fit, fit)[["Pr(>F)"]][2]
  }, error = function(e) NA_real_)

  list(
    mesor = unname(co["(Intercept)"]),
    amplitude = sqrt(b1^2 + b2^2),
    acrophase = atan2(b2, b1),
    pvalue = pvalue
  )
}


#' Test differential rhythmicity between two groups via lm() interaction
#'
#' Mirrors the `cosinor::cosinor.lm(Y ~ time(time) + group + amp.acro(group))`
#' vs `Y ~ time(time) + group` comparison, but built directly on `lm()`:
#' reduced model shares amplitude/acrophase across groups, full model lets
#' the cosine/sine coefficients differ by group. The F-test on the added
#' interaction terms is the differential-rhythmicity p-value.
#'
#' @keywords internal
#' @noRd
.cosinor_diff_manual <- function(y, time, group, period) {

  w <- 2 * pi / period
  df <- data.frame(
    y = y,
    cc = cos(w * time),
    ss = sin(w * time),
    group = factor(group)
  )

  tryCatch({
    fit_reduced <- stats::lm(y ~ cc + ss + group, data = df)
    fit_full <- stats::lm(y ~ (cc + ss) * group, data = df)
    a <- stats::anova(fit_reduced, fit_full)
    a[["Pr(>F)"]][2]
  }, error = function(e) NA_real_)
}


# ======================================================= #== #
# RUN SINGLE-GENE COSINOR
# ======================================================= #== #
#' Run single-component cosinor regression per gene
#'
#' Uses `cosinor::cosinor.lm()` + `cosinor2::cosinor.detect()` when
#' `use_pkg = TRUE`, otherwise falls back to a manual `lm()`-based fit
#' (see `.cosinor_fit_manual()`) so the analysis works without the optional
#' `cosinor`/`cosinor2` packages installed.
#'
#' @keywords internal
#' @noRd
.run_cosinor <- function(
    mat,
    metadata,
    time_col,
    cycle_length,
    log,
    use_pkg = TRUE
) {

  time_vec <- metadata[[time_col]]

  rows <- lapply(rownames(mat), function(g) {

    y <- as.numeric(mat[g, ])

    if (isTRUE(use_pkg)) {

      df <- data.frame(Y = y, time = time_vec)

      fit <- tryCatch(
        cosinor::cosinor.lm(Y ~ time(time), period = cycle_length, data = df),
        error = function(e) NULL
      )

      if (!is.null(fit)) {

        tt <- summary(fit)$transformed.table
        detect <- tryCatch(cosinor2::cosinor.detect(fit), error = function(e) NULL)

        acr <- tt["acr", "estimate"]
        phase_h <- ((acr %% (2 * pi)) / (2 * pi)) * cycle_length

        return(
          data.frame(
            gene = g,
            cosinor_mesor = tt["(Intercept)", "estimate"],
            cosinor_amplitude = tt["amp", "estimate"],
            cosinor_acrophase = acr,
            cosinor_phase = phase_h,
            cosinor_pvalue = if (!is.null(detect)) detect[1, "p"] else NA_real_,
            cosinor_engine = "cosinor2",
            stringsAsFactors = FALSE
          )
        )
      }
      # package call failed for this gene: fall through to manual fit
    }

    fit_m <- .cosinor_fit_manual(y, time_vec, cycle_length)
    phase_h <- ((fit_m$acrophase %% (2 * pi)) / (2 * pi)) * cycle_length

    data.frame(
      gene = g,
      cosinor_mesor = fit_m$mesor,
      cosinor_amplitude = fit_m$amplitude,
      cosinor_acrophase = fit_m$acrophase,
      cosinor_phase = phase_h,
      cosinor_pvalue = fit_m$pvalue,
      cosinor_engine = "manual_lm",
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$cosinor_padj <- stats::p.adjust(out$cosinor_pvalue, method = "BH")
  rownames(out) <- NULL
  out
}


# ======================================================= #== #
# RUN COSINOR DIFFERENTIAL RHYTHMICITY (GROUP INTERACTION)
# ======================================================= #== #
#' Test differential rhythmicity between two groups via cosinor interaction
#'
#' Uses `cosinor::cosinor.lm()` with `amp.acro(group)` when `use_pkg = TRUE`,
#' otherwise falls back to the manual `lm()` interaction test (see
#' `.cosinor_diff_manual()`).
#'
#' @keywords internal
#' @noRd
.run_cosinor_diff_group <- function(
    mat,
    metadata,
    time_col,
    group_col,
    cycle_length,
    log,
    use_pkg = TRUE
) {

  time_vec <- metadata[[time_col]]
  group_vec <- factor(metadata[[group_col]])

  rows <- lapply(rownames(mat), function(g) {

    y <- as.numeric(mat[g, ])

    pvalue <- NA_real_

    if (isTRUE(use_pkg)) {

      df <- data.frame(Y = y, time = time_vec, group = group_vec)

      pvalue <- tryCatch({
        fit_full <- cosinor::cosinor.lm(
          Y ~ time(time) + group + amp.acro(group),
          period = cycle_length,
          data = df
        )
        fit_reduced <- cosinor::cosinor.lm(
          Y ~ time(time) + group,
          period = cycle_length,
          data = df
        )
        a <- stats::anova(fit_reduced$fit, fit_full$fit)
        a[["Pr(>F)"]][2]
      }, error = function(e) NA_real_)
    }

    if (!isTRUE(use_pkg) || is.na(pvalue)) {
      pvalue <- .cosinor_diff_manual(y, time_vec, group_vec, cycle_length)
    }

    data.frame(gene = g, diff_rhythm_pvalue = pvalue, stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  out$diff_rhythm_padj <- stats::p.adjust(out$diff_rhythm_pvalue, method = "BH")
  rownames(out) <- NULL
  out
}


# ======================================================= #== #
# MERGE RHYTHMICITY RESULTS
# ======================================================= #== #
#' Merge MetaCycle, cosinor and differential-rhythmicity tables
#' @keywords internal
#' @noRd
.merge_rhythm_results <- function(
    meta_res,
    cos_res,
    diff_res,
    padj_cutoff = 0.05
) {

  tabs <- list(meta_res, cos_res, diff_res)
  tabs <- tabs[!vapply(tabs, is.null, logical(1))]

  if (!length(tabs)) {
    return(data.frame())
  }

  merged <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), tabs)

  merged$rhythmic_metacycle <- if ("metacycle_padj" %in% colnames(merged)) {
    merged$metacycle_padj < padj_cutoff
  } else {
    NA
  }

  merged$rhythmic_cosinor <- if ("cosinor_padj" %in% colnames(merged)) {
    merged$cosinor_padj < padj_cutoff
  } else {
    NA
  }

  merged$rhythmic_by <- mapply(
    function(m, c) {
      if (is.na(m) && is.na(c)) return(NA_character_)
      m <- isTRUE(m)
      c <- isTRUE(c)
      if (m && c) "both"
      else if (m) "metacycle_only"
      else if (c) "cosinor_only"
      else "none"
    },
    merged$rhythmic_metacycle,
    merged$rhythmic_cosinor
  )

  order_col <- intersect(c("metacycle_padj", "cosinor_padj"), colnames(merged))

  if (length(order_col)) {
    merged <- merged[order(merged[[order_col[1]]], na.last = TRUE), , drop = FALSE]
  }

  rownames(merged) <- NULL
  merged
}




# ======================================================= #== #
# RHYTHMICITY DIAGNOSTIC PLOTS
# ======================================================= #== #

#' Resolve readable gene labels
#'
#' Uses a supplied annotation table or symbol columns already present in the
#' rhythmicity summary. Ensembl IDs are retained as the internal identifiers.
#'
#' @param gene_ids Character vector of gene identifiers.
#' @param summary Rhythmicity summary table.
#' @param gene_annotation Optional annotation data frame containing gene IDs
#'   and gene symbols.
#' @param gene_id_col Gene identifier column in `gene_annotation`.
#' @param gene_symbol_col Optional gene-symbol column. When `NULL`, DEGgo
#'   searches common symbol-column names automatically.
#'
#' @return A named character vector with gene IDs as names and display labels
#'   as values.
#'
#' @keywords internal
#' @noRd
.deggo_rhythm_gene_labels <- function(
    gene_ids,
    summary = NULL,
    gene_annotation = NULL,
    gene_id_col = "gene_id",
    gene_symbol_col = NULL
) {

  gene_ids <- as.character(gene_ids)

  labels <- stats::setNames(
    gene_ids,
    gene_ids
  )

  # ======================================================= #
  # 1. Search the rhythmicity summary
  # ======================================================= #

  if (!is.null(summary)) {

    summary <- as.data.frame(
      summary,
      stringsAsFactors = FALSE
    )

    if (nrow(summary) > 0L) {

      summary_id_candidates <- unique(
        c(
          "gene",
          gene_id_col,
          "gene_id",
          "matrix_gene_id",
          "feature",
          "EnsemblID",
          "ENSEMBL"
        )
      )

      summary_symbol_candidates <- unique(
        c(
          gene_symbol_col,
          "gene_symbol",
          "symbol",
          "Gene",
          "SYMBOL",
          "gene_name",
          "external_gene_name",
          "feature_name"
        )
      )

      summary_id_candidates <- summary_id_candidates[
        !is.na(summary_id_candidates) &
          nzchar(summary_id_candidates)
      ]

      summary_symbol_candidates <- summary_symbol_candidates[
        !is.na(summary_symbol_candidates) &
          nzchar(summary_symbol_candidates)
      ]

      summary_id_col <- intersect(
        summary_id_candidates,
        colnames(summary)
      )

      summary_symbol_col <- intersect(
        summary_symbol_candidates,
        colnames(summary)
      )

      if (
        length(summary_id_col) > 0L &&
        length(summary_symbol_col) > 0L
      ) {

        summary_id_col <- summary_id_col[1L]
        summary_symbol_col <- summary_symbol_col[1L]

        matched_symbols <- as.character(
          summary[[summary_symbol_col]][
            match(
              gene_ids,
              as.character(summary[[summary_id_col]])
            )
          ]
        )

        valid <- !is.na(matched_symbols) &
          nzchar(trimws(matched_symbols)) &
          matched_symbols != "NA"

        labels[valid] <- matched_symbols[valid]
      }
    }
  }

  # ======================================================= #
  # 2. Search the external annotation table
  # ======================================================= #

  unresolved <- names(labels)[
    is.na(labels) |
      !nzchar(labels) |
      labels == names(labels)
  ]

  if (
    length(unresolved) > 0L &&
    !is.null(gene_annotation)
  ) {

    gene_annotation <- as.data.frame(
      gene_annotation,
      stringsAsFactors = FALSE
    )

    if (nrow(gene_annotation) > 0L) {

      annotation_id_candidates <- unique(
        c(
          gene_id_col,
          "matrix_gene_id",
          "gene_id",
          "gene",
          "feature",
          "EnsemblID",
          "ENSEMBL"
        )
      )

      annotation_symbol_candidates <- unique(
        c(
          gene_symbol_col,
          "gene_symbol",
          "symbol",
          "Gene",
          "SYMBOL",
          "gene_name",
          "external_gene_name",
          "feature_name"
        )
      )

      annotation_id_candidates <- annotation_id_candidates[
        !is.na(annotation_id_candidates) &
          nzchar(annotation_id_candidates)
      ]

      annotation_symbol_candidates <- annotation_symbol_candidates[
        !is.na(annotation_symbol_candidates) &
          nzchar(annotation_symbol_candidates)
      ]

      annotation_id_col <- intersect(
        annotation_id_candidates,
        colnames(gene_annotation)
      )

      annotation_symbol_col <- intersect(
        annotation_symbol_candidates,
        colnames(gene_annotation)
      )

      if (
        length(annotation_id_col) > 0L &&
        length(annotation_symbol_col) > 0L
      ) {

        annotation_id_col <- annotation_id_col[1L]
        annotation_symbol_col <- annotation_symbol_col[1L]

        matched_symbols <- as.character(
          gene_annotation[[annotation_symbol_col]][
            match(
              unresolved,
              as.character(
                gene_annotation[[annotation_id_col]]
              )
            )
          ]
        )

        valid <- !is.na(matched_symbols) &
          nzchar(trimws(matched_symbols)) &
          matched_symbols != "NA"

        labels[unresolved[valid]] <- matched_symbols[valid]
      }
    }
  }

  labels
}


#' Format rhythmicity p-values
#'
#' @keywords internal
#' @noRd
.deggo_rhythm_format_p <- function(x) {

  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }

  if (x < 0.001) {
    return(formatC(x, format = "e", digits = 2))
  }

  formatC(x, format = "f", digits = 3)
}



#' Add gene symbols to a rhythmicity result table
#'
#' @keywords internal
#' @noRd
.deggo_add_rhythm_annotation <- function(
    tab,
    gene_annotation = NULL,
    gene_id_col = "gene_id",
    gene_symbol_col = NULL
) {

  if (
    is.null(tab) ||
    !is.data.frame(tab) ||
    nrow(tab) == 0L ||
    is.null(gene_annotation) ||
    is.null(gene_symbol_col)
  ) {
    return(tab)
  }

  if (
    !gene_id_col %in% colnames(gene_annotation) ||
    !gene_symbol_col %in% colnames(gene_annotation)
  ) {
    return(tab)
  }

  annotation_map <- gene_annotation[
    ,
    c(
      gene_id_col,
      gene_symbol_col
    ),
    drop = FALSE
  ]

  colnames(annotation_map) <- c(
    "gene",
    "gene_symbol"
  )

  annotation_map$gene <- as.character(
    annotation_map$gene
  )

  annotation_map$gene_symbol <- as.character(
    annotation_map$gene_symbol
  )

  annotation_map <- annotation_map[
    !is.na(annotation_map$gene) &
      nzchar(trimws(annotation_map$gene)),
    ,
    drop = FALSE
  ]

  annotation_map <- annotation_map[
    !duplicated(annotation_map$gene),
    ,
    drop = FALSE
  ]

  original_gene_order <- as.character(
    tab$gene
  )

  out <- merge(
    tab,
    annotation_map,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )

  out <- out[
    match(
      original_gene_order,
      out$gene
    ),
    ,
    drop = FALSE
  ]

  rownames(out) <- NULL

  out
}




# ======================================================= #== #
# RUN DEGGO RHYTHMICITY (main entry point)
# ======================================================= #== #

#' Run DEGgo rhythmicity (circadian) analysis
#'
#' Standalone rhythmicity detection workflow using MetaCycle (`meta2d`) and
#' single-component cosinor regression, with an optional cosinor-based test
#' for differential rhythmicity between two groups.
#'
#' This function is independent from [run_deggo()] and operates on
#' time-course expression data supplied as either a numeric expression matrix
#' or a DESeq2 `DESeqDataSet`.
#'
#' @param expr A genes-by-samples numeric matrix/data frame of expression
#'   values, or a DESeq2 `DESeqDataSet` object.
#' @param metadata Sample metadata data frame. Must contain a sample
#'   identifier column defined by `sample_col` and a numeric time column
#'   defined by `time_col`.
#' @param sample_col Column in `metadata` containing sample identifiers that
#'   match the column names of `expr`.
#' @param time_col Column in `metadata` containing numeric time values, such
#'   as hours or Zeitgeber Time.
#' @param group_col Optional column in `metadata` defining exactly two groups
#'   for differential-rhythmicity testing.
#' @param assay Assay to extract when `expr` is a DESeq2 `DESeqDataSet`. One
#'   of `"vst"`, `"normalized"`, `"log2_normalized"`, or `"raw"`.
#' @param methods Rhythmicity-detection methods to run. One or both of
#'   `"meta2d"` and `"cosinor"`.
#' @param period_range Numeric vector of length two defining MetaCycle's period
#'   search range in the same units as `time_col`.
#' @param cycle_length Assumed period for the single-component cosinor model.
#' @param cycMethod Character vector of MetaCycle methods to combine. Allowed
#'   values include `"ARS"`, `"JTK"`, and `"LS"`.
#' @param padj_cutoff Adjusted p-value cutoff used to classify rhythmic genes.
#' @param cosinor_engine Cosinor fitting engine. `"auto"` uses the
#'   `cosinor` and `cosinor2` packages when installed and otherwise falls back
#'   to the manual `lm()` implementation. `"package"` requires the optional
#'   packages. `"manual"` always uses the internal `lm()` implementation.
#' @param gene_annotation Optional gene-annotation data frame used to map
#'   expression-matrix row identifiers to readable gene symbols.
#' @param gene_id_col Gene-identifier column in `gene_annotation`.
#' @param gene_symbol_col Optional gene-symbol column in `gene_annotation`.
#'   When `NULL`, common symbol-column names are detected automatically.
#' @param show_gene_id Logical. Include the original gene identifier in
#'   individual rhythmicity plot titles.
#' @param output_dir Output directory.
#' @param project_name Optional project name.
#' @param generate_plots Logical. Generate rhythmicity diagnostic plots.
#' @param n_top_plots Number of top rhythmic genes to plot individually.
#' @param txtsize Base text size.
#' @param seed Random seed.
#' @param verbose Logical. Print progress messages.
#'
#' @return A `deggo_rhythm_results` object containing the combined summary,
#'   method-specific result tables, matched metadata, expression matrix,
#'   annotation settings, parameters, output paths, and generated plots.
#'
#' @export
run_deggo_rhythmicity <- function(
    expr,
    metadata,
    sample_col = "sample",
    time_col = "time",
    group_col = NULL,
    assay = c(
      "vst",
      "raw",
      "normalized",
      "log2_normalized"
    ),
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
    padj_cutoff = 0.05,
    cosinor_engine = c(
      "auto",
      "package",
      "manual"
    ),
    output_dir = "DEGgo_rhythmicity",
    project_name = "DEGgo rhythmicity analysis",
    generate_plots = TRUE,
    n_top_plots = 20,
    txtsize = 12,
    seed = 4173,
    gene_annotation = NULL,
    gene_id_col = "gene_id",
    gene_symbol_col = NULL,
    show_gene_id = TRUE,
    verbose = TRUE
) {

  set.seed(seed)

  assay <- match.arg(assay)

  cosinor_engine <- match.arg(
    cosinor_engine
  )

  methods <- match.arg(
    methods,
    choices = c(
      "meta2d",
      "cosinor"
    ),
    several.ok = TRUE
  )

  methods <- unique(methods)

  log <- .deggo_msg(
    verbose = verbose,
    prefix = "DEGgo"
  )

  log(
    "Starting DEGgo rhythmicity analysis",
    type = "header"
  )

  t_start <- Sys.time()

  # ======================================================= #
  # 1. Validate arguments
  # ======================================================= #

  if (
    !is.numeric(period_range) ||
    length(period_range) != 2L ||
    anyNA(period_range) ||
    any(!is.finite(period_range)) ||
    period_range[1L] <= 0 ||
    period_range[2L] <= period_range[1L]
  ) {
    stop(
      "'period_range' must contain two finite values with ",
      "0 < minimum < maximum.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(cycle_length) ||
    length(cycle_length) != 1L ||
    is.na(cycle_length) ||
    !is.finite(cycle_length) ||
    cycle_length <= 0
  ) {
    stop(
      "'cycle_length' must be one positive finite value.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(padj_cutoff) ||
    length(padj_cutoff) != 1L ||
    is.na(padj_cutoff) ||
    padj_cutoff <= 0 ||
    padj_cutoff >= 1
  ) {
    stop(
      "'padj_cutoff' must be between 0 and 1.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(n_top_plots) ||
    length(n_top_plots) != 1L ||
    is.na(n_top_plots) ||
    n_top_plots < 0
  ) {
    stop(
      "'n_top_plots' must be a non-negative number.",
      call. = FALSE
    )
  }

  n_top_plots <- as.integer(
    n_top_plots
  )

  allowed_cyc_methods <- c(
    "ARS",
    "JTK",
    "LS"
  )

  if (
    !is.character(cycMethod) ||
    !length(cycMethod) ||
    any(!cycMethod %in% allowed_cyc_methods)
  ) {
    stop(
      "'cycMethod' must contain one or more of: ",
      paste(
        allowed_cyc_methods,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  cycMethod <- unique(cycMethod)

  # ======================================================= #
  # 2. Prepare expression and metadata once
  # ======================================================= #

  prepared <- .prepare_rhythm_input(
    expr = expr,
    metadata = metadata,
    sample_col = sample_col,
    time_col = time_col,
    group_col = group_col,
    assay = assay
  )

  mat <- prepared$mat
  md <- prepared$metadata

  if (anyDuplicated(rownames(mat))) {
    stop(
      "The prepared expression matrix contains duplicated gene IDs.",
      call. = FALSE
    )
  }

  log(
    paste0(
      nrow(mat),
      " genes x ",
      ncol(mat),
      " samples prepared for rhythmicity analysis."
    ),
    type = "info"
  )

  time_values <- sort(
    unique(md[[time_col]])
  )

  log(
    paste0(
      "Time points: ",
      paste(
        time_values,
        collapse = ", "
      )
    ),
    type = "info"
  )

  # ======================================================= #
  # 3. Validate gene annotation
  # ======================================================= #

  detected_symbol_col <- gene_symbol_col

  if (!is.null(gene_annotation)) {

    gene_annotation <- as.data.frame(
      gene_annotation,
      stringsAsFactors = FALSE
    )

    if (!gene_id_col %in% colnames(gene_annotation)) {
      stop(
        "gene_annotation must contain ID column '",
        gene_id_col,
        "'.",
        call. = FALSE
      )
    }

    gene_annotation[[gene_id_col]] <- as.character(
      gene_annotation[[gene_id_col]]
    )

    gene_annotation <- gene_annotation[
      !is.na(gene_annotation[[gene_id_col]]) &
        nzchar(
          trimws(
            gene_annotation[[gene_id_col]]
          )
        ),
      ,
      drop = FALSE
    ]

    if (is.null(detected_symbol_col)) {

      symbol_candidates <- c(
        "symbol",
        "gene_symbol",
        "Gene",
        "SYMBOL",
        "gene_name",
        "external_gene_name",
        "feature"
      )

      available_symbol_cols <- intersect(
        symbol_candidates,
        setdiff(
          colnames(gene_annotation),
          gene_id_col
        )
      )

      if (length(available_symbol_cols) > 0L) {
        detected_symbol_col <- available_symbol_cols[1L]
      }
    }

    if (
      !is.null(detected_symbol_col) &&
      !detected_symbol_col %in%
      colnames(gene_annotation)
    ) {
      stop(
        "gene_annotation does not contain symbol column '",
        detected_symbol_col,
        "'.",
        call. = FALSE
      )
    }

    if (!is.null(detected_symbol_col)) {
      gene_annotation[[detected_symbol_col]] <- as.character(
        gene_annotation[[detected_symbol_col]]
      )
    }

    if (!is.null(detected_symbol_col)) {

      valid_symbol <- !is.na(
        gene_annotation[[detected_symbol_col]]
      ) &
        nzchar(
          trimws(
            gene_annotation[[detected_symbol_col]]
          )
        )

      annotation_order <- order(
        !valid_symbol
      )

      gene_annotation <- gene_annotation[
        annotation_order,
        ,
        drop = FALSE
      ]
    }

    gene_annotation <- gene_annotation[
      !duplicated(
        gene_annotation[[gene_id_col]]
      ),
      ,
      drop = FALSE
    ]
  }

  # ======================================================= #
  # 4. Dependencies and cosinor engine
  # ======================================================= #

  if (
    "meta2d" %in% methods &&
    !requireNamespace(
      "MetaCycle",
      quietly = TRUE
    )
  ) {
    stop(
      "Package 'MetaCycle' is required for MetaCycle analysis.",
      call. = FALSE
    )
  }

  use_cosinor_pkg <- FALSE

  if ("cosinor" %in% methods) {

    package_engine_available <-
      requireNamespace(
        "cosinor",
        quietly = TRUE
      ) &&
      requireNamespace(
        "cosinor2",
        quietly = TRUE
      )

    use_cosinor_pkg <- switch(
      cosinor_engine,

      package = {
        if (!package_engine_available) {
          stop(
            "cosinor_engine = 'package' requires packages ",
            "'cosinor' and 'cosinor2'.",
            call. = FALSE
          )
        }

        TRUE
      },

      manual = FALSE,

      auto = package_engine_available
    )

    if (isTRUE(use_cosinor_pkg)) {
      log(
        "Cosinor engine: 'cosinor'/'cosinor2' packages.",
        type = "info"
      )
    } else {
      log(
        "Cosinor engine: manual lm()-based fit.",
        type = "info"
      )
    }
  }

  # ======================================================= #
  # 5. Output directories
  # ======================================================= #

  dirs <- .deggo_rhythm_dirs(
    output_dir
  )

  # ======================================================= #
  # 6. MetaCycle
  # ======================================================= #

  meta_res <- NULL

  if ("meta2d" %in% methods) {

    log(
      "Running MetaCycle (meta2d)...",
      type = "step"
    )

    meta_res <- .run_metacycle(
      mat = mat,
      metadata = md,
      time_col = time_col,
      minper = period_range[1L],
      maxper = period_range[2L],
      cycMethod = cycMethod,
      log = log
    )

    if (is.null(meta_res)) {
      log(
        "MetaCycle returned no usable results.",
        type = "warn"
      )
    }
  }

  # ======================================================= #
  # 7. Cosinor
  # ======================================================= #

  cos_res <- NULL
  diff_res <- NULL

  if ("cosinor" %in% methods) {

    log(
      "Running single-component cosinor regression...",
      type = "step"
    )

    cos_res <- .run_cosinor(
      mat = mat,
      metadata = md,
      time_col = time_col,
      cycle_length = cycle_length,
      log = log,
      use_pkg = use_cosinor_pkg
    )

    if (!is.null(group_col)) {

      group_values <- as.character(
        md[[group_col]]
      )

      group_values <- group_values[
        !is.na(group_values) &
          nzchar(trimws(group_values))
      ]

      group_levels <- unique(
        group_values
      )

      if (length(group_levels) == 2L) {

        log(
          paste0(
            "Testing differential rhythmicity: ",
            paste(
              group_levels,
              collapse = " vs "
            )
          ),
          type = "step"
        )

        diff_res <- .run_cosinor_diff_group(
          mat = mat,
          metadata = md,
          time_col = time_col,
          group_col = group_col,
          cycle_length = cycle_length,
          log = log,
          use_pkg = use_cosinor_pkg
        )

      } else {

        log(
          paste0(
            "Differential rhythmicity skipped: '",
            group_col,
            "' contains ",
            length(group_levels),
            " levels; exactly two are required."
          ),
          type = "warn"
        )
      }
    }
  }

  # ======================================================= #
  # 8. Merge method-specific results
  # ======================================================= #

  merged <- .merge_rhythm_results(
    meta_res = meta_res,
    cos_res = cos_res,
    diff_res = diff_res,
    padj_cutoff = padj_cutoff
  )

  if (!is.data.frame(merged)) {
    merged <- as.data.frame(
      merged,
      stringsAsFactors = FALSE
    )
  }

  # ======================================================= #
  # 9. Add gene symbols
  # ======================================================= #

  meta_res <- .deggo_add_rhythm_annotation(
    tab = meta_res,
    gene_annotation = gene_annotation,
    gene_id_col = gene_id_col,
    gene_symbol_col = detected_symbol_col
  )

  cos_res <- .deggo_add_rhythm_annotation(
    tab = cos_res,
    gene_annotation = gene_annotation,
    gene_id_col = gene_id_col,
    gene_symbol_col = detected_symbol_col
  )

  diff_res <- .deggo_add_rhythm_annotation(
    tab = diff_res,
    gene_annotation = gene_annotation,
    gene_id_col = gene_id_col,
    gene_symbol_col = detected_symbol_col
  )

  merged <- .deggo_add_rhythm_annotation(
    tab = merged,
    gene_annotation = gene_annotation,
    gene_id_col = gene_id_col,
    gene_symbol_col = detected_symbol_col
  )

  # ======================================================= #
  # 10. Write result tables
  # ======================================================= #

  write_result <- function(
    object,
    filename
  ) {

    if (
      !is.null(object) &&
      is.data.frame(object) &&
      nrow(object) > 0L
    ) {
      utils::write.table(
        object,
        file.path(
          dirs$results,
          filename
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  }

  write_result(
    meta_res,
    "metacycle_results.tsv"
  )

  write_result(
    cos_res,
    "cosinor_results.tsv"
  )

  write_result(
    diff_res,
    "cosinor_differential_rhythmicity.tsv"
  )

  write_result(
    merged,
    "rhythmicity_summary.tsv"
  )

  # ======================================================= #
  # 11. Generate plots
  # ======================================================= #

  plots <- list()

  expression_label <- switch(
    assay,
    vst = "VST expression",
    raw = "log2(count + 1)",
    normalized = "log2 normalized expression",
    log2_normalized = "log2 normalized expression"
  )

  if (
    isTRUE(generate_plots) &&
    nrow(merged) > 0L
  ) {

    log(
      paste0(
        "Generating rhythmicity plots in: ",
        normalizePath(
          dirs$plots,
          winslash = "/",
          mustWork = FALSE
        )
      ),
      type = "info"
    )

    log(
      paste0(
        "Plot summary dimensions: ",
        nrow(merged),
        " genes x ",
        ncol(merged),
        " columns."
      ),
      type = "info"
    )

    plots <- .deggo_rhythm_plots(
      mat = mat,
      metadata = md,
      summary = merged,
      time_col = time_col,
      sample_col = sample_col,
      group_col = group_col,
      cycle_length = cycle_length,
      output_dir = dirs$plots,
      n_top = n_top_plots,
      txtsize = txtsize,
      expression_label = expression_label,
      padj_cutoff = padj_cutoff,
      gene_annotation = gene_annotation,
      gene_id_col = gene_id_col,
      gene_symbol_col = detected_symbol_col,
      show_gene_id = show_gene_id,
      fit_full_cycle = TRUE,
      duplicate_cycle = FALSE,
      shade_phases = TRUE,
      lights_on = 0,
      lights_off = 12,
      phase_prefix = "ZT",
      palette = "nature",
      theme_style = "classic",
      mode = "light",
      point_size = 2.2,
      point_alpha = 0.85,
      line_width = 0.8,
      width = 6.5,
      height = 5,
      dpi = 300,
      export_pdf = TRUE,
      rank_by = "consensus",
      significant_only = TRUE,
      log = log
    )

    n_gene_plots <- if (
      !is.null(plots$gene_curves)
    ) {
      length(plots$gene_curves)
    } else {
      0L
    }

    log(
      paste0(
        "Gene-curve plots generated: ",
        n_gene_plots
      ),
      type = "info"
    )

  } else if (
    isTRUE(generate_plots) &&
    nrow(merged) == 0L
  ) {

    log(
      "Plots were skipped because the merged result table is empty.",
      type = "warn"
    )
  }

  # ======================================================= #
  # 12. Analysis statistics
  # ======================================================= #

  n_meta <- if (
    nrow(merged) > 0L &&
    "rhythmic_metacycle" %in%
    colnames(merged)
  ) {
    sum(
      merged$rhythmic_metacycle %in% TRUE
    )
  } else {
    NA_integer_
  }

  n_cosinor <- if (
    nrow(merged) > 0L &&
    "rhythmic_cosinor" %in%
    colnames(merged)
  ) {
    sum(
      merged$rhythmic_cosinor %in% TRUE
    )
  } else {
    NA_integer_
  }

  n_both <- if (
    nrow(merged) > 0L &&
    "rhythmic_by" %in%
    colnames(merged)
  ) {
    sum(
      merged$rhythmic_by == "both",
      na.rm = TRUE
    )
  } else {
    NA_integer_
  }

  n_diff <- if (
    nrow(merged) > 0L &&
    "diff_rhythm_padj" %in%
    colnames(merged)
  ) {
    sum(
      merged$diff_rhythm_padj <
        padj_cutoff,
      na.rm = TRUE
    )
  } else {
    NA_integer_
  }

  duration <- as.numeric(
    difftime(
      Sys.time(),
      t_start,
      units = "secs"
    )
  )

  log(
    paste0(
      "Rhythmic genes - MetaCycle: ",
      ifelse(
        is.na(n_meta),
        "not run",
        n_meta
      ),
      "; cosinor: ",
      ifelse(
        is.na(n_cosinor),
        "not run",
        n_cosinor
      ),
      "; consensus: ",
      ifelse(
        is.na(n_both),
        "not available",
        n_both
      ),
      "."
    ),
    type = "info"
  )

  if (!is.na(n_diff)) {
    log(
      paste0(
        "Differentially rhythmic genes: ",
        n_diff,
        " at FDR < ",
        padj_cutoff,
        "."
      ),
      type = "info"
    )
  }

  log(
    "DEGgo rhythmicity analysis complete.",
    type = "done",
    duration = duration
  )

  # ======================================================= #
  # 13. Return
  # ======================================================= #

  structure(
    list(
      summary = merged,
      metacycle = meta_res,
      cosinor = cos_res,
      cosinor_diff_group = diff_res,
      metadata = md,
      expr = mat,
      gene_annotation = gene_annotation,
      stats = list(
        n_genes_tested = nrow(mat),
        n_samples = ncol(mat),
        n_timepoints = length(time_values),
        rhythmic_metacycle = n_meta,
        rhythmic_cosinor = n_cosinor,
        rhythmic_consensus = n_both,
        differential_rhythmic = n_diff
      ),
      params = list(
        methods = methods,
        period_range = period_range,
        cycle_length = cycle_length,
        cycMethod = cycMethod,
        padj_cutoff = padj_cutoff,
        cosinor_engine = if (
          "cosinor" %in% methods
        ) {
          if (isTRUE(use_cosinor_pkg)) {
            "package"
          } else {
            "manual"
          }
        } else {
          NA_character_
        },
        sample_col = sample_col,
        time_col = time_col,
        group_col = group_col,
        assay = assay,
        gene_id_col = gene_id_col,
        gene_symbol_col = detected_symbol_col,
        show_gene_id = show_gene_id,
        project_name = project_name,
        seed = seed
      ),
      output_dir = output_dir,
      output_dirs = dirs,
      plots = plots,
      version = "0.1.0"
    ),
    class = "deggo_rhythm_results"
  )
}



# ======================================================= #== #
# RUN_DEGGO INTEGRATION: OPTIONAL RHYTHMICITY STEP
# ======================================================= #== #
#' Run rhythmicity analysis as part of the run_deggo() pipeline
#'
#' Internal helper called from [run_deggo()] when `rhythmicity_analysis =
#' TRUE`. Reuses the dds object already fitted by the DE pipeline and the
#' project's own metadata, and stores results at `de_results$rhythmicity`.
#' Failures are caught and logged as warnings so a rhythmicity error never
#' aborts the main DE/GO workflow.
#'
#' @keywords internal
#' @noRd
.deggo_make_rhythmicity <- function(
    de_results,
    output_dir,
    analysis_mode,
    time_col,
    group_col = NULL,
    assay = "vst",
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
    padj_cutoff = 0.05,
    cosinor_engine = "auto",
    generate_plots = TRUE,
    n_top_plots = 20,
    txtsize = 12,
    gene_annotation = NULL,
    gene_id_col = "gene_id",
    gene_symbol_col = NULL,
    show_gene_id = TRUE,
    project_name = "DEGgo rhythmicity analysis",
    log
) {

  metadata <- de_results$metadata
  dds <- de_results$dds

  if (is.null(dds)) {
    log(
      "Rhythmicity analysis skipped: no dds object available.",
      type = "warn"
    )

    return(de_results)
  }

  if (
    is.null(metadata) ||
    !time_col %in% colnames(metadata)
  ) {
    log(
      paste0(
        "Rhythmicity analysis skipped: metadata column '",
        time_col,
        "' was not found."
      ),
      type = "warn"
    )

    return(de_results)
  }

  if (!is.numeric(metadata[[time_col]])) {
    log(
      paste0(
        "Rhythmicity analysis skipped: metadata column '",
        time_col,
        "' is not numeric."
      ),
      type = "warn"
    )

    return(de_results)
  }

  log(
    "[RHYTHM] Running rhythmicity analysis",
    type = "step"
  )

  rhythm_dir <- file.path(
    output_dir,
    paste0(
      analysis_mode,
      "_rhythmicity"
    )
  )

  rhythm_results <- tryCatch(
    run_deggo_rhythmicity(
      expr = dds,
      metadata = metadata,
      sample_col = if (
        "sample" %in% colnames(metadata)
      ) {
        "sample"
      } else {
        colnames(metadata)[1L]
      },
      time_col = time_col,
      group_col = group_col,
      assay = assay,
      methods = methods,
      period_range = period_range,
      cycle_length = cycle_length,
      cycMethod = cycMethod,
      padj_cutoff = padj_cutoff,
      cosinor_engine = cosinor_engine,
      output_dir = rhythm_dir,
      project_name = project_name,
      generate_plots = generate_plots,
      n_top_plots = n_top_plots,
      txtsize = txtsize,
      gene_annotation = gene_annotation,
      gene_id_col = gene_id_col,
      gene_symbol_col = gene_symbol_col,
      show_gene_id = show_gene_id,
      verbose = TRUE
    ),
    error = function(e) {

      log(
        paste0(
          "Rhythmicity analysis failed: ",
          conditionMessage(e)
        ),
        type = "error"
      )

      NULL
    }
  )

  de_results$rhythmicity <- rhythm_results

  de_results
}
