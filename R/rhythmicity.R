# ============================================================ #
# DEGgo rhythmicity analysis (MetaCycle + cosinor)
# ============================================================ #
#
# Standalone module for detecting rhythmic (e.g. circadian) genes in
# time-course expression data using MetaCycle (meta2d) and single-component
# cosinor regression (cosinor2), plus an optional cosinor-based test for
# differential rhythmicity between two groups.
#
# This module is intentionally independent from run_deggo(): it can be run
# on any expression matrix / DESeq2 dds object with a numeric time column in
# metadata.


# ========================================================= #
# CREATE RHYTHMICITY DIRECTORIES
# ========================================================= #
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


# ========================================================= #
# PREPARE RHYTHMICITY INPUT
# ========================================================= #
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
    assay = "vst"
) {

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (!sample_col %in% colnames(metadata)) {
    stop("metadata must contain a '", sample_col, "' column.", call. = FALSE)
  }

  if (!time_col %in% colnames(metadata)) {
    stop("metadata must contain a '", time_col, "' column.", call. = FALSE)
  }

  if (!is.numeric(metadata[[time_col]])) {
    stop("metadata[['", time_col, "']] must be numeric (time, e.g. hours).", call. = FALSE)
  }

  if (!is.null(group_col) && !group_col %in% colnames(metadata)) {
    stop("metadata must contain a '", group_col, "' column.", call. = FALSE)
  }

  is_dds <- inherits(expr, "DESeqDataSet")

  if (is_dds) {

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("Package 'DESeq2' is required to extract expression from a dds object.", call. = FALSE)
    }

    mat <- switch(
      assay,
      raw = DESeq2::counts(expr, normalized = FALSE),
      normalized = DESeq2::counts(expr, normalized = TRUE),
      log2_normalized = log2(DESeq2::counts(expr, normalized = TRUE) + 1),
      vst = SummarizedExperiment::assay(DESeq2::vst(expr, blind = FALSE)),
      stop("Unsupported assay: ", assay, call. = FALSE)
    )

  } else {
    mat <- as.matrix(expr)
  }

  if (is.null(colnames(mat))) {
    stop("Expression input must have column names matching sample IDs.", call. = FALSE)
  }

  common <- intersect(colnames(mat), metadata[[sample_col]])

  if (length(common) < 4) {
    stop(
      "At least 4 matched samples with time points are required ",
      "for rhythmicity analysis (found ", length(common), ").",
      call. = FALSE
    )
  }

  metadata <- metadata[metadata[[sample_col]] %in% common, , drop = FALSE]
  metadata <- metadata[order(metadata[[time_col]]), , drop = FALSE]
  rownames(metadata) <- NULL

  mat <- mat[, metadata[[sample_col]], drop = FALSE]
  storage.mode(mat) <- "double"

  list(mat = mat, metadata = metadata)
}


# ========================================================= #
# RUN METACYCLE
# ========================================================= #
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


# ========================================================= #
# MANUAL LM()-BASED COSINOR FALLBACK
# ========================================================= #
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


# ========================================================= #
# RUN SINGLE-GENE COSINOR
# ========================================================= #
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


# ========================================================= #
# RUN COSINOR DIFFERENTIAL RHYTHMICITY (GROUP INTERACTION)
# ========================================================= #
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


# ========================================================= #
# MERGE RHYTHMICITY RESULTS
# ========================================================= #
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




# ========================================================= #
# RHYTHMICITY DIAGNOSTIC PLOTS
# ========================================================= #

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

  labels <- stats::setNames(gene_ids, gene_ids)

  symbol_candidates <- c(
    "Gene",
    "gene",
    "feature",
    "symbol",
    "SYMBOL",
    "gene_symbol",
    "gene_name",
    "external_gene_name"
  )

  # -------------------------------------------------------
  # 1. Prefer symbol columns already present in summary
  # -------------------------------------------------------

  if (!is.null(summary) && nrow(summary)) {

    summary_symbol_cols <- intersect(
      symbol_candidates,
      setdiff(colnames(summary), "gene")
    )

    if (length(summary_symbol_cols)) {

      symbol_col <- summary_symbol_cols[1]

      idx <- match(gene_ids, as.character(summary$gene))
      values <- as.character(summary[[symbol_col]][idx])

      valid <- !is.na(values) &
        nzchar(trimws(values)) &
        values != "NA"

      labels[valid] <- values[valid]
    }
  }

  # -------------------------------------------------------
  # 2. Use external annotation table
  # -------------------------------------------------------

  if (!is.null(gene_annotation)) {

    gene_annotation <- as.data.frame(
      gene_annotation,
      stringsAsFactors = FALSE
    )

    if (!gene_id_col %in% colnames(gene_annotation)) {
      stop(
        "gene_annotation must contain the gene ID column '",
        gene_id_col,
        "'.",
        call. = FALSE
      )
    }

    if (is.null(gene_symbol_col)) {

      available <- intersect(
        symbol_candidates,
        setdiff(colnames(gene_annotation), gene_id_col)
      )

      if (length(available)) {
        gene_symbol_col <- available[1]
      }
    }

    if (
      !is.null(gene_symbol_col) &&
      gene_symbol_col %in% colnames(gene_annotation)
    ) {

      annotation <- gene_annotation[
        !duplicated(gene_annotation[[gene_id_col]]),
        c(gene_id_col, gene_symbol_col),
        drop = FALSE
      ]

      idx <- match(
        gene_ids,
        as.character(annotation[[gene_id_col]])
      )

      values <- as.character(annotation[[gene_symbol_col]][idx])

      valid <- !is.na(values) &
        nzchar(trimws(values)) &
        values != "NA"

      labels[valid] <- values[valid]
    }
  }

  # Clean labels
  labels <- trimws(labels)

  invalid <- is.na(labels) | !nzchar(labels) | labels == "NA"
  labels[invalid] <- names(labels)[invalid]

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


#' Generate rhythmicity diagnostic plots
#'
#' @param mat Numeric genes-by-samples expression matrix.
#' @param metadata Sample metadata matched to the columns of `mat`.
#' @param summary Combined rhythmicity results table.
#' @param time_col Numeric time column in `metadata`.
#' @param group_col Optional grouping column.
#' @param cycle_length Assumed rhythm period in hours.
#' @param output_dir Plot output directory.
#' @param n_top Number of top-ranked genes to plot.
#' @param txtsize Base text size.
#' @param gene_annotation Optional gene annotation table.
#' @param gene_id_col Gene-ID column in `gene_annotation`.
#' @param gene_symbol_col Optional symbol column in `gene_annotation`.
#' @param show_gene_id Logical. Include the Ensembl ID below the gene symbol.
#' @param fit_full_cycle Logical. Draw fitted curves across a complete cycle.
#'
#' @return A named list of generated ggplot objects.
#'
#' @keywords internal
#' @noRd
.deggo_rhythm_plots <- function(
    mat,
    metadata,
    summary,
    time_col,
    group_col = NULL,
    cycle_length = 24,
    output_dir,
    n_top = 20,
    txtsize = 12,
    gene_annotation = NULL,
    gene_id_col = "gene_id",
    gene_symbol_col = NULL,
    show_gene_id = TRUE,
    fit_full_cycle = TRUE
) {

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  plots <- list()

  if (is.null(summary) || !nrow(summary)) {
    return(plots)
  }

  mat <- as.matrix(mat)
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (is.null(rownames(mat))) {
    stop("mat must have gene identifiers as row names.", call. = FALSE)
  }

  if (!"gene" %in% colnames(summary)) {
    stop("summary must contain a 'gene' column.", call. = FALSE)
  }

  if (!time_col %in% colnames(metadata)) {
    stop(
      "metadata does not contain time column '",
      time_col,
      "'.",
      call. = FALSE
    )
  }

  if (!is.numeric(metadata[[time_col]])) {
    stop(
      "metadata[['",
      time_col,
      "']] must be numeric.",
      call. = FALSE
    )
  }

  if (
    !is.null(group_col) &&
    !group_col %in% colnames(metadata)
  ) {
    stop(
      "metadata does not contain group column '",
      group_col,
      "'.",
      call. = FALSE
    )
  }

  # =======================================================
  # Safe -log10 transformation
  # =======================================================

  safe_neglog10 <- function(x) {
    x <- as.numeric(x)
    x[x <= 0] <- .Machine$double.xmin
    -log10(x)
  }

  # =======================================================
  # 1. MetaCycle period distribution
  # =======================================================

  if ("metacycle_period" %in% colnames(summary)) {

    period_data <- summary[
      is.finite(summary$metacycle_period),
      ,
      drop = FALSE
    ]

    if (nrow(period_data)) {

      p <- ggplot2::ggplot(
        period_data,
        ggplot2::aes(x = .data[["metacycle_period"]])
      ) +
        ggplot2::geom_histogram(
          binwidth = 1,
          boundary = 0,
          na.rm = TRUE
        ) +
        ggplot2::labs(
          title = "MetaCycle period distribution",
          x = "Estimated period (hours)",
          y = "Number of genes"
        ) +
        .deggo_theme(txtsize = txtsize)

      ggplot2::ggsave(
        filename = file.path(
          output_dir,
          "period_histogram.png"
        ),
        plot = p,
        width = 6,
        height = 4,
        dpi = 300,
        bg = "white"
      )

      plots$period_histogram <- p
    }
  }

  # =======================================================
  # 2. Phase distribution
  # =======================================================

  phase_col <- intersect(
    c("metacycle_phase", "cosinor_phase"),
    colnames(summary)
  )

  if (length(phase_col)) {

    phase_col <- phase_col[1]

    phase_data <- summary[
      is.finite(summary[[phase_col]]),
      ,
      drop = FALSE
    ]

    if (nrow(phase_data)) {

      p <- ggplot2::ggplot(
        phase_data,
        ggplot2::aes(x = .data[[phase_col]])
      ) +
        ggplot2::geom_histogram(
          binwidth = 2,
          boundary = 0,
          na.rm = TRUE
        ) +
        ggplot2::scale_x_continuous(
          limits = c(0, cycle_length),
          breaks = seq(0, cycle_length, by = 4)
        ) +
        ggplot2::labs(
          title = "Rhythmic phase distribution",
          subtitle = paste("Phase source:", phase_col),
          x = "Peak phase (hours)",
          y = "Number of genes"
        ) +
        .deggo_theme(txtsize = txtsize)

      ggplot2::ggsave(
        filename = file.path(
          output_dir,
          "phase_histogram.png"
        ),
        plot = p,
        width = 6,
        height = 4,
        dpi = 300,
        bg = "white"
      )

      plots$phase_histogram <- p
    }
  }

  # =======================================================
  # 3. MetaCycle versus cosinor agreement
  # =======================================================

  if (
    all(
      c("metacycle_padj", "cosinor_padj") %in%
      colnames(summary)
    )
  ) {

    agreement_df <- summary

    agreement_df$metacycle_score <- safe_neglog10(
      agreement_df$metacycle_padj
    )

    agreement_df$cosinor_score <- safe_neglog10(
      agreement_df$cosinor_padj
    )

    agreement_df <- agreement_df[
      is.finite(agreement_df$metacycle_score) &
        is.finite(agreement_df$cosinor_score),
      ,
      drop = FALSE
    ]

    if (nrow(agreement_df)) {

      threshold <- -log10(0.05)

      p <- ggplot2::ggplot(
        agreement_df,
        ggplot2::aes(
          x = .data[["metacycle_score"]],
          y = .data[["cosinor_score"]]
        )
      ) +
        ggplot2::geom_point(
          alpha = 0.5,
          na.rm = TRUE
        ) +
        ggplot2::geom_hline(
          yintercept = threshold,
          linetype = "dashed"
        ) +
        ggplot2::geom_vline(
          xintercept = threshold,
          linetype = "dashed"
        ) +
        ggplot2::labs(
          title = "MetaCycle versus cosinor agreement",
          x = expression(-log[10]("MetaCycle adjusted p-value")),
          y = expression(-log[10]("Cosinor adjusted p-value"))
        ) +
        .deggo_theme(txtsize = txtsize)

      ggplot2::ggsave(
        filename = file.path(
          output_dir,
          "metacycle_vs_cosinor.png"
        ),
        plot = p,
        width = 6,
        height = 6,
        dpi = 300,
        bg = "white"
      )

      plots$agreement <- p
    }
  }

  # =======================================================
  # 4. Select top rhythmic genes
  # =======================================================

  rank_candidates <- intersect(
    c(
      "metacycle_padj",
      "cosinor_padj",
      "metacycle_pvalue",
      "cosinor_pvalue"
    ),
    colnames(summary)
  )

  ranked_summary <- summary

  if (length(rank_candidates)) {

    rank_col <- rank_candidates[1]

    ranked_summary <- ranked_summary[
      order(
        ranked_summary[[rank_col]],
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
  }

  top_genes <- unique(
    utils::head(
      as.character(ranked_summary$gene),
      n_top
    )
  )

  top_genes <- intersect(top_genes, rownames(mat))

  if (!length(top_genes)) {
    return(plots)
  }

  # Resolve symbol labels
  gene_labels <- .deggo_rhythm_gene_labels(
    gene_ids = top_genes,
    summary = summary,
    gene_annotation = gene_annotation,
    gene_id_col = gene_id_col,
    gene_symbol_col = gene_symbol_col
  )

  # =======================================================
  # 5. Expression and fitted-curve plots
  # =======================================================

  time_vec <- metadata[[time_col]]

  if (isTRUE(fit_full_cycle)) {
    time_grid <- seq(
      0,
      cycle_length,
      length.out = 300
    )
  } else {
    time_grid <- seq(
      min(time_vec, na.rm = TRUE),
      max(time_vec, na.rm = TRUE),
      length.out = 300
    )
  }

  curve_plots <- list()

  for (g in top_genes) {

    symbol <- unname(gene_labels[g])

    if (is.na(symbol) || !nzchar(symbol)) {
      symbol <- g
    }

    df <- data.frame(
      sample = colnames(mat),
      time = time_vec,
      expression = as.numeric(mat[g, ]),
      stringsAsFactors = FALSE
    )

    if (!is.null(group_col)) {
      df$group <- factor(metadata[[group_col]])
    }

    result_row <- summary[
      as.character(summary$gene) == g,
      ,
      drop = FALSE
    ]

    fit_df <- NULL

    required_fit_cols <- c(
      "cosinor_mesor",
      "cosinor_amplitude",
      "cosinor_acrophase"
    )

    can_fit <- nrow(result_row) &&
      all(required_fit_cols %in% colnames(result_row)) &&
      all(
        is.finite(
          as.numeric(
            result_row[1, required_fit_cols]
          )
        )
      )

    if (can_fit) {

      mesor <- as.numeric(result_row$cosinor_mesor[1])
      amplitude <- as.numeric(result_row$cosinor_amplitude[1])
      acrophase <- as.numeric(result_row$cosinor_acrophase[1])

      fit_df <- data.frame(
        time = time_grid,
        fitted = mesor +
          amplitude *
          cos(
            2 * pi * time_grid / cycle_length -
              acrophase
          )
      )
    }

    # -----------------------------------------------------
    # Plot title and statistical subtitle
    # -----------------------------------------------------

    plot_title <- symbol

    if (
      isTRUE(show_gene_id) &&
      !identical(symbol, g)
    ) {
      plot_title <- paste0(symbol, " (", g, ")")
    }

    subtitle_parts <- character()

    if (
      nrow(result_row) &&
      "cosinor_padj" %in% colnames(result_row)
    ) {
      subtitle_parts <- c(
        subtitle_parts,
        paste0(
          "Cosinor FDR = ",
          .deggo_rhythm_format_p(
            result_row$cosinor_padj[1]
          )
        )
      )
    }

    if (
      nrow(result_row) &&
      "metacycle_padj" %in% colnames(result_row)
    ) {
      subtitle_parts <- c(
        subtitle_parts,
        paste0(
          "MetaCycle FDR = ",
          .deggo_rhythm_format_p(
            result_row$metacycle_padj[1]
          )
        )
      )
    }

    if (
      nrow(result_row) &&
      "cosinor_phase" %in% colnames(result_row) &&
      is.finite(result_row$cosinor_phase[1])
    ) {
      subtitle_parts <- c(
        subtitle_parts,
        paste0(
          "Peak = ZT",
          formatC(
            result_row$cosinor_phase[1],
            format = "f",
            digits = 1
          )
        )
      )
    }

    plot_subtitle <- if (length(subtitle_parts)) {
      paste(subtitle_parts, collapse = " | ")
    } else {
      NULL
    }

    # -----------------------------------------------------
    # Build plot
    # -----------------------------------------------------

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[["time"]],
        y = .data[["expression"]]
      )
    )

    if (!is.null(group_col)) {
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(
            color = .data[["group"]]
          ),
          size = 2
        )
    } else {
      p <- p +
        ggplot2::geom_point(size = 2)
    }

    if (!is.null(fit_df)) {
      p <- p +
        ggplot2::geom_line(
          data = fit_df,
          mapping = ggplot2::aes(
            x = .data[["time"]],
            y = .data[["fitted"]]
          ),
          inherit.aes = FALSE,
          linewidth = 0.7
        )
    }

    p <- p +
      ggplot2::scale_x_continuous(
        limits = c(0, cycle_length),
        breaks = seq(0, cycle_length, by = 4),
        minor_breaks = NULL
      ) +
      ggplot2::labs(
        title = plot_title,
        subtitle = plot_subtitle,
        x = "Zeitgeber time (hours)",
        y = "VST expression",
        color = if (!is.null(group_col)) group_col else NULL
      ) +
      .deggo_theme(
        txtsize = txtsize,
        x.ang = 0
      )

    safe_symbol <- make.names(symbol)

    filename <- paste0(
      "gene_",
      safe_symbol,
      "_",
      make.names(g),
      "_fit.png"
    )

    ggplot2::ggsave(
      filename = file.path(output_dir, filename),
      plot = p,
      width = 6.5,
      height = 4.5,
      dpi = 300,
      bg = "white"
    )

    curve_plots[[symbol]] <- p
  }

  plots$gene_curves <- curve_plots
  plots
}



# ========================================================= #
# RUN DEGGO RHYTHMICITY (main entry point)
# ========================================================= #

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
    assay = c("vst", "normalized", "log2_normalized", "raw"),
    methods = c("meta2d", "cosinor"),
    period_range = c(20, 28),
    cycle_length = 24,
    cycMethod = c("ARS", "JTK", "LS"),
    padj_cutoff = 0.05,
    cosinor_engine = c("auto", "package", "manual"),
    gene_annotation = NULL,
    gene_id_col = "gene_id",
    gene_symbol_col = NULL,
    show_gene_id = FALSE,
    output_dir = "DEGgo_rhythmicity_out",
    project_name = NULL,
    generate_plots = TRUE,
    n_top_plots = 20,
    seed = 4173,
    verbose = TRUE
) {

  # ------------------------------------------------------- #
  # Argument validation
  # ------------------------------------------------------- #

  set.seed(seed)

  assay <- match.arg(assay)

  methods <- match.arg(
    methods,
    choices = c("meta2d", "cosinor"),
    several.ok = TRUE
  )

  methods <- unique(methods)

  cosinor_engine <- match.arg(cosinor_engine)

  if (
    !is.numeric(period_range) ||
    length(period_range) != 2L ||
    anyNA(period_range) ||
    any(!is.finite(period_range)) ||
    period_range[1] <= 0 ||
    period_range[2] <= period_range[1]
  ) {
    stop(
      "'period_range' must be a numeric vector of length 2 with ",
      "0 < period_range[1] < period_range[2].",
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
      "'cycle_length' must be one positive finite numeric value.",
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
      "'padj_cutoff' must be a single numeric value between 0 and 1.",
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
      "'n_top_plots' must be a non-negative numeric value.",
      call. = FALSE
    )
  }

  n_top_plots <- as.integer(n_top_plots)

  allowed_cyc_methods <- c("ARS", "JTK", "LS")

  if (
    !is.character(cycMethod) ||
    !length(cycMethod) ||
    any(!cycMethod %in% allowed_cyc_methods)
  ) {
    stop(
      "'cycMethod' must contain one or more of: ",
      paste(allowed_cyc_methods, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  cycMethod <- unique(cycMethod)

  log <- .deggo_msg(
    verbose = verbose,
    prefix = "DEGgo"
  )

  log(
    "Starting DEGgo rhythmicity analysis",
    type = "header"
  )

  t_start <- Sys.time()

  # ------------------------------------------------------- #
  # Validate optional annotation
  # ------------------------------------------------------- #

  if (!is.null(gene_annotation)) {

    gene_annotation <- as.data.frame(
      gene_annotation,
      stringsAsFactors = FALSE
    )

    if (!gene_id_col %in% colnames(gene_annotation)) {
      stop(
        "gene_annotation must contain the gene ID column '",
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
        nzchar(trimws(gene_annotation[[gene_id_col]])),
      ,
      drop = FALSE
    ]

    gene_annotation <- gene_annotation[
      !duplicated(gene_annotation[[gene_id_col]]),
      ,
      drop = FALSE
    ]

    if (!is.null(gene_symbol_col)) {

      if (!gene_symbol_col %in% colnames(gene_annotation)) {
        stop(
          "gene_annotation does not contain gene-symbol column '",
          gene_symbol_col,
          "'.",
          call. = FALSE
        )
      }

      gene_annotation[[gene_symbol_col]] <- as.character(
        gene_annotation[[gene_symbol_col]]
      )
    }
  }

  # ------------------------------------------------------- #
  # Dependency checks
  # ------------------------------------------------------- #

  if (
    "meta2d" %in% methods &&
    !requireNamespace("MetaCycle", quietly = TRUE)
  ) {
    stop(
      "Package 'MetaCycle' is required for meta2d-based rhythmicity ",
      "analysis. Install it with install.packages('MetaCycle').",
      call. = FALSE
    )
  }

  use_cosinor_pkg <- FALSE

  if ("cosinor" %in% methods) {

    pkg_ok <- requireNamespace("cosinor", quietly = TRUE) &&
      requireNamespace("cosinor2", quietly = TRUE)

    use_cosinor_pkg <- switch(
      cosinor_engine,

      package = {
        if (!pkg_ok) {
          stop(
            "cosinor_engine = 'package' requires both 'cosinor' and ",
            "'cosinor2'. Install them or use cosinor_engine = 'manual'.",
            call. = FALSE
          )
        }

        TRUE
      },

      manual = FALSE,

      auto = pkg_ok
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

  # ------------------------------------------------------- #
  # Prepare output directories and inputs
  # ------------------------------------------------------- #

  dirs <- .deggo_rhythm_dirs(output_dir)

  prepped <- .prepare_rhythm_input(
    expr = expr,
    metadata = metadata,
    sample_col = sample_col,
    time_col = time_col,
    group_col = group_col,
    assay = assay
  )

  mat <- prepped$mat
  md <- prepped$metadata

  if (is.null(rownames(mat))) {
    stop(
      "The prepared expression matrix must contain gene identifiers as row names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(rownames(mat))) {
    stop(
      "The prepared expression matrix contains duplicated gene identifiers.",
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

  time_values <- sort(unique(md[[time_col]]))

  log(
    paste0(
      "Time points: ",
      paste(time_values, collapse = ", ")
    ),
    type = "info"
  )

  # ------------------------------------------------------- #
  # Run MetaCycle
  # ------------------------------------------------------- #

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
      minper = period_range[1],
      maxper = period_range[2],
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

  # ------------------------------------------------------- #
  # Run cosinor analyses
  # ------------------------------------------------------- #

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

      group_values <- as.character(md[[group_col]])
      group_values <- group_values[
        !is.na(group_values) &
          nzchar(trimws(group_values))
      ]

      grp_levels <- unique(group_values)

      if (length(grp_levels) == 2L) {

        log(
          paste0(
            "Testing differential rhythmicity: ",
            paste(grp_levels, collapse = " vs ")
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
            "Skipping differential rhythmicity test: '",
            group_col,
            "' has ",
            length(grp_levels),
            " non-missing levels; exactly 2 are required."
          ),
          type = "warn"
        )
      }
    }
  }

  # ------------------------------------------------------- #
  # Merge results
  # ------------------------------------------------------- #

  merged <- .merge_rhythm_results(
    meta_res = meta_res,
    cos_res = cos_res,
    diff_res = diff_res,
    padj_cutoff = padj_cutoff
  )

  # ------------------------------------------------------- #
  # Add gene symbols to result tables
  # ------------------------------------------------------- #

  detected_symbol_col <- gene_symbol_col

  if (
    !is.null(gene_annotation) &&
    is.null(detected_symbol_col)
  ) {

    symbol_candidates <- c(
      "Gene",
      "feature",
      "symbol",
      "SYMBOL",
      "gene_symbol",
      "gene_name",
      "external_gene_name"
    )

    available_symbol_cols <- intersect(
      symbol_candidates,
      setdiff(colnames(gene_annotation), gene_id_col)
    )

    if (length(available_symbol_cols)) {
      detected_symbol_col <- available_symbol_cols[1]
    }
  }

  add_annotation <- function(tab) {

    if (
      is.null(tab) ||
      !nrow(tab) ||
      is.null(gene_annotation) ||
      is.null(detected_symbol_col)
    ) {
      return(tab)
    }

    annotation_map <- gene_annotation[
      ,
      c(gene_id_col, detected_symbol_col),
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
      !duplicated(annotation_map$gene),
      ,
      drop = FALSE
    ]

    out <- merge(
      annotation_map,
      tab,
      by = "gene",
      all.y = TRUE,
      sort = FALSE
    )

    if ("gene" %in% colnames(tab)) {
      original_order <- match(
        as.character(tab$gene),
        as.character(out$gene)
      )

      if (all(!is.na(original_order))) {
        out <- out[original_order, , drop = FALSE]
      }
    }

    rownames(out) <- NULL
    out
  }

  meta_res <- add_annotation(meta_res)
  cos_res <- add_annotation(cos_res)
  diff_res <- add_annotation(diff_res)
  merged <- add_annotation(merged)

  # ------------------------------------------------------- #
  # Write result tables
  # ------------------------------------------------------- #

  if (!is.null(meta_res) && nrow(meta_res)) {
    utils::write.table(
      meta_res,
      file.path(
        dirs$results,
        "metacycle_results.tsv"
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  if (!is.null(cos_res) && nrow(cos_res)) {
    utils::write.table(
      cos_res,
      file.path(
        dirs$results,
        "cosinor_results.tsv"
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  if (!is.null(diff_res) && nrow(diff_res)) {
    utils::write.table(
      diff_res,
      file.path(
        dirs$results,
        "cosinor_differential_rhythmicity.tsv"
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  if (nrow(merged)) {
    utils::write.table(
      merged,
      file.path(
        dirs$results,
        "rhythmicity_summary.tsv"
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  # ------------------------------------------------------- #
  # Generate plots
  # ------------------------------------------------------- #

  plots <- NULL

  if (
    isTRUE(generate_plots) &&
    nrow(merged) &&
    n_top_plots > 0L
  ) {

    log(
      "Generating rhythmicity diagnostic plots...",
      type = "step"
    )

    plots <- tryCatch(
      .deggo_rhythm_plots(
        mat = mat,
        metadata = md,
        summary = merged,
        time_col = time_col,
        group_col = group_col,
        cycle_length = cycle_length,
        output_dir = dirs$plots,
        n_top = n_top_plots,
        gene_annotation = gene_annotation,
        gene_id_col = gene_id_col,
        gene_symbol_col = detected_symbol_col,
        show_gene_id = show_gene_id
      ),
      error = function(e) {
        log(
          paste(
            "Plot generation failed:",
            conditionMessage(e)
          ),
          type = "warn"
        )

        NULL
      }
    )
  }

  # ------------------------------------------------------- #
  # Analysis summary
  # ------------------------------------------------------- #

  n_meta <- if (
    nrow(merged) &&
    "rhythmic_metacycle" %in% colnames(merged)
  ) {
    sum(merged$rhythmic_metacycle %in% TRUE)
  } else {
    NA_integer_
  }

  n_cosinor <- if (
    nrow(merged) &&
    "rhythmic_cosinor" %in% colnames(merged)
  ) {
    sum(merged$rhythmic_cosinor %in% TRUE)
  } else {
    NA_integer_
  }

  n_both <- if (
    nrow(merged) &&
    "rhythmic_by" %in% colnames(merged)
  ) {
    sum(merged$rhythmic_by == "both", na.rm = TRUE)
  } else {
    NA_integer_
  }

  n_diff <- if (
    nrow(merged) &&
    "diff_rhythm_padj" %in% colnames(merged)
  ) {
    sum(
      merged$diff_rhythm_padj < padj_cutoff,
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
      ifelse(is.na(n_meta), "not run", n_meta),
      "; cosinor: ",
      ifelse(is.na(n_cosinor), "not run", n_cosinor),
      "; consensus: ",
      ifelse(is.na(n_both), "not available", n_both),
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

  # ------------------------------------------------------- #
  # Return object
  # ------------------------------------------------------- #

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
        cosinor_engine = if ("cosinor" %in% methods) {
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



# ========================================================= #
# RUN_DEGGO INTEGRATION: OPTIONAL RHYTHMICITY STEP
# ========================================================= #
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
    group_col,
    assay,
    methods,
    period_range,
    cycle_length,
    cycMethod,
    padj_cutoff,
    cosinor_engine,
    generate_plots,
    n_top_plots,
    project_name,
    log
) {

  metadata <- de_results$metadata
  dds <- de_results$dds

  if (is.null(dds)) {
    log("Rhythmicity analysis skipped: no dds object available.", type = "warn")
    return(de_results)
  }

  if (is.null(metadata) || !time_col %in% colnames(metadata)) {
    log(
      paste0(
        "Rhythmicity analysis skipped: metadata column '", time_col,
        "' not found. Set rhythmicity_time_col to a numeric time/ZT column ",
        "to enable this step."
      ),
      type = "warn"
    )
    return(de_results)
  }

  if (!is.numeric(metadata[[time_col]])) {
    log(
      paste0(
        "Rhythmicity analysis skipped: metadata column '", time_col,
        "' is not numeric."
      ),
      type = "warn"
    )
    return(de_results)
  }

  log("[RHYTHM] Running rhythmicity analysis (MetaCycle + cosinor)", type = "step")

  rhythm_dir <- file.path(output_dir, paste0(analysis_mode, "_rhythmicity"))

  rhythm_results <- tryCatch(
    run_deggo_rhythmicity(
      expr = dds,
      metadata = metadata,
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
      verbose = TRUE
    ),
    error = function(e) {
      log(paste("Rhythmicity analysis failed:", conditionMessage(e)), type = "error")
      NULL
    }
  )

  de_results$rhythmicity <- rhythm_results
  de_results
}
