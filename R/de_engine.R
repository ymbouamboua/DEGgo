# =========================================================
# RUN DIFFERENTIAL EXPRESSION
# =========================================================

#' Run differential expression analysis
#'
#' Runs differential expression analysis using DESeq2, edgeR, or limma-voom.
#'
#' @param counts Count matrix with genes as rows and samples as columns.
#' @param metadata Sample metadata with samples as rows.
#' @param method Differential expression method.
#'   One of \code{"DESeq2"}, \code{"edgeR"}, or \code{"limma"}.
#' @param design_formula Design formula used for modeling.
#' @param contrast Optional contrast. For DESeq2, use a character vector such
#'   as \code{c("condition", "treated", "control")}. For edgeR/limma, use
#'   either a coefficient name/index or a limma-style contrast string.
#'
#' @return A list containing the model object, result table, design matrix,
#'   method, design formula, and contrast.
#'
#' @keywords internal
#' @noRd
run_de <- function(
    counts,
    metadata,
    method = c("DESeq2", "edgeR", "limma"),
    design_formula = ~ condition,
    contrast = NULL
) {

  method <- match.arg(method)

  log <- .deggo_msg(verbose = TRUE, prefix = "DEGgo")
  t0 <- Sys.time()

  counts <- round(as.matrix(counts))
  storage.mode(counts) <- "integer"

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (is.null(rownames(metadata))) {
    stop("metadata must have sample names as rownames.", call. = FALSE)
  }

  if (!identical(colnames(counts), rownames(metadata))) {
    stop(
      "Sample order mismatch between counts and metadata. ",
      "Run prepare_counts_metadata() first.",
      call. = FALSE
    )
  }

  design_vars <- all.vars(design_formula)
  missing_design_vars <- setdiff(design_vars, colnames(metadata))

  if (length(missing_design_vars) > 0) {
    stop(
      "Variables missing from metadata: ",
      paste(missing_design_vars, collapse = ", "),
      call. = FALSE
    )
  }

  for (v in design_vars) {
    if (is.character(metadata[[v]])) {
      metadata[[v]] <- factor(metadata[[v]])
    }
  }

  if ("condition" %in% colnames(metadata)) {
    metadata$condition <- factor(metadata$condition)

    if (nlevels(metadata$condition) < 2) {
      stop(
        "metadata$condition must contain at least two groups.",
        call. = FALSE
      )
    }

    if ("WT" %in% levels(metadata$condition)) {
      metadata$condition <- stats::relevel(metadata$condition, ref = "WT")
    }
  }

  design <- stats::model.matrix(design_formula, data = metadata)

  contrast_label <- if (is.null(contrast)) {
    "default"
  } else {
    paste(contrast, collapse = "_")
  }

  # -------------------------------------------------------
  # DESeq2
  # -------------------------------------------------------

  if (method == "DESeq2") {

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("Package 'DESeq2' is required.", call. = FALSE)
    }

    log("Running DESeq2 analysis...", type = "info")

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = metadata,
      design = design_formula
    )

    dds <- DESeq2::DESeq(dds)

    if (is.null(contrast)) {
      res <- DESeq2::results(dds)
    } else {
      res <- DESeq2::results(dds, contrast = contrast)
    }

    res_df <- as.data.frame(res)
    res_df$ENSEMBL <- rownames(res_df)
    res_df$contrast <- contrast_label

    rownames(res_df) <- NULL

    sig_n <- sum(
      !is.na(res_df$padj) &
        res_df$padj <= 0.05 &
        abs(res_df$log2FoldChange) >= 0.25,
      na.rm = TRUE
    )

    log("Significant genes:", sig_n, type = "info")
    log(
      "DESeq2 analysis completed.",
      type = "done",
      duration = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )

    return(list(
      method = method,
      design = design_formula,
      design_matrix = design,
      contrast = contrast,
      contrast_label = contrast_label,
      dds = dds,
      res_df = res_df
    ))
  }

  # -------------------------------------------------------
  # edgeR
  # -------------------------------------------------------

  if (method == "edgeR") {

    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("Package 'edgeR' is required.", call. = FALSE)
    }

    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("Package 'limma' is required for edgeR contrasts.", call. = FALSE)
    }

    log("Running edgeR analysis...", type = "info")

    y <- edgeR::DGEList(counts = counts)
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, design)

    fit <- edgeR::glmQLFit(y, design)

    if (is.null(contrast)) {

      qlf <- edgeR::glmQLFTest(fit, coef = ncol(design))

    } else if (is.character(contrast) && length(contrast) == 1) {

      if (contrast %in% colnames(design)) {
        qlf <- edgeR::glmQLFTest(fit, coef = contrast)
      } else {
        contrast_matrix <- limma::makeContrasts(
          contrasts = contrast,
          levels = design
        )
        qlf <- edgeR::glmQLFTest(fit, contrast = contrast_matrix[, 1])
      }

    } else if (is.numeric(contrast)) {

      qlf <- edgeR::glmQLFTest(fit, coef = contrast)

    } else {

      qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
    }

    tab <- edgeR::topTags(qlf, n = Inf)$table

    res_df <- data.frame(
      baseMean = NA_real_,
      log2FoldChange = tab$logFC,
      lfcSE = NA_real_,
      stat = tab$F,
      pvalue = tab$PValue,
      padj = tab$FDR,
      ENSEMBL = rownames(tab),
      contrast = contrast_label,
      stringsAsFactors = FALSE
    )

    rownames(res_df) <- NULL

    sig_n <- sum(
      !is.na(res_df$padj) &
        res_df$padj <= 0.05 &
        abs(res_df$log2FoldChange) >= 0.25,
      na.rm = TRUE
    )

    log("Significant genes:", sig_n, type = "info")
    log(
      "edgeR analysis completed.",
      type = "done",
      duration = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )

    return(list(
      method = method,
      design = design_formula,
      design_matrix = design,
      contrast = contrast,
      contrast_label = contrast_label,
      dge = y,
      fit = fit,
      res_df = res_df
    ))
  }

  # -------------------------------------------------------
  # limma-voom
  # -------------------------------------------------------

  if (method == "limma") {

    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("Package 'edgeR' is required for limma-voom.", call. = FALSE)
    }

    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("Package 'limma' is required.", call. = FALSE)
    }

    log("Running limma-voom analysis...", type = "info")

    y <- edgeR::DGEList(counts = counts)
    y <- edgeR::calcNormFactors(y)

    v <- limma::voom(y, design, plot = FALSE)
    fit <- limma::lmFit(v, design)

    if (is.null(contrast)) {

      fit2 <- limma::eBayes(fit)
      tab <- limma::topTable(
        fit2,
        coef = ncol(design),
        number = Inf,
        sort.by = "P"
      )

    } else if (is.character(contrast) && length(contrast) == 1) {

      if (contrast %in% colnames(design)) {

        fit2 <- limma::eBayes(fit)
        tab <- limma::topTable(
          fit2,
          coef = contrast,
          number = Inf,
          sort.by = "P"
        )

      } else {

        contrast_matrix <- limma::makeContrasts(
          contrasts = contrast,
          levels = design
        )

        fit2 <- limma::contrasts.fit(fit, contrast_matrix)
        fit2 <- limma::eBayes(fit2)

        tab <- limma::topTable(
          fit2,
          coef = 1,
          number = Inf,
          sort.by = "P"
        )
      }

    } else if (is.numeric(contrast)) {

      fit2 <- limma::eBayes(fit)
      tab <- limma::topTable(
        fit2,
        coef = contrast,
        number = Inf,
        sort.by = "P"
      )

    } else {

      stop(
        "For limma, contrast must be NULL, a coefficient name/index, ",
        "or a limma-style contrast string.",
        call. = FALSE
      )
    }

    res_df <- data.frame(
      baseMean = NA_real_,
      log2FoldChange = tab$logFC,
      lfcSE = NA_real_,
      stat = tab$t,
      pvalue = tab$P.Value,
      padj = tab$adj.P.Val,
      ENSEMBL = rownames(tab),
      contrast = contrast_label,
      stringsAsFactors = FALSE
    )

    rownames(res_df) <- NULL

    sig_n <- sum(
      !is.na(res_df$padj) &
        res_df$padj <= 0.05 &
        abs(res_df$log2FoldChange) >= 0.25,
      na.rm = TRUE
    )

    log("Significant genes:", sig_n, type = "info")
    log(
      "limma-voom analysis completed.",
      type = "done",
      duration = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )

    return(list(
      method = method,
      design = design_formula,
      design_matrix = design,
      contrast = contrast,
      contrast_label = contrast_label,
      voom = v,
      fit = fit2,
      res_df = res_df
    ))
  }
}


# =========================================================
# PAIRWISE CONTRAST GENERATION
# =========================================================

#' Generate pairwise contrasts automatically
#'
#' Creates DESeq2-compatible contrast definitions from one or multiple
#' metadata grouping variables.
#'
#' @param metadata Sample metadata data frame.
#' @param group_cols Character vector of metadata columns used to define groups.
#' @param contrast_col Name of the generated grouping column.
#' @param mode Contrast generation mode: \code{"all"},
#'   \code{"within_first"}, or \code{"within_second"}.
#'
#' @return A named list of DESeq2 contrast vectors.
#'
#' @keywords internal
#' @noRd
make_pairwise_contrasts <- function(
    metadata,
    group_cols,
    contrast_col = "comparison_group",
    mode = c("all", "within_first", "within_second")
) {

  mode <- match.arg(mode)

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  missing_cols <- setdiff(group_cols, colnames(metadata))

  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(group_cols) < 1) {
    stop("'group_cols' must contain at least one column.", call. = FALSE)
  }

  if (mode %in% c("within_first", "within_second") &&
      length(group_cols) < 2) {
    stop(
      "Modes 'within_first' and 'within_second' require at least two group columns.",
      call. = FALSE
    )
  }

  metadata[[contrast_col]] <- apply(
    metadata[, group_cols, drop = FALSE],
    1,
    paste,
    collapse = "_"
  )

  groups <- sort(unique(metadata[[contrast_col]]))

  if (length(groups) < 2) {
    stop(
      "At least two comparison groups are required.",
      call. = FALSE
    )
  }

  if (mode == "all") {

    pairs <- utils::combn(groups, 2, simplify = FALSE)

  } else if (mode == "within_first") {

    split_groups <- split(metadata, metadata[[group_cols[1]]])

    pairs <- unlist(
      lapply(split_groups, function(df) {
        x <- sort(unique(df[[contrast_col]]))
        if (length(x) < 2) return(list())
        utils::combn(x, 2, simplify = FALSE)
      }),
      recursive = FALSE
    )

  } else {

    split_groups <- split(metadata, metadata[[group_cols[2]]])

    pairs <- unlist(
      lapply(split_groups, function(df) {
        x <- sort(unique(df[[contrast_col]]))
        if (length(x) < 2) return(list())
        utils::combn(x, 2, simplify = FALSE)
      }),
      recursive = FALSE
    )
  }

  if (length(pairs) == 0) {
    stop(
      "No valid pairwise contrasts could be generated.",
      call. = FALSE
    )
  }

  out <- lapply(pairs, function(x) {
    c(contrast_col, x[2], x[1])
  })

  names(out) <- vapply(
    pairs,
    function(x) {
      paste0(
        make.names(x[2]),
        "_vs_",
        make.names(x[1])
      )
    },
    character(1)
  )

  out
}


# =========================================================
# RUN MULTIPLE DESEQ2 PAIRWISE CONTRASTS
# =========================================================

#' Run multiple DESeq2 pairwise comparisons
#'
#' Performs differential expression analysis across multiple automatically
#' generated or user-defined pairwise contrasts.
#'
#' @param counts Raw count matrix with genes as rows and samples as columns.
#' @param metadata Sample metadata data frame.
#' @param group_cols Character vector of metadata columns used to define groups.
#' @param contrast_col Name of the generated comparison column.
#' @param contrast_list Optional named list of DESeq2 contrasts.
#' @param design_formula Optional DESeq2 design formula.
#' @param pairwise_mode Contrast generation mode.
#'
#' @return A list containing the DESeq2 object, result tables, summary,
#'   metadata, contrast column, contrast levels, and contrast list.
#'
#' @keywords internal
#' @noRd
run_deseq2_pairwise <- function(
    counts,
    metadata,
    group_cols,
    contrast_col = "comparison_group",
    contrast_list = NULL,
    design_formula = NULL,
    pairwise_mode = c("all", "within_first", "within_second")
) {

  pairwise_mode <- match.arg(pairwise_mode)
  log <- .deggo_msg(verbose = TRUE, prefix = "DEGgo")
  t0 <- Sys.time()

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required.", call. = FALSE)
  }

  counts <- round(as.matrix(counts))
  storage.mode(counts) <- "integer"

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if ("sample" %in% colnames(metadata)) {
    rownames(metadata) <- metadata$sample
  }

  if (missing(group_cols) || !length(group_cols)) {
    stop("'group_cols' must contain at least one metadata column.", call. = FALSE)
  }

  missing_cols <- setdiff(group_cols, colnames(metadata))
  if (length(missing_cols)) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  sample_ids <- colnames(counts)
  missing_meta <- setdiff(sample_ids, rownames(metadata))

  if (length(missing_meta)) {
    stop(
      "Samples missing in metadata rownames: ",
      paste(missing_meta, collapse = ", "),
      call. = FALSE
    )
  }

  metadata <- metadata[match(sample_ids, rownames(metadata)), , drop = FALSE]

  if (!identical(colnames(counts), rownames(metadata))) {
    stop("Sample order mismatch after matching.", call. = FALSE)
  }

  metadata[[contrast_col]] <- apply(
    metadata[, group_cols, drop = FALSE],
    1,
    paste,
    collapse = "_"
  )

  metadata[[contrast_col]] <- factor(metadata[[contrast_col]])

  if (is.null(design_formula)) {
    design_formula <- stats::as.formula(paste("~", contrast_col))
  }

  missing_design <- setdiff(all.vars(design_formula), colnames(metadata))
  if (length(missing_design)) {
    stop(
      "Variables missing from metadata: ",
      paste(missing_design, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(contrast_list)) {
    contrast_list <- make_pairwise_contrasts(
      metadata = metadata,
      group_cols = group_cols,
      contrast_col = contrast_col,
      mode = pairwise_mode
    )
  }

  if (is.null(names(contrast_list)) || any(names(contrast_list) == "")) {
    names(contrast_list) <- paste0("contrast_", seq_along(contrast_list))
  }

  available_levels <- levels(metadata[[contrast_col]])

  needed_levels <- unique(
    unlist(
      lapply(contrast_list, function(x) x[2:3]),
      use.names = FALSE
    )
  )

  missing_levels <- setdiff(needed_levels, available_levels)

  if (length(missing_levels)) {
    stop(
      "Missing contrast levels in metadata$",
      contrast_col,
      ": ",
      paste(missing_levels, collapse = ", "),
      "\nAvailable levels: ",
      paste(available_levels, collapse = ", "),
      call. = FALSE
    )
  }

  log(
    paste(
      "Pairwise contrast levels:",
      paste(available_levels, collapse = ", ")
    ),
    type = "info"
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = design_formula
  )

  log("Running DESeq2 model for pairwise analysis...", type = "info")
  dds <- DESeq2::DESeq(dds)

  out <- vector("list", length(contrast_list))
  samples <- vector("list", length(contrast_list))

  names(out) <- names(contrast_list)
  names(samples) <- names(contrast_list)

  for (nm in names(contrast_list)) {

    contrast <- contrast_list[[nm]]
    levels_use <- contrast[2:3]

    samples[[nm]] <- rownames(
      metadata[
        metadata[[contrast_col]] %in% levels_use,
        ,
        drop = FALSE
      ]
    )

    res <- DESeq2::results(
      dds,
      contrast = contrast
    )

    res_df <- as.data.frame(res)
    res_df$ENSEMBL <- rownames(res_df)
    res_df$comparison <- nm
    res_df$contrast <- paste(contrast, collapse = "_")

    res_df <- res_df[
      order(res_df$padj, -abs(res_df$log2FoldChange)),
      ,
      drop = FALSE
    ]

    rownames(res_df) <- NULL
    out[[nm]] <- res_df
  }

  summary <- do.call(
    rbind,
    lapply(names(out), function(nm) {

      x <- out[[nm]]

      data.frame(
        comparison = nm,
        total_genes = nrow(x),
        significant = sum(
          !is.na(x$padj) &
            x$padj <= 0.05 &
            abs(x$log2FoldChange) >= 0.25,
          na.rm = TRUE
        ),
        up = sum(
          !is.na(x$padj) &
            x$padj <= 0.05 &
            x$log2FoldChange >= 0.25,
          na.rm = TRUE
        ),
        down = sum(
          !is.na(x$padj) &
            x$padj <= 0.05 &
            x$log2FoldChange <= -0.25,
          na.rm = TRUE
        ),
        n_samples = length(samples[[nm]]),
        stringsAsFactors = FALSE
      )
    })
  )

  log(
    "Pairwise DESeq2 analysis completed.",
    type = "done",
    duration = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )

  list(
    dds = dds,
    results = out,
    samples = samples,
    summary = summary,
    metadata = metadata,
    contrast_col = contrast_col,
    contrast_levels = available_levels,
    contrast_list = contrast_list,
    design = design_formula
  )
}
