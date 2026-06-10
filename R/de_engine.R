# =========================================================
# RUN DIFFERENTIAL EXPRESSION
# =========================================================

run_de <- function(
    counts,
    metadata,
    method = c("DESeq2", "edgeR", "limma"),
    design_formula = ~ condition,
    contrast = NULL
) {
  
  method <- match.arg(method)
  
  if (!"condition" %in% colnames(metadata)) {
    stop("metadata must contain 'condition'.", call. = FALSE)
  }
  
  metadata$condition <- factor(metadata$condition)
  
  if ("WT" %in% levels(metadata$condition)) {
    metadata$condition <- stats::relevel(metadata$condition, ref = "WT")
  }
  
  # -------------------------------------------------------
  # DESeq2
  # -------------------------------------------------------
  
  if (method == "DESeq2") {
    
    log_msg("Running DESeq2 analysis...")
    
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
    
    return(list(
      method = method,
      object = dds,
      dds = dds,
      res_df = res_df
    ))
  }
  
  # -------------------------------------------------------
  # edgeR
  # -------------------------------------------------------
  
  if (method == "edgeR") {
    
    log_msg("Running edgeR analysis...")
    
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("Package 'edgeR' is required.", call. = FALSE)
    }
    
    y <- edgeR::DGEList(counts = counts)
    y <- edgeR::calcNormFactors(y)
    
    design <- stats::model.matrix(design_formula, data = metadata)
    
    y <- edgeR::estimateDisp(y, design)
    
    fit <- edgeR::glmQLFit(y, design)
    
    coef_use <- if (is.null(contrast)) {
      ncol(design)
    } else {
      contrast
    }
    
    if (is.numeric(coef_use) || is.character(coef_use)) {
      qlf <- edgeR::glmQLFTest(fit, coef = coef_use)
    } else {
      qlf <- edgeR::glmQLFTest(fit, contrast = coef_use)
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
      stringsAsFactors = FALSE
    )
    
    return(list(
      method = method,
      object = y,
      fit = fit,
      res_df = res_df
    ))
  }
  
  # -------------------------------------------------------
  # limma-voom
  # -------------------------------------------------------
  
  if (method == "limma") {
    
    log_msg("Running limma-voom analysis...")
    
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("Package 'edgeR' is required.", call. = FALSE)
    }
    
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("Package 'limma' is required.", call. = FALSE)
    }
    
    y <- edgeR::DGEList(counts = counts)
    y <- edgeR::calcNormFactors(y)
    
    design <- stats::model.matrix(design_formula, data = metadata)
    
    v <- limma::voom(y, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
    
    coef_use <- if (is.null(contrast)) {
      ncol(design)
    } else {
      contrast
    }
    
    fit2 <- limma::eBayes(fit)
    
    tab <- limma::topTable(
      fit2,
      coef = coef_use,
      number = Inf,
      sort.by = "P"
    )
    
    res_df <- data.frame(
      baseMean = NA_real_,
      log2FoldChange = tab$logFC,
      lfcSE = NA_real_,
      stat = tab$t,
      pvalue = tab$P.Value,
      padj = tab$adj.P.Val,
      ENSEMBL = rownames(tab),
      stringsAsFactors = FALSE
    )
    
    return(list(
      method = method,
      object = v,
      fit = fit2,
      res_df = res_df
    ))
  }
}



#' Generate pairwise contrasts automatically
#'
#' Creates DESeq2-compatible contrast definitions from one or
#' multiple metadata grouping variables.
#'
#' This helper is used internally by
#' \code{run_deseq2_pairwise()} to automatically generate
#' pairwise comparisons between experimental groups.
#'
#' @param metadata Sample metadata data frame.
#' @param group_cols Character vector of metadata columns used
#'   to define groups.
#' @param contrast_col Name of the generated grouping column.
#' @param mode Character specifying how contrasts are generated:
#'   \itemize{
#'   \item \code{"all"}: all pairwise combinations.
#'   \item \code{"within_first"}: compare levels within the first grouping variable.
#'   \item \code{"within_second"}: compare levels within the second grouping variable.
#'   }
#'
#' @return A named list of DESeq2 contrast vectors.
#'
#' @details
#' For example, using:
#'
#' \preformatted{
#' group_cols = c("condition", "sex")
#' }
#'
#' creates combined groups such as:
#'
#' \preformatted{
#' WT_M
#' WT_F
#' TAU22_M
#' TAU22_F
#' }
#'
#' and automatically generates pairwise contrasts according
#' to the selected mode.
#'
#' @keywords internal
#' @noRd
#' 
make_pairwise_contrasts <- function(
    metadata,
    group_cols,
    contrast_col = "comparison_group",
    mode = c("all", "within_first", "within_second")
) {
  
  mode <- match.arg(mode)
  
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  
  metadata[[contrast_col]] <- apply(
    metadata[, group_cols, drop = FALSE],
    1,
    paste,
    collapse = "_"
  )
  
  groups <- sort(unique(metadata[[contrast_col]]))
  
  if (mode == "all") {
    pairs <- utils::combn(groups, 2, simplify = FALSE)
  }
  
  if (mode == "within_first") {
    # Example group_cols = c("condition", "sex")
    # compares sex levels within each condition:
    # WT_F vs WT_M, TAU22_F vs TAU22_M
    split_groups <- split(metadata, metadata[[group_cols[1]]])
    
    pairs <- unlist(lapply(split_groups, function(df) {
      sort(unique(df[[contrast_col]])) |>
        utils::combn(2, simplify = FALSE)
    }), recursive = FALSE)
  }
  
  if (mode == "within_second") {
    # Example group_cols = c("condition", "sex")
    # compares condition levels within each sex:
    # WT_M vs TAU22_M, WT_F vs TAU22_F
    split_groups <- split(metadata, metadata[[group_cols[2]]])
    
    pairs <- unlist(lapply(split_groups, function(df) {
      sort(unique(df[[contrast_col]])) |>
        utils::combn(2, simplify = FALSE)
    }), recursive = FALSE)
  }
  
  out <- lapply(pairs, function(x) {
    c(contrast_col, x[2], x[1])
  })
  
  names(out) <- vapply(pairs, function(x) {
    paste0(x[2], "_vs_", x[1])
  }, character(1))
  
  out
}


#' Run multiple DESeq2 pairwise comparisons
#'
#' Performs differential expression analysis across multiple
#' automatically generated or user-defined pairwise contrasts.
#'
#' Groups are constructed from one or more metadata columns,
#' allowing flexible comparisons such as:
#'
#' \itemize{
#'   \item genotype comparisons
#'   \item sex comparisons
#'   \item treatment comparisons
#'   \item genotype-by-sex comparisons
#' }
#'
#' @param counts Raw count matrix with genes as rows and
#'   samples as columns.
#' @param metadata Sample metadata data frame.
#' @param group_cols Character vector of metadata columns used
#'   to define comparison groups.
#' @param contrast_col Name of the generated comparison column.
#' @param contrast_list Optional named list of DESeq2 contrasts.
#'   If \code{NULL}, contrasts are generated automatically.
#' @param design_formula DESeq2 design formula. If \code{NULL},
#'   a formula based on \code{contrast_col} is created.
#' @param pairwise_mode Contrast generation mode:
#'   \itemize{
#'   \item \code{"all"}
#'   \item \code{"within_first"}
#'   \item \code{"within_second"}
#'   }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{dds}: DESeq2 dataset.
#'   \item \code{results}: list of DESeq2 result tables.
#'   \item \code{summary}: summary statistics for each contrast.
#'   \item \code{metadata}: processed metadata.
#'   \item \code{contrast_col}: grouping column used.
#'   \item \code{contrast_levels}: available comparison levels.
#' }
#'
#' @details
#' When \code{contrast_list = NULL},
#' \code{make_pairwise_contrasts()} is used to generate all
#' requested pairwise comparisons automatically.
#'
#' Example:
#'
#' \preformatted{
#' group_cols = c("condition", "sex")
#' }
#'
#' creates:
#'
#' \preformatted{
#' WT_M
#' WT_F
#' TAU22_M
#' TAU22_F
#' }
#'
#' and generates contrasts according to
#' \code{pairwise_mode}.
#'
#' @keywords internal
#' @noRd
#' 
run_deseq2_pairwise <- function(
    counts,
    metadata,
    group_cols,
    contrast_col = "comparison_group",
    contrast_list,
    design_formula = NULL,
    pairwise_mode = c("all", "within_first", "within_second")
) {
  
  pairwise_mode <- match.arg(pairwise_mode)
  
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required.", call. = FALSE)
  }
  
  counts <- round(as.matrix(counts))
  mode(counts) <- "integer"
  
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  
  if ("sample" %in% colnames(metadata)) {
    rownames(metadata) <- metadata$sample
  }
  
  if (missing(group_cols) || length(group_cols) < 1) {
    stop("'group_cols' must contain at least one metadata column.", call. = FALSE)
  }
  
  if (is.null(contrast_list)) {
    contrast_list <- make_pairwise_contrasts(
      metadata = metadata,
      group_cols = group_cols,
      contrast_col = contrast_col,
      mode = pairwise_mode
    )
  }
  
  missing_cols <- setdiff(group_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  sample_ids <- colnames(counts)
  
  missing_meta <- setdiff(sample_ids, rownames(metadata))
  if (length(missing_meta) > 0) {
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
    design_formula <- stats::as.formula(
      paste("~", contrast_col)
    )
  }
  
  available_levels <- levels(metadata[[contrast_col]])
  
  needed_levels <- unique(unlist(lapply(contrast_list, function(x) x[2:3])))
  missing_levels <- setdiff(needed_levels, available_levels)
  
  if (length(missing_levels) > 0) {
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
  
  message("Pairwise contrast levels: ", paste(available_levels, collapse = ", "))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = design_formula
  )
  
  dds <- DESeq2::DESeq(dds)
  
  out <- lapply(names(contrast_list), function(nm) {
    
    contrast <- contrast_list[[nm]]
    
    res <- DESeq2::results(
      dds,
      contrast = contrast
    )
    
    res_df <- as.data.frame(res)
    res_df$ENSEMBL <- rownames(res_df)
    res_df$comparison <- nm
    
    res_df <- res_df[
      order(res_df$padj, -abs(res_df$log2FoldChange)),
      ,
      drop = FALSE
    ]
    
    rownames(res_df) <- NULL
    res_df
  })
  
  names(out) <- names(contrast_list)
  
  summary <- do.call(
    rbind,
    lapply(names(out), function(nm) {
      x <- out[[nm]]
      data.frame(
        comparison = nm,
        total_genes = nrow(x),
        significant = sum(
          x$padj < 0.05 & abs(x$log2FoldChange) >= 0.25,
          na.rm = TRUE
        ),
        up = sum(
          x$padj < 0.05 & x$log2FoldChange >= 0.25,
          na.rm = TRUE
        ),
        down = sum(
          x$padj < 0.05 & x$log2FoldChange <= -0.25,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    })
  )
  
  list(
    dds = dds,
    results = out,
    summary = summary,
    metadata = metadata,
    contrast_col = contrast_col,
    contrast_levels = available_levels
  )
}
