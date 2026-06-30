# ========================================================= #
# RUN DIFFERENTIAL EXPRESSION
# ========================================================= #

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
.run_de <- function(
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

  # ------------------------------------------------------- #
  # DESeq2
  # ------------------------------------------------------- #

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

  # ------------------------------------------------------- #
  # edgeR
  # ------------------------------------------------------- #

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

  # ------------------------------------------------------- #
  # limma-voom
  # ------------------------------------------------------- #

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




# ============================================================ #
# Single mode
# ============================================================ #
#' Run DEGgo single-comparison analysis
#'
#' Internal helper used by `run_deggo()` to perform a standard
#' differential expression analysis workflow for a single contrast.
#'
#' The workflow includes:
#' \itemize{
#'   \item Differential expression analysis.
#'   \item Gene annotation.
#'   \item Significant DEG extraction.
#'   \item Volcano plot generation.
#'   \item PCA visualization (DESeq2 only).
#'   \item Heatmap generation.
#'   \item Gene Ontology enrichment analysis.
#'   \item Result export and summary generation.
#' }
#'
#' @param counts Integer count matrix.
#' @param metadata Sample metadata data frame.
#' @param method Differential expression method.
#' @param design_formula Design formula passed to DESeq2, edgeR or limma.
#' @param contrast Contrast vector for differential expression analysis.
#' @param orgdb Annotation database object.
#' @param dirs Output directory list generated by `.deggo_dirs()`.
#' @param ontology GO ontology ("BP", "MF", or "CC").
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param logfc_cutoff Absolute log2 fold-change cutoff.
#' @param top_n_heatmap Number of top genes displayed in heatmaps.
#' @param top_n_labels Number of labeled genes in volcano plots.
#' @param txtsize Base text size used for plots.
#' @param log DEGgo logging function.
#'
#' @return A list containing differential expression results,
#' significant genes, plots, GO enrichment results and summary statistics.
#'
#' @keywords internal
#'
.run_deggo_single <- function(
    counts, metadata, method, design_formula, contrast,
    orgdb, dirs, ontology,
    padj_cutoff, logfc_cutoff,
    top_n_heatmap, top_n_labels,
    txtsize, log
) {

  log("[7/11] Running differential expression analysis", type = "step")

  de_results <- .run_de(
    counts = counts,
    metadata = metadata,
    method = method,
    design_formula = design_formula,
    contrast = contrast
  )

  dds <- de_results$dds %||% de_results$object

  if (method == "DESeq2") {
    dds <- .annotate_dds(dds, orgdb = orgdb)
  }

  res_df <- de_results$res_df

  log("[8/11] Annotating DE results", type = "step")

  res_df$SYMBOL <- map_ensembl_to_feature(
    ensembl_vec = res_df$ENSEMBL,
    orgdb = orgdb
  )

  res_df <- map_entrez_ids(
    res_df = res_df,
    orgdb = orgdb
  )

  processed <- process_deg_results(
    res_df = res_df,
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff
  )

  res_df <- processed$res_df
  sig_deg <- processed$sig_deg

  log("[9/11] Generating volcano, PCA and heatmap", type = "step")

  volcano_plot <- plot_volcano(
    res_df = res_df,
    top_n_labels = top_n_labels,
    txtsize = txtsize,
    output_dir = dirs$volcano,
    filename = "single_comparison_Volcano_Plot",
    title = "single_comparison",
    logfc_cutoff = logfc_cutoff,
    padj_cutoff = padj_cutoff
  )

  pca_plot <- NULL
  heatmap_matrix <- NULL

  if (method == "DESeq2") {

    pca_results <- plot_pca(
      dds = dds,
      metadata = metadata,
      output_dir = dirs$pca
    )

    pca_plot <- pca_results$pca_plot

    heatmap_matrix <- plot_heatmap(
      vsd = pca_results$vsd,
      res_df = res_df,
      metadata = metadata,
      top_n_heatmap = top_n_heatmap,
      padj_cutoff = padj_cutoff,
      output_dir = dirs$heatmaps,
      main = "single_comparison",
      filename = "single_comparison_Heatmap",
      annotation_cols = intersect(
        c("condition", "treatment", "sex", "tissue"),
        colnames(metadata)
      ),
      order_by = intersect(
        c("condition", "tissue", "sex", "treatment"),
        colnames(metadata)
      ),
      fontsize_row = txtsize,
      fontsize_col = max(5, txtsize - 1)
    )
  }

  log("[10/11] Running GO enrichment", type = "step")

  go_results <- run_go_enrichment(
    sig_deg = sig_deg,
    comparison = "single_comparison",
    ontology = ontology,
    output_dir = dirs$go,
    orgdb = orgdb
  )

  go_plot <- NULL
  go_df <- go_results$go_results

  if (!is.null(go_df) && nrow(go_df) > 0) {

    go_plot <- plot_go_terms(
      go_df = go_df,
      comparison = "single_comparison GO enrichment",
      top_n = 10,
      style = "bw",
      txtsize = txtsize
    )

    ggplot2::ggsave(
      file.path(dirs$go, "single_comparison_GO_terms.png"),
      go_plot,
      width = 8,
      height = 6,
      dpi = 300
    )

    ggplot2::ggsave(
      file.path(dirs$go, "single_comparison_GO_terms.pdf"),
      go_plot,
      width = 8,
      height = 6
    )
  }

  log("[11/11] Exporting final results", type = "step")

  utils::write.table(
    res_df,
    file.path(dirs$results, "single_comparison.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  utils::write.table(
    sig_deg,
    file.path(dirs$significant, "single_comparison_significant.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  de_results$dds <- dds
  de_results$res_df <- res_df

  de_results$results <- list(
    single_comparison = res_df
  )

  de_results$sig_deg <- list(
    single_comparison = sig_deg
  )

  de_results$volcano_plot <- volcano_plot
  de_results$volcano_plots <- list(
    single_comparison = volcano_plot
  )

  de_results$pca_plot <- pca_plot
  de_results$pca <- list(
    single_comparison = pca_plot
  )

  de_results$heatmap_matrix <- heatmap_matrix
  de_results$heatmaps <- list(
    single_comparison = heatmap_matrix
  )

  de_results$go_results <- list(
    single_comparison = go_results
  )

  de_results$go_plot <- go_plot
  de_results$go_plots <- list(
    single_comparison = go_plot
  )

  de_results$summary <- data.frame(
    comparison = "single_comparison",
    total_genes = nrow(res_df),
    significant = nrow(sig_deg),
    up = sum(sig_deg$log2FoldChange > logfc_cutoff, na.rm = TRUE),
    down = sum(sig_deg$log2FoldChange < -logfc_cutoff, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (method == "DESeq2") {
    de_results$gene_annotation <- data.frame(
      ENSEMBL = rownames(dds),
      SYMBOL = SummarizedExperiment::rowData(dds)$SYMBOL,
      ENTREZID = SummarizedExperiment::rowData(dds)$ENTREZID,
      GENENAME = SummarizedExperiment::rowData(dds)$GENENAME,
      stringsAsFactors = FALSE
    )
  }

  de_results
}


# ============================================================ #
# Pairwise mode
# ============================================================ #
#' Run DEGgo pairwise differential expression analysis
#'
#' Internal helper used by `run_deggo()` to perform automated
#' pairwise differential expression analyses using DESeq2.
#'
#' The workflow includes:
#' \itemize{
#'   \item Pairwise contrast generation and testing.
#'   \item PCA analysis.
#'   \item DEG annotation.
#'   \item Volcano plot generation.
#'   \item Heatmap generation.
#'   \item GO enrichment analysis.
#'   \item Summary table generation.
#' }
#'
#' @param counts Integer count matrix.
#' @param metadata Sample metadata data frame.
#' @param orgdb Annotation database object.
#' @param dirs Output directory list generated by `.deggo_dirs()`.
#' @param ontology GO ontology ("BP", "MF", or "CC").
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param logfc_cutoff Absolute log2 fold-change cutoff.
#' @param top_n_heatmap Number of top genes displayed in heatmaps.
#' @param top_n_labels Number of labeled genes in volcano plots.
#' @param txtsize Base text size used for plots.
#' @param pairwise_group_cols Metadata columns defining pairwise groups.
#' @param pairwise_contrast_col Column storing generated contrasts.
#' @param pairwise_contrasts Optional named list of pairwise contrasts.
#' @param pairwise_mode Pairwise comparison mode.
#' @param log DEGgo logging function.
#'
#' @return A list containing pairwise differential expression results,
#' significant genes, plots, GO enrichment results and summary statistics.
#'
#' @keywords internal
#'
.run_deggo_pairwise <- function(
    counts, metadata, orgdb, dirs, ontology,
    padj_cutoff, logfc_cutoff,
    top_n_heatmap, top_n_labels,
    txtsize,
    pairwise_group_cols,
    pairwise_contrast_col,
    pairwise_contrasts,
    pairwise_mode,
    log
) {

  log("[7/11] Running pairwise DESeq2 contrasts", type = "step")

  de_results <- .run_deseq2_pairwise(
    counts = counts,
    metadata = metadata,
    group_cols = pairwise_group_cols,
    contrast_col = pairwise_contrast_col,
    contrast_list = pairwise_contrasts,
    pairwise_mode = pairwise_mode
  )

  de_results$dds <- .annotate_dds(de_results$dds, orgdb = orgdb)

  log("[8/11] PCA analysis", type = "step")

  de_results$pca <- .make_pca_list(
    dds = de_results$dds,
    md = de_results$metadata,
    pca_dir = dirs$pca,
    txtsize = txtsize
  )

  log("[9/11] Annotating DE tables", type = "step")

  de_results$results <- lapply(
    de_results$results,
    annotate_de_results,
    orgdb = orgdb,
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff
  )

  de_results$sig_deg <- list()
  de_results$volcano_plots <- list()
  de_results$heatmaps <- list()
  de_results$go_results <- list()
  de_results$go_plots <- list()

  log("[10/11] Exporting plots, heatmaps and GO", type = "step")

  vsd_pairwise <- DESeq2::vst(de_results$dds, blind = FALSE)

  for (nm in names(de_results$results)) {

    log(paste("Processing comparison:", nm), type = "info")

    res_df <- de_results$results[[nm]]

    sig_df <- res_df[
      !is.na(res_df$padj) &
        res_df$padj < padj_cutoff &
        abs(res_df$log2FoldChange) > logfc_cutoff,
      ,
      drop = FALSE
    ]

    de_results$sig_deg[[nm]] <- sig_df

    utils::write.table(
      res_df,
      file.path(dirs$results, paste0(nm, ".tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    utils::write.table(
      sig_df,
      file.path(dirs$significant, paste0(nm, "_significant.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    de_results$volcano_plots[[nm]] <- plot_volcano(
      res_df = res_df,
      top_n_labels = top_n_labels,
      txtsize = txtsize,
      output_dir = dirs$volcano,
      filename = paste0(nm, "_Volcano_Plot"),
      title = nm,
      logfc_cutoff = logfc_cutoff,
      padj_cutoff = padj_cutoff
    )

    de_results$heatmaps[[nm]] <- plot_heatmap(
      vsd = vsd_pairwise,
      res_df = res_df,
      metadata = de_results$metadata,
      sample_subset = de_results$samples[[nm]],
      top_n_heatmap = top_n_heatmap,
      padj_cutoff = padj_cutoff,
      output_dir = dirs$heatmaps,
      main = nm,
      filename = paste0(nm, "_Heatmap"),
      annotation_cols = intersect(
        c("condition", "treatment", "sex", "tissue"),
        colnames(de_results$metadata)
      ),
      order_by = intersect(
        c("condition", "tissue", "sex", "treatment"),
        colnames(de_results$metadata)
      ),
      fontsize_row = txtsize,
      fontsize_col = max(5, txtsize - 1)
    )

    de_results$go_results[[nm]] <- run_go_enrichment(
      sig_deg = sig_df,
      comparison = nm,
      ontology = ontology,
      output_dir = dirs$go,
      orgdb = orgdb
    )

    go_df <- de_results$go_results[[nm]]$go_results

    if (!is.null(go_df) && nrow(go_df) > 0) {

      de_results$go_plots[[nm]] <- plot_go_terms(
        go_df = go_df,
        comparison = paste0(nm, " GO enrichment"),
        top_n = 10,
        style = "bw",
        txtsize = txtsize
      )

      ggplot2::ggsave(
        file.path(dirs$go, paste0(nm, "_GO_terms.png")),
        de_results$go_plots[[nm]],
        width = 8,
        height = 6,
        dpi = 300
      )

      ggplot2::ggsave(
        file.path(dirs$go, paste0(nm, "_GO_terms.pdf")),
        de_results$go_plots[[nm]],
        width = 8,
        height = 6
      )
    }
  }

  de_results$go_merged <- .merge_pairwise_go(
    go_results = de_results$go_results,
    go_dir = dirs$go,
    ontology = ontology
  )

  log("[11/11] Summarizing pairwise results", type = "step")

  de_results$summary <- do.call(
    rbind,
    lapply(names(de_results$results), function(nm) {

      x <- de_results$results[[nm]]

      data.frame(
        comparison = nm,
        total_genes = nrow(x),
        significant = sum(
          !is.na(x$padj) &
            x$padj < padj_cutoff &
            abs(x$log2FoldChange) > logfc_cutoff,
          na.rm = TRUE
        ),
        up = sum(
          !is.na(x$padj) &
            x$padj < padj_cutoff &
            x$log2FoldChange > logfc_cutoff,
          na.rm = TRUE
        ),
        down = sum(
          !is.na(x$padj) &
            x$padj < padj_cutoff &
            x$log2FoldChange < -logfc_cutoff,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    })
  )

  de_results
}



# ========================================================= #
# PAIRWISE CONTRAST GENERATION
# ========================================================= #

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
.make_pairwise_contrasts <- function(
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


# ========================================================= #
# RUN MULTIPLE DESEQ2 PAIRWISE CONTRASTS
# ========================================================= #

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
.run_deseq2_pairwise <- function(
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
    contrast_list <- .make_pairwise_contrasts(
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
