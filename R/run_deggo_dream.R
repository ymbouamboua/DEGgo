# ============================================================ #
# dream differential expression engine
# ============================================================ #

.run_deggo_dream <- function(
    counts,
    metadata,
    design_formula,
    contrast,
    orgdb = NULL,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25,
    ddf = c(
      "adaptive",
      "Satterthwaite",
      "Kenward-Roger"
    ),
    n_cores = 1L,
    log = NULL
) {
  ddf <- match.arg(ddf)

  required_packages <- c(
    "edgeR",
    "limma",
    "variancePartition",
    "BiocParallel"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_packages)) {
    stop(
      "The dream method requires the following package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Install them with:\n",
      "BiocManager::install(c(",
      paste(
        sprintf('"%s"', missing_packages),
        collapse = ", "
      ),
      "))",
      call. = FALSE
    )
  }

  if (!inherits(design_formula, "formula")) {
    stop(
      "'design_formula' must be a formula.",
      call. = FALSE
    )
  }

  if (
    is.null(contrast) ||
    length(contrast) != 3L
  ) {
    stop(
      paste(
        "For dream, 'contrast' must have the form",
        "c(variable, numerator, denominator)."
      ),
      call. = FALSE
    )
  }

  contrast <- as.character(contrast)

  contrast_var <- contrast[1]
  numerator <- contrast[2]
  denominator <- contrast[3]

  if (!contrast_var %in% colnames(metadata)) {
    stop(
      "Contrast variable '",
      contrast_var,
      "' was not found in metadata.",
      call. = FALSE
    )
  }

  formula_variables <- all.vars(design_formula)

  if (!contrast_var %in% formula_variables) {
    stop(
      "Contrast variable '",
      contrast_var,
      "' is not present in 'design_formula'.",
      call. = FALSE
    )
  }

  metadata <- as.data.frame(
    metadata,
    stringsAsFactors = FALSE
  )

  contrast_values <- as.character(
    metadata[[contrast_var]]
  )

  missing_levels <- setdiff(
    c(numerator, denominator),
    unique(contrast_values)
  )

  if (length(missing_levels)) {
    stop(
      "Contrast level(s) missing from metadata column '",
      contrast_var,
      "': ",
      paste(missing_levels, collapse = ", "),
      call. = FALSE
    )
  }

  # The denominator becomes the reference level. Consequently,
  # the numerator coefficient represents numerator versus denominator.
  metadata[[contrast_var]] <- stats::relevel(
    factor(contrast_values),
    ref = denominator
  )

  # Convert other character variables used in the model to factors.
  for (variable in formula_variables) {
    if (
      variable %in% colnames(metadata) &&
      is.character(metadata[[variable]])
    ) {
      metadata[[variable]] <- factor(
        metadata[[variable]]
      )
    }
  }

  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"

  if (is.null(rownames(counts))) {
    stop(
      "'counts' must have gene identifiers as row names.",
      call. = FALSE
    )
  }

  if (is.null(colnames(counts))) {
    stop(
      "'counts' must have sample identifiers as column names.",
      call. = FALSE
    )
  }

  if (anyNA(counts)) {
    stop(
      "'counts' contains missing values.",
      call. = FALSE
    )
  }

  if (any(!is.finite(counts))) {
    stop(
      "'counts' contains non-finite values.",
      call. = FALSE
    )
  }

  if (any(counts < 0)) {
    stop(
      "'counts' contains negative values.",
      call. = FALSE
    )
  }

  if (any(abs(counts - round(counts)) > 1e-8)) {
    stop(
      "'counts' must contain raw integer-like counts for dream.",
      call. = FALSE
    )
  }

  counts <- round(counts)
  storage.mode(counts) <- "integer"

  if (
    is.null(rownames(metadata)) ||
    any(!nzchar(rownames(metadata)))
  ) {
    stop(
      "'metadata' must have sample identifiers as row names.",
      call. = FALSE
    )
  }

  missing_metadata <- setdiff(
    colnames(counts),
    rownames(metadata)
  )

  if (length(missing_metadata)) {
    stop(
      "Sample(s) missing from metadata: ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }

  metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE
  ]

  stopifnot(
    identical(
      colnames(counts),
      rownames(metadata)
    )
  )

  if (is.function(log)) {
    log(
      "[DE] Running dream mixed-effects analysis",
      type = "step"
    )
  }

  # ---------------------------------------------------------- #
  # edgeR normalization
  # ---------------------------------------------------------- #

  dge <- edgeR::DGEList(
    counts = counts
  )

  keep_library <- dge$samples$lib.size > 0

  if (!all(keep_library)) {
    stop(
      "One or more samples have a library size of zero.",
      call. = FALSE
    )
  }

  dge <- edgeR::calcNormFactors(
    dge,
    method = "TMM"
  )

  # Filtering should normally already have been performed by DEGgo.
  # This secondary guard removes genes with no information.
  keep_gene <- rowSums(dge$counts) > 0

  dge <- dge[
    keep_gene,
    ,
    keep.lib.sizes = FALSE
  ]

  if (!nrow(dge)) {
    stop(
      "No expressed genes remain for dream analysis.",
      call. = FALSE
    )
  }

  dge <- edgeR::calcNormFactors(
    dge,
    method = "TMM"
  )

  # ---------------------------------------------------------- #
  # Parallel backend
  # ---------------------------------------------------------- #

  n_cores <- max(
    1L,
    as.integer(n_cores)[1]
  )

  bp_param <- if (n_cores == 1L) {
    BiocParallel::SerialParam(
      progressbar = FALSE
    )
  } else {
    BiocParallel::SnowParam(
      workers = n_cores,
      type = "SOCK",
      progressbar = FALSE
    )
  }

  # ---------------------------------------------------------- #
  # voomWithDreamWeights + dream
  # ---------------------------------------------------------- #

  voom_object <- variancePartition::voomWithDreamWeights(
    counts = dge,
    formula = design_formula,
    data = metadata,
    plot = FALSE,
    BPPARAM = bp_param
  )

  fit <- variancePartition::dream(
    exprObj = voom_object,
    formula = design_formula,
    data = metadata,
    ddf = ddf,
    BPPARAM = bp_param
  )

  fit <- variancePartition::eBayes(
    fit
  )

  coefficient_names <- colnames(
    fit$coefficients
  )

  expected_coefficient <- paste0(
    contrast_var,
    make.names(numerator)
  )

  # Depending on factor names and versions, model.matrix() can apply
  # slightly different syntactic transformations.
  coefficient_candidates <- unique(c(
    expected_coefficient,
    make.names(
      paste0(
        contrast_var,
        numerator
      )
    ),
    grep(
      paste0(
        "^",
        make.names(contrast_var),
        ".*",
        make.names(numerator),
        "$"
      ),
      coefficient_names,
      value = TRUE
    )
  ))

  coefficient_name <- coefficient_candidates[
    coefficient_candidates %in% coefficient_names
  ][1]

  if (
    length(coefficient_name) == 0L ||
    is.na(coefficient_name)
  ) {
    stop(
      paste0(
        "Could not identify the dream coefficient for ",
        numerator,
        " versus ",
        denominator,
        ". Available coefficients: ",
        paste(coefficient_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  comparison_name <- paste0(
    make.names(numerator),
    "_vs_",
    make.names(denominator)
  )

  result_table <- variancePartition::topTable(
    fit,
    coef = coefficient_name,
    number = Inf,
    sort.by = "P"
  )

  result_table <- as.data.frame(
    result_table,
    stringsAsFactors = FALSE
  )

  result_table$gene_id <- rownames(
    result_table
  )

  rownames(result_table) <- NULL

  # ---------------------------------------------------------- #
  # Standard DEGgo column names
  # ---------------------------------------------------------- #

  result_table$baseMean <- result_table$AveExpr
  result_table$log2FoldChange <- result_table$logFC
  result_table$lfcSE <- if (
    "t" %in% colnames(result_table)
  ) {
    abs(
      result_table$logFC /
        result_table$t
    )
  } else {
    NA_real_
  }

  result_table$stat <- result_table$t
  result_table$pvalue <- result_table$P.Value
  result_table$padj <- result_table$adj.P.Val

  result_table$significant <- with(
    result_table,
    !is.na(padj) &
      padj < padj_cutoff &
      abs(log2FoldChange) >= logfc_cutoff
  )

  result_table$regulation <- ifelse(
    result_table$significant &
      result_table$log2FoldChange > 0,
    "Up",
    ifelse(
      result_table$significant &
        result_table$log2FoldChange < 0,
      "Down",
      "NS"
    )
  )

  result_table$comparison <- comparison_name
  result_table$contrast_variable <- contrast_var
  result_table$contrast_numerator <- numerator
  result_table$contrast_denominator <- denominator
  result_table$method <- "dream"

  first_columns <- c(
    "gene_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "significant",
    "regulation",
    "comparison"
  )

  result_table <- result_table[
    ,
    c(
      intersect(first_columns, colnames(result_table)),
      setdiff(colnames(result_table), first_columns)
    ),
    drop = FALSE
  ]

  significant_table <- result_table[
    result_table$significant %in% TRUE,
    ,
    drop = FALSE
  ]

  # log2 CPM matrix for PCA, heatmaps and downstream plotting
  normalized_counts <- edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = 2
  )

  output <- list(
    results = stats::setNames(
      list(result_table),
      comparison_name
    ),
    sig_deg = stats::setNames(
      list(significant_table),
      comparison_name
    ),
    metadata = metadata,
    method = "dream",
    analysis_mode = "single",
    contrast = contrast,
    comparison = comparison_name,
    coefficient = coefficient_name,
    fit = fit,
    dream_fit = fit,
    voom = voom_object,
    voom_object = voom_object,
    dge = dge,
    normalized_counts = normalized_counts,
    transformed_counts = normalized_counts,

    # Compatibility fields. There is no DESeqDataSet for dream.
    dds = NULL,
    vst = NULL
  )

  if (is.function(log)) {
    log(
      paste0(
        "[DE] dream completed: ",
        nrow(result_table),
        " genes tested; ",
        nrow(significant_table),
        " significant"
      ),
      type = "done"
    )
  }

  output
}
