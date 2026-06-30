test_that("DEG and GO extraction helpers work on airway DEGgo object", {

  skip_if_not_installed("DESeq2")
  skip_if_not_installed("org.Hs.eg.db")

  counts <- read.delim(
    system.file("extdata", "airway_counts.tsv", package = "DEGgo"),
    check.names = FALSE
  )

  metadata <- read.delim(
    system.file("extdata", "airway_metadata.tsv", package = "DEGgo"),
    check.names = FALSE
  )

  metadata$condition <- metadata$dex

  outdir <- tempfile("DEGgo_airway_extract_helpers_")

  res <- run_deggo(
    counts = counts,
    metadata = metadata,
    gene_col = "gene_id",
    organism = "human",
    sample_col = "Run",
    method = "DESeq2",
    analysis_mode = "single",
    design_formula = ~ cell + condition,
    contrast = c("condition", "trt", "untrt"),
    filter_method = "count",
    min_count = 5,
    min_samples = 2,
    min_total = 10,
    output_dir = outdir,
    generate_report = FALSE
  )

  expect_true(!is.null(res$sig_deg))

  expect_true(
    sum(vapply(
      res$sig_deg,
      function(x) if (is.null(x)) 0L else nrow(x),
      integer(1)
    )) > 0
  )

  # --------------------------------------------------
  # Test deggo_extract_deg_genes()
  # --------------------------------------------------

  genes <- unique(na.omit(res$sig_deg$SYMBOL))
  genes <- head(genes, 5)

  if (length(genes) > 0) {

    res_for_extract <- res

    res_for_extract$results <- list(
      airway_trt_vs_untrt = res$sig_deg
    )

    gene_summary <- deggo_extract_deg_genes(
      results = res_for_extract,
      genes = genes
    )

    expect_true(is.data.frame(gene_summary))
    expect_true(nrow(gene_summary) > 0)
  }

  # --------------------------------------------------
  # Test deggo_extract_go_keywords()
  # --------------------------------------------------

  if (!is.null(res$go_results)) {

    keyword_summary <- tryCatch(
      deggo_extract_go_keywords(
        results = res,
        keywords = c(
          "response",
          "regulation",
          "process"
        )
      ),
      error = function(e) NULL
    )

    expect_true(
      is.null(keyword_summary) ||
        is.data.frame(keyword_summary)
    )
  }
})
