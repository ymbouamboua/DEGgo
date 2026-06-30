test_that("airway run provides reusable DEGgo objects", {

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

  outdir <- tempfile("DEGgo_airway_derived_")

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
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25,
    output_dir = outdir,
    generate_report = FALSE
  )

  expect_true(is.list(res))
  expect_true(!is.null(res$dds))
  expect_true(!is.null(res$metadata))
  expect_true(!is.null(res$sig_deg))
  expect_true(sum(vapply(res$sig_deg, nrow, numeric(1))) > 0)
})
