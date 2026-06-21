test_that("DEGgo runs on airway example data", {
  
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
  
  outdir <- tempfile("DEGgo_airway_test_")
  
  results <- run_deggo(
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
  
  expect_true(is.list(results))
  expect_true(!is.null(results$summary))
  expect_true(nrow(results$summary) > 0)
  expect_true(!is.null(results$sig_deg))
  expect_true(nrow(results$sig_deg) > 0)
  expect_true(file.exists(outdir))
})