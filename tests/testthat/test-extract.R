# tests/testthat/test-extract.R

test_that("extract_expression works on airway DEGgo results", {

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
    generate_report = FALSE
  )

  expr <- extract_expression(
    dds = res$dds,
    metadata = res$metadata,
    genes = head(rownames(res$dds), 3),
    assay = "vst",
    gene_col = "ENSEMBL"
  )

  expect_true(is.data.frame(expr))
  expect_true(nrow(expr) > 0)
})
