# tests/testthat/test-deggo-derived-objects.R

test_that("DEGgo downstream functions work from run_deggo output", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("ggplot2")
  
  set.seed(1)
  
  counts <- data.frame(
    gene_id = paste0("ENSMUSG", sprintf("%011d", 1:40)),
    gene_name = paste0("Gene", 1:40),
    Ctrl_1 = rnbinom(40, mu = 50, size = 1),
    Ctrl_2 = rnbinom(40, mu = 55, size = 1),
    Treat_1 = rnbinom(40, mu = 80, size = 1),
    Treat_2 = rnbinom(40, mu = 85, size = 1),
    check.names = FALSE
  )
  
  metadata <- data.frame(
    sample = c("Ctrl_1", "Ctrl_2", "Treat_1", "Treat_2"),
    condition = factor(c("Ctrl", "Ctrl", "Treat", "Treat")),
    treatment = factor(c("Ctrl", "Ctrl", "Treat", "Treat")),
    stringsAsFactors = FALSE
  )
  
  res <- run_deggo(
    counts = counts,
    metadata = metadata,
    organism = "mouse",
    method = "DESeq2",
    analysis_mode = "single",
    contrast = c("condition", "Treat", "Ctrl"),
    design_formula = ~ condition,
    filter_method = "none",
    generate_report = FALSE,
    output_dir = tempdir()
  )
  
  expect_true(is.list(res))
  expect_true("summary" %in% names(res))
  expect_true(is.data.frame(res$summary))
  
  if (!is.null(res$results)) {
    p <- plot_volcano(res$results, title = "Test volcano")
    expect_s3_class(p, "ggplot")
  }
  
  if (!is.null(res$dds)) {
    expr <- extract_expression(
      dds = res$dds,
      metadata = res$metadata,
      genes = head(SummarizedExperiment::rowData(res$dds)$SYMBOL, 2),
      assay = "raw",
      gene_col = "SYMBOL"
    )
    
    expect_true(is.data.frame(expr))
  }
  
  res_extract <- res
  res_extract$results <- data.frame(
    comparison = "single_comparison",
    SYMBOL = c("Gene1", "Gene2", "Gene3"),
    log2FoldChange = c(1.2, -1.3, 0.1),
    padj = c(0.01, 0.02, 0.9),
    stringsAsFactors = FALSE
  )
  
  out <- deggo_extract_deg_genes(
    results = res_extract,
    genes = c("Gene1", "Gene2")
  )
  
  expect_true(is.data.frame(out))
})