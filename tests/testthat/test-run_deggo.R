test_that("run_deggo works on a minimal DESeq2 dataset", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("org.Mm.eg.db")

  set.seed(1)

  counts <- data.frame(
    gene_id = paste0("ENSMUSG", sprintf("%011d", 1:30)),
    gene_name = paste0("Gene", 1:30),
    Ctrl_1 = rnbinom(30, mu = 50, size = 1),
    Ctrl_2 = rnbinom(30, mu = 55, size = 1),
    Treat_1 = rnbinom(30, mu = 80, size = 1),
    Treat_2 = rnbinom(30, mu = 85, size = 1),
    check.names = FALSE
  )

  metadata <- data.frame(
    sample = c("Ctrl_1", "Ctrl_2", "Treat_1", "Treat_2"),
    condition = c("Ctrl", "Ctrl", "Treat", "Treat"),
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
  expect_true("res_df" %in% names(res) || "deg_results" %in% names(res))
  expect_true(is.data.frame(res$summary))
})
