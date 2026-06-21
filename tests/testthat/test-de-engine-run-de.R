# tests/testthat/test-de-engine-run-de.R

test_that("run_de works with DESeq2 on minimal data", {
  skip_if_not_installed("DESeq2")
  
  set.seed(1)
  
  counts <- matrix(
    rnbinom(80, mu = 50, size = 1),
    nrow = 20,
    ncol = 4
  )
  
  rownames(counts) <- paste0("Gene", 1:20)
  colnames(counts) <- c("Ctrl_1", "Ctrl_2", "Treat_1", "Treat_2")
  
  metadata <- data.frame(
    sample = colnames(counts),
    condition = factor(c("Ctrl", "Ctrl", "Treat", "Treat")),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )
  
  out <- run_de(
    counts = counts,
    metadata = metadata,
    method = "DESeq2",
    design_formula = ~ condition,
    contrast = c("condition", "Treat", "Ctrl")
  )
  
  expect_true(is.list(out))
  expect_true(any(c("res_df", "results", "dds", "object") %in% names(out)))
})