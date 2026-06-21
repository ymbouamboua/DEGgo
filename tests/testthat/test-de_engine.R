test_that("minimal DE engine input objects are valid", {
  set.seed(1)

  counts <- matrix(
    rpois(40, 20),
    nrow = 10,
    ncol = 4
  )

  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- c("S1", "S2", "S3", "S4")

  metadata <- data.frame(
    sample = colnames(counts),
    condition = factor(c("A", "A", "B", "B")),
    stringsAsFactors = FALSE
  )

  expect_true(is.matrix(counts))
  expect_equal(ncol(counts), nrow(metadata))
  expect_equal(colnames(counts), metadata$sample)
})
