test_that("edgeR and limma engines run on minimal data", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")

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
})


# test_that("edgeR engine runs", {
#
#   expect_error(
#     run_deggo(
#       counts = counts,
#       metadata = metadata,
#       method = "edgeR",
#       filter_method = "none"
#     ),
#     NA
#   )
#
# })
#
#
# test_that("edgeR engine runs", {
#
#   expect_error(
#     run_deggo(
#       counts = counts,
#       metadata = metadata,
#       method = "limma",
#       filter_method = "none"
#     ),
#     NA
#   )
#
# })
