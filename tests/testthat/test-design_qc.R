test_that("design_qc detects repeated cage structure", {
  skip_if_not_installed("DESeq2")

  set.seed(1)
  counts <- matrix(
    rnbinom(100 * 8, mu = 100, size = 10),
    nrow = 100,
    dimnames = list(paste0("g", 1:100), paste0("s", 1:8))
  )
  metadata <- data.frame(
    sample = paste0("s", 1:8),
    Group = rep(c("CT", "ABX"), each = 4),
    Cage = rep(c("C1", "C2", "C3", "C4"), each = 2)
  )

  res <- design_qc(
    counts,
    metadata,
    sample_col = "sample",
    group_col = "Group",
    unit_col = "Cage",
    transform = "vst",
    run_permanova = FALSE,
    verbose = FALSE
  )

  expect_s3_class(res, "deggo_design_qc")
  expect_true(res$structure$has_repeated_units)
  expect_false(res$structure$unit_crosses_group)
  expect_equal(res$recommendation$method, "dream")
})

test_that("aggregation sums counts within unit", {
  counts <- matrix(
    1:24,
    nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("s", 1:8))
  )
  metadata <- data.frame(
    sample = paste0("s", 1:8),
    Group = rep(c("CT", "ABX"), each = 4),
    Cage = rep(c("C1", "C2", "C3", "C4"), each = 2)
  )

  agg <- aggregate_counts_by_unit(
    counts,
    metadata,
    sample_col = "sample",
    unit_col = "Cage",
    group_col = "Group"
  )

  expect_equal(ncol(agg$counts), 4)
  expect_equal(agg$counts[, "C1"], rowSums(counts[, 1:2, drop = FALSE]))
  expect_equal(as.character(agg$metadata["C1", "Group"]), "CT")
})
