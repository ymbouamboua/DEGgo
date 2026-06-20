library(testthat)
library(DEGgo)

# ---------------------------------------------------------
# Minimal dataset
# ---------------------------------------------------------

counts <- data.frame(
  gene_id = paste0("ENSMUSG", 1:10),
  gene_name = paste0("Gene", 1:10),
  Sample1 = c(100,50,20,10,5,0,1,2,3,4),
  Sample2 = c(120,60,25,15,8,0,1,2,3,5),
  check.names = FALSE
)

metadata <- data.frame(
  sample = c("Sample1","Sample2"),
  tissue = c("BAT","BAT"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------
# check_raw_counts()
# ---------------------------------------------------------

test_that("check_raw_counts runs without error", {
  expect_no_error(
    check_raw_counts(
      counts = counts,
      gene_col = "gene_id",
      feature_col = "gene_name"
    )
  )
})

# ---------------------------------------------------------
# remove_flagged_samples()
# ---------------------------------------------------------

test_that("remove_flagged_samples removes failed samples", {

  qc_table <- data.frame(
    sample = c("Sample1","Sample2"),
    recommend_remove = c(FALSE, TRUE)
  )

  res <- remove_flagged_samples(
    counts = counts,
    metadata = metadata,
    qc_table = qc_table,
    sample_col = "sample",
    remove_col = "recommend_remove",
    gene_cols = c("gene_id","gene_name")
  )

  expect_true(is.list(res))
  expect_true("counts" %in% names(res))
  expect_true("metadata" %in% names(res))

  expect_equal(nrow(res$metadata), 1)
  expect_true("Sample1" %in% res$metadata$sample)

})

# ---------------------------------------------------------
# marker_score_check()
# ---------------------------------------------------------

test_that("marker_score_check returns expected outputs", {

  counts2 <- data.frame(
    gene_id = c("ENSMUSG1","ENSMUSG2"),
    gene_name = c("Ucp1","Adipoq"),
    Sample1 = c(100,5),
    Sample2 = c(80,3),
    check.names = FALSE
  )

  marker_sets <- list(
    BAT = c("Ucp1"),
    WAT = c("Adipoq")
  )

  out <- marker_score_check(
    counts = counts2,
    metadata = metadata,
    marker_sets = marker_sets,
    sample_col = "sample",
    group_col = "tissue",
    feature_col = "gene_name"
  )

  expect_true(is.list(out))

  expect_true("scores" %in% names(out))
  expect_true("swaps" %in% names(out))

  expect_true(is.data.frame(out$scores))

})

# ---------------------------------------------------------
# run_sample_qc()
# ---------------------------------------------------------

test_that("run_sample_qc errors when fewer than two samples are matched", {
  counts_qc <- data.frame(
    gene_id = paste0("ENSMUSG", 1:20),
    gene_name = paste0("Gene", 1:20),
    Sample1 = rpois(20, 50),
    Sample2 = rpois(20, 60),
    Sample3 = rpois(20, 55),
    check.names = FALSE
  )

  metadata_qc <- data.frame(
    sample = c("Sample1", "Sample2", "Sample3"),
    condition = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    run_sample_qc(
      counts = counts_qc,
      metadata = metadata_qc,
      output_dir = tempdir()
    ),
    "At least two matched samples are required"
  )
})

