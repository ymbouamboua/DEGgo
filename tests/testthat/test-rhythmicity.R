test_that("run_deggo_rhythmicity detects rhythmic vs non-rhythmic genes", {
  skip_if_not_installed("MetaCycle")
  skip_if_not_installed("cosinor2")
  skip_if_not_installed("cosinor")

  tc <- make_toy_timecourse()
  out_dir <- tempfile("deggo_rhythm_")

  res <- run_deggo_rhythmicity(
    expr = tc$mat,
    metadata = tc$metadata,
    time_col = "time",
    group_col = "group",
    methods = c("meta2d", "cosinor"),
    output_dir = out_dir,
    generate_plots = FALSE,
    verbose = FALSE
  )

  expect_s3_class(res, "deggo_rhythm_results")
  expect_true(is.data.frame(res$summary))
  expect_setequal(res$summary$gene, rownames(tc$mat))

  rhy_shared <- res$summary[res$summary$gene == "RHY_SHARED", ]
  flat <- res$summary[res$summary$gene == "FLAT", ]

  expect_true(rhy_shared$cosinor_padj < 0.05)
  expect_true(is.na(flat$cosinor_padj) || flat$cosinor_padj > rhy_shared$cosinor_padj)

  expect_true(file.exists(file.path(out_dir, "rhythmicity_results", "rhythmicity_summary.tsv")))
})


test_that("run_deggo_rhythmicity flags differential rhythmicity between groups", {
  skip_if_not_installed("MetaCycle")
  skip_if_not_installed("cosinor2")
  skip_if_not_installed("cosinor")

  tc <- make_toy_timecourse()
  out_dir <- tempfile("deggo_rhythm_")

  res <- run_deggo_rhythmicity(
    expr = tc$mat,
    metadata = tc$metadata,
    time_col = "time",
    group_col = "group",
    methods = "cosinor",
    output_dir = out_dir,
    generate_plots = FALSE,
    verbose = FALSE
  )

  expect_true(is.data.frame(res$cosinor_diff_group))

  diff_gene <- res$cosinor_diff_group[res$cosinor_diff_group$gene == "RHY_DIFF", ]
  shared_gene <- res$cosinor_diff_group[res$cosinor_diff_group$gene == "RHY_SHARED", ]

  expect_true(diff_gene$diff_rhythm_pvalue < shared_gene$diff_rhythm_pvalue)
  expect_true(diff_gene$diff_rhythm_pvalue < 0.05)
})


test_that("run_deggo_rhythmicity skips differential test with >2 groups", {
  skip_if_not_installed("MetaCycle")
  skip_if_not_installed("cosinor2")
  skip_if_not_installed("cosinor")

  tc <- make_toy_timecourse()
  tc$metadata$group <- rep(c("A", "B", "C", "D"), length.out = nrow(tc$metadata))
  out_dir <- tempfile("deggo_rhythm_")

  res <- run_deggo_rhythmicity(
    expr = tc$mat,
    metadata = tc$metadata,
    time_col = "time",
    group_col = "group",
    methods = "cosinor",
    output_dir = out_dir,
    generate_plots = FALSE,
    verbose = FALSE
  )

  expect_null(res$cosinor_diff_group)
})


test_that("deggo_extract_rhythmic_genes filters by gene", {
  skip_if_not_installed("MetaCycle")
  skip_if_not_installed("cosinor2")
  skip_if_not_installed("cosinor")

  tc <- make_toy_timecourse()
  out_dir <- tempfile("deggo_rhythm_")

  res <- run_deggo_rhythmicity(
    expr = tc$mat,
    metadata = tc$metadata,
    time_col = "time",
    methods = "cosinor",
    output_dir = out_dir,
    generate_plots = FALSE,
    verbose = FALSE
  )

  out <- deggo_extract_rhythmic_genes(res, genes = "RHY_SHARED")
  expect_equal(nrow(out), 1)
  expect_equal(out$gene, "RHY_SHARED")
})


test_that("manual lm() cosinor fallback agrees with the package engine", {
  skip_if_not_installed("MetaCycle")
  skip_if_not_installed("cosinor2")
  skip_if_not_installed("cosinor")

  tc <- make_toy_timecourse()

  res_pkg <- run_deggo_rhythmicity(
    expr = tc$mat, metadata = tc$metadata, time_col = "time", group_col = "group",
    methods = "cosinor", cosinor_engine = "package",
    output_dir = tempfile("deggo_rhythm_pkg_"), generate_plots = FALSE, verbose = FALSE
  )

  res_man <- run_deggo_rhythmicity(
    expr = tc$mat, metadata = tc$metadata, time_col = "time", group_col = "group",
    methods = "cosinor", cosinor_engine = "manual",
    output_dir = tempfile("deggo_rhythm_man_"), generate_plots = FALSE, verbose = FALSE
  )

  expect_equal(res_pkg$params$cosinor_engine, "package")
  expect_equal(res_man$params$cosinor_engine, "manual")
  expect_true(all(res_pkg$cosinor$cosinor_engine == "cosinor2"))
  expect_true(all(res_man$cosinor$cosinor_engine == "manual_lm"))

  pkg <- res_pkg$cosinor[order(res_pkg$cosinor$gene), ]
  man <- res_man$cosinor[order(res_man$cosinor$gene), ]

  # Mesor, amplitude magnitude and p-values should match exactly (same
  # underlying lm() fit); acrophase may differ by pi (atan vs atan2), so we
  # compare cos(acrophase) which is invariant to that ambiguity only when
  # amplitude sign is also reconciled -- instead just check numeric equality
  # allowing for the documented +/- pi flip.
  expect_equal(pkg$cosinor_mesor, man$cosinor_mesor, tolerance = 1e-8)
  expect_equal(pkg$cosinor_amplitude, man$cosinor_amplitude, tolerance = 1e-8)
  expect_equal(pkg$cosinor_pvalue, man$cosinor_pvalue, tolerance = 1e-8)

  phase_diff <- abs(pkg$cosinor_acrophase - man$cosinor_acrophase) %% (2 * pi)
  expect_true(all(phase_diff < 1e-6 | abs(phase_diff - pi) < 1e-6))

  # Differential rhythmicity p-values should match exactly between engines
  diff_pkg <- res_pkg$cosinor_diff_group[order(res_pkg$cosinor_diff_group$gene), ]
  diff_man <- res_man$cosinor_diff_group[order(res_man$cosinor_diff_group$gene), ]
  expect_equal(diff_pkg$diff_rhythm_pvalue, diff_man$diff_rhythm_pvalue, tolerance = 1e-8)
})
