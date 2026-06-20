test_that("generate_deggo_report errors for unsupported format", {
  results <- list(
    output_dir = tempdir(),
    summary = data.frame(
      comparison = "test",
      total_genes = 10,
      significant = 0,
      up = 0,
      down = 0,
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(
      sample = c("S1", "S2"),
      condition = c("Ctrl", "Treat"),
      stringsAsFactors = FALSE
    ),
    run_params = list(
      analysis_mode = "single",
      method = "DESeq2",
      organism = "mouse"
    ),
    version = "0.1.0"
  )

  expect_error(
    generate_deggo_report(
      results = results,
      output_dir = tempdir(),
      formats = "docx"
    )
  )
})


test_that("generate_deggo_report accepts html format", {

  expect_error(
    generate_deggo_report(
      results = list(),
      output_dir = tempdir(),
      formats = "html"
    ),
    NA
  )

})


test_that("generate_deggo_report accepts pdf format", {

  expect_error(
    generate_deggo_report(
      results = list(),
      output_dir = tempdir(),
      formats = "pdf"
    ),
    NA
  )

})


test_that("generate_deggo_report accepts multiple formats", {

  expect_error(
    generate_deggo_report(
      results = list(),
      output_dir = tempdir(),
      formats = c("html","pdf")
    ),
    NA
  )

})

