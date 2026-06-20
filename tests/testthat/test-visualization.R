test_that("plot_gene_expression returns ggplot", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    gene = rep("Lep", 4),
    sample = paste0("S", 1:4),
    expression = c(1, 2, 3, 4),
    treatment = c("A", "A", "B", "B"),
    assay = rep("vst", 4),
    stringsAsFactors = FALSE
  )

  p <- plot_gene_expression(
    df,
    gene = "Lep",
    x = "treatment",
    color = "treatment"
  )

  expect_s3_class(p, "ggplot")
})
