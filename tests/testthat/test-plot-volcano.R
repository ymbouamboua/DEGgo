test_that("plot_volcano returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  
  res_df <- data.frame(
    SYMBOL = paste0("Gene", 1:20),
    log2FoldChange = c(rnorm(10, 1), rnorm(10, -1)),
    padj = c(rep(0.001, 5), rep(0.2, 15)),
    stringsAsFactors = FALSE
  )
  
  p <- plot_volcano(
    res_df,
    title = "Test volcano",
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25
  )
  
  expect_s3_class(p, "ggplot")
})