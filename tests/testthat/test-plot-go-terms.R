test_that("plot_go_terms returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  
  go_df <- data.frame(
    ID = paste0("GO:", sprintf("%07d", 1:10)),
    Description = paste("Biological process", 1:10),
    GeneRatio = paste0(1:10, "/", 100),
    BgRatio = paste0(10:19, "/", 1000),
    FoldEnrichment = seq(1.5, 4, length.out = 10),
    p.adjust = seq(0.001, 0.05, length.out = 10),
    qvalue = seq(0.001, 0.05, length.out = 10),
    Count = 1:10,
    Regulation = rep(c("Up", "Down"), 5),
    stringsAsFactors = FALSE
  )
  
  p <- plot_go_terms(
    go_df = go_df,
    comparison = "minimal_test",
    top_n = 5
  )
  
  expect_s3_class(p, "ggplot")
})