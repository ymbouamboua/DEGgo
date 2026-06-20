test_that("run_go_enrichment handles too few genes without error", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Mm.eg.db")
  
  sig_deg <- data.frame(
    SYMBOL = c("Sox9", "Amh", "Star", "Cyp11a1", "Cyp17a1"),
    ENTREZID = c("20682", "11705", "20845", "13070", "13074"),
    log2FoldChange = c(1.2, 1.1, -1.3, 1.5, -1.1),
    padj = c(0.01, 0.02, 0.03, 0.04, 0.05),
    stringsAsFactors = FALSE
  )
  
  expect_no_error(
    go <- run_go_enrichment(
      sig_deg = sig_deg,
      comparison = "minimal_test",
      ontology = "BP",
      orgdb = org.Mm.eg.db::org.Mm.eg.db,
      output_dir = tempdir()
    )
  )
  
  expect_null(go)
})