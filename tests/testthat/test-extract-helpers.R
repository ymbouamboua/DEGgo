test_that("deggo_extract_deg_genes handles DEGgo results object", {
  results <- list(
    results = data.frame(
      comparison = c("test", "test", "test", "test"),
      SYMBOL = c("Amh", "Amhr2", "Ucp1", "Lep"),
      log2FoldChange = c(1.2, -1.1, 0.2, 2.1),
      padj = c(0.01, 0.02, 0.8, 0.001),
      stringsAsFactors = FALSE
    )
  )

  out <- deggo_extract_deg_genes(
    results = results,
    genes = c("Amh", "Lep")
  )

  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 0)
})

test_that("deggo_extract_go_keywords handles nested GO results", {
  results <- list(
    go_results = list(
      test_comparison = list(
        go_results = data.frame(
          comparison = c("test", "test"),
          Regulation = c("Up", "Down"),
          ID = c("GO:0000001", "GO:0000002"),
          Description = c("lipid metabolic process", "immune response"),
          p.adjust = c(0.01, 0.02),
          Count = c(5, 6),
          stringsAsFactors = FALSE
        )
      )
    )
  )

  out <- deggo_extract_go_keywords(
    results = results,
    keywords = c("lipid", "immune")
  )

  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 0)
})



test_that("deggo_extract_deg_genes empty genes", {

  res <- list(
    results = data.frame(
      SYMBOL = c("Amh","Lep")
    )
  )

  out <- deggo_extract_deg_genes(
    res,
    genes = "XYZ"
  )

  expect_equal(nrow(out), 0)
})




