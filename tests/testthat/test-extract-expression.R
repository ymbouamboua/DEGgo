test_that("extract_expression works with a minimal DESeq2 object", {
  skip_if_not_installed("DESeq2")
  
  set.seed(1)
  
  count_mat <- matrix(
    rnbinom(40, mu = 50, size = 1),
    nrow = 10,
    ncol = 4
  )
  
  rownames(count_mat) <- paste0("Gene", 1:10)
  colnames(count_mat) <- c("Ctrl_1", "Ctrl_2", "Treat_1", "Treat_2")
  
  metadata <- data.frame(
    sample = colnames(count_mat),
    condition = factor(c("Ctrl", "Ctrl", "Treat", "Treat")),
    row.names = colnames(count_mat),
    stringsAsFactors = FALSE
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = metadata,
    design = ~ condition
  )
  
  SummarizedExperiment::rowData(dds)$SYMBOL <- rownames(count_mat)
  
  expr <- extract_expression(
    dds = dds,
    metadata = metadata,
    genes = c("Gene1", "Gene2"),
    assay = "raw",
    gene_col = "SYMBOL"
  )
  
  expect_true(is.data.frame(expr))
  expect_true(all(c("gene", "sample", "expression") %in% colnames(expr)))
  expect_true(all(c("Gene1", "Gene2") %in% expr$gene))
})