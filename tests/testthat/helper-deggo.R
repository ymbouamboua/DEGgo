make_toy_counts <- function() {
  counts <- data.frame(
    gene_id = c(
      "ENSG00000152583.1",
      "ENSG00000120129.2",
      "ENSG00000171862.3",
      "ENSG00000141510.4",
      "ENSG00000155657.5",
      "ENSG00000157764.6"
    ),
    gene_name = c("SPARCL1", "DUSP1", "PTEN", "TP53", "TTN", "BRAF"),
    S1 = c(100, 90, 80, 40, 20, 10),
    S2 = c(110, 95, 78, 45, 25, 12),
    S3 = c(500, 400, 75, 42, 22, 11),
    S4 = c(520, 420, 74, 44, 23, 13),
    check.names = FALSE
  )
  counts
}

make_toy_metadata <- function() {
  data.frame(
    sample = c("S1", "S2", "S3", "S4"),
    condition = c("ctrl", "ctrl", "trt", "trt"),
    batch = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
}