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
make_toy_timecourse <- function(seed = 8542) {
  set.seed(seed)

  time <- rep(seq(0, 44, by = 4), 2)
  group <- rep(c("A", "B"), each = length(time) / 2)
  n <- length(time)

  # Rhythmic gene, same phase/amplitude in both groups
  rhy_shared <- 10 + 3 * cos(2 * pi * time / 24 - 1) + rnorm(n, sd = 0.3)

  # Rhythmic gene with a phase shift between groups (differential rhythmicity)
  acro <- ifelse(group == "A", -pi / 2, pi / 2)
  rhy_diff <- 10 + 3 * cos(2 * pi * time / 24 + acro) + rnorm(n, sd = 0.3)

  # Flat, non-rhythmic gene
  flat <- 10 + rnorm(n, sd = 0.3)

  mat <- rbind(rhy_shared, rhy_diff, flat)
  rownames(mat) <- c("RHY_SHARED", "RHY_DIFF", "FLAT")
  colnames(mat) <- paste0("S", seq_len(n))

  metadata <- data.frame(
    sample = colnames(mat),
    time = time,
    group = group,
    stringsAsFactors = FALSE
  )

  list(mat = mat, metadata = metadata)
}
