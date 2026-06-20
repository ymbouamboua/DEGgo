fake_go <- data.frame(
  Description = c(
    "lipid metabolism",
    "immune response",
    "thermogenesis"
  ),
  p.adjust = c(0.01,0.02,0.03)
)

fake_results <- list(
  go_results = list(
    comparison1 = list(
      go_results = fake_go
    )
  )
)