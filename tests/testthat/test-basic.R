# test-basic

test_that("DEGgo package loads", {
  expect_true(requireNamespace("DEGgo", quietly = TRUE))
})

test_that("run_deggo function exists", {
  expect_true(exists("run_deggo", where = asNamespace("DEGgo")))
})

test_that(".deggo_theme function exists", {
  expect_true(exists(".deggo_theme", where = asNamespace("DEGgo")))
})