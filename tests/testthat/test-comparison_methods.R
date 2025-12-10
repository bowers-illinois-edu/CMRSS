# Tests for comparison methods (method_caughey, method_chen_li, method_berger_boos)

test_that("method_caughey returns valid confidence intervals", {
  skip_if_not(
    requireNamespace("RIQITE", quietly = TRUE),
    "RIQITE package not available (install from GitHub: li-xinran/RIQITE)"
  )

  set.seed(12345)

  n <- 40
  m <- 20
  Z <- sample(c(rep(1, m), rep(0, n - m)))
  Y <- rnorm(n) + 0.3 * Z

  k_vec <- floor(c(0.7, 0.8, 0.9) * n)

  # Test with Wilcoxon
  result_wil <- method_caughey(Z, Y, k_vec,
                               method.list = list(name = "Wilcoxon"),
                               nperm = 500,
                               alpha = 0.10)

  expect_s3_class(result_wil, "data.frame")
  expect_true("lower" %in% names(result_wil))
  expect_true("upper" %in% names(result_wil))
  expect_length(result_wil$lower, length(k_vec))

  # Test with Stephenson
  result_ste <- method_caughey(Z, Y, k_vec,
                               method.list = list(name = "Stephenson", s = 3),
                               nperm = 500,
                               alpha = 0.10)

  expect_s3_class(result_ste, "data.frame")
  expect_length(result_ste$lower, length(k_vec))
})


test_that("method_chen_li returns valid confidence intervals", {
  skip_if_not(
    requireNamespace("RIQITE", quietly = TRUE),
    "RIQITE package not available (install from GitHub: li-xinran/RIQITE)"
  )

  skip_if_not(
    requireNamespace("extraDistr", quietly = TRUE),
    "extraDistr package not available"
  )

  set.seed(23456)

  n <- 30
  m <- 15
  Z <- sample(c(rep(1, m), rep(0, n - m)))
  Y <- rnorm(n) + 0.5 * Z

  k_vec <- floor(c(0.6, 0.8) * n)

  # Test sample inference (N = n)
  result <- method_chen_li(Z, Y,
                           N = n,
                           k_vec = k_vec,
                           treat.method.list = list(name = "Stephenson", s = 3),
                           control.method.list = list(name = "Stephenson", s = 3),
                           nperm = 500,
                           ndraw = 100,
                           alpha = 0.10)

  expect_s3_class(result, "data.frame")
  expect_true("k" %in% names(result))
  expect_true("lower" %in% names(result))
  expect_length(result$lower, length(k_vec))
})


test_that("method_berger_boos returns valid confidence intervals", {
  skip_if_not(
    requireNamespace("RIQITE", quietly = TRUE),
    "RIQITE package not available (install from GitHub: li-xinran/RIQITE)"
  )

  skip_if_not(
    requireNamespace("extraDistr", quietly = TRUE),
    "extraDistr package not available"
  )

  set.seed(34567)

  n <- 30
  m <- 15
  Z <- sample(c(rep(1, m), rep(0, n - m)))
  Y <- rnorm(n) + 0.4 * Z

  k_vec <- floor(c(0.7, 0.9) * n)

  result <- method_berger_boos(Z, Y,
                                N = n,
                                k_vec = k_vec,
                                treat.method.list = list(name = "Stephenson", s = 3),
                                control.method.list = list(name = "Stephenson", s = 3),
                                nperm = 500,
                                ndraw = 100,
                                alpha = 0.10)

  expect_s3_class(result, "data.frame")
  expect_true("k" %in% names(result))
  expect_true("lower" %in% names(result))
  expect_length(result$lower, length(k_vec))
})


test_that("comparison methods work with electric_teachers data", {
  skip_if_not(
    requireNamespace("RIQITE", quietly = TRUE),
    "RIQITE package not available (install from GitHub: li-xinran/RIQITE)"
  )

  data(electric_teachers)

  Z <- electric_teachers$TxAny
  Y <- electric_teachers$gain
  n <- length(Z)

  k_vec <- floor(c(0.7, 0.8, 0.9) * n)

  # Test method_caughey with real data
  result <- method_caughey(Z, Y, k_vec,
                           method.list = list(name = "Stephenson", s = 6),
                           nperm = 500,
                           alpha = 0.10)

  expect_s3_class(result, "data.frame")
  expect_length(result$lower, length(k_vec))
  expect_true(all(result$lower <= result$upper))
})


test_that("rank_score_riqite produces max-normalized scores", {
  n <- 50

  # Wilcoxon should be 1:n divided by n
  result_wil <- rank_score_riqite(n, list(name = "Wilcoxon"))
  expected_wil <- c(1:n) / n
  expect_equal(result_wil, expected_wil)
  expect_equal(max(result_wil), 1)

  # Stephenson should be choose divided by max
  s <- 4
  result_ste <- rank_score_riqite(n, list(name = "Stephenson", s = s))
  raw <- choose(c(1:n) - 1, s - 1)
  expected_ste <- raw / max(raw)
  expect_equal(result_ste, expected_ste)
  expect_equal(max(result_ste), 1)
})
