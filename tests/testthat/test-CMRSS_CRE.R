# Tests for CRE (Completely Randomized Experiment) functions

test_that("comb_p_val_cre returns valid p-value", {
  set.seed(12345)

  n <- 30
  m <- 15
  Z <- sample(c(rep(1, m), rep(0, n - m)))
  Y <- rnorm(n) + 0.5 * Z

  k <- floor(0.8 * n)
  c <- 0

  # Define methods list
  methods.list <- list(
    list(name = "Wilcoxon", scale = FALSE),
    list(name = "Stephenson", s = 3, scale = FALSE)
  )

  pval <- comb_p_val_cre(Z, Y, k, c, methods.list, nperm = 500)

  expect_true(is.numeric(pval))
  expect_true(pval >= 0 && pval <= 1)
})


test_that("com_conf_quant_larger_cre returns correct structure for set='treat'", {
  set.seed(23456)

  n <- 40
  m <- 20
  Z <- sample(c(rep(1, m), rep(0, n - m)))
  Y <- rnorm(n) + 0.3 * Z

  methods.list <- list(
    list(name = "Stephenson", s = 2, std = TRUE, scale = TRUE),
    list(name = "Stephenson", s = 6, std = TRUE, scale = TRUE)
  )

  ci <- com_conf_quant_larger_cre(Z, Y,
                                  methods.list = methods.list,
                                  nperm = 500,
                                  set = "treat",
                                  alpha = 0.10,
                                  tol = 0.1)

  # Should return m values (one for each treated unit's quantile)
  expect_length(ci, m)
  expect_true(all(is.numeric(ci)))
  # Lower bounds should be non-increasing (or -Inf)
  finite_ci <- ci[is.finite(ci)]
  if (length(finite_ci) > 1) {
    expect_true(all(diff(finite_ci) >= -0.001))  # Allow small numerical tolerance
  }
})


test_that("com_conf_quant_larger_cre returns correct structure for set='control'", {
  set.seed(34567)

  n <- 40
  m <- 20
  Z <- sample(c(rep(1, m), rep(0, n - m)))
  Y <- rnorm(n) + 0.3 * Z

  methods.list <- list(
    list(name = "Stephenson", s = 2, std = TRUE, scale = TRUE),
    list(name = "Stephenson", s = 6, std = TRUE, scale = TRUE)
  )

  ci <- com_conf_quant_larger_cre(Z, Y,
                                  methods.list = methods.list,
                                  nperm = 500,
                                  set = "control",
                                  alpha = 0.10,
                                  tol = 0.1)

  # Should return n-m values (one for each control unit's quantile)
  expect_length(ci, n - m)
  expect_true(all(is.numeric(ci)))
})


test_that("com_conf_quant_larger_cre returns correct structure for set='all'", {
  set.seed(45678)

  n <- 30
  m <- 15
  Z <- sample(c(rep(1, m), rep(0, n - m)))
  Y <- rnorm(n) + 0.3 * Z

  methods.list <- list(
    list(name = "Stephenson", s = 2, std = TRUE, scale = TRUE),
    list(name = "Stephenson", s = 6, std = TRUE, scale = TRUE)
  )

  ci <- com_conf_quant_larger_cre(Z, Y,
                                  methods.list = methods.list,
                                  nperm = 500,
                                  set = "all",
                                  alpha = 0.10,
                                  tol = 0.1)

  # Should return n values (one for each unit's quantile)
  expect_length(ci, n)
  expect_true(all(is.numeric(ci)))
})


test_that("CRE functions work with electric_teachers data", {
  data(electric_teachers)

  Z <- electric_teachers$TxAny
  Y <- electric_teachers$gain

  # Define Stephenson methods
  methods.list <- list(
    list(name = "Stephenson", s = 2, std = TRUE, scale = TRUE),
    list(name = "Stephenson", s = 6, std = TRUE, scale = TRUE)
  )

  # Test comb_p_val_cre
  n <- length(Z)
  k <- floor(0.9 * n)

  pval <- comb_p_val_cre(Z, Y, k, c = 0, methods.list, nperm = 500)
  expect_true(pval >= 0 && pval <= 1)
})
