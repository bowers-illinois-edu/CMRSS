# Test CMRSS stratified randomized experiment (SCRE) functionality
#
# This test verifies that pval_comb_block works correctly for stratified
# randomized experiments with multiple blocks.

test_that("pval_comb_block works for stratified experiments", {
  # Need a solver for this test
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  set.seed(12345) # Set seed for reproducibility

  s <- 3    # number of strata
  n <- 10   # units per stratum
  m <- 5    # treated per stratum
  N <- s * n
  k <- floor(0.8 * N)
  c <- 0

  # Create block assignments and treatment vector
  Z_block <- factor(rep(1:s, each = n))
  Z <- rep(0, N)
  for (i in 1:s) {
    block_idx <- which(Z_block == i)
    Z[sample(block_idx, m)] <- 1
  }

  # Generate outcome data
  Y0 <- rnorm(N)
  Y1 <- Y0 + 0.5  # Small treatment effect
  Y <- Z * Y1 + (1 - Z) * Y0

  # Set up methods list for CMRSS (one method per block)
  methods.list.all <- list()
  for (i in 1:s) {
    methods.list.all[[i]] <- list(name = "Wilcoxon", scale = FALSE)
  }

  # Run pval_comb_block
  result <- pval_comb_block(
    Z = Z, Y = Y, k = k, c = c,
    block = Z_block,
    methods.list.all = list(methods.list.all),
    statistic = TRUE
  )

  # Check that result has expected structure

  expect_true("p.value" %in% names(result))
  expect_true("test.stat" %in% names(result))

  # Check that p-value is valid (between 0 and 1)
  expect_true(result["p.value"] >= 0 && result["p.value"] <= 1)
})

test_that("SCRE with single block matches CRE behavior", {
  # Skip if RIQITE is not available
  skip_if_not(
    requireNamespace("RIQITE", quietly = TRUE),
    "RIQITE package not available (install from GitHub: li-xinran/RIQITE)"
  )

  # Need a solver for this test
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  set.seed(99999)

  n <- 40
  m <- 20
  k <- floor(0.8 * n)
  c <- 0

  Z <- c(rep(1, m), rep(0, n - m))
  Z <- sample(Z)
  Y <- rnorm(n)
  block <- factor(rep(1, n))

  # RIQITE for CRE
  riqite_result <- RIQITE::pval_quantile(
    Z = Z, Y = Y, k = k, c = c,
    method.list = list(name = "Wilcoxon")
  )

  # CMRSS with single block (should behave like CRE)
  cmrss_result <- pval_comb_block(
    Z = Z, Y = Y, k = k, c = c,
    block = block,
    methods.list.all = list(list(list(name = "Wilcoxon", scale = FALSE))),
    statistic = FALSE
  )

  # Results should be similar
  expect_equal(riqite_result, cmrss_result, tolerance = 0.02)
})
