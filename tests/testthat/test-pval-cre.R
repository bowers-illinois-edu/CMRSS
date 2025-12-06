# Test comparing CMRSS CRE results with RIQITE package
#
# This test requires the RIQITE package for comparison functions.
# Install from GitHub: li-xinran/RIQITE
#
# The test is skipped if RIQITE is not available.

test_that("pvalue in one stratum matches RIQITE", {
  # Skip if RIQITE is not available
  skip_if_not(
    requireNamespace("RIQITE", quietly = TRUE),
    "RIQITE package not available (install from GitHub: li-xinran/RIQITE)"
  )

  # Also need a solver for CMRSS
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  set.seed(54321) # Set seed for reproducibility

  n <- 50  # Fixed size for reproducibility
  m <- 25  # Number of treated units
  Z <- c(rep(1, m), rep(0, n - m))
  Z <- sample(Z)  # Shuffle

  Y <- rnorm(n)
  k <- floor(0.8 * n)
  c <- 0

  r <- 3  # Stephenson parameter
  block <- factor(rep(1, n))

  # RIQITE method list format
  method.list.riqite.wil <- list(name = "Wilcoxon")
  method.list.riqite.ste <- list(name = "Stephenson", s = r)

  # CMRSS method list format (needs scale parameter, one list per block)
  methods.list.cmrss.wil <- list(list(list(name = "Wilcoxon", scale = FALSE)))
  methods.list.cmrss.ste <- list(list(list(name = "Stephenson", scale = FALSE, s = r)))

  # Compare RIQITE::pval_quantile with CMRSS::pval_comb_block for Wilcoxon
  riqite_wil <- RIQITE::pval_quantile(
    Z = Z, Y = Y, k = k, c = c,
    method.list = method.list.riqite.wil
  )

  cmrss_wil <- pval_comb_block(
    Z = Z, Y = Y, k = k, c = c,
    block = block,
    methods.list.all = methods.list.cmrss.wil,
    statistic = FALSE
  )

  expect_equal(riqite_wil, cmrss_wil, tolerance = 0.01)

  # Compare RIQITE::pval_quantile with CMRSS::pval_comb_block for Stephenson
  riqite_ste <- RIQITE::pval_quantile(
    Z = Z, Y = Y, k = k, c = c,
    method.list = method.list.riqite.ste
  )

  cmrss_ste <- pval_comb_block(
    Z = Z, Y = Y, k = k, c = c,
    block = block,
    methods.list.all = methods.list.cmrss.ste,
    statistic = FALSE
  )

  expect_equal(riqite_ste, cmrss_ste, tolerance = 0.01)
})
