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
  k <- floor(0.8 * s * m)
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

# CMRSS tests H_{k,c}^treat (treated-only) with k in 1..sum(Z); RIQITE tests
# H_{k,c}     (all-units)    with k in 1..n.  These are different hypotheses,
# but their LP minima for the rank-sum statistic coincide under the translation
#
#     k_R (RIQITE, all-units) = k_C (CMRSS, treated-only) + (n - m).
#
# Reason: in the LP, xi[i] only enters the statistic through Z[i] * xi[i], so
# controls' xi is unconstrained.  Assigning all (n - m) controls to xi = c
# absorbs (n - m) of RIQITE's all-units constraint slack, leaving exactly
# k_C = k_R - (n - m) treated units constrained.  With a shared Z.perm matrix
# the two packages return identical p-values.

test_that("SCRE single block matches RIQITE CRE under shared null and k translation", {
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
  k_C   <- floor(0.8 * m)        # treated-only k for CMRSS, = 16
  k_R   <- k_C + (n - m)         # equivalent all-units k for RIQITE, = 36
  c_val <- 0

  Z <- c(rep(1, m), rep(0, n - m))
  Z <- sample(Z)
  Y <- rnorm(n)
  block <- factor(rep(1, n))

  nperm <- 1000
  set.seed(11111)
  Z.perm <- RIQITE::assign_CRE(n, m, nperm)

  method.list.riqite <- list(name = "Wilcoxon")
  methods.list.cmrss <- list(list(list(name = "Wilcoxon", scale = FALSE)))

  block.sum <- summary_block(Z, block)
  weight    <- weight_scheme(block.sum, "asymp.opt")
  scores    <- list(score_all_blocks(block.sum$nb, methods.list.cmrss[[1]]))

  riq_null <- RIQITE::null_dist(n, m, nperm = nperm, Z.perm = Z.perm,
                                method.list = method.list.riqite)
  cms_null <- com_null_dist_block(Z, block, methods.list.cmrss, scores,
                                  null.max = nperm, weight, block.sum,
                                  Z.perm = Z.perm)

  riqite_result <- RIQITE::pval_quantile(
    Z = Z, Y = Y, k = k_R, c = c_val,
    method.list = method.list.riqite,
    stat.null = riq_null, switch = FALSE
  )

  cmrss_result <- pval_comb_block(
    Z = Z, Y = Y, k = k_C, c = c_val,
    block = block,
    methods.list.all = methods.list.cmrss,
    stat.null = cms_null, statistic = FALSE
  )

  expect_equal(as.numeric(riqite_result), as.numeric(cmrss_result),
               tolerance = 1e-6)
})
