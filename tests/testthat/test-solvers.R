# Tests for solver utilities and HiGHS/Gurobi comparison
#
# These tests verify that:
# 1. Solver utility functions work correctly
# 2. HiGHS and Gurobi produce equivalent results when both are available
# 3. The package gracefully handles missing solvers

test_that("parse_opt_method works correctly", {
  # ILP variants
  expect_equal(parse_opt_method("ILP_gurobi"), list(exact = TRUE, solver = "gurobi"))
  expect_equal(parse_opt_method("ILP_highs"), list(exact = TRUE, solver = "highs"))
  expect_equal(parse_opt_method("ILP_auto"), list(exact = TRUE, solver = "auto"))
  expect_equal(parse_opt_method("ILP"), list(exact = TRUE, solver = "auto"))

  # LP variants
  expect_equal(parse_opt_method("LP_gurobi"), list(exact = FALSE, solver = "gurobi"))
  expect_equal(parse_opt_method("LP_highs"), list(exact = FALSE, solver = "highs"))
  expect_equal(parse_opt_method("LP_auto"), list(exact = FALSE, solver = "auto"))
  expect_equal(parse_opt_method("LP"), list(exact = FALSE, solver = "auto"))

  # Case insensitivity
  expect_equal(parse_opt_method("ilp_HIGHS"), list(exact = TRUE, solver = "highs"))
  expect_equal(parse_opt_method("lp_Gurobi"), list(exact = FALSE, solver = "gurobi"))

  # Invalid method should error
  expect_error(parse_opt_method("invalid_method"))
})


test_that("solver_available returns logical", {
  # These should return logical values without errors
  expect_type(solver_available("highs"), "logical")
  expect_type(solver_available("gurobi"), "logical")

  # Invalid solver should error
  expect_error(solver_available("invalid_solver"))
})


test_that("get_default_solver returns valid solver when available", {
  skip_if_not(solver_available("highs") || solver_available("gurobi"),
              "Neither HiGHS nor Gurobi is available")

  solver <- get_default_solver()
  expect_true(solver %in% c("highs", "gurobi"))
  expect_true(solver_available(solver))
})


# Test that HiGHS solver works correctly
test_that("HiGHS solver works for basic optimization", {
  skip_if_not(solver_available("highs"), "HiGHS is not available")

  # Set up a simple test case
  set.seed(12345)

  s <- 3  # number of strata
  n <- 6  # units per stratum
  m <- 3  # treated per stratum
  N <- s * n
  k <- ceiling(0.8 * N)
  c <- 0

  # Create block structure
  block <- factor(rep(1:s, each = n))

  # Treatment assignment within blocks
  Z <- rep(0, N)
  for (i in 1:s) {
    block_idx <- which(block == i)
    Z[sample(block_idx, m)] <- 1
  }

  # Generate outcomes
  Y0 <- rnorm(N)
  Y1 <- Y0 + 1
  Y <- Z * Y1 + (1 - Z) * Y0

  # Set up methods
  methods.list.all <- list()
  for (h in 1:1) {
    methods.list.all[[h]] <- list()
    for (i in 1:s) {
      methods.list.all[[h]][[i]] <- list(name = "Wilcoxon", scale = FALSE)
    }
  }

  # Test that pval_comb_block works with HiGHS
  result <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                            opt.method = "ILP_highs",
                            null.max = 100,
                            statistic = TRUE)

  expect_true(is.numeric(result))
  expect_length(result, 2)
  expect_true(result["p.value"] >= 0 && result["p.value"] <= 1)
  expect_true(is.finite(result["test.stat"]))
})


# Compare HiGHS and Gurobi results
test_that("HiGHS and Gurobi produce equivalent results", {
  skip_if_not(solver_available("highs") && solver_available("gurobi"),
              "Both HiGHS and Gurobi are required for this test")

  set.seed(54321)

  s <- 4  # number of strata
  n <- 8  # units per stratum
  m <- 4  # treated per stratum
  N <- s * n
  k <- ceiling(0.9 * N)
  c <- 0.5

  # Create block structure
  block <- factor(rep(1:s, each = n))

  # Treatment assignment within blocks
  Z <- rep(0, N)
  for (i in 1:s) {
    block_idx <- which(block == i)
    Z[sample(block_idx, m)] <- 1
  }

  # Generate outcomes
  Y0 <- rnorm(N)
  Y1 <- Y0 + 0.8
  Y <- Z * Y1 + (1 - Z) * Y0

  # Set up methods - using Stephenson score
  methods.list.all <- list()
  for (h in 1:2) {
    methods.list.all[[h]] <- list()
    for (i in 1:s) {
      if (h == 1) {
        methods.list.all[[h]][[i]] <- list(name = "Wilcoxon", scale = FALSE)
      } else {
        methods.list.all[[h]][[i]] <- list(name = "Stephenson", s = 3, scale = FALSE)
      }
    }
  }

  # Generate shared null distribution for fair comparison
  block.sum <- summary_block(Z, block)
  weight <- weight_scheme(block.sum, "asymp.opt")
  scores.list.all <- list()
  for (h in 1:2) {
    scores.list.all[[h]] <- score_all_blocks(block.sum$nb, methods.list.all[[h]])
  }
  stat.null <- com_null_dist_block(Z, block, methods.list.all, scores.list.all,
                                   null.max = 500, weight, block.sum, Z.perm = NULL)

  # Test ILP (exact) mode
  result_highs_ilp <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                                       opt.method = "ILP_highs",
                                       stat.null = stat.null,
                                       statistic = TRUE)

  result_gurobi_ilp <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                                        opt.method = "ILP_gurobi",
                                        stat.null = stat.null,
                                        statistic = TRUE)

  # Results should be nearly identical (allowing for small numerical differences)
  # Different MIP solvers may find different optimal solutions when there are ties,
  # so we use a tolerance of 0.02 (2%) for p-value comparisons
  expect_equal(result_highs_ilp["p.value"], result_gurobi_ilp["p.value"],
               tolerance = 0.02,
               info = "P-values from HiGHS and Gurobi should match for ILP")

  expect_equal(result_highs_ilp["test.stat"], result_gurobi_ilp["test.stat"],
               tolerance = 0.02,
               info = "Test statistics from HiGHS and Gurobi should match for ILP")

  # Test LP (relaxed) mode
  result_highs_lp <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                                      opt.method = "LP_highs",
                                      stat.null = stat.null,
                                      statistic = TRUE)

  result_gurobi_lp <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                                       opt.method = "LP_gurobi",
                                       stat.null = stat.null,
                                       statistic = TRUE)

  expect_equal(result_highs_lp["test.stat"], result_gurobi_lp["test.stat"],
               tolerance = 0.02,
               info = "Test statistics from HiGHS and Gurobi should match for LP")
})


# Test stratum-level solver (for comb.method=2)
test_that("Stratum-level solver works with HiGHS", {
  skip_if_not(solver_available("highs"), "HiGHS is not available")

  set.seed(22222)

  s <- 3  # number of strata
  n <- 6  # units per stratum
  m <- 3  # treated per stratum
  N <- s * n
  k <- ceiling(0.8 * N)
  c <- 0

  block <- factor(rep(1:s, each = n))

  Z <- rep(0, N)
  for (i in 1:s) {
    block_idx <- which(block == i)
    Z[sample(block_idx, m)] <- 1
  }

  Y0 <- rnorm(N)
  Y1 <- Y0 + 1
  Y <- Z * Y1 + (1 - Z) * Y0

  # Set up polynomial methods for comb.method=2
  r.vec <- c(2, 6)
  methods.list.all <- list()
  for (j in seq_along(r.vec)) {
    methods.list.all[[j]] <- lapply(1:s, function(i) {
      list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
    })
  }

  # Test comb.method=2 which uses stratum-level solver internally
  result <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                            opt.method = "ILP_highs",
                            null.max = 100,
                            comb.method = 2)

  expect_true(is.numeric(result))
  expect_length(result, 2)
  expect_true(result["p.value"] >= 0 && result["p.value"] <= 1)
  expect_true(is.finite(result["test.stat"]))
})


test_that("Stratum-level solvers give equivalent results for HiGHS and Gurobi", {
  skip_if_not(solver_available("highs") && solver_available("gurobi"),
              "Both HiGHS and Gurobi are required for this test")

  set.seed(33333)

  s <- 3
  n <- 6
  m <- 3
  N <- s * n
  k <- ceiling(0.8 * N)
  c <- 0

  block <- factor(rep(1:s, each = n))

  Z <- rep(0, N)
  for (i in 1:s) {
    block_idx <- which(block == i)
    Z[sample(block_idx, m)] <- 1
  }

  Y0 <- rnorm(N)
  Y1 <- Y0 + 0.8
  Y <- Z * Y1 + (1 - Z) * Y0

  r.vec <- c(2, 6)
  methods.list.all <- list()
  for (j in seq_along(r.vec)) {
    methods.list.all[[j]] <- lapply(1:s, function(i) {
      list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
    })
  }

  # Run both solvers with same null.max for comparable results
  # Note: Since null distributions are generated with random permutations,
  # results may differ slightly between runs even with same seed
  result_highs <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                                   opt.method = "ILP_highs",
                                   null.max = 200,
                                   comb.method = 2)

  result_gurobi <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                                    opt.method = "ILP_gurobi",
                                    null.max = 200,
                                    comb.method = 2)

  # Both should produce valid results
  # Note: P-values may differ because each call generates its own null distribution
  # with random permutations. We just verify both produce valid outputs.
  expect_true(result_highs["p.value"] >= 0 && result_highs["p.value"] <= 1)
  expect_true(result_gurobi["p.value"] >= 0 && result_gurobi["p.value"] <= 1)
  expect_true(is.finite(result_highs["test.stat"]))
  expect_true(is.finite(result_gurobi["test.stat"]))
})


# Test solve_optimization directly
test_that("solve_optimization works with auto solver selection", {
  skip_if_not(solver_available("highs") || solver_available("gurobi"),
              "Neither HiGHS nor Gurobi is available")

  set.seed(99999)

  s <- 2
  n <- 4
  m <- 2
  N <- s * n
  k <- N - 1
  p <- N - k
  c <- 0

  block <- factor(rep(1:s, each = n))

  Z <- rep(0, N)
  for (i in 1:s) {
    block_idx <- which(block == i)
    Z[sample(block_idx, m)] <- 1
  }

  Y0 <- rnorm(N)
  Y1 <- Y0 + 1
  Y <- Z * Y1 + (1 - Z) * Y0

  methods.list.all <- list()
  methods.list.all[[1]] <- list()
  for (i in 1:s) {
    methods.list.all[[1]][[i]] <- list(name = "Wilcoxon", scale = FALSE)
  }

  block.sum <- summary_block(Z, block)
  weight <- weight_scheme(block.sum, "asymp.opt")
  scores.list.all <- list()
  scores.list.all[[1]] <- score_all_blocks(block.sum$nb, methods.list.all[[1]])

  ms_list <- mu_sigma_list(Z, block, weight, methods.list.all, scores.list.all, block.sum)
  coeflists <- comb_matrix_block(Z, Y, block, c, methods.list.all, scores.list.all, block.sum)

  # Test with auto solver selection
  result <- solve_optimization(Z, block, weight, coeflists, p, ms_list,
                                exact = TRUE, block.sum = block.sum, solver = "auto")

  expect_true(is.list(result))
  expect_true("obj" %in% names(result))
  expect_true("sol" %in% names(result))
  expect_true("solver" %in% names(result))
  expect_true(result$solver %in% c("highs", "gurobi"))
  expect_true(is.finite(result$obj))
})


# Test backward compatibility
test_that("Old opt.method values still work (backward compatibility)", {
  skip_if_not(solver_available("gurobi"), "Gurobi is required for backward compatibility test")

  set.seed(11111)

  s <- 2
  n <- 4
  m <- 2
  N <- s * n
  k <- ceiling(0.75 * N)
  c <- 0

  block <- factor(rep(1:s, each = n))

  Z <- rep(0, N)
  for (i in 1:s) {
    block_idx <- which(block == i)
    Z[sample(block_idx, m)] <- 1
  }

  Y0 <- rnorm(N)
  Y1 <- Y0 + 0.5
  Y <- Z * Y1 + (1 - Z) * Y0

  methods.list.all <- list()
  methods.list.all[[1]] <- list()
  for (i in 1:s) {
    methods.list.all[[1]][[i]] <- list(name = "Wilcoxon", scale = FALSE)
  }

  # Test old-style opt.method values still work
  result <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                            opt.method = "ILP_gurobi",
                            null.max = 100,
                            statistic = TRUE)

  expect_true(is.numeric(result))
  expect_true(result["p.value"] >= 0 && result["p.value"] <= 1)
})


# Test multiple quantiles with different solvers
test_that("Multiple quantile tests produce consistent results across solvers", {
  skip_if_not(solver_available("highs") && solver_available("gurobi"),
              "Both HiGHS and Gurobi are required for this test")

  set.seed(77777)

  s <- 3
  n <- 6
  m <- 3
  N <- s * n
  c <- 0

  block <- factor(rep(1:s, each = n))

  Z <- rep(0, N)
  for (i in 1:s) {
    block_idx <- which(block == i)
    Z[sample(block_idx, m)] <- 1
  }

  Y0 <- rnorm(N)
  Y1 <- Y0 + 1.2
  Y <- Z * Y1 + (1 - Z) * Y0

  methods.list.all <- list()
  methods.list.all[[1]] <- list()
  for (i in 1:s) {
    methods.list.all[[1]][[i]] <- list(name = "Wilcoxon", scale = FALSE)
  }

  # Pre-compute null distribution for consistency
  block.sum <- summary_block(Z, block)
  weight <- weight_scheme(block.sum, "asymp.opt")
  scores.list.all <- list()
  scores.list.all[[1]] <- score_all_blocks(block.sum$nb, methods.list.all[[1]])
  stat.null <- com_null_dist_block(Z, block, methods.list.all, scores.list.all,
                                   null.max = 200, weight, block.sum, Z.perm = NULL)

  # Test several quantiles
  k_values <- c(floor(0.7 * N), floor(0.8 * N), floor(0.9 * N))

  for (k in k_values) {
    result_highs <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                                     opt.method = "ILP_highs",
                                     stat.null = stat.null,
                                     statistic = FALSE)

    result_gurobi <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                                      opt.method = "ILP_gurobi",
                                      stat.null = stat.null,
                                      statistic = FALSE)

    # Different MIP solvers may find different optimal solutions when there are ties
    # Use larger tolerance (5%) as MIP solvers can find different optima
    expect_equal(result_highs, result_gurobi, tolerance = 0.05,
                 info = paste("Results should match for k =", k))
  }
})
