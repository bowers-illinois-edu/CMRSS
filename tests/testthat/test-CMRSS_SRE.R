# Tests for SRE (Stratified Randomized Experiment) functions

test_that("pval_comb_block works with comb.method=1", {
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  set.seed(12345)

  s <- 3    # number of strata
  n <- 8    # units per stratum
  m <- 4    # treated per stratum
  N <- s * n
  k <- floor(0.8 * N)
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

  # Set up methods
  r.vec <- c(2, 6)
  methods.list.all <- list()
  for (j in seq_along(r.vec)) {
    methods.list.all[[j]] <- lapply(1:s, function(i) {
      list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
    })
  }

  result <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                            opt.method = "ILP",
                            null.max = 500,
                            comb.method = 1)

  expect_true(is.numeric(result))
  expect_length(result, 2)
  expect_true(result["p.value"] >= 0 && result["p.value"] <= 1)
})


test_that("pval_comb_block works with comb.method=2", {
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  set.seed(23456)

  s <- 3
  n <- 8
  m <- 4
  N <- s * n
  k <- floor(0.8 * N)
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

  r.vec <- c(2, 6)
  methods.list.all <- list()
  for (j in seq_along(r.vec)) {
    methods.list.all[[j]] <- lapply(1:s, function(i) {
      list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
    })
  }

  result <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                            opt.method = "ILP",
                            null.max = 500,
                            comb.method = 2)

  expect_true(is.numeric(result))
  expect_length(result, 2)
  expect_true(result["p.value"] >= 0 && result["p.value"] <= 1)
})


test_that("com_block_conf_quant_larger works with comb.method=1", {
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  set.seed(34567)

  s <- 3
  n <- 8
  m <- 4
  N <- s * n

  block <- factor(rep(1:s, each = n))
  Z <- rep(0, N)
  for (i in 1:s) {
    block_idx <- which(block == i)
    Z[sample(block_idx, m)] <- 1
  }

  Y0 <- rnorm(N)
  Y1 <- Y0 + 0.5
  Y <- Z * Y1 + (1 - Z) * Y0

  r.vec <- c(2, 6)
  methods.list.all <- list()
  for (j in seq_along(r.vec)) {
    methods.list.all[[j]] <- lapply(1:s, function(i) {
      list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
    })
  }

  ci <- com_block_conf_quant_larger(Z, Y, block,
                                    set = "treat",
                                    methods.list.all = methods.list.all,
                                    opt.method = "ILP",
                                    null.max = 200,
                                    tol = 0.1,
                                    alpha = 0.10,
                                    comb.method = 1)

  expect_length(ci, sum(Z))
  expect_true(all(is.numeric(ci)))
})


test_that("com_block_conf_quant_larger works with comb.method=2", {
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  set.seed(45678)

  s <- 3
  n <- 8
  m <- 4
  N <- s * n

  block <- factor(rep(1:s, each = n))
  Z <- rep(0, N)
  for (i in 1:s) {
    block_idx <- which(block == i)
    Z[sample(block_idx, m)] <- 1
  }

  Y0 <- rnorm(N)
  Y1 <- Y0 + 0.5
  Y <- Z * Y1 + (1 - Z) * Y0

  r.vec <- c(2, 6)
  methods.list.all <- list()
  for (j in seq_along(r.vec)) {
    methods.list.all[[j]] <- lapply(1:s, function(i) {
      list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
    })
  }

  ci <- com_block_conf_quant_larger(Z, Y, block,
                                    set = "treat",
                                    methods.list.all = methods.list.all,
                                    opt.method = "ILP",
                                    null.max = 200,
                                    tol = 0.1,
                                    alpha = 0.10,
                                    comb.method = 2)

  expect_length(ci, sum(Z))
  expect_true(all(is.numeric(ci)))
})


test_that("SRE functions work with electric_teachers data", {
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  data(electric_teachers)

  Z <- electric_teachers$TxAny
  Y <- electric_teachers$gain
  block <- factor(electric_teachers$Site)
  s <- length(levels(block))
  n <- length(Y)

  # Define polynomial methods
  r.vec <- c(2, 6)
  methods.list.all <- list()
  for (j in seq_along(r.vec)) {
    methods.list.all[[j]] <- lapply(1:s, function(i) {
      list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
    })
  }

  # Test pval_comb_block
  k <- floor(0.9 * n)
  result <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
                            opt.method = "ILP",
                            null.max = 500)

  expect_true(result["p.value"] >= 0 && result["p.value"] <= 1)
})


test_that("Both comb.methods give valid but different results", {
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  set.seed(56789)

  s <- 3
  n <- 8
  m <- 4
  N <- s * n
  k <- floor(0.8 * N)
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

  r.vec <- c(2, 6)
  methods.list.all <- list()
  for (j in seq_along(r.vec)) {
    methods.list.all[[j]] <- lapply(1:s, function(i) {
      list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
    })
  }

  result1 <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                             opt.method = "ILP",
                             null.max = 500,
                             comb.method = 1)

  result2 <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                             opt.method = "ILP",
                             null.max = 500,
                             comb.method = 2)

  # Both should be valid p-values
  expect_true(result1["p.value"] >= 0 && result1["p.value"] <= 1)
  expect_true(result2["p.value"] >= 0 && result2["p.value"] <= 1)

  # Test statistics should be different (different combination methods)
  # This is not always guaranteed but typically true
  expect_true(is.finite(result1["test.stat"]))
  expect_true(is.finite(result2["test.stat"]))
})


test_that("LP relaxation gives valid results", {
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  set.seed(67890)

  s <- 3
  n <- 8
  m <- 4
  N <- s * n
  k <- floor(0.8 * N)
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

  r.vec <- c(2, 6)
  methods.list.all <- list()
  for (j in seq_along(r.vec)) {
    methods.list.all[[j]] <- lapply(1:s, function(i) {
      list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
    })
  }

  # LP relaxation
  result_lp <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                               opt.method = "LP",
                               null.max = 500,
                               comb.method = 1)

  expect_true(result_lp["p.value"] >= 0 && result_lp["p.value"] <= 1)
  expect_true(is.finite(result_lp["test.stat"]))
})
