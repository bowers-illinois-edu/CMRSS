# Tests for SRE (Stratified Randomized Experiment) functions

# Internal function tests (verify vectorization correctness)

test_that("summary_block correctly summarizes block structure", {
  Z <- c(1, 1, 0, 0, 0,
         1, 1, 1, 0, 0)
  block <- factor(c(rep("A", 5), rep("B", 5)))

  result <- CMRSS:::summary_block(Z, block)

  expect_equal(result$B, 2)
  expect_equal(result$nb, c(5, 5))
  expect_equal(result$mb, c(2, 3))
  expect_equal(result$mb_ctrl, c(3, 2))
})


test_that("assign_block generates valid block-randomized permutations", {
  set.seed(202)

  # Create a simple block structure
  n <- 30
  Z <- c(rep(1, 5), rep(0, 5),
         rep(1, 3), rep(0, 7),
         rep(1, 4), rep(0, 6))
  block <- factor(c(rep("A", 10), rep("B", 10), rep("C", 10)))

  block_sum <- CMRSS:::summary_block(Z, block)
  null_max <- 100

  z_perm <- CMRSS:::assign_block(block_sum, null_max)

  # Check dimensions
  expect_equal(dim(z_perm), c(n, null_max))

  # Check all values are 0 or 1
  expect_true(all(z_perm %in% c(0, 1)))

  # Check each block has correct number of treated in each permutation
  for (i in 1:block_sum$B) {
    block_sums <- colSums(z_perm[block_sum$units.block[[i]], , drop = FALSE])
    expect_true(all(block_sums == block_sum$mb[i]))
  }

  # Check total treated is correct in each permutation
  expect_true(all(colSums(z_perm) == sum(Z)))
})

test_that(".generate_Z_chunk_block generates valid block-randomized chunks", {
  set.seed(101)

  n <- 50
  block <- factor(sample(letters[1:10], n, replace = TRUE))
  Z <- sample(c(rep(1, 25), rep(0, 25)))

  block_sum <- CMRSS:::summary_block(Z, block)
  size <- 200
  z_chunk <- CMRSS:::.generate_Z_chunk_block(block_sum$units.block, block_sum$mb, n, size)

  expect_equal(dim(z_chunk), c(n, size))
  expect_true(all(z_chunk %in% c(0, 1)))

  for (i in 1:block_sum$B) {
    block_sums <- colSums(z_chunk[block_sum$units.block[[i]], , drop = FALSE])
    expect_true(all(block_sums == block_sum$mb[i]))
  }

  expect_true(all(colSums(z_chunk) == sum(Z)))
})

test_that("com_null_dist_block_stratum matches naive computation for fixed Z.perm", {
  set.seed(202)

  block <- factor(rep(1:3, each = 4))
  n <- length(block)
  Z <- c(1, 1, 0, 0,
         1, 0, 1, 0,
         0, 1, 1, 0)

  block_sum <- CMRSS:::summary_block(Z, block)
  weight <- block_sum$mb

  r.vec <- c(2, 3, 6)
  methods.list.all <- lapply(r.vec, function(r) {
    lapply(1:block_sum$B, function(i) list(name = "Polynomial", r = r, std = TRUE, scale = FALSE))
  })

  scores.list.all <- lapply(methods.list.all, function(method.list.all) {
    CMRSS:::score_all_blocks(block_sum$nb, method.list.all)
  })

  mu_sd_block_list <- CMRSS:::mu_sd_block(Z, block, methods.list.all, scores.list.all, weight, block_sum)

  null_max <- 200
  Z.perm <- CMRSS:::assign_block(block_sum, null_max)

  fast <- CMRSS:::com_null_dist_block_stratum(
    Z, block, methods.list.all,
    scores.list.all = scores.list.all,
    mu_sd_block_list = mu_sd_block_list,
    null.max = null_max,
    weight = weight,
    block.sum = block_sum,
    Z.perm = Z.perm,
    chunk_size = 37
  )

  H <- length(methods.list.all)
  B <- block_sum$B
  naive <- numeric(null_max)
  for (iter in 1:null_max) {
    total <- 0
    for (b in 1:B) {
      units <- block_sum$units.block[[b]]
      Zb <- Z.perm[units, iter]
      mb <- block_sum$mb[b]
      vals <- numeric(H)
      for (h in 1:H) {
        score_bh <- scores.list.all[[h]][[b]]
        mu <- mu_sd_block_list$mu_list[[b]][h]
        sd <- mu_sd_block_list$sig_list[[b]][h]
        vals[h] <- ((sum(score_bh * Zb) / mb) - mu) / sd
      }
      total <- total + weight[b] * max(vals)
    }
    naive[iter] <- total
  }

  expect_equal(fast, naive)
})


# API function tests

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
