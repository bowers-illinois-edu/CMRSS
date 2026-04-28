# Tests for the `p` convention in `pval_comb_block` (R/CMRSS_SRE.R).
#
# Why this file exists:
#
# `pval_comb_block` tests the treated-only hypothesis H_{k,c}^treat. The LP
# coverage constraint at R/CMRSS_SRE.R:1034 sets `p <- m - k`, so the valid
# input range is k in 1..m_total (where m_total = sum(Z)). Outside that range
# `p` becomes negative, the LP is infeasible, the test statistic returns Inf,
# and the function silently reports `p.value = 0` -- a false rejection. These
# tests pin (i) that the function fails visibly outside the valid range, and
# (ii) that valid inputs produce well-formed output.
#
# Note: the inversion sibling `com_block_conf_quant_larger_trt`
# (R/CMRSS_SRE.R:1171, 1234) uses `p <- n - k` (all-units convention). It
# tests a different hypothesis than `pval_comb_block` despite the `_trt`
# suffix. No cross-function comparison test is included here; whether the two
# functions should be reconciled is a separate decision.

test_that("pval_comb_block does not silently return p.value = 0 when k > m_total", {
  skip_if_not(solver_available("highs"), "HiGHS not available")

  set.seed(12345)
  s <- 3; n_per <- 8; m_per <- 4
  N <- s * n_per
  m_total <- s * m_per
  k <- m_total + 1  # out of valid treated-only range

  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- Z * (rnorm(N) + 0.5) + (1 - Z) * rnorm(N)

  methods.list.all <- lapply(c(2, 6), function(rj) {
    lapply(1:s, function(i) list(name = "Polynomial", r = rj, std = TRUE, scale = FALSE))
  })

  # If the function errors on out-of-domain k, treat that as the desired
  # behavior (visible failure). If it returns, inspect the return value.
  # Suppress warnings so a stray HiGHS warning does not mask the assertion;
  # the bug is that the function returns p.value = 0 + test.stat = Inf even
  # when the LP is infeasible.
  res <- suppressWarnings(tryCatch(
    pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
                    opt.method = "ILP_highs", null.max = 200, comb.method = 1),
    error = function(e) c(p.value = NA_real_, test.stat = NA_real_)
  ))

  bad <- isTRUE(unname(res["p.value"]) == 0 &&
                is.infinite(unname(res["test.stat"])))
  expect_false(bad,
               info = "silent false rejection on k > m_total: p.value=0, test.stat=Inf")
})


test_that("pval_comb_block returns finite stat and p.value in [0, 1] for valid k", {
  skip_if_not(solver_available("highs"), "HiGHS not available")

  set.seed(12345)
  s <- 3; n_per <- 8; m_per <- 4
  N <- s * n_per
  m_total <- s * m_per

  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- Z * (rnorm(N) + 0.5) + (1 - Z) * rnorm(N)

  methods.list.all <- lapply(c(2, 6), function(rj) {
    lapply(1:s, function(i) list(name = "Polynomial", r = rj, std = TRUE, scale = FALSE))
  })

  for (k in c(1, floor(m_total / 2), m_total - 1)) {
    res <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
                           opt.method = "ILP_highs", null.max = 200, comb.method = 1)
    expect_true(is.finite(unname(res["test.stat"])),
                info = sprintf("test.stat not finite at k = %d", k))
    expect_true(unname(res["p.value"]) >= 0 && unname(res["p.value"]) <= 1,
                info = sprintf("p.value out of [0, 1] at k = %d", k))
  }
})


test_that("pval_comb_block is feasible at boundary k = m_total (p = 0)", {
  skip_if_not(solver_available("highs"), "HiGHS not available")

  # k = m_total is the most extreme treated quantile (largest treated
  # effect). p = m - k = 0, which is the tightest valid LP constraint.
  # The LP must remain feasible at this boundary.

  set.seed(12345)
  s <- 3; n_per <- 8; m_per <- 4
  N <- s * n_per
  m_total <- s * m_per
  k <- m_total

  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- Z * (rnorm(N) + 0.5) + (1 - Z) * rnorm(N)

  methods.list.all <- lapply(c(2, 6), function(rj) {
    lapply(1:s, function(i) list(name = "Polynomial", r = rj, std = TRUE, scale = FALSE))
  })

  res <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
                         opt.method = "ILP_highs", null.max = 200, comb.method = 1)
  expect_true(is.finite(unname(res["test.stat"])))
  expect_true(unname(res["p.value"]) >= 0 && unname(res["p.value"]) <= 1)
})
