# Tests for the `p` convention in `pval_comb_block` (R/CMRSS_SRE.R).
#
# Why this file exists:
#
# Commit 762c4d08 (David Kim, 2025-12-18) changed line 1034 of CMRSS_SRE.R from
#   p <- N - k
# to
#   p <- m - k
# and left the original commented out at line 1033. The function's docstring
# still says "k between 1 and n", its example uses k = floor(0.9 * n), and the
# inversion sibling com_block_conf_quant_larger_trt continues to use
# p <- n - k. Under the new formula, any call with k > sum(Z) makes the LP
# infeasible; HiGHS warns, the test statistic returns Inf, and the function
# silently reports p.value = 0 -- a false rejection.
#
# These tests pin the desired behavior on the assumption that p = n - k is the
# correct convention (matching the paper's reference code at
# combined_stephenson_tests/code/codes_20251026.R line 1426, the docstring,
# the docstring example, the SRE test inputs, and the inversion path).
# They are skipped until that convention is confirmed by David Kim. When
# confirmed, remove the skip() calls and apply the fix at R/CMRSS_SRE.R:1034.

test_that("pval_comb_block does not silently return p.value = 0 on infeasible LP", {
  skip("Pending: commit 762c4d08 review with David Kim. See file header.")

  set.seed(12345)
  s <- 3; n_per <- 8; m_per <- 4
  N <- s * n_per
  k <- floor(0.8 * N)  # 19, larger than total treated (12)

  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- Z * (rnorm(N) + 0.5) + (1 - Z) * rnorm(N)

  methods.list.all <- lapply(c(2, 6), function(rj) {
    lapply(1:s, function(i) list(name = "Polynomial", r = rj, std = TRUE, scale = FALSE))
  })

  res <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
                         opt.method = "ILP_highs", null.max = 200, comb.method = 1)

  # Either the LP must be feasible (p = n - k = 5, valid allocation exists),
  # or the function must signal an error / NA / non-zero p-value. What it
  # must NOT do is return p.value = 0 + test.stat = Inf.
  expect_false(isTRUE(res["p.value"] == 0 && is.infinite(res["test.stat"])))
})


test_that("pval_comb_block test statistic matches the inversion path at the same (k, c)", {
  skip("Pending: commit 762c4d08 review with David Kim. See file header.")

  # When pval_comb_block(k, c) and com_block_conf_quant_larger_trt are
  # evaluated on the same data at the same (k, c), the underlying LP that
  # produces the test statistic should be identical. This pins the two
  # functions to share a single `p` convention. Under the current code they
  # do not, and pval_comb_block's LP is infeasible whenever k > sum(Z).

  set.seed(12345)
  s <- 3; n_per <- 8; m_per <- 4
  N <- s * n_per
  k <- floor(0.8 * N)
  c_val <- 0

  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- Z * (rnorm(N) + 0.5) + (1 - Z) * rnorm(N)

  methods.list.all <- lapply(c(2, 6), function(rj) {
    lapply(1:s, function(i) list(name = "Polynomial", r = rj, std = TRUE, scale = FALSE))
  })

  # Reconstruct the LP that pval_comb_block calls, but with p = n - k
  # (the convention used by com_block_conf_quant_larger_trt at line 1171).
  opt_params <- CMRSS:::parse_opt_method("ILP_highs")
  block.sum <- CMRSS:::summary_block(Z, block)
  block_f <- block.sum$block
  weight <- CMRSS:::weight_scheme(block.sum, "asymp.opt")
  scores.list.all <- lapply(methods.list.all, function(ml) {
    CMRSS:::score_all_blocks(block.sum$nb, ml)
  })
  ms_list <- CMRSS:::mu_sigma_list(Z, block_f, weight, methods.list.all,
                                   scores.list.all, block.sum)
  coeflists <- CMRSS:::comb_matrix_block(Z, Y, block_f, c_val, methods.list.all,
                                         scores.list.all, block.sum)
  stat_min_inversion <- CMRSS:::solve_optimization(
    Z, block_f, weight, coeflists, p = N - k, ms_list,
    opt_params$exact, block.sum, opt_params$solver)$obj

  res <- pval_comb_block(Z, Y, k, c = c_val, block, methods.list.all,
                         opt.method = "ILP_highs", null.max = 200, comb.method = 1)

  expect_equal(unname(res["test.stat"]), stat_min_inversion, tolerance = 1e-6)
})


test_that("pval_comb_block at k = n - m + 1 gives test stat consistent with treated-quantile reading", {
  skip("Pending: commit 762c4d08 review with David Kim. See file header.")

  # Per ZL24 (cited at main.tex line 518): a p-value for H_{k,c} (all-units
  # quantile) is also valid for H_{k - n_control, c}^treat. So at k = n - m + 1
  # (treated-quantile-1, the largest treated effect), p = n - k = m - 1, which
  # is feasible whether you read k as all-units or as 'just above the
  # treated-only range'. This is a sanity check that should pass under either
  # convention -- if it fails, the discrepancy is more serious than a single
  # constant offset.

  set.seed(12345)
  s <- 3; n_per <- 8; m_per <- 4
  N <- s * n_per
  m_total <- s * m_per
  k <- N - m_total + 1  # = 13, just above the treated range

  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- Z * (rnorm(N) + 0.5) + (1 - Z) * rnorm(N)

  methods.list.all <- lapply(c(2, 6), function(rj) {
    lapply(1:s, function(i) list(name = "Polynomial", r = rj, std = TRUE, scale = FALSE))
  })

  res <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
                         opt.method = "ILP_highs", null.max = 200, comb.method = 1)
  expect_true(is.finite(res["test.stat"]))
  expect_true(res["p.value"] >= 0 && res["p.value"] <= 1)
})
