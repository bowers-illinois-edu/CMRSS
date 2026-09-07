# Tests for PLAN item 1B: the column range of comb_matrix_block_stratum.
#
# `comb_matrix_block` enumerates exempt counts 0..m_b, one per treated unit
# in the block plus zero.  `comb_matrix_block_stratum` enumerates 0..n_b,
# one per unit in the block.  The index means the number of treated units in
# the block whose effect exceeds the threshold, so it cannot exceed m_b, and
# the wider range looked like it might be a bug.
#
# It is not, and these tests say why.  Inside `min_stat` the exempt count is
# `min(m, n - k)`, so once the requested count reaches m_b every further
# column exempts the same units and the statistic stops moving.  The extra
# columns are duplicates of the last real one.
#
# What matters is that they are duplicates carrying a larger index.  The
# stratum solver picks exactly one column per block (`sum x = 1`) and
# charges the chosen column's index against a shared budget
# (`sum index * x <= p`, R/solvers.R:38-46 of HiGHS_sol_stratum_com).  A
# duplicate that costs more budget for the same objective value can never
# be part of an optimal solution when the budget binds, and gives the same
# objective when it does not.  So dropping the columns above m_b leaves
# every answer alone and removes work that grows with n_b - m_b.
#
# These tests pin all three steps: the saturation, the domination, and the
# p-values that come out the far end.

# One stratified design, fixed.
stratum_design <- function() {
  Z <- c(1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0)
  Y <- c(0.916, 0.152, -0.361, -0.438, -0.561, 0.293,
         2.273, 1.581, -0.369, 1.010, 0.107, 0.721)
  block <- factor(rep(1:2, each = 6))
  methods.list.all <- lapply(c(2, 6), function(r)
    lapply(1:2, function(i)
      list(name = "Polynomial", r = r, std = TRUE, scale = FALSE)))
  list(Z = Z, Y = Y, block = block, methods.list.all = methods.list.all)
}


test_that("min_stat saturates once the exempt count reaches the treated count", {
  # The fact the whole item turns on.  R/CMRSS_CRE.R:275 exempts
  # min(m, n - k) units, so asking for more than m changes nothing.
  set.seed(31)
  nb <- 10
  mb <- 4
  Zb <- sample(rep(0:1, times = c(nb - mb, mb)))
  Yb <- round(rnorm(nb), 3)
  score <- rank_score(nb, list(name = "Polynomial", r = 3,
                               std = TRUE, scale = FALSE))

  vals <- vapply(0:nb, function(ii)
    min_stat(Zb, Yb, nb - ii, 0.2, score = score), numeric(1))

  # Columns 0..mb move; columns mb..nb are all the same number.
  expect_equal(length(unique(vals[(mb + 1):(nb + 1)])), 1L)
  expect_gt(length(unique(vals[1:(mb + 1)])), 1L)
})


test_that("the stratum matrix stops at the treated count", {
  d <- stratum_design()
  block.sum <- summary_block(d$Z, d$block)
  weight <- weight_scheme(block.sum, "asymp.opt")

  cl <- max_comb_matrix_block_stratum(d$Z, d$Y, d$block, 0.2,
                                      d$methods.list.all, NULL, block.sum,
                                      weight = weight)

  for (i in seq_along(cl)) {
    mb <- block.sum$mb[i]
    nb <- block.sum$nb[i]
    # One column per exempt count from 0 to m_b, where the old range ran to
    # n_b.  On this design that is 4 columns rather than 7 per block.
    expect_equal(ncol(cl[[i]]), mb + 1L)
    expect_equal(cl[[i]][1, ], 0:mb)
    expect_lt(ncol(cl[[i]]), nb + 1L)
  }
})


test_that("restoring the dropped columns would not change the solver's answer", {
  skip_if_not(solver_available("highs"), "HiGHS not available")

  # Why dropping them is safe, checked against the code as it now stands.
  # Rebuild the wider matrix the old range produced --- the extra columns
  # repeat the last statistic while carrying indices m_b + 1 up to n_b ---
  # and confirm the solver returns the same objective at every budget.  A
  # duplicate that costs more budget for the same value can never improve a
  # minimum.
  d <- stratum_design()
  block.sum <- summary_block(d$Z, d$block)
  weight <- weight_scheme(block.sum, "asymp.opt")

  widen <- function(mat, nb) {
    mb <- ncol(mat) - 1L
    extra_idx <- seq.int(mb + 1L, nb)
    extra <- rbind(extra_idx, rep(mat[2, mb + 1L], length(extra_idx)))
    cbind(mat, extra)
  }

  for (cc in c(-0.5, 0, 0.2, 1)) {
    narrow <- max_comb_matrix_block_stratum(d$Z, d$Y, d$block, cc,
                                            d$methods.list.all, NULL,
                                            block.sum, weight = weight)
    wide <- lapply(seq_along(narrow),
                   function(i) widen(narrow[[i]], block.sum$nb[i]))

    for (p in 0:sum(block.sum$mb)) {
      expect_equal(
        solve_stratum_optimization(narrow, p, exact = TRUE,
                                   solver = "highs")$obj,
        solve_stratum_optimization(wide, p, exact = TRUE,
                                   solver = "highs")$obj,
        info = sprintf("c = %g, p = %d", cc, p)
      )
    }
  }
})


test_that("comb.method = 2 p-values are unchanged", {
  skip_if_not(solver_available("highs"), "HiGHS not available")

  # The end-to-end pin.  Captured 2026-09-07 at commit 91bd44b, before the
  # column range was narrowed.
  d <- stratum_design()

  set.seed(808)
  res <- pval_comb_block(d$Z, d$Y, k = 5, c = 0, block = d$block,
                         methods.list.all = d$methods.list.all,
                         comb.method = 2, opt.method = "ILP_highs",
                         null.max = 200)

  expect_equal(unname(res[2]), 0.8895742683, tolerance = 1e-8)
})
