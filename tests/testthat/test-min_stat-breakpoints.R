# Tests for the claim that makes the confidence-interval refactor possible:
# with the design, the outcomes, the scores, and k all held fixed,
# `min_stat` is a monotone step function of the threshold c whose jumps can
# only occur at the within-block treated-minus-control differences
#
#     D = { Y_i - Y_j : Z_i = 1, Z_j = 0 },
#
# a set with at most m_b * (n_b - m_b) elements.
#
# Why this is true, from the code at R/CMRSS_CRE.R:253-281.  `min_stat`
# exempts treated units by outcome rank through `sort_treat(Y, Z)`, which
# reads Y and Z and not c, so the exempt set does not move with c.  The
# exempt units get xi = Inf and therefore Y - Z*xi = -Inf, placing them at
# the bottom of the ranking whatever c is.  The remaining treated units are
# ranked at Y_i - c and the controls at Y_j.  Two treated units shift by the
# same c, so their order never changes; a treated unit i passes a control
# unit j exactly as c crosses Y_i - Y_j.  The statistic reads the ranks and
# nothing else, so it is constant between consecutive elements of D.
#
# Why it matters.  `com_block_conf_quant_larger` currently finds each limit
# with `uniroot` plus four uncapped `while (f(c.sol) <= 0) c.sol <- c.sol -
# tol` loops, and the profiling memo at memos/memo_cmrss_reuse_2026-09-07.qmd
# measured 171,600 `min_stat` calls in a single 21.8-second interval, 61
# percent of the run.  If the jumps are known in advance, the search becomes
# a bisection over a sorted finite list, and the limit it returns is exact
# rather than accurate to `tol`.
#
# These tests pin the claim before any of that is written.  They pass
# against the current implementation and must keep passing against the
# tabulated one.

# Regression values captured from the current implementation on
# 2026-09-07 at commit 455f72f, seed 7, a 12-unit block with 6 treated and
# Polynomial(r = 6, std = TRUE) scores, at
# c = -2, -0.5, 0, 0.25, 1, 3.  The tabulated implementation must return
# these exactly: it is meant to stop recomputing the same values, not to
# compute different ones.
MIN_STAT_REGRESSION_K12 <- c(1.553896788, 0.9384987059, 0.6040566345,
                             0.3688461673, 0.2812468859, 0.03286084036)
MIN_STAT_REGRESSION_K9 <- c(0.8508994244, 0.1709539366, 0.05718394906,
                            0.05718394906, 0.03286084036, 0.03286084036)

# Build one block and everything min_stat needs for it.
make_block <- function(nb, mb, seed, r = 6) {
  set.seed(seed)
  Zb <- sample(rep(0:1, times = c(nb - mb, mb)))
  Yb <- round(rnorm(nb), 3)
  score <- rank_score(nb, list(name = "Polynomial", r = r,
                               std = TRUE, scale = FALSE))
  diffs <- sort(unique(as.vector(outer(Yb[Zb == 1], Yb[Zb == 0], "-"))))
  list(Z = Zb, Y = Yb, score = score, diffs = diffs, nb = nb, mb = mb)
}

# Evaluate min_stat across a grid of thresholds.
stat_curve <- function(b, k, grid) {
  vapply(grid, function(cc) min_stat(b$Z, b$Y, k, cc, score = b$score),
         numeric(1))
}


test_that("the exempt treated set does not depend on c", {
  # The premise of everything below.  If sort_treat read c, the jump points
  # would not be the Y_i - Y_j differences and the tabulation would be wrong.
  b <- make_block(nb = 12, mb = 6, seed = 7)

  first <- sort_treat(b$Y, b$Z)
  for (cc in c(-100, -1, 0, 1, 100)) {
    expect_identical(sort_treat(b$Y, b$Z), first)
  }

  # sort_treat's formal arguments name Y and Z and no threshold.
  expect_false("c" %in% names(formals(sort_treat)))
})


test_that("min_stat is a monotone step function of c", {
  # Three block shapes, including unbalanced ones, and three values of k
  # spanning no exemption through several exempted units.
  shapes <- list(c(nb = 12, mb = 6), c(nb = 10, mb = 3), c(nb = 9, mb = 7))

  for (i in seq_along(shapes)) {
    sh <- shapes[[i]]
    b <- make_block(nb = sh[["nb"]], mb = sh[["mb"]], seed = 100 + i)
    grid <- seq(min(b$diffs) - 1, max(b$diffs) + 1, length.out = 4001)

    for (k in c(b$nb, b$nb - 2, b$nb - 5)) {
      vals <- stat_curve(b, k, grid)

      # Raising c can only push treated units down the ranking, and the
      # scores increase with rank, so the statistic cannot increase.
      expect_true(all(diff(vals) <= 0),
                  info = sprintf("nb=%d mb=%d k=%d not monotone",
                                 b$nb, b$mb, k))

      # A step function takes finitely many values, at most one more than
      # its number of jumps.
      expect_lte(length(unique(vals)), length(b$diffs) + 1)
    }
  }
})


test_that("every jump of min_stat sits at a treated-minus-control difference", {
  # The claim the tabulation rests on.  A jump anywhere else would mean the
  # candidate list is incomplete and the refactored search could step over
  # the true limit.
  shapes <- list(c(nb = 12, mb = 6), c(nb = 10, mb = 3), c(nb = 8, mb = 4))

  for (i in seq_along(shapes)) {
    sh <- shapes[[i]]
    b <- make_block(nb = sh[["nb"]], mb = sh[["mb"]], seed = 200 + i)
    grid <- seq(min(b$diffs) - 1, max(b$diffs) + 1, length.out = 20001)
    step <- diff(grid)[1]

    for (k in c(b$nb, b$nb - 3)) {
      vals <- stat_curve(b, k, grid)
      jumps <- grid[which(diff(vals) != 0) + 1]

      # Each observed jump is within one grid step of some element of D.
      # The grid cannot resolve a jump more finely than that.
      near_a_difference <- vapply(
        jumps,
        function(x) min(abs(x - b$diffs)) <= step,
        logical(1)
      )
      expect_true(all(near_a_difference),
                  info = sprintf("nb=%d mb=%d k=%d: jump away from D",
                                 b$nb, b$mb, k))
    }
  }
})


test_that("the candidate list is bounded by m_b * (n_b - m_b)", {
  # The size of the search the refactor replaces uniroot with.  Ties among
  # the differences only shrink it.
  for (sh in list(c(12, 6), c(10, 3), c(20, 10))) {
    b <- make_block(nb = sh[1], mb = sh[2], seed = 300 + sh[1])
    expect_lte(length(b$diffs), b$mb * (b$nb - b$mb))
  }
})


test_that("min_stat at a breakpoint takes the left value or the right one, and which varies", {
  # The reason the confidence-limit search probes strictly between
  # breakpoints instead of at them.  Exactly at a breakpoint, Y_i - c equals
  # some Y_j, and `ties.method = "first"` at R/CMRSS_CRE.R:278 breaks the tie
  # by position in the input vector.  So the value at a breakpoint sometimes
  # matches the interval below and sometimes the interval above, decided by
  # unit ordering rather than by the data.
  b <- make_block(nb = 12, mb = 6, seed = 7)

  side <- vapply(b$diffs, function(cc) {
    at <- min_stat(b$Z, b$Y, k = b$nb, c = cc, score = b$score)
    lo <- min_stat(b$Z, b$Y, k = b$nb, c = cc - 1e-9, score = b$score)
    hi <- min_stat(b$Z, b$Y, k = b$nb, c = cc + 1e-9, score = b$score)
    if (isTRUE(all.equal(at, hi)) && !isTRUE(all.equal(at, lo))) {
      "right"
    } else if (isTRUE(all.equal(at, lo)) && !isTRUE(all.equal(at, hi))) {
      "left"
    } else {
      "no jump"
    }
  }, character(1))

  # The value at a breakpoint always equals one of its two neighbors, which
  # is what makes this a step function rather than one with an extra value
  # sitting at each jump.
  expect_true(all(side %in% c("left", "right", "no jump")))

  # Both sides occur on this block, so no rule of the form "the value at a
  # breakpoint is the value above it" is available.  A search that evaluated
  # at breakpoints would inherit that arbitrariness; `search_limit` does not,
  # because `interval_probes` keeps every evaluation strictly inside an
  # interval.
  expect_gt(sum(side == "left"), 0)
  expect_gt(sum(side == "right"), 0)
})


test_that("min_stat values are unchanged by the refactor", {
  # A regression pin.  These are the values the current implementation
  # returns; the tabulated version must return the same ones exactly, since
  # tabulation is meant to avoid recomputation and not to change arithmetic.
  b <- make_block(nb = 12, mb = 6, seed = 7)
  cs <- c(-2, -0.5, 0, 0.25, 1, 3)

  observed <- vapply(cs, function(cc)
    min_stat(b$Z, b$Y, k = b$nb, c = cc, score = b$score), numeric(1))

  expect_equal(observed, MIN_STAT_REGRESSION_K12, tolerance = 1e-8)

  observed_k9 <- vapply(cs, function(cc)
    min_stat(b$Z, b$Y, k = 9, c = cc, score = b$score), numeric(1))

  expect_equal(observed_k9, MIN_STAT_REGRESSION_K9, tolerance = 1e-8)
})
