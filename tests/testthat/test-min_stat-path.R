# Tests for computing a whole block's min_stat sequence from one ranking.
#
# `comb_matrix_block` (R/CMRSS_SRE.R) needs, for each block and each rank
# score function, the value of `min_stat` at every exempt count
# ii = 0, 1, ..., m_b.  It gets them with m_b + 1 separate calls, and each
# call runs `sort_treat` and a full `rank` over the block.  The profiling in
# memos/memo_cmrss_reuse_2026-09-07.qmd found `min_stat` at 61 percent of a
# confidence-interval run and `rank`, called from inside it, at 48.
#
# Those m_b + 1 values are one ranking's worth of information.  Write r0 for
# the ranks of the ii = 0 vector, where treated units sit at Y_i - c and
# controls at Y_j.  The null exempts treated units in order of decreasing
# outcome (`sort_treat` sorts ascending and `min_stat` takes the tail), and
# the treated units' order in r0 is their outcome order, since every treated
# unit shifts by the same c.  So the ii exempted units are exactly the ii
# treated units with the largest r0.  Two consequences:
#
#   - No remaining unit has an exempted unit below it, so every remaining
#     rank is its r0 shifted up by ii.
#   - The ii exempted units take ranks 1..ii in some order, and all of them
#     are treated, so they contribute sum(score[1:ii]) whatever that order
#     is.  Their internal order does not have to be worked out.
#
# Which gives
#
#     min_stat(ii) = sum(score[1:ii]) + sum over remaining treated of
#                    score[ii + r0(unit)].
#
# These tests check that identity against `min_stat` itself rather than
# assuming it, on random designs with ties and at exact breakpoints, since
# those are where a rank identity is most likely to break.

skip_unless_path <- function() {
  skip_if_not(
    exists("min_stat_path", where = asNamespace("CMRSS"), inherits = FALSE),
    "min_stat_path() not implemented yet"
  )
}

# The loop the identity is meant to replace, written out plainly.
min_stat_loop <- function(Zb, Yb, c, score) {
  nb <- length(Yb)
  mb <- sum(Zb)
  vapply(0:mb, function(ii) min_stat(Zb, Yb, nb - ii, c, score = score),
         numeric(1))
}


test_that("min_stat_path reproduces the loop it replaces", {
  skip_unless_path()

  set.seed(11)
  worst <- 0
  cases <- 0

  for (rep in 1:200) {
    nb <- sample(4:16, 1)
    mb <- sample(1:(nb - 1), 1)
    Zb <- sample(rep(0:1, times = c(nb - mb, mb)))
    # Rounding to two digits forces ties, both among outcomes and between a
    # shifted treated outcome and a control outcome.
    Yb <- round(rnorm(nb), 2)
    r <- sample(c(1, 2, 3, 6, 10), 1)
    score <- rank_score(nb, list(name = "Polynomial", r = r,
                                 std = TRUE, scale = FALSE))

    # Thresholds include exact breakpoints, where Y_i - c equals some Y_j
    # and ties.method = "first" decides the ranking by input order.
    d <- as.vector(outer(Yb[Zb == 1], Yb[Zb == 0], "-"))
    thresholds <- c(sample(d, min(3, length(d))), 0, rnorm(2))

    for (cc in thresholds) {
      worst <- max(worst, max(abs(
        min_stat_path(Zb, Yb, cc, score) - min_stat_loop(Zb, Yb, cc, score)
      )))
      cases <- cases + 1
    }
  }

  expect_gt(cases, 1000)
  expect_lt(worst, 1e-10)
})


test_that("min_stat_path handles blocks with one treated unit and with one control", {
  skip_unless_path()

  # The two ends of the range, where the exempt set is almost empty or
  # almost everything.
  for (shape in list(c(nb = 6, mb = 1), c(nb = 6, mb = 5))) {
    nb <- shape[["nb"]]
    mb <- shape[["mb"]]
    set.seed(nb * 100 + mb)
    Zb <- sample(rep(0:1, times = c(nb - mb, mb)))
    Yb <- round(rnorm(nb), 3)
    score <- rank_score(nb, list(name = "Polynomial", r = 3,
                                 std = TRUE, scale = FALSE))

    expect_equal(min_stat_path(Zb, Yb, 0.25, score),
                 min_stat_loop(Zb, Yb, 0.25, score))
    expect_length(min_stat_path(Zb, Yb, 0.25, score), mb + 1)
  }
})


test_that("comb_matrix_block returns the values it returned before", {
  # A regression pin on the consumer rather than the helper. Unlike the
  # binary-search change, this one is meant to leave every number alone: it
  # computes the same quantities from one ranking instead of m_b + 1.
  Z <- c(0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0)
  Y <- c(-1.971, 0.821, -0.831, -0.194, -1.366, -1.023,
         0.487, 1.211, 0.638, -0.560, 0.740, -1.815)
  block <- factor(rep(1:2, each = 6))
  methods.list.all <- lapply(c(2, 6), function(r)
    lapply(1:2, function(i)
      list(name = "Polynomial", r = r, std = TRUE, scale = FALSE)))

  block.sum <- summary_block(Z, block)
  scores.list.all <- lapply(methods.list.all,
                            function(m) score_all_blocks(block.sum$nb, m))

  cf <- comb_matrix_block(Z, Y, block, 0.3, methods.list.all,
                          scores.list.all, block.sum)

  # Captured 2026-09-07 at commit cb69196.
  expected <- list(
    list(c(1.571428571, 1.142857143, 1, 0.8571428571),
         c(1.571428571, 1.142857143, 1, 0.8571428571)),
    list(c(0.4790265961, 0.07544475516, 0.06289046231, 0.01642172904),
         c(0.4790265961, 0.07544475516, 0.06289046231, 0.01642172904))
  )

  for (h in 1:2) {
    for (b in 1:2) {
      # Row 1 is the exempt count, row 2 the statistic.
      expect_equal(cf[[h]][[b]][1, ], 0:3)
      expect_equal(cf[[h]][[b]][2, ], expected[[h]][[b]], tolerance = 1e-8)
    }
  }
})
