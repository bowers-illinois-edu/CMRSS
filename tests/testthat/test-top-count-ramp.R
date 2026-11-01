# Two rank scores that look at the top j ranks alone (CMRSS 0.2.23).
#
# Why they exist. Every nondecreasing score is a nonnegative mixture of
# indicators "the unit's rank is above n - j" (its increments are the
# weights), so a polynomial score of degree r is an average of top-j
# counts over many j. The two scores here are the ends of that
# construction: the top-j count itself, which is the indicator of the
# top j ranks, and a ramp that is zero below rank n - j and rises in
# equal steps to 1 at rank n, so that it still orders the units within
# the top j. Both are specified like the other scores, as a method list
# with a name and a parameter, so that every function of the package
# that builds scores by name (null_dist, min_stat, comb_p_val_cre,
# score_all_blocks for the stratified test) takes them without change.
#
#   list(name = "TopCount", j = 3)   score 1 at ranks n - 2, n - 1, n; 0 below
#   list(name = "Ramp", j = 4)       score 1/4, 2/4, 3/4, 1 at the top four ranks; 0 below
#
# j must be a whole number from 1 to n. The sum of the top-j count over
# the treated units is the number of treated units among the top j
# ranks, whose distribution under complete randomization is
# hypergeometric, which is the check that ties the score to something
# known.

test_that("the top-j count is the indicator of the top j ranks and the ramp rises through them", {
  n <- 10
  expect_equal(rank_score(n, list(name = "TopCount", j = 3)), c(rep(0, 7), 1, 1, 1))
  expect_equal(rank_score(n, list(name = "Ramp", j = 4)), c(rep(0, 6), 0.25, 0.5, 0.75, 1))
  # j = n: the count is the constant 1 and the ramp is the standardized
  # Wilcoxon score i / n
  expect_equal(rank_score(n, list(name = "TopCount", j = n)), rep(1, n))
  expect_equal(rank_score(n, list(name = "Ramp", j = n)), (1:n) / n)
  # j = 1: both are the indicator of the top unit, which is also
  # Stephenson's score with s = n divided by its top value
  expect_equal(rank_score(n, list(name = "TopCount", j = 1)), c(rep(0, n - 1), 1))
  expect_equal(rank_score(n, list(name = "Ramp", j = 1)), c(rep(0, n - 1), 1))
  expect_equal(rank_score(n, list(name = "TopCount", j = 1)), rank_score(n, list(name = "Stephenson", s = n)))
})

test_that("both scores are nondecreasing in the rank for every j, which the worst-case construction needs", {
  n <- 25
  for (j in 1:n) {
    expect_true(all(diff(rank_score(n, list(name = "TopCount", j = j))) >= 0))
    expect_true(all(diff(rank_score(n, list(name = "Ramp", j = j))) >= 0))
  }
})

test_that("scale = TRUE standardizes both, and the scores are plain numeric vectors", {
  n <- 40
  for (nm in c("TopCount", "Ramp")) {
    s <- rank_score(n, list(name = nm, j = 5, scale = TRUE))
    expect_type(s, "double")
    expect_length(s, n)
    expect_equal(mean(s), 0, tolerance = 1e-10)
    expect_equal(sd(s), 1, tolerance = 1e-10)
  }
})

test_that("a j outside 1..n, a missing j, and an unknown score name are refused with a message that names the fault", {
  expect_error(rank_score(10, list(name = "TopCount", j = 0)), "j")
  expect_error(rank_score(10, list(name = "TopCount", j = 11)), "j")
  expect_error(rank_score(10, list(name = "Ramp", j = 2.5)), "j")
  expect_error(rank_score(10, list(name = "Ramp")), "j")
  expect_error(rank_score(10, list(name = "Hinge", j = 2)), "Hinge")
})

test_that("the treated sum of the top-j count is the number of treated units among the top j ranks", {
  set.seed(1)
  n <- 30; m <- 12; j <- 5
  Z <- sample(rep(c(1, 0), c(m, n - m)))
  Y <- rnorm(n)
  score <- rank_score(n, list(name = "TopCount", j = j))
  r <- rank(Y, ties.method = "first")
  expect_equal(sum(score[r[Z == 1]]), sum(r[Z == 1] > n - j))
  # and at k = n (no exempt unit) min_stat() returns the same number
  expect_equal(min_stat(Z, Y, k = n, c = 0, method.list = list(name = "TopCount", j = j)), sum(r[Z == 1] > n - j))
})

test_that("under complete randomization the top-j count is hypergeometric: its exact null distribution over every assignment matches dhyper", {
  n <- 8; m <- 3; j <- 2
  # every assignment of 3 treated among 8, as the columns of Z.perm
  combs <- combn(n, m)
  Z.perm <- apply(combs, 2, function(idx) { z <- numeric(n); z[idx] <- 1; z })
  stat <- null_dist(n, m, method.list = list(name = "TopCount", j = j), Z.perm = Z.perm)
  expect_length(stat, ncol(combs))
  observed <- as.numeric(table(factor(stat, levels = 0:j))) / ncol(combs)
  expect_equal(observed, dhyper(0:j, m = j, n = n - j, k = m))
})

test_that("the worst case over the exempt units sends the treated units with the largest outcomes to the bottom, for both scores", {
  # Ten units, five treated, the treated units holding ranks 6 to 10.
  # Under "at most two of the five treated gained" (k = n - 2 in
  # min_stat's convention, the two largest treated outcomes exempt), the
  # two exempt units go below every other unit, the remaining treated
  # units take ranks 8, 9, 10 (the controls move down two places in the
  # ranking: ranks 3 to 7), and the statistic is the sum of the scores at
  # those three ranks.
  n <- 10
  Y <- c(1:5, 6:10)
  Z <- c(rep(0, 5), rep(1, 5))
  for (ml in list(list(name = "TopCount", j = 4), list(name = "Ramp", j = 4))) {
    score <- rank_score(n, ml)
    expect_equal(min_stat(Z, Y, k = n, c = 0, method.list = ml), sum(score[6:10]))
    expect_equal(min_stat(Z, Y, k = n - 2, c = 0, method.list = ml), sum(score[c(1, 2, 8, 9, 10)]))
  }
})

test_that("comb_p_val_cre takes a set of top-j counts and ramps and returns a p-value in [0, 1], small when the treated hold the top ranks", {
  set.seed(2)
  n <- 40; m <- 20
  Z <- sample(rep(c(1, 0), c(m, n - m)))
  Y <- rnorm(n) + 3 * Z          # every treated unit gained, so the treated hold the top ranks
  Z.perm <- assign_CRE(n, m, 2000)
  methods <- list(list(name = "TopCount", j = 5), list(name = "TopCount", j = 10), list(name = "Ramp", j = 10))
  p <- suppressWarnings(comb_p_val_cre(Z, Y, k = n, c = 0, methods.list = methods, Z.perm = Z.perm, nperm = 2000))
  expect_true(p >= 0 && p <= 1)
  expect_lt(p, 0.05)
  # with no effect the same test keeps its level over repeated experiments
  rej <- replicate(200, {
    Zr <- sample(rep(c(1, 0), c(m, n - m)))
    suppressWarnings(comb_p_val_cre(Zr, rnorm(n), k = n, c = 0, methods.list = methods, Z.perm = Z.perm, nperm = 2000)) <= 0.10
  })
  expect_lt(mean(rej), 0.17)   # 0.10 plus about three binomial standard errors of 0.021
})

test_that("pval_comb_block takes a top-j count with its own j in each stratum", {
  skip_if_not(solver_available("highs"), "HiGHS is not installed")
  set.seed(3)
  nb <- c(10, 20); mb <- c(4, 8)
  block <- factor(rep(1:2, nb))
  Z <- unlist(lapply(1:2, function(b) sample(rep(c(1, 0), c(mb[b], nb[b] - mb[b])))))
  Y <- rnorm(sum(nb)) + 2 * Z
  # statistic 1 counts the top 2 ranks of the stratum of 10 and the top 4
  # of the stratum of 20; statistic 2 is the ramp over the same ranks
  methods <- list(list(list(name = "TopCount", j = 2), list(name = "TopCount", j = 4)),
                  list(list(name = "Ramp", j = 2), list(name = "Ramp", j = 4)))
  r <- suppressWarnings(pval_comb_block(Z, Y, k = sum(mb), c = 0, block = block, methods.list.all = methods,
                                        null.max = 1000, statistic = TRUE, opt.method = "ILP_highs"))
  expect_true(r[["p.value"]] >= 0 && r[["p.value"]] <= 1)
  expect_lt(r[["p.value"]], 0.10)
  # and below k = m the p-value cannot be smaller than at k = m
  r2 <- suppressWarnings(pval_comb_block(Z, Y, k = sum(mb) - 2, c = 0, block = block, methods.list.all = methods,
                                         null.max = 1000, statistic = TRUE, opt.method = "ILP_highs"))
  expect_gte(r2[["p.value"]], r[["p.value"]])
})
