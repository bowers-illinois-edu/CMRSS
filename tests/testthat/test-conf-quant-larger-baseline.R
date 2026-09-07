# Baseline for the confidence limits that `com_block_conf_quant_larger`
# returns today, recorded before the min_stat tabulation changes how the
# search for each limit is done.
#
# What the refactor is allowed to change, and what it is not.  Tabulating
# min_stat does not change min_stat's arithmetic, so the statistic at any
# given c must come out identical; that is pinned in
# test-min_stat-breakpoints.R.  What does change is where the search stops.
# Today `com_block_conf_quant_larger_trt` runs `uniroot` and then four
# uncapped `while (f(c.sol) <= 0) c.sol <- c.sol - tol` loops, so a limit is
# located only to within `tol`.  The evidence is visible in the numbers
# below: with tol = 0.05 every finite limit is a multiple of 0.05.  After
# the refactor each limit is an element of the finite candidate set
# { Y_i - Y_j : Z_i = 1, Z_j = 0 } and is exact.
#
# So the correct test is not that the limits are unchanged.  It is that they
# move by no more than `tol`, which is the accuracy the old search ever
# claimed.  A move larger than that is a bug in the refactor; a move smaller
# than that is the improvement being bought.
#
# Note on the existing suite: the assertions in
# test-com_block_conf_quant_larger-sets.R compare the wrapper's output
# against a manual call of the same underlying function, so both sides move
# together and those tests survive the refactor untouched. Nothing in the
# suite pinned a limit against a stored number before this file.
#
# Note on the paper: reproducing the numbers in the submitted manuscript is
# a separate exercise. Those came from CMRSS 0.2.5 at commit 1371f315 under
# Gurobi, driven by the scripts in ~/repos/combined_stephenson_tests/code/,
# and they are not checked here. Under the R&R constraint in HANDOFF.md that
# replication run happens before any change that moves a published number.

# The design, fixed so the baseline is reproducible. Small on purpose: this
# runs in about half a second, where the profiling memo's n = 200 design
# takes 22.
baseline_design <- function() {
  set.seed(20260907)
  B <- 3
  per <- 8
  N <- B * per
  block <- factor(rep(seq_len(B), each = per))
  Z <- unlist(lapply(split(seq_len(N), block),
                     function(ix) sample(rep(0:1, each = per / 2))))
  Y <- round(rnorm(N) + 0.5 * Z, 4)
  methods.list.all <- lapply(c(2, 6), function(r)
    lapply(seq_len(B), function(i)
      list(name = "Polynomial", r = r, std = TRUE, scale = FALSE)))
  list(Z = Z, Y = Y, block = block, methods.list.all = methods.list.all,
       tol = 0.05, alpha = 0.1, null.max = 500)
}

# Captured 2026-09-07 from the implementation at commit 455f72f.
CI_TREAT_BASELINE <- c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,
                       -1.400000, -0.700000, -0.450000, 0.150000)


test_that("the treated-unit confidence limits stay within tol of the baseline", {
  skip_if_not(solver_available("highs"), "HiGHS not available")

  d <- baseline_design()
  set.seed(20260907)
  ci <- com_block_conf_quant_larger(
    d$Z, d$Y, d$block, set = "treat",
    methods.list.all = d$methods.list.all,
    opt.method = "ILP_highs",
    null.max = d$null.max, tol = d$tol, alpha = d$alpha
  )

  expect_length(ci, sum(d$Z))

  # The infinite limits mark quantiles the data cannot bound below. Which
  # ones those are is a property of the design and the alpha, not of the
  # search, so the refactor must not change the pattern.
  expect_identical(is.infinite(as.numeric(ci)),
                   is.infinite(CI_TREAT_BASELINE))

  # Each finite limit moves by at most tol.
  finite <- is.finite(CI_TREAT_BASELINE)
  expect_lte(max(abs(as.numeric(ci)[finite] - CI_TREAT_BASELINE[finite])),
             d$tol)
})


test_that("today's finite limits sit on the tol grid rather than on data values", {
  # The imprecision the refactor removes, stated as a fact about the
  # baseline rather than as a description of the algorithm. Every finite
  # limit is a multiple of tol = 0.05, which is what a search that steps by
  # tol produces and what a search over the differences Y_i - Y_j would not.
  finite <- CI_TREAT_BASELINE[is.finite(CI_TREAT_BASELINE)]
  on_grid <- abs(finite / 0.05 - round(finite / 0.05)) < 1e-8
  expect_true(all(on_grid))
})


test_that("after the refactor each finite limit is one of the candidate values", {
  # The check that only means something once the tabulation exists. It is
  # the positive statement of what the refactor buys: the limit is not near
  # a breakpoint, it is a breakpoint.
  skip_if_not(
    exists("min_stat_breakpoints", where = asNamespace("CMRSS"),
           inherits = FALSE),
    "min_stat_breakpoints() not implemented yet"
  )
  skip_if_not(solver_available("highs"), "HiGHS not available")

  d <- baseline_design()
  set.seed(20260907)
  ci <- com_block_conf_quant_larger(
    d$Z, d$Y, d$block, set = "treat",
    methods.list.all = d$methods.list.all,
    opt.method = "ILP_highs",
    null.max = d$null.max, tol = d$tol, alpha = d$alpha
  )

  # Pool the candidate sets across blocks: a limit for the whole design has
  # to be a breakpoint of some block.
  candidates <- unlist(lapply(split(seq_along(d$Z), d$block), function(ix) {
    as.vector(outer(d$Y[ix][d$Z[ix] == 1], d$Y[ix][d$Z[ix] == 0], "-"))
  }))

  finite_limits <- as.numeric(ci)[is.finite(as.numeric(ci))]
  distance_to_nearest <- vapply(
    finite_limits,
    function(x) min(abs(x - candidates)),
    numeric(1)
  )
  expect_lt(max(distance_to_nearest), 1e-8)
})
