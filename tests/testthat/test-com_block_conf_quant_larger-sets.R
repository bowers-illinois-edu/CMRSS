# Structural tests for the three set values of com_block_conf_quant_larger.
#
# Per the paper (Kim, Bowers, Li, main.tex), the wrapper supports:
# - set = "treat":   m = sum(Z) lower CIs for tau_{treat(k)}, k = 1..m,
#                    obtained by inverting H_{k,c}^treat (eq:H_kc_t).
# - set = "control": (n - m) lower CIs for tau_{control(k)}, k = 1..(n - m),
#                    obtained by swapping Z <- 1 - Z and Y <- -Y (main.tex
#                    lines 479-485 and 940).  Under the swap, individual
#                    effects are preserved via the relabeling
#                    tilde_tau_{treat(k)} = tau_{control(k)}, so no further
#                    sign or index transformation is needed on the output.
# - set = "all":     n CIs for tau_{(k)}, pooled from the treat and control
#                    branches each run at level alpha/2 (Bonferroni, main.tex
#                    lines 502 and 948-950).  The user's alpha is the joint
#                    (1 - alpha) simultaneous coverage level.
#
# These tests pin the wrapper's structural behavior by replicating its body
# with set.seed pre-set so the RNG consumed for the null permutations is
# identical between the wrapper call and the hand-written replication.  They
# do not verify finite-sample coverage; that would require a Monte Carlo
# coverage study and is out of scope here.

test_that("set values return vectors of expected length and numeric type", {
  skip_if_not(solver_available("highs") || solver_available("gurobi"),
              "Neither HiGHS nor Gurobi solver is available")

  set.seed(42)
  s <- 3; n_per <- 8; m_per <- 4
  N <- s * n_per
  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- rnorm(N) + 0.5 * Z

  methods.list.all <- list(
    lapply(1:s, function(i) list(name = "Wilcoxon", scale = FALSE))
  )

  ci_treat <- com_block_conf_quant_larger(
    Z, Y, block, set = "treat",
    methods.list.all = methods.list.all,
    opt.method = "ILP_highs", null.max = 200, tol = 0.1, alpha = 0.10
  )
  ci_control <- com_block_conf_quant_larger(
    Z, Y, block, set = "control",
    methods.list.all = methods.list.all,
    opt.method = "ILP_highs", null.max = 200, tol = 0.1, alpha = 0.10
  )
  ci_all <- com_block_conf_quant_larger(
    Z, Y, block, set = "all",
    methods.list.all = methods.list.all,
    opt.method = "ILP_highs", null.max = 200, tol = 0.1, alpha = 0.10
  )

  expect_length(ci_treat,   sum(Z))
  expect_length(ci_control, N - sum(Z))
  expect_length(ci_all,     N)
  expect_true(all(is.numeric(ci_treat)))
  expect_true(all(is.numeric(ci_control)))
  expect_true(all(is.numeric(ci_all)))
})


test_that("set = 'control' equals com_block_conf_quant_larger_trt(1 - Z, -Y) with paired slicing", {
  skip_if_not(solver_available("highs") || solver_available("gurobi"),
              "Neither HiGHS nor Gurobi solver is available")

  set.seed(42)
  s <- 3; n_per <- 8; m_per <- 4
  N <- s * n_per
  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- rnorm(N) + 0.5 * Z

  methods.list.all <- list(
    lapply(1:s, function(i) list(name = "Wilcoxon", scale = FALSE))
  )
  alpha_lvl <- 0.10

  set.seed(99)
  ci_wrap <- com_block_conf_quant_larger(
    Z, Y, block, set = "control",
    methods.list.all = methods.list.all,
    opt.method = "ILP_highs", null.max = 200, tol = 0.1, alpha = alpha_lvl
  )

  set.seed(99)
  Z_swap <- 1 - Z
  Y_swap <- -Y
  ci_full <- com_block_conf_quant_larger_trt(
    Z_swap, Y_swap, block, k.vec = NULL,
    methods.list.all = methods.list.all,
    opt.method = "ILP_highs", null.max = 200, tol = 0.1, alpha = alpha_lvl
  )
  ci_manual <- ci_full[(N - sum(Z_swap) + 1):N]

  expect_equal(as.numeric(ci_wrap), as.numeric(ci_manual), tolerance = 1e-6)
})


test_that("set = 'all' equals the sorted Bonferroni pool of treat (alpha/2) and control (alpha/2)", {
  skip_if_not(solver_available("highs") || solver_available("gurobi"),
              "Neither HiGHS nor Gurobi solver is available")

  set.seed(42)
  s <- 3; n_per <- 8; m_per <- 4
  N <- s * n_per
  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- rnorm(N) + 0.5 * Z

  methods.list.all <- list(
    lapply(1:s, function(i) list(name = "Wilcoxon", scale = FALSE))
  )
  alpha_lvl <- 0.10

  set.seed(123)
  ci_all <- com_block_conf_quant_larger(
    Z, Y, block, set = "all",
    methods.list.all = methods.list.all,
    opt.method = "ILP_highs", null.max = 200, tol = 0.1, alpha = alpha_lvl
  )

  # Replicate the wrapper's set = "all" body: two _trt calls in order, each at
  # alpha / 2, then sort.  RNG consumption order matches because set.seed is
  # reset to the same value before the manual replication.
  set.seed(123)
  ci_t_full <- com_block_conf_quant_larger_trt(
    Z, Y, block, k.vec = NULL,
    methods.list.all = methods.list.all,
    opt.method = "ILP_highs", null.max = 200, tol = 0.1, alpha = alpha_lvl / 2
  )
  ci_t <- ci_t_full[(N - sum(Z) + 1):N]

  Z_swap <- 1 - Z
  Y_swap <- -Y
  ci_c_full <- com_block_conf_quant_larger_trt(
    Z_swap, Y_swap, block, k.vec = NULL,
    methods.list.all = methods.list.all,
    opt.method = "ILP_highs", null.max = 200, tol = 0.1, alpha = alpha_lvl / 2
  )
  ci_c <- ci_c_full[(N - sum(Z_swap) + 1):N]

  ci_manual <- sort(c(ci_t, ci_c))

  expect_equal(as.numeric(ci_all), as.numeric(ci_manual), tolerance = 1e-6)
})


test_that("set = 'all' uses alpha/2 per branch (treat-half is tighter than set = 'treat' at alpha)", {
  skip_if_not(solver_available("highs") || solver_available("gurobi"),
              "Neither HiGHS nor Gurobi solver is available")

  # If the wrapper passes alpha/2 to each branch (Bonferroni), the treat-half
  # of set = "all" at alpha = 0.10 should equal set = "treat" at alpha = 0.05.
  # Both calls invoke the same _trt(Z, Y, ...) code path.

  set.seed(42)
  s <- 3; n_per <- 8; m_per <- 4
  N <- s * n_per
  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- rnorm(N) + 0.5 * Z

  methods.list.all <- list(
    lapply(1:s, function(i) list(name = "Wilcoxon", scale = FALSE))
  )

  set.seed(777)
  ci_all <- com_block_conf_quant_larger(
    Z, Y, block, set = "all",
    methods.list.all = methods.list.all,
    opt.method = "ILP_highs", null.max = 200, tol = 0.1, alpha = 0.10
  )

  set.seed(777)
  ci_treat_half <- com_block_conf_quant_larger(
    Z, Y, block, set = "treat",
    methods.list.all = methods.list.all,
    opt.method = "ILP_highs", null.max = 200, tol = 0.1, alpha = 0.05
  )

  # The m largest entries of ci_all should be the treat-half (since the
  # control-half values are at most equal to those for symmetric problems
  # and the merge is sorted ascending).  Equivalently: sort(ci_all) contains
  # the treat-half values at known positions.  The cleanest assertion is that
  # ci_treat_half is a subset of ci_all.
  m_treat <- sum(Z)
  expect_true(all(ci_treat_half %in% ci_all),
              info = "treat-half at alpha = 0.05 not found in set = 'all' at alpha = 0.10")
})
