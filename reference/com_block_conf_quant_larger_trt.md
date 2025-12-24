# Helper function for Simultaneous Inference for multiple quantiles on SRE

Output is a lower limit of prediction intervals for prespecified
quantiles.

## Usage

``` r
com_block_conf_quant_larger_trt(
  Z,
  Y,
  block,
  k.vec = NULL,
  methods.list.all,
  weight.name = "asymp.opt",
  opt.method = "ILP_auto",
  comb.method = 1,
  stat.null = NULL,
  null.max = 10^4,
  Z.perm = NULL,
  tol = 0.01,
  alpha = 0.1
)
```

## Arguments

- Z:

  An \\n\\ dimensional treatment assignment vector.

- Y:

  An \\n\\ dimensional observed outcome vector.

- block:

  An \\n\\ dimensional vector specifying block of each units.

- k.vec:

  Optional vector of specific quantile indices to compute. If NULL,
  computes for all quantiles.

- methods.list.all:

  A list of lists of lists. Corresponds to the method for each stratum
  for each different stratified rank sum statistic.

- weight.name:

  Weighting method to be implemented. If "asymp.opt", asymptotically
  optimal scheme under a class of local alternatives is adjusted, if
  "dist.free", design-free scheme is adjusted.

- opt.method:

  Optimization method. Options include:

  - "ILP_gurobi": Integer linear programming with Gurobi solver

  - "LP_gurobi": Linear programming relaxation with Gurobi solver

  - "ILP_highs": Integer linear programming with HiGHS solver
    (open-source)

  - "LP_highs": Linear programming relaxation with HiGHS solver

  - "ILP" or "ILP_auto": Integer LP with auto-selected solver

  - "LP" or "LP_auto": Linear programming with auto-selected solver

  HiGHS is recommended for users without a Gurobi license. Both solvers
  produce equivalent results.

- comb.method:

  Integer specifying the combination method:

  - 1: Combine statistics across strata first, then take maximum across
    methods (default)

  - 2: Take maximum across methods within each stratum first, then
    combine across strata

- stat.null:

  An vector whose empirical distribution approximates the randomization
  distribution of the combined stratified rank sum statistic.

- null.max:

  A positive integer representing the number of permutations for
  approximating the randomization distribution of the rank sum
  statistic.

- tol:

  Tolerance for root-finding algorithm.

- alpha:

  Significance level for confidence intervals.

## Value

Vector of lower confidence limits for the specified quantiles.
