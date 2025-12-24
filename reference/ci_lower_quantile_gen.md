# Generalized confidence intervals with hypergeometric correction

Internal function for computing confidence intervals with correction for
generalizing to a larger population.

## Usage

``` r
ci_lower_quantile_gen(
  Z,
  Y,
  N = 2 * length(Z),
  K.vec = ceiling(c(0.6, 0.7, 0.8, 0.9) * N),
  alpha = 0.05,
  gamma = 0.5,
  ndraw = 10^4,
  treat.method.list = list(name = "Stephenson", s = 6),
  control.method.list = list(name = "Stephenson", s = 6),
  score = NULL,
  stat.null = NULL,
  nperm = 10^4,
  Z.perm = NULL,
  set = "all",
  alpha.ratio.treat = 0.5,
  tol = 10^(-3)
)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector (1 = treated, 0 =
  control).

- Y:

  An n-dimensional observed outcome vector.

- N:

  Population size.

- K.vec:

  Vector of quantile ranks in the population.

- alpha:

  Significance level.

- gamma:

  Proportion of alpha for hypergeometric correction.

- ndraw:

  Number of Monte Carlo draws.

- treat.method.list:

  Method specification for treated units.

- control.method.list:

  Method specification for control units.

- score:

  Optional pre-computed score vector.

- stat.null:

  Optional pre-computed null distribution.

- nperm:

  Number of permutations for null distribution.

- Z.perm:

  Optional permutation matrix.

- set:

  Set of quantiles: "treat", "control", or "all".

- alpha.ratio.treat:

  For set="all", proportion of alpha allocated to treated.

- tol:

  Tolerance for root-finding.

## Value

A data frame with columns k, lower, and upper.
