# Confidence intervals for effect quantiles (experimental)

Computes (1-alpha) simultaneous confidence/prediction intervals (lower
bounds) for effect quantiles among treated, control, or all units using
the RIQITE package.

## Usage

``` r
ci_lower_quantile_exp(
  Z,
  Y,
  treat.method.list = list(name = "Stephenson", s = 6),
  control.method.list = list(name = "Stephenson", s = 6),
  score = NULL,
  stat.null = NULL,
  nperm = 10^4,
  Z.perm = NULL,
  alpha = 0.05,
  set = "treat",
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

- alpha:

  Significance level.

- set:

  Set of quantiles: "treat", "control", or "all".

- alpha.ratio.treat:

  For set="all", proportion of alpha allocated to treated.

- tol:

  Tolerance for root-finding.

## Value

A numeric vector of lower confidence limits.
