# Simultaneous inference for multiple quantiles on CRE (treated units)

Helper function for computing lower limits of prediction intervals for
quantiles across the treated units using combined p-value.

## Usage

``` r
com_conf_quant_larger_trt(
  Z,
  Y,
  methods.list = NULL,
  nperm = 10^4,
  k.vec = NULL,
  Z.perm = NULL,
  alpha = 0.05,
  tol = 10^(-3),
  ind.sort.treat = NULL
)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector.

- Y:

  An n-dimensional observed outcome vector.

- methods.list:

  A list of method specifications.

- nperm:

  Number of permutations for null distribution.

- k.vec:

  Vector of quantile indices to compute intervals for.

- Z.perm:

  Optional permutation matrix.

- alpha:

  Significance level.

- tol:

  Tolerance for root-finding.

- ind.sort.treat:

  Optional sorted indices of treated units.

## Value

A numeric vector of lower confidence limits.
