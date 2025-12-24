# Calculate minimum p-value from multiple rank sum statistics

Computes the minimum of p-values calibrated by the distribution of
minimum of tail probabilities.

## Usage

``` r
min_p_multiple_rank_sum(
  Z,
  Y,
  k,
  c,
  methods.list,
  Z.perm = NULL,
  nperm,
  stat.null.mult = NULL
)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector.

- Y:

  An n-dimensional observed outcome vector.

- k:

  Quantile index for the hypothesis.

- c:

  Threshold for the null hypothesis.

- methods.list:

  A list of method specifications.

- Z.perm:

  Optional permutation matrix.

- nperm:

  Number of permutations for null distribution.

- stat.null.mult:

  Optional pre-computed null distribution matrix.

## Value

A numeric minimum p-value.
