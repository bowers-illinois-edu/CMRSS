# Calculate p-value for single method

Computes the p-value assuming at most n-k treated units have effects \>
c, using a single rank sum statistic.

## Usage

``` r
pval_cre(
  Z,
  Y,
  k,
  c,
  method.list,
  score = NULL,
  stat.null = NULL,
  nperm = 10^3,
  Z.perm = NULL,
  ind.sort.treat = NULL
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

- method.list:

  A list specifying the rank score method.

- score:

  Optional pre-computed score vector.

- stat.null:

  Optional pre-computed null distribution.

- nperm:

  Number of permutations for null distribution.

- Z.perm:

  Optional permutation matrix.

- ind.sort.treat:

  Optional sorted indices of treated units.

## Value

A numeric p-value.
