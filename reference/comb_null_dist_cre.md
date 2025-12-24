# Null distribution for minimum p-values

Computes the distribution of the minimum p-value, particularly Monte
Carlo samples from the combined null distribution F in Theorem 1 of the
paper.

## Usage

``` r
comb_null_dist_cre(
  n,
  m,
  methods.list,
  Z.perm = NULL,
  nperm = 10^4,
  stat.null.mult = NULL
)
```

## Arguments

- n:

  Total number of units.

- m:

  Number of treated units.

- methods.list:

  A list of method specifications.

- Z.perm:

  Optional permutation matrix.

- nperm:

  Number of permutations for null distribution.

- stat.null.mult:

  Optional pre-computed null distribution matrix.

## Value

A numeric vector of minimum p-values under the null.
