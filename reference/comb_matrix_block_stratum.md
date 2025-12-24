# Compute stratum-wise standardized coefficient matrices

Calculates the stratum-wise standardized test statistic values for the
stratum-level combination method (comb.method = 2).

## Usage

``` r
comb_matrix_block_stratum(
  Z,
  Y,
  block,
  c,
  methods.list.all,
  weight,
  scores.list.all = NULL,
  block.sum = NULL
)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector.

- Y:

  An n-dimensional observed outcome vector.

- block:

  An n-dimensional vector specifying block membership.

- c:

  Threshold for the null hypothesis.

- methods.list.all:

  A list of H method specifications.

- weight:

  A B-dimensional vector of block weights.

- scores.list.all:

  Optional pre-computed scores.

- block.sum:

  Optional pre-computed block summary.

## Value

A list with H elements, each containing B matrices with stratum-wise
standardized test statistics.
