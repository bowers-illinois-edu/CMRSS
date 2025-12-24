# Maximum of stratum-wise standardized statistics

Computes the maximum of standardized statistics among each stratum
across multiple methods.

## Usage

``` r
max_comb_matrix_block_stratum(
  Z,
  Y,
  block,
  c,
  methods.list.all,
  scores.list.all = NULL,
  block.sum = NULL,
  weight
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

- scores.list.all:

  Optional pre-computed scores.

- block.sum:

  Optional pre-computed block summary.

- weight:

  A B-dimensional vector of block weights.

## Value

A list of B matrices, one per block, containing the maximum standardized
test statistic values.
