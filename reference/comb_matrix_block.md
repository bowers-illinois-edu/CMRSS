# Compute coefficient matrices for optimization

Calculates the coefficient matrices needed for the optimization problem
in the combined stratified rank sum test.

## Usage

``` r
comb_matrix_block(
  Z,
  Y,
  block,
  c,
  methods.list.all,
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

  A list of lists of method specifications.

- scores.list.all:

  Optional pre-computed scores.

- block.sum:

  Optional pre-computed block summary.

## Value

A list with H elements, each containing a list of B matrices. Each
matrix is 2 x (mb+1) containing the number of units and minimum
statistic values for different coverage scenarios.
