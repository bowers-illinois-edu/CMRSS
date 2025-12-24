# Compute weighted stratified rank sum statistic

Calculates a single weighted stratified rank sum statistic combining
block-specific statistics.

## Usage

``` r
single_weight_strat_rank_sum_stat(
  Z,
  Y,
  block,
  method.list.all,
  score.list.all = NULL,
  weight,
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

- method.list.all:

  A list of method specifications, one per block.

- score.list.all:

  Optional pre-computed scores.

- weight:

  A B-dimensional vector of block weights.

- block.sum:

  Optional pre-computed block summary.

## Value

A scalar test statistic value.
