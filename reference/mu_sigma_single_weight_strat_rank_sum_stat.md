# Compute mean and standard deviation of stratified statistic

Calculates the mean and standard deviation of a weighted stratified rank
sum statistic under the null hypothesis.

## Usage

``` r
mu_sigma_single_weight_strat_rank_sum_stat(
  Z,
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

A list with components:

- mu:

  Overall weighted mean

- sigma:

  Overall weighted standard deviation

- mu_weight:

  Block-wise weighted means (mean_block \* weight)

- mu_block:

  Block-wise unweighted means

- sd_block:

  Block-wise weighted standard deviations
