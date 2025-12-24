# Compute means and standard deviations for multiple statistics

Calculates the mean and standard deviation for multiple weighted
stratified rank sum statistics.

## Usage

``` r
mu_sigma_list(
  Z,
  block,
  weight,
  methods.list.all,
  scores.list.all = NULL,
  block.sum = NULL
)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector.

- block:

  An n-dimensional vector specifying block membership.

- weight:

  A B-dimensional vector of block weights.

- methods.list.all:

  A list of lists of method specifications.

- scores.list.all:

  Optional pre-computed scores.

- block.sum:

  Optional pre-computed block summary.

## Value

A list with components `mu` and `sigma`, each an H-dimensional vector
for H statistics.
