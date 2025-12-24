# Block-wise means and standard deviations for multiple methods

Computes the mean and standard deviation for each block under multiple
rank score methods.

## Usage

``` r
mu_sd_block(
  Z,
  block,
  methods.list.all,
  scores.list.all = NULL,
  weight,
  block.sum = NULL
)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector.

- block:

  An n-dimensional vector specifying block membership.

- methods.list.all:

  A list of H method specifications, each containing B block-specific
  methods.

- scores.list.all:

  Optional pre-computed scores.

- weight:

  A B-dimensional vector of block weights.

- block.sum:

  Optional pre-computed block summary.

## Value

A list with:

- mu_list:

  List of B vectors, each of length H with block-wise means

- sig_list:

  List of B vectors, each of length H with block-wise SDs
