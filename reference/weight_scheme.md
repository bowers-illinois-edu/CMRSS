# Compute block weights

Computes weights for combining block-specific statistics according to a
specified weighting scheme.

## Usage

``` r
weight_scheme(block.sum, weight.name = "asymp.opt")
```

## Arguments

- block.sum:

  A block summary object from `summary_block`.

- weight.name:

  Weighting scheme: "asymp.opt" for asymptotically optimal weights, or
  "dis.free" for distribution-free weights.

## Value

A B-dimensional vector of block weights.
