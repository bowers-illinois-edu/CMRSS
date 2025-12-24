# Generate block-randomized treatment assignments

Generates permuted treatment assignments that respect the block
structure of a stratified randomized experiment.

## Usage

``` r
assign_block(block.sum, null.max = 10^4)
```

## Arguments

- block.sum:

  A block summary object from `summary_block`.

- null.max:

  Number of permutations to generate.

## Value

An n x null.max matrix of permuted treatment assignments.
