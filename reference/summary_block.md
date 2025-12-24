# Summarize block structure

Computes summary statistics for a stratified randomized experiment,
including the number of blocks, units per block, and treated units per
block.

## Usage

``` r
summary_block(Z, block)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector.

- block:

  An n-dimensional vector or factor specifying block membership.

## Value

A list with components:

- block:

  The block factor

- B:

  Number of blocks

- nb:

  Vector of block sizes

- mb:

  Vector of treated units per block

- mb_ctrl:

  Vector of control units per block

- units.block:

  List of unit indices per block

- block.levels:

  Levels of the block factor
