# Generate null distribution of combined stratified statistics

Generates the randomization null distribution of the maximum among
weighted and normalized stratified rank sum statistics.

## Usage

``` r
com_null_dist_block(
  Z,
  block,
  methods.list.all,
  scores.list.all = NULL,
  null.max = 10^4,
  weight,
  block.sum = NULL,
  Z.perm = NULL,
  chunk_size = NULL
)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector.

- block:

  An n-dimensional vector specifying block membership.

- methods.list.all:

  A list of lists of method specifications.

- scores.list.all:

  Optional pre-computed scores.

- null.max:

  Number of permutations.

- weight:

  A B-dimensional vector of block weights.

- block.sum:

  Optional pre-computed block summary.

- Z.perm:

  Optional pre-computed permutation matrix.

- chunk_size:

  Optional number of permutations per chunk when `Z.perm` is `NULL`. If
  `NULL`, chosen adaptively based on `n`.

## Value

A numeric vector of the null distribution (length null.max).
