# Generate complete randomization assignments

Generates a matrix of permuted treatment assignments for completely
randomized experiments.

## Usage

``` r
assign_CRE(n, m, nperm)
```

## Arguments

- n:

  Total number of units.

- m:

  Number of treated units.

- nperm:

  Number of permutations. If Inf, generates all possible assignments.

## Value

An n x nperm matrix where each column is a permuted treatment
assignment.
