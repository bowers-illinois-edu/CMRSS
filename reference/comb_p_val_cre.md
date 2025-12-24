# Combined p-value for quantile treatment effects in CRE

Tests the null hypothesis \\H_0: \tau\_{(k)} \leq c\\, where
\\\tau\_{(k)}\\ denotes the individual treatment effect at rank k. The
test combines multiple rank sum statistics to improve power.

## Usage

``` r
comb_p_val_cre(
  Z,
  Y,
  k,
  c,
  methods.list,
  Z.perm = NULL,
  nperm,
  stat.null.mult = NULL
)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector (1 = treated, 0 =
  control).

- Y:

  An n-dimensional observed outcome vector.

- k:

  An integer between 1 and n specifying which quantile of individual
  effect is of interest.

- c:

  A numeric value specifying the threshold for the null hypothesis.

- methods.list:

  A list of method specifications for the rank sum statistics. Each
  element should be a list with:

  - `name`: "Wilcoxon", "Stephenson", or "Polynomial"

  - `s`: (for Stephenson) parameter s

  - `r`: (for Polynomial) power parameter

  - `std`: (for Polynomial) logical, use Puri normalization

  - `scale`: logical, standardize scores

- Z.perm:

  An n x nperm matrix of permuted treatment assignments for
  approximating the null distribution. If NULL, generated automatically.

- nperm:

  Number of permutations for approximating the null distribution.

- stat.null.mult:

  A matrix whose empirical distribution approximates the randomization
  distribution of multiple rank statistics. If NULL, computed from
  Z.perm.

## Value

A numeric p-value for testing the specified null hypothesis.

## Details

Calculate a valid p-value, based on multiple rank sum statistics, for
testing the null hypothesis about quantiles of individual treatment
effects in completely randomized experiments (CRE).

## See also

[`pval_comb_block`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md)
for stratified experiments

## Examples

``` r
# Simple example with Wilcoxon and Stephenson statistics
set.seed(123)
n <- 30
Z <- sample(c(rep(1, 15), rep(0, 15)))
Y <- rnorm(n) + 0.5 * Z  # Treatment effect of 0.5

# Define methods: Wilcoxon and Stephenson with s=3
methods.list <- list(
  list(name = "Wilcoxon", scale = FALSE),
  list(name = "Stephenson", s = 3, scale = FALSE)
)

# Test if the 80th percentile effect is <= 0
k <- floor(0.8 * n)
pval <- comb_p_val_cre(Z, Y, k, c = 0, methods.list, nperm = 1000)
```
