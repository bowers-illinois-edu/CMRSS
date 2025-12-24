# Multivariate hypergeometric correction

Computes the correction term based on the multivariate hypergeometric
distribution for generalizing inference to larger populations.

## Usage

``` r
correct_hypergeom(
  N = 100,
  n = 50,
  K.vec = c(60, 80),
  kk.vec = c(30, 40),
  ndraw = 10^4,
  hg.draw = NULL
)
```

## Arguments

- N:

  Population size.

- n:

  Sample size.

- K.vec:

  Vector of quantile ranks in the population.

- kk.vec:

  Vector of corresponding sample quantile ranks.

- ndraw:

  Number of draws for Monte Carlo approximation.

- hg.draw:

  Optional pre-computed hypergeometric draws.

## Value

Probability value for the correction.
