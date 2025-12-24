# Threshold for hypergeometric correction

Calculates the correction term Delta using the proposed choice of kk.vec
for generalizing inference to larger populations.

## Usage

``` r
threshold_correct_hypergeom(
  N = 100,
  n = 50,
  K.vec = c(60, 80),
  alpha = 0.05,
  ndraw = 10^4,
  hg.draw = NULL,
  tol = 10^(-2)
)
```

## Arguments

- N:

  Population size.

- n:

  Sample size.

- K.vec:

  Vector of quantile ranks in the population.

- alpha:

  Significance level.

- ndraw:

  Number of draws for Monte Carlo approximation.

- hg.draw:

  Optional pre-computed hypergeometric draws.

- tol:

  Tolerance for bisection search.

## Value

A list with:

- kappa:

  The correction factor

- prob:

  The resulting probability

- kk.vec:

  The adjusted sample quantile ranks
