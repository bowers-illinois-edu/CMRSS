# Berger and Boos (1994) Method

Method using the approach from Berger and Boos (1994) for confidence
intervals that combine treated and control inference.

## Usage

``` r
method_berger_boos(
  Z,
  Y,
  N,
  k_vec,
  treat.method.list = list(name = "Stephenson", s = 6),
  control.method.list = list(name = "Stephenson", s = 6),
  simul = TRUE,
  stat.null = NULL,
  score = NULL,
  Z.perm = NULL,
  alpha = 0.05,
  gamma = 0.5,
  alpha.ratio.treat = 0.5,
  ndraw = 10^5,
  nperm = 10^4,
  tol = 10^(-3)
)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector (1 = treated, 0 =
  control).

- Y:

  An n-dimensional observed outcome vector.

- N:

  Population size for generalization. Set to length(Z) for sample
  inference.

- k_vec:

  Vector of quantile ranks to compute intervals for.

- treat.method.list:

  Method specification for treated units.

- control.method.list:

  Method specification for control units.

- simul:

  Logical; if TRUE (default), compute simultaneous intervals.

- stat.null:

  Optional pre-computed null distribution.

- score:

  Optional pre-computed score vector.

- Z.perm:

  Optional permutation matrix.

- alpha:

  Significance level.

- gamma:

  Proportion of alpha for hypergeometric correction.

- alpha.ratio.treat:

  Proportion of alpha allocated to treated units.

- ndraw:

  Number of Monte Carlo draws for hypergeometric correction.

- nperm:

  Number of permutations for null distribution.

- tol:

  Tolerance for root-finding.

## Value

A data frame with columns k, lower, and upper, containing the
element-wise maximum of treated and control confidence limits.

## Details

This method separately computes confidence intervals from the treated
and control perspectives, then takes the element-wise maximum. This
provides a more conservative but potentially more robust interval.

## References

Berger, R. L. and Boos, D. D. (1994). P Values Maximized Over a
Confidence Set for the Nuisance Parameter. *Journal of the American
Statistical Association*, 89(427), 1012-1016.

## See also

[`method_chen_li`](https://bowers-illinois-edu.github.io/CMRSS/reference/method_chen_li.md)
for the Chen and Li combination method

## Examples

``` r
if (FALSE) { # \dontrun{
data(electric_teachers)
Z <- electric_teachers$TxAny
Y <- electric_teachers$gain
n <- length(Z)

# Using Berger-Boos approach
result <- method_berger_boos(Z, Y,
                             N = n,
                             k_vec = floor(c(0.6, 0.7, 0.8, 0.9) * n),
                             treat.method.list = list(name = "Stephenson", s = 6),
                             control.method.list = list(name = "Stephenson", s = 6),
                             nperm = 10000,
                             alpha = 0.05)
print(result)
} # }
```
