# Chen and Li Combination Method

Method combining treated and control inference using Bonferroni
correction, based on Chen and Li (2024).

## Usage

``` r
method_chen_li(
  Z,
  Y,
  N,
  k_vec,
  control.method.list = list(name = "Stephenson", s = 6),
  treat.method.list = list(name = "Stephenson", s = 6),
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

- control.method.list:

  Method specification for control units.

- treat.method.list:

  Method specification for treated units.

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

A data frame with columns k, lower, and upper.

## Details

This method combines inference from both treated and control
perspectives using a Bonferroni-style correction. It provides confidence
intervals for all effect quantiles by separately computing prediction
intervals for treated and control units, then combining them.

When N \> n (the sample size), the method incorporates a hypergeometric
correction to generalize inference to the larger population.

## References

Chen, H. and Li, X. (2024). Randomization Inference for Quantile
Treatment Effects with Varying Numbers of Treated and Control Units.

## See also

[`method_caughey`](https://bowers-illinois-edu.github.io/CMRSS/reference/method_caughey.md)
for single-direction inference,
[`com_conf_quant_larger_cre`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_conf_quant_larger_cre.md)
for the CMRSS combined method

## Examples

``` r
if (FALSE) { # \dontrun{
data(electric_teachers)
Z <- electric_teachers$TxAny
Y <- electric_teachers$gain
n <- length(Z)

# Confidence intervals for all effect quantiles
result <- method_chen_li(Z, Y,
                         N = n,  # sample inference
                         k_vec = 1:n,
                         treat.method.list = list(name = "Stephenson", s = 6),
                         control.method.list = list(name = "Stephenson", s = 6),
                         nperm = 10000,
                         alpha = 0.10)
print(head(result))
} # }
```
