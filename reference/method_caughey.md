# Caughey et al. (2023) Method

Wrapper for the method from Caughey et al. (2023) using the RIQITE
package. This provides confidence intervals for effect quantiles using a
single rank statistic.

## Usage

``` r
method_caughey(
  Z,
  Y,
  k_vec,
  method.list = list(name = "Wilcoxon"),
  nperm = 10^4,
  alpha = 0.05
)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector (1 = treated, 0 =
  control).

- Y:

  An n-dimensional observed outcome vector.

- k_vec:

  Vector of quantile ranks to compute intervals for.

- method.list:

  A list specifying the rank statistic method:

  - `name`: "Wilcoxon" or "Stephenson"

  - `s`: (for Stephenson) the parameter s

- nperm:

  Number of permutations for the null distribution.

- alpha:

  Significance level.

## Value

A data frame with columns k, lower, and upper.

## Details

This method is based on Caughey, Dafoe, Li, and Miratrix (2023) and uses
the RIQITE package implementation. It provides simultaneous confidence
intervals for specified quantiles using a single rank statistic.

## References

Caughey, D., Dafoe, A., Li, X., and Miratrix, L. (2023). Randomization
Inference for Treatment Effect Quantiles. *Journal of the American
Statistical Association*.

## See also

[`method_chen_li`](https://bowers-illinois-edu.github.io/CMRSS/reference/method_chen_li.md)
for the Chen and Li combination method,
[`com_conf_quant_larger_cre`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_conf_quant_larger_cre.md)
for the CMRSS combined method

## Examples

``` r
if (FALSE) { # \dontrun{
data(electric_teachers)
Z <- electric_teachers$TxAny
Y <- electric_teachers$gain

# Using Stephenson statistic with s = 6
result <- method_caughey(Z, Y,
                         k_vec = floor(c(0.7, 0.8, 0.9) * length(Z)),
                         method.list = list(name = "Stephenson", s = 6),
                         nperm = 10000,
                         alpha = 0.05)
print(result)
} # }
```
