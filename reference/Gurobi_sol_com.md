# Optimization function using Gurobi solver (original implementation)

Solves the integer/linear programming problem for computing the minimum
test statistic in stratified randomized experiments using the Gurobi
solver.

## Usage

``` r
Gurobi_sol_com(
  Z,
  block,
  weight,
  coeflists,
  p,
  ms_list,
  exact = TRUE,
  block.sum = NULL
)
```

## Arguments

- Z:

  An n-dimensional treatment assignment vector.

- block:

  An n-dimensional vector specifying block of each unit.

- weight:

  A B-dimensional vector of block weights.

- coeflists:

  A list of H elements, each containing B matrices with block-specific
  test statistics for different values of k.

- p:

  Upper bound on the number of units with effects greater than c.

- ms_list:

  A list containing mu (means) and sigma (standard deviations) for each
  test statistic.

- exact:

  Logical; if TRUE, solve as integer linear program (ILP). If FALSE,
  solve as linear program (LP relaxation).

- block.sum:

  Optional pre-computed block summary from summary_block().

## Value

A list with:

- sol:

  Solution vector

- obj:

  Optimal objective value

## Details

This is the original Gurobi-based implementation. For users without a
Gurobi license,
[`HiGHS_sol_com`](https://bowers-illinois-edu.github.io/CMRSS/reference/HiGHS_sol_com.md)
provides equivalent functionality using the open-source HiGHS solver.

## See also

[`HiGHS_sol_com`](https://bowers-illinois-edu.github.io/CMRSS/reference/HiGHS_sol_com.md)
for the HiGHS implementation,
[`solve_optimization`](https://bowers-illinois-edu.github.io/CMRSS/reference/solve_optimization.md)
for the unified wrapper

## Examples

``` r
if (FALSE) { # \dontrun{
# This is typically called internally by pval_comb_block()
result <- Gurobi_sol_com(Z, block, weight, coeflists, p, ms_list, exact = TRUE)
} # }
```
