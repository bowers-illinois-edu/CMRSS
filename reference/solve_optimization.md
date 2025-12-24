# Unified solver wrapper

A unified interface for solving the optimization problem using either
Gurobi or HiGHS solver.

## Usage

``` r
solve_optimization(
  Z,
  block,
  weight,
  coeflists,
  p,
  ms_list,
  exact = TRUE,
  block.sum = NULL,
  solver = "auto"
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

- solver:

  Character string specifying the solver to use. Options are:

  - "auto": Automatically select available solver (default)

  - "highs": Use HiGHS solver (open-source)

  - "gurobi": Use Gurobi solver (commercial, requires license)

## Value

A list with:

- sol:

  Solution vector

- obj:

  Optimal objective value

- solver:

  The solver that was used

## Details

This function provides a unified interface for both solvers. When solver
= "auto", it will prefer Gurobi if available (for backward
compatibility), otherwise fall back to HiGHS.

Both solvers produce equivalent results for the same optimization
problem, though there may be minor numerical differences (typically \<
1e-6) due to different internal algorithms.

## See also

[`HiGHS_sol_com`](https://bowers-illinois-edu.github.io/CMRSS/reference/HiGHS_sol_com.md),
[`Gurobi_sol_com`](https://bowers-illinois-edu.github.io/CMRSS/reference/Gurobi_sol_com.md),
[`solver_available`](https://bowers-illinois-edu.github.io/CMRSS/reference/solver_available.md),
[`get_default_solver`](https://bowers-illinois-edu.github.io/CMRSS/reference/get_default_solver.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Automatically select solver
result <- solve_optimization(Z, block, weight, coeflists, p, ms_list)

# Explicitly use HiGHS
result <- solve_optimization(Z, block, weight, coeflists, p, ms_list, solver = "highs")

# Explicitly use Gurobi
result <- solve_optimization(Z, block, weight, coeflists, p, ms_list, solver = "gurobi")
} # }
```
