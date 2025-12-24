# Unified solver wrapper for stratum-level optimization

A unified interface for solving the stratum-level optimization problem
using either Gurobi or HiGHS solver.

## Usage

``` r
solve_stratum_optimization(coeflist, p, exact = TRUE, solver = "auto")
```

## Arguments

- coeflist:

  A list of B matrices containing the stratum-wise standardized test
  statistics for different coverage scenarios.

- p:

  Upper bound on the number of units with effects greater than c.

- exact:

  Logical; if TRUE, solve as integer linear program (ILP). If FALSE,
  solve as linear program (LP relaxation).

- solver:

  Character string specifying the solver to use: "auto", "highs", or
  "gurobi".

## Value

A list with:

- sol:

  Solution vector

- obj:

  Optimal objective value

- solver:

  The solver that was used

## See also

[`solve_optimization`](https://bowers-illinois-edu.github.io/CMRSS/reference/solve_optimization.md)
for the standard optimization
