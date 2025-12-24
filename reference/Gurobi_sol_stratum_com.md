# Stratum-level optimization using Gurobi solver

Solves the optimization problem for the stratum-level combination method
(comb.method = 2) using Gurobi.

## Usage

``` r
Gurobi_sol_stratum_com(coeflist, p, exact = TRUE)
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

## Value

A list with:

- sol:

  Solution vector

- obj:

  Optimal objective value

## See also

[`HiGHS_sol_stratum_com`](https://bowers-illinois-edu.github.io/CMRSS/reference/HiGHS_sol_stratum_com.md)
for the HiGHS implementation
