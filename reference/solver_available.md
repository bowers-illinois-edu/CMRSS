# Solver Utilities for CMRSS Package

The package supports two optimization solvers:

- **Gurobi**: Commercial solver (requires license). Fast and robust.

- **HiGHS**: Open-source solver. No license required, good performance.

Both solvers produce equivalent results for the same optimization
problem. HiGHS is recommended for users without a Gurobi license. Check
if a solver is available

## Usage

``` r
solver_available(solver)
```

## Arguments

- solver:

  Character string: "gurobi" or "highs"

## Value

Logical indicating if the solver is available

## Details

This file contains optimization solver functions for solving the
integer/linear programming problems in stratified randomized
experiments. It supports both Gurobi and HiGHS solvers.

## Examples

``` r
solver_available("highs")
solver_available("gurobi")
```
