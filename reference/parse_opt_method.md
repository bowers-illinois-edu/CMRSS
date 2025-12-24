# Parse optimization method string

Parses the opt.method parameter to determine solver and exactness.

## Usage

``` r
parse_opt_method(opt.method)
```

## Arguments

- opt.method:

  Character string specifying optimization method:

  - "ILP_gurobi": Integer linear programming with Gurobi

  - "LP_gurobi": Linear programming relaxation with Gurobi

  - "ILP_highs": Integer linear programming with HiGHS

  - "LP_highs": Linear programming relaxation with HiGHS

  - "ILP_auto" or "ILP": Integer linear programming with auto-selected
    solver

  - "LP_auto" or "LP": Linear programming with auto-selected solver

## Value

A list with:

- exact:

  Logical indicating if ILP (TRUE) or LP (FALSE)

- solver:

  Character string: "gurobi", "highs", or "auto"

## Examples

``` r
parse_opt_method("ILP_gurobi")  # list(exact = TRUE, solver = "gurobi")
#> $exact
#> [1] TRUE
#> 
#> $solver
#> [1] "gurobi"
#> 
parse_opt_method("LP_highs")    # list(exact = FALSE, solver = "highs")
#> $exact
#> [1] FALSE
#> 
#> $solver
#> [1] "highs"
#> 
parse_opt_method("ILP")         # list(exact = TRUE, solver = "auto")
#> $exact
#> [1] TRUE
#> 
#> $solver
#> [1] "auto"
#> 
```
