# Get default solver

Returns the default solver based on availability. Prefers Gurobi if
available (to maintain backward compatibility with earlier
versions/workflows), and falls back to HiGHS.

## Usage

``` r
get_default_solver()
```

## Value

Character string indicating the available solver

## Examples

``` r
get_default_solver()
#> [1] "highs"
```
