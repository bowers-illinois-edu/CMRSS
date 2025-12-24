# Sort treated units by outcome rank

Returns the indices of treated units sorted by their outcome values in
increasing order.

## Usage

``` r
sort_treat(Y, Z)
```

## Arguments

- Y:

  An n-dimensional observed outcome vector.

- Z:

  An n-dimensional binary treatment assignment vector.

## Value

Integer vector of indices of treated units, sorted by increasing
outcome.
