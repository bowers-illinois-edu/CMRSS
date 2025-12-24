# Rank score function for comparison methods

Computes rank scores for Wilcoxon or Stephenson statistics. This version
normalizes scores by the maximum value, matching the RIQITE package
convention.

## Usage

``` r
rank_score_riqite(n, method.list = list(name = "Wilcoxon"))
```

## Arguments

- n:

  Number of units.

- method.list:

  A list specifying the scoring method:

  - `name`: "Wilcoxon" or "Stephenson"

  - `s`: (for Stephenson) the parameter s

## Value

A numeric vector of length n containing the normalized rank scores.
