# Stratum-specific weighted maximum rank sum statistic

Calculates the weighted maximum rank sum statistic for a single stratum
across multiple methods.

## Usage

``` r
single_weight_max_rank_sum_stat_stratum(
  Z,
  Y,
  methods.list,
  score.list = NULL,
  mu_vec,
  sig_vec,
  weight.b
)
```

## Arguments

- Z:

  Treatment assignment vector for the stratum.

- Y:

  Outcome vector for the stratum.

- methods.list:

  A list of H method specifications.

- score.list:

  Optional pre-computed scores.

- mu_vec:

  Vector of H means.

- sig_vec:

  Vector of H standard deviations.

- weight.b:

  Weight for this stratum.

## Value

Weighted maximum test statistic value.
