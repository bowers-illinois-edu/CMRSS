# Electric Circuits Professional Development Study Data

Data from a randomized experiment studying the effects of professional
development on elementary teachers' science teaching, specifically
related to teaching about electric circuits.

## Usage

``` r
electric_teachers
```

## Format

A data frame with 233 rows and 7 variables:

- Site:

  Character. Study site identifier (S2-S8). The experiment was conducted
  across 7 sites.

- T.ID:

  Numeric. Unique teacher identifier.

- Tx:

  Character. Treatment group assignment with four levels:

  - A: Treatment condition A

  - B: Treatment condition B

  - C: Treatment condition C

  - D: Control condition

- TxAny:

  Numeric. Binary indicator for any treatment (1) vs control (0).
  Teachers in groups A, B, or C have TxAny=1; teachers in group D have
  TxAny=0.

- Per.1:

  Numeric. Pre-test score (percentage correct, 0-100).

- Per.2:

  Numeric. Post-test score (percentage correct, 0-100).

- gain:

  Numeric. Gain score calculated as Per.2 - Per.1.

## Source

Heller, J. L., Shinohara, M., Miratrix, L., Hesketh, S. R., and Daehler,
K. R. (2010). Learning science for teaching: Effects of professional
development on elementary teachers, classrooms, and students.
*Proceedings from Society for Research on Educational Effectiveness*.

## Details

This dataset comes from a study of professional development programs for
elementary science teachers. Teachers were randomly assigned to
different professional development conditions within each site
(stratified randomized experiment). The outcome measures assess
teachers' content knowledge about electric circuits before and after the
intervention.

The dataset can be used to demonstrate methods for inference about
quantiles of individual treatment effects in stratified randomized
experiments using the CMRSS package.

## Examples

``` r
data(electric_teachers)

# Basic exploration
head(electric_teachers)
#> # A tibble: 6 × 7
#>   Site   T.ID Tx    TxAny Per.1 Per.2  gain
#>   <chr> <dbl> <chr> <dbl> <dbl> <dbl> <dbl>
#> 1 S2     2083 C         1  53.3  66.7  13.3
#> 2 S2     2085 C         1  36.7  80    43.3
#> 3 S2     2086 C         1  30    83.3  53.3
#> 4 S2     2088 C         1  46.7  70    23.3
#> 5 S2     2090 C         1  40    60    20  
#> 6 S2     2093 C         1  46.7  73.3  26.7
table(electric_teachers$Site, electric_teachers$Tx)
#>     
#>       A  B  C  D
#>   S2  6  7 12 11
#>   S3 23 11  0  8
#>   S4  5  7  0  6
#>   S5  0  0  9  6
#>   S6 15  9  9 10
#>   S7  0  9  6 10
#>   S8  9 19  8 18

# Create treatment indicator for any treatment vs control
Z <- electric_teachers$TxAny
Y <- electric_teachers$gain
block <- factor(electric_teachers$Site)

# Summary by treatment status
tapply(Y, Z, mean)
#>         0         1 
#>  2.463333 22.235366 

# Number of treated and control units per site
table(block, Z)
#>      Z
#> block  0  1
#>    S2 11 25
#>    S3  8 34
#>    S4  6 12
#>    S5  6  9
#>    S6 10 33
#>    S7 10 15
#>    S8 18 36

if (FALSE) { # \dontrun{
##### Example 1: CRE Analysis (ignoring site stratification) #####

# Define Stephenson statistics with different s values
s.vec <- c(2, 6, 10, 30)
methods.list.cre <- lapply(s.vec, function(s) {
  list(name = "Stephenson", s = s, std = TRUE, scale = TRUE)
})

# Confidence intervals for effect quantiles among treated
ci.treat <- com_conf_quant_larger_cre(Z, Y,
                                      methods.list = methods.list.cre,
                                      nperm = 10000,
                                      set = "treat",
                                      alpha = 0.05)

##### Example 2: SRE Analysis (accounting for site stratification) #####

s <- length(levels(block))  # 7 sites

# Define polynomial rank statistics for each stratum
r.vec <- c(2, 6, 10)
methods.list.sre <- list()
for (j in seq_along(r.vec)) {
  methods.list.sre[[j]] <- lapply(1:s, function(i) {
    list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
  })
}

# Confidence intervals accounting for stratification
ci.treat.sre <- com_block_conf_quant_larger(Z, Y, block,
                                            set = "treat",
                                            methods.list.all = methods.list.sre,
                                            weight.name = "asymp.opt",
                                            opt.method = "ILP",
                                            alpha = 0.05)
} # }
```
