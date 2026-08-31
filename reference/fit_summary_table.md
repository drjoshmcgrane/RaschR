# Test-of-fit summary as a table

Returns the model, estimation method, trait chi-square, calibration and
person fit-residual moments, fit-location correlations, the approximate
asymptotic chi-square flag count, and disordered-threshold count as a
two-column table. MFRM and EFRM summaries label their fitted
item-by-facet or item-by-frame columns as response cells. A wholly
unavailable probability family is labelled as unavailable; a partial
family reports the tested and unavailable counts.

## Usage

``` r
fit_summary_table(fit)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  or a paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) (which
  reports its own headline set: convergence, pairwise chi-square, object
  separation, thresholds structure, and dependence effects when
  estimated).

## Value

A data frame with columns `statistic` and `value`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
fit_summary_table(rasch(X))
#>                                             statistic                   value
#> 1                                               Model                     PCM
#> 2                                          Estimation pairwise conditional ML
#> 3                                           Converged                     yes
#> 4                                          Iterations                       3
#> 5  Approximate asymptotic total item-trait chi-square                  25.067
#> 6                                  Degrees of freedom                      24
#> 7       Approximate asymptotic item-trait probability                   0.402
#> 8                                     Class intervals                       5
#> 9                              Item fit residual mean                   -0.21
#> 10                               Item fit residual SD                    0.61
#> 11                         Item fit residual skewness                    0.64
#> 12                         Item fit residual kurtosis                   -0.78
#> 13                           Person fit residual mean                   -0.24
#> 14                             Person fit residual SD                    0.66
#> 15                       Person fit residual skewness                    0.78
#> 16                       Person fit residual kurtosis                    0.02
#> 17                   Fit-location correlation (items)                  -0.182
#> 18                 Fit-location correlation (persons)                   0.004
#> 19     Items with approximate asymptotic Holm p < .05                  0 of 6
#> 20                              Disordered thresholds                    none
```
