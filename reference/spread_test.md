# Spread-parameter test for dependence within subtests

Andrich's (1985) least-upper-bound screen: the spread component
\\\lambda\\ of a polytomous item (half the distance between successive
thresholds in the principal-components parameterisation, estimated here
by
[`pcml_pc`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md))
cannot fall below the value implied by the binomial distribution when
the item is a subtest of equally difficult, independent dichotomous
items; different difficulties only raise it. Sampling uncertainty
matters when the fitted spread lies near the bound. The function
therefore reports whether the point estimate is below the bound
separately from a one-sided test of \\H_0: \lambda \geq \lambda_0\\
against \\H_1: \lambda \< \lambda_0\\. Evidence of dependence requires
the adjusted one-sided probability to be below `alpha`; a point estimate
below the bound alone is not treated as a verdict. Applied to the
superitems recorded by
[`combine_items`](https://drjoshmcgrane.github.io/rasch/reference/combine_items.md).
The binomial bound applies only when every component was dichotomous; a
composite containing a polytomous item is shown but its bound and
verdict are withheld. The input calibration and the principal-components
refit must both converge.

## Usage

``` r
spread_test(fit, maxit = 60, tol = 1e-08, alpha = 0.05, p_adjust = "holm")
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- maxit, tol:

  Passed to the
  [`pcml_pc`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  refit.

- alpha:

  Significance level for the one-sided dependence screen.

- p_adjust:

  Multiplicity adjustment across the eligible superitems; one of
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html). An
  eligible superitem remains in the family if its probability is
  unavailable.

## Value

A data frame with one row per recorded superitem: `item`, `m`, whether
the binomial bound is `eligible`, the `spread` estimate and its `se`,
the bound `lub` (available for dichotomous-component subtests with
maximum scores 2 to 8), `z` = (spread - lub)/se, the one-sided `p` and
adjusted `p_adj`, `below_bound` for the point-estimate comparison, and
`dependent` for adjusted evidence at `alpha`. Items not formed by
[`combine_items()`](https://drjoshmcgrane.github.io/rasch/reference/combine_items.md)
are omitted. The result retains `alpha` and `p_adjust` as attributes.

## References

Andrich, D. (1985). An elaboration of Guttman scaling with Rasch models
for measurement. In N. B. Tuma (Ed.), Sociological Methodology 1985 (pp.
33–80). Jossey-Bass.

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

## Examples

``` r
set.seed(1); N <- 600
d0 <- seq(-1.5, 1.5, length.out = 8)
X <- matrix(rbinom(N * 8, 1, plogis(outer(rnorm(N), d0, "-"))), N, 8)
X[, 5] <- ifelse(runif(N) < .85, X[, 4], 1 - X[, 4])
X[, 6] <- ifelse(runif(N) < .75, X[, 4], 1 - X[, 4]) # a dependent triple
colnames(X) <- paste0("I", 1:8)
fit2 <- combine_items(rasch(X), list(c("I4", "I5", "I6"), c("I1", "I2", "I3")))
spread_test(fit2)
#> Spread-parameter screen (Andrich 1985): one-sided evidence below the binomial bound (holm adjustment; alpha 0.050)
#>      item m eligible spread    se bound       z       p   p_adj below_bound
#>  I4+I5+I6 3        * -0.118 0.051 0.550 -13.110 < 0.001 < 0.001           *
#>  I1+I2+I3 3        *  0.501 0.081 0.550  -0.610   0.271   0.271           *
#>  dependent
#>          *
#>           
```
