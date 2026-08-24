# Residual correlations for local dependence (Yen's Q3)

The pairwise correlations of the standardised response residuals are
Yen's (1984) Q3 statistics. Under unidimensionality and local
independence the off-diagonal values sit near `-1/(L-1)`; large positive
values flag local dependence between item pairs. Following Christensen,
Makransky and Horton (2017), each Q3 is also reported relative to the
average off-diagonal value (`q3_star`). There is no universal
adjusted-Q3 critical value: it depends on sample size, test length,
category structure, and missingness. The default therefore reports the
statistics without a binary flag. A user-supplied `flag` is an
explicitly heuristic screening threshold, not a calibrated significance
test.

## Usage

``` r
residual_correlations(fit, flag = NULL)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- flag:

  Optional heuristic excess above the average off-diagonal Q3 at which a
  pair is flagged. The default `NULL` withholds binary flags.

## Value

A list with the Q3 `matrix`, the adjusted-Q3 `star_matrix` (each Q3 less
the average off-diagonal value, diagonal empty), the `average`
off-diagonal value, `pairs` (every item pair with `q3`, `q3_star` and a
`flagged` indicator, sorted by `q3`), and the subset of `flagged` pairs.

## References

Yen, W. M. (1984). Effects of local item dependence on the fit and
equating performance of the three-parameter logistic model. *Applied
Psychological Measurement*, 8(2), 125-145.

Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
for Yen's Q3: identification of local dependence in the Rasch model
using residual correlations. *Applied Psychological Measurement*, 41(3),
178-194.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 8)
X <- matrix(rbinom(300 * 8, 1, plogis(outer(rnorm(300), d, "-"))), 300, 8)
colnames(X) <- paste0("I", 1:8)
residual_correlations(rasch(X))$average
#> [1] -0.1297187
```
