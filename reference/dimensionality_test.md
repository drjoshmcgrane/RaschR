# Residual-component test of unidimensionality

Estimates each person separately on two item subsets and compares the
two estimates with a per-person t-test (Smith 2002). By default the
subsets are defined by the sign of a residual-component loading (the
first by default; any leading component may be chosen); they can also be
nominated manually (for example, by content). Under unidimensionality
and local independence the two subset estimates are independent given
the person location, so
`t = (theta_A - theta_B) / sqrt(se_A^2 + se_B^2)` is approximately
standard normal and about `alpha` of the tests should reach
significance. Persons with an extreme score on either subset are
excluded (their weighted-likelihood estimates are most biased there).
The proportion of significant tests is reported with an exact
(Clopper-Pearson) binomial confidence interval; a lower bound above
`alpha` signals multidimensionality.

## Usage

``` r
dimensionality_test(
  fit,
  alpha = 0.05,
  items_positive = NULL,
  items_negative = NULL,
  component = 1,
  min_score_points = 15L
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- alpha:

  Nominal significance level for the per-person t-tests.

- items_positive, items_negative:

  Optional character vectors naming the two item subsets; both must be
  given (disjoint, at least two items each), otherwise the sign of a
  residual component defines the split.

- component:

  Which residual principal component's loading sign defines the default
  split (ignored when subsets are named). Default the first component.

- min_score_points:

  Score-point threshold below which the verdict carries a caution.
  Andrich and Marais (2019) recommend subtests of roughly 15 score
  points for stable subtest estimates; shorter subsets (the norm for
  ordinary dichotomous tests) still receive a verdict, with a `caution`
  field noting the reduced stability. A quiet verdict under caution is
  inconclusive, not clean: with a four-item subtest the test lacks power
  where nonparametric alternatives still flag. The procedure is
  deliberately conservative – in cross-package comparison it held an
  exact null (no false flags) while flagging a balanced planted second
  dimension in about two-thirds of replicates where the DETECT index
  flagged all; quasi-exact matrix-sampling tests showed elevated null
  rates on the same data.

## Value

A list with the proportion of significant tests, its exact confidence
interval, the sample sizes (`n` used, `n_excluded_extreme`), the item
split and its source, a `multidimensional` verdict, a `caution` note
when the subtests fall short of `min_score_points`, and `paired_t`, the
paired t-test of the two subset means (the group-level comparison, which
requires pairing because both estimates come from the same persons).
When the comparison itself is unavailable (undefined split, degenerate
subsets, too few persons) the list carries a `note` explaining why and
`multidimensional = NA`.

## References

Smith, E. V. Jr. (2002). Detecting and evaluating the impact of
multidimensionality using item fit statistics and principal component
analysis of residuals. Journal of Applied Measurement, 3(2), 205–231.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 8)
X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
colnames(X) <- paste0("I", 1:8)
dimensionality_test(rasch(X))$multidimensional
#> [1] FALSE
```
