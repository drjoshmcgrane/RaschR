# Fit a many-facet Rasch model

Fits an additive many-facet Rasch model (Linacre 1989) to scored
responses indexed by person, item, and one or more facets such as rater,
task, or occasion. Facet severities, item thresholds, person locations,
and fit statistics are reported on a common logit scale.

## Usage

``` r
rasch_mfrm(
  data,
  person,
  item = NULL,
  score = NULL,
  facets,
  items = NULL,
  n_groups = NULL,
  adjust_N = NA,
  na_codes = -1,
  interaction = NULL,
  factors = NULL,
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- data:

  Long-format data frame.

- person:

  Name of the person identifier column.

- item:

  Name of the item column.

- score:

  Name of the integer score column (categories from 0; gaps are
  collapsed per item with a note).

- facets:

  Character vector naming one or more facet columns (for example a rater
  column).

- items:

  Optional character vector of item score columns for data in wide
  format: one row per person-by-facet combination (for example one row
  per script per rater) with one column per item or criterion. The long
  form (`item` + `score`) remains available for data where the facet
  varies within items.

- n_groups:

  Number of class intervals for the item-trait chi-square; `NULL` (the
  default) applies the class-interval rule of Andrich and Marais (2019,
  ch. 15) (at least 50 non-extreme persons per interval, at most 10
  intervals, at least 2).

- adjust_N:

  Optional reference sample size for the chi-square.

- na_codes:

  Score values to read as missing (default `-1`); any negative score is
  also treated as missing.

- interaction:

  Optional name of one facet to interact with the items (interactive
  facet mode). See Details.

- factors:

  Optional person factors for DIF analysis: a character vector naming
  columns constant within person, or a data frame with one row per data
  row or unique person. Facets belong in `facets`, not here.

- maxit, tol:

  Newton-Raphson iteration cap and convergence tolerance.

## Value

An object of classes `"rasch_mfrm"` and `"rasch"`. Model-specific
components are `facet_effects`, `item_effects`, `item_thresholds`, and
`facet_spec`. Interactive fits also contain `interaction_test` and
`interaction_effects`. `fit_resid` averages virtual-item residuals
within a margin; `fit_resid_pooled` is the response-weighted pooled
statistic, with degrees of freedom in `df_fit`.

## Details

For person \\n\\, item \\i\\, and facet levels \\f_1,\ldots,f_Q\\, the
additive model is \$\$P(X\_{ni\mathbf{f}}=x)=\frac{\exp\\x\theta_n-
\sum\_{k=1}^{x}\[\delta\_{ik}+\sum\_{q=1}^{Q}\rho\_{qf_q}\]\\}
{\sum\_{y=0}^{m_i}\exp\\y\theta_n-
\sum\_{k=1}^{y}\[\delta\_{ik}+\sum\_{q=1}^{Q}\rho\_{qf_q}\]\\}.\$\$
Positive facet values therefore denote greater severity. The item
thresholds have a common sum-zero origin and the levels of each facet
sum to zero. If `interaction` is requested, an item-by-level term is
added with both its item and facet margins constrained to sum to zero.

Estimation represents each observed item-by-facet combination as a
virtual item and imposes the additive structure in the pairwise
conditional likelihood. The person parameter cancels before calibration.
The covariance of the structural parameters is the transformed Godambe
sandwich covariance.

Facet levels must be connected through common persons and items. A facet
nested within an item or a person-disjoint block can be confounded with
the item location. The function checks the structural rank and response
graph before fitting the model.

An item-by-facet interaction retains equal discrimination but allows
facet differences to vary by item. The omnibus Wald test in
`interaction_test` is the primary test; cell tests are Holm-adjusted
follow-ups.

## References

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Linacre, J. M. (1989). Many-Facet Rasch Measurement. Chicago: MESA
Press.

## See also

[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md),
[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md),
and
[`simulate_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/simulate_mfrm.md).

## Examples

``` r
set.seed(1)
simP <- function(th, tau) {
  x <- 0:length(tau)
  p <- exp(x * th - c(0, cumsum(tau)))
  p / sum(p)
}
persons <- sprintf("P%03d", 1:120); raters <- paste0("R", 1:4)
th <- setNames(rnorm(120, 0, 1.3), persons)
rho <- setNames(c(-0.6, -0.2, 0.2, 0.6), raters)
tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
d <- expand.grid(person = persons, item = names(tau), rater = raters,
                 stringsAsFactors = FALSE)
d$score <- mapply(function(p, i, r)
  sample(0:2, 1, prob = simP(th[p], tau[[i]] + rho[r])),
  d$person, d$item, d$rater)
fit <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                  facets = "rater")
fit$facet_effects$rater
#>    level   severity         se   n  infit_ms outfit_ms    fit_resid
#> R1    R1 -0.6741794 0.08866656 354 1.0748808 1.0555789  0.352198069
#> R2    R2 -0.1984195 0.07726631 354 1.0559605 1.0389336  0.337676916
#> R3    R3  0.1850839 0.07784122 354 1.0009494 0.9972743 -0.007416974
#> R4    R4  0.6875150 0.07293532 354 0.9666351 0.9568306 -0.336987303
#>    fit_resid_pooled df_fit
#> R1      0.736766065  322.5
#> R2      0.548246161  322.5
#> R3      0.008326603  322.5
#> R4     -0.508911000  322.5
```
