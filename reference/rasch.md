# Fit a Rasch model

Fits the partial credit model (PCM) or rating scale model (RSM) by
pairwise conditional maximum likelihood. Person locations are Warm
weighted likelihood estimates. The fitted object contains item and
person fit, targeting, reliability, threshold diagnostics, residuals,
and a score-to-measure table.

## Usage

``` r
rasch(
  data,
  model = c("PCM", "RSM"),
  id = NULL,
  factors = NULL,
  items = NULL,
  n_groups = NULL,
  adjust_N = NA,
  anchors = NULL,
  na_codes = -1,
  key = NULL,
  pc_components = NULL,
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- data:

  Persons-by-items integer score matrix (categories from 0), or a data
  frame also containing ID and person-factor columns. Missing values are
  allowed subject to the identification and ignorability conditions
  described above.

- model:

  Either `"PCM"` (partial credit) or `"RSM"` (rating scale).

- id:

  Optional name of an ID column in `data`, or a vector of IDs; carried
  through to the person estimates.

- factors:

  Optional character vector of person-factor column names in `data` (for
  DIF analysis), a data frame of factors, or one grouping vector with
  one entry per data row.

- items:

  Optional character vector naming the item columns; by default every
  column not named in `id` or `factors`.

- n_groups:

  Number of class intervals for the item-trait chi-square and ANOVA item
  fit. The default `NULL` applies the rule of Andrich and Marais (2019,
  ch. 15): as many intervals of at least 50 non-extreme persons as the
  sample allows, at most 10, at least 2. The resolved value is stored in
  `fit$n_groups`.

- adjust_N:

  Optional reference sample size used to rescale the item-trait
  chi-squares. See Details.

- anchors:

  Optional anchor table for equating: a data frame with columns `item`,
  `k`, and `tau`; see
  [`pcml`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md).
  Anchors determine the scale origin.

- na_codes:

  Values to read as missing. Defaults to `-1`, the conventional
  missing-response code; any negative score is also treated as missing,
  since valid category scores start at zero.

- key:

  Optional multiple-choice key: a named item-to-option vector, an
  item/key table, or an item/option/score table. See Details.

- pc_components:

  `NULL` (the default) estimates all PCM thresholds freely. Values from
  1 to 4 use the principal-components form in
  [`pcml_pc`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md):
  location, then spread, skewness, and kurtosis. This can stabilise
  sparse categories. Component estimates are stored in the estimation
  details. Available for PCM fits without anchors.

- maxit, tol:

  Newton-Raphson iteration cap and convergence tolerance of the pairwise
  conditional estimation.

## Value

An object of class `"rasch"`. Its principal components are the item
summary, threshold table, person table, score table, residuals,
reliability, targeting, item-trait statistics, threshold diagnostics,
and estimation details. The component `summary_stats` contains the
distribution summaries, fit-location correlations, and the cell
degrees-of-freedom factor.

## Details

For scores \\x=0,\ldots,m_i\\, the PCM is
\$\$P(X\_{ni}=x)=\frac{\exp\\x\theta_n-\sum\_{k=1}^{x}\delta\_{ik}\\}
{\sum\_{y=0}^{m_i}\exp\\y\theta_n-\sum\_{k=1}^{y}\delta\_{ik}\\}.\$\$
The RSM constrains \\\delta\_{ik}=\beta_i+\tau_k\\, where \\\beta_i\\ is
the item location and \\\tau_k\\ is common across items. Dichotomous
items are the one-threshold case of the PCM.

Pairwise conditioning removes \\\theta_n\\ from the item likelihood.
Missing responses are omitted from pairwise contributions, and person
measures are estimated within each observed item pattern. The observed
item-pair graph must identify a common scale. This covers planned linked
designs and ignorable missingness; informative missingness can still
bias the estimates.

The fit residual is the log-of-mean-square statistic described by
Andrich and Marais (2019, ch. 23). It is approximately standard normal
under fit; positive values indicate under-discrimination and negative
values indicate over-discrimination. The item-trait chi-square and
class-interval F tests are large-sample diagnostic approximations and
should be considered with the residual statistics, effect sizes, and
item content.

Multiple-choice responses may be scored from a named item-to-key vector,
an item/key table, or an item/option/score table. A slash separates
alternative correct options. The third form assigns integer category
scores to nominated options and fits the resulting item as polytomous;
unlisted options score zero. Raw responses are retained in `fit$mc` for
distractor analysis.

If `adjust_N` is supplied, each item-trait chi-square is multiplied by
the reference sample size divided by the number of classified persons.
The scaling is global: an item answered by a subset retains its
proportionally smaller share of the reference sample.

## References

Rasch, G. (1960). Probabilistic Models for Some Intelligence and
Attainment Tests. Copenhagen: Danish Institute for Educational Research.
(Expanded edition, 1980, Chicago: University of Chicago Press.)

Rasch, G. (1961). On general laws and the meaning of measurement in
psychology. In Proceedings of the Fourth Berkeley Symposium on
Mathematical Statistics and Probability (Vol. 4, pp. 321–333). Berkeley:
University of California Press.

Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
Rasch model for ordered response categories using principal components.
Journal of Applied Measurement, 4(3), 205–221.

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. Psychometrika, 54(3), 427–450.

## See also

[`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md),
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md),
[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md),
[`test_information`](https://drjoshmcgrane.github.io/rasch/reference/test_information.md),
and
[`run_app`](https://drjoshmcgrane.github.io/rasch/reference/run_app.md).

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 8)
X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(X, model = "PCM")
fit$items
#>  item max location    se fit_resid  df_fit natural_resid infit_ms outfit_ms
#>    I1   1   -2.160 0.144     0.167 423.500         0.170    1.027     0.988
#>    I2   1   -1.352 0.113    -0.576 423.500        -0.555    1.030     0.920
#>    I3   1   -0.894 0.106    -1.427 423.500        -1.331    1.005     0.876
#>    I4   1   -0.253 0.098    -0.684 423.500        -0.668    1.061     0.975
#>    I5   1    0.391 0.098     0.868 423.500         0.896    1.160     1.094
#>    I6   1    0.890 0.100    -0.563 423.500        -0.549    1.074     0.963
#>    I7   1    1.468 0.113     0.186 423.500         0.189    1.105     1.018
#>    I8   1    1.910 0.128     0.229 423.500         0.233    0.997     1.012
#>  infit_z outfit_z  chisq df     p p_adj p_bonf F_anova p_anova
#>    0.332   -0.007  5.367  6 0.498 0.569  1.000   0.550   0.770
#>    0.514   -0.694  8.361  6 0.213 0.569  1.000   1.387   0.218
#>    0.107   -1.473 13.544  6 0.035 0.281  0.281   2.771   0.012
#>    1.424   -0.373  7.579  6 0.271 0.569  1.000   1.493   0.179
#>    3.694    1.432  2.837  6 0.829 0.829  1.000   0.499   0.809
#>    1.566   -0.429  5.526  6 0.478 0.569  1.000   0.941   0.465
#>    1.767    0.196  5.761  6 0.451 0.569  1.000   0.708   0.643
#>   -0.018    0.130  6.873  6 0.333 0.569  1.000   0.956   0.455
fit$psi$PSI
#> [1] 0.4907417
```
