# Fit and diagnose a Rasch model by pairwise conditional estimation

Runs a complete Rasch analysis: Andrich and Luo pairwise conditional
maximum likelihood item estimation (see
[`pcml`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)), Warm
weighted likelihood person estimates per missing-data pattern, item and
person fit residuals (the log-of-mean-square statistic of Andrich and
Marais 2019, ch. 23, with its untransformed natural form and degrees of
freedom), infit and outfit, the item-trait interaction chi-square and
the class-interval ANOVA item-fit F, the person separation index with
and without extremes, Cronbach's alpha, targeting, threshold
diagnostics, and the score-to-measure table.

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

  Optional reference sample size; if supplied, item-trait chi-squares
  are rescaled to this size (a sample-size adjustment for the
  sensitivity of the chi-square to large samples). The scaling is
  proportional and global – every item's chi-square is multiplied by
  `adjust_N` over the number of classified persons – so an item answered
  by a subset of persons keeps its proportionally smaller share of the
  notional sample rather than being inflated to the full `adjust_N`.

- anchors:

  Optional anchor table for equating: a data frame with columns `item`,
  `k`, and `tau` fixing nominated thresholds at known values; see
  [`pcml`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md).
  With anchors in place the scale origin comes from the anchors, so
  person measures are directly comparable across separately analysed
  datasets.

- na_codes:

  Values to read as missing. Defaults to `-1`, the conventional
  missing-response code; any negative score is also treated as missing,
  since valid category scores start at zero.

- key:

  Optional multiple-choice scoring key, in any of three forms. (1) A
  named vector or data frame with columns `item` and `key` naming each
  item's correct option: scored 0/1 (case-insensitive after trimming;
  blanks become missing). (2) Double keying: several correct options
  separated by `"/"` (for example `"A/C"`), all scoring 1. (3)
  Polytomous option scoring (Andrich and Styles 2011): a data frame with
  columns `item`, `option`, and `score` assigning an integer score to
  every credited option (unlisted options score 0), so informative
  distractors receive partial credit and the item is fitted as
  polytomous; see
  [`distractor_rescore`](https://drjoshmcgrane.github.io/rasch/reference/distractor_rescore.md)
  for an evidence-based proposal. Raw responses are retained in `fit$mc`
  for
  [`distractor_analysis`](https://drjoshmcgrane.github.io/rasch/reference/distractor_analysis.md)
  and
  [`plot_distractors`](https://drjoshmcgrane.github.io/rasch/reference/plot_distractors.md).

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

The fit residual follows Andrich and Marais (2019, ch. 23) exactly:
standardised residuals are squared and summed over each item's
(person's) observed cells among non-extreme persons, compared with the
summed cell degrees of freedom (the model-testing degrees of freedom,
cells minus estimated parameters, apportioned equally over cells), and
symmetrised by the log-of-mean-square transform \\f (\ln Y^2 - \ln
f)/\sqrt{V\[Y^2\]}\\ with model-based variance \\V\[Y^2\] = \sum
(C_4/V^2 - 1)\\. Values are approximately N(0,1) under fit; the
conventional flagging value is 2.5 (Andrich and Marais 2019, ch. 15).
Negative values indicate over-discrimination (Guttman-like responses),
positive values under-discrimination.

A calibration note that applies to the whole test-of-fit suite: the
item-trait chi-square evaluates unconditional residuals at estimated
person measures and refers the sum of correlated item statistics to a
chi-square with summed degrees of freedom, following the convention of
Andrich and Marais (2019); the class-interval ANOVA F and the
Wilson-Hilferty-style transformations are approximations of the same
kind. Null simulation with this package shows rejection rates near but
not exactly at the nominal level (conservative in the settings
examined), and the direction can vary with design. These statistics are
therefore approximate, convention-faithful diagnostics for ordering and
flagging misfit – not exactly calibrated hypothesis tests; where exact
calibration matters, use the simulation tools
([`sim_replicate`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md))
to build parametric-bootstrap reference distributions for the observed
design.

Missing responses are omitted from the pairwise contributions and person
estimates are computed for each observed item pattern. This supports
planned linked designs and ignorable response missingness when the
co-observation graph identifies one scale. It does not make informative
missingness harmless: if response availability depends on an unmodelled
response or person process, item and person estimates can be biased. The
missingness mechanism and design connectedness therefore remain part of
the substantive analysis.

## References

Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
Rasch model for ordered response categories using principal components.
Journal of Applied Measurement, 4(3), 205–221.

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. Psychometrika, 54(3), 427–450.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 8)
X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(X, model = "PCM")
fit$items
#>   item max   location         se  fit_resid df_fit natural_resid  infit_ms
#> 1   I1   1 -2.1604709 0.14430225  0.1666750  423.5     0.1695416 1.0268996
#> 2   I2   1 -1.3518390 0.11268892 -0.5756998  423.5    -0.5548927 1.0304617
#> 3   I3   1 -0.8944690 0.10569675 -1.4271840  423.5    -1.3314617 1.0045536
#> 4   I4   1 -0.2525755 0.09842604 -0.6843318  423.5    -0.6675975 1.0605330
#> 5   I5   1  0.3908159 0.09808533  0.8675845  423.5     0.8959149 1.1604892
#> 6   I6   1  0.8900794 0.10034984 -0.5631943  423.5    -0.5486780 1.0739741
#> 7   I7   1  1.4679700 0.11257588  0.1864529  423.5     0.1887381 1.1051993
#> 8   I8   1  1.9104891 0.12794128  0.2287863  423.5     0.2332503 0.9970331
#>   outfit_ms     infit_z     outfit_z     chisq df          p     p_adj
#> 1 0.9882788  0.33225795 -0.006659097  5.367357  6 0.49763112 0.5687213
#> 2 0.9200216  0.51421357 -0.694204299  8.360771  6 0.21284560 0.5687213
#> 3 0.8764090  0.10655603 -1.472784349 13.543617  6 0.03517114 0.2813691
#> 4 0.9752165  1.42357529 -0.372765909  7.579412  6 0.27056336 0.5687213
#> 5 1.0935539  3.69355737  1.432422372  2.836605  6 0.82905758 0.8290576
#> 6 0.9633468  1.56569054 -0.429193066  5.525811  6 0.47834263 0.5687213
#> 7 1.0181545  1.76661209  0.196435036  5.760869  6 0.45050380 0.5687213
#> 8 1.0119034 -0.01774722  0.129649879  6.873008  6 0.33275124 0.5687213
#>      p_bonf   F_anova    p_anova
#> 1 1.0000000 0.5500185 0.77003278
#> 2 1.0000000 1.3872863 0.21785564
#> 3 0.2813691 2.7706872 0.01171615
#> 4 1.0000000 1.4928791 0.17862223
#> 5 1.0000000 0.4990124 0.80919806
#> 6 1.0000000 0.9406784 0.46536016
#> 7 1.0000000 0.7079669 0.64333086
#> 8 1.0000000 0.9560240 0.45459628
fit$psi$PSI
#> [1] 0.4907417
```
