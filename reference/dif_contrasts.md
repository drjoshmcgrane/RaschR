# Planned DIF contrasts

Tests a specified family of one-degree-of-freedom DIF contrasts. By
default, contrasts are derived from the factor structure: differences
for two-level factors, polynomial trends for ordered factors, and
pairwise or level-against-rest comparisons for nominal factors. Leading
contrasts are crossed to form two-factor interactions. User-supplied
cell weights are also accepted.

## Usage

``` r
dif_contrasts(
  fit,
  factors = NULL,
  items = NULL,
  within = NULL,
  id = NULL,
  contrasts = "auto",
  p_adjust = "holm",
  alpha = 0.05,
  flag_logits = 0.5,
  min_n = 20
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- factors:

  A data frame of person factors, a character vector naming factors
  nominated in the fit, or a single grouping vector. Defaults to every
  factor stored in the fit.

- items:

  Item names or indices to test; all items by default.

- within:

  Names of factors that vary within person (for example time). Detected
  automatically when `id` is supplied and a factor varies within an id.

- id:

  Person identifier with one entry per row, or the name of a nominated
  factor holding it; required for stacked designs where the same person
  occupies several rows.

- contrasts:

  `"auto"` (derive the family from the factor structure) or a named list
  of numeric cell-weight vectors, each named by the design-cell labels
  (factor levels joined by `":"`). Weights are rescaled so the positive
  and negative parts each sum to one.

- p_adjust:

  Familywise adjustment across items and contrasts. The default is
  `"holm"`.

- alpha:

  Significance level for the adjusted probabilities.

- flag_logits:

  Absolute estimate flagged as practically significant.

- min_n:

  Cells with fewer responders to an item are dropped from that item's
  resolution, with a note.

## Value

A list of class `"rasch_dif_contrasts"`: `table` (one row per item and
contrast: estimate in logits, SE, statistic, df where a t test was used,
raw and adjusted p, 95 per cent interval, `significant`, `practical`,
`within`), `family` (the derived questions with their cell weights), the
settings, and any `notes`.

## Details

Each logit contrast is calculated from resolved item locations. Weights
are scaled so their positive and negative parts each sum to one. With
repeated persons, inference uses person-level residual contrast scores:
within-person contrasts are tested against zero, between-person
contrasts use person means, and mixed interactions compare within-person
contrast scores between groups. The resolved logit estimate is retained,
but its calibration-based standard error is withheld because it does not
include repeated-person dependence.

For independent rows, a contrast with weights \\\mathbf{c}\\ is
\$\$\Delta_i=\mathbf{c}^{\mathsf T}\delta_i,\qquad
\operatorname{SE}(\Delta_i)= \sqrt{\mathbf{c}^{\mathsf
T}\mathbf{V}\_i\mathbf{c}}.\$\$ In a repeated-measures design, a
within-person contrast is formed from the standardised residuals,
\$\$s_p=\sum_l c_l z\_{pl},\$\$ and tested over persons. Between-person
contrasts use Welch tests of person means; mixed interactions use Welch
tests of the within-person contrast scores. The sign of each residual
test is aligned with the resolved logit contrast.

## References

Maxwell, S. E. and Delaney, H. D. (2004). Designing Experiments and
Analyzing Data (2nd ed.). Erlbaum.

Andrich, D. and Hagquist, C. (2015). Real and artificial differential
item functioning in polytomous items. Educational and Psychological
Measurement, 75(2), 185–207.

Hagquist, C. and Andrich, D. (2017). Recent advances in analysis of
differential item functioning in health research using the Rasch model.
Health and Quality of Life Outcomes, 15, 181.

## See also

[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
and
[`dif_size`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md).

## Examples

``` r
set.seed(1); n <- 600
d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(data.frame(X, grp = g), factors = "grp")
dif_contrasts(fit, items = c("I3", "I5"))
#> Planned DIF contrasts (1 questions x 2 items; holm over the family)
#>   grp: b - a
#> 
#>  item   contrast estimate    se statistic   p_adj significant practical
#>    I3 grp: b - a    0.907 0.204     4.441 < 0.001           *         *
#>    I5 grp: b - a   -0.570 0.200    -2.846   0.004           *         *
```
