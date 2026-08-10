# DIF magnitude in logits with pairwise comparisons

Quantifies differential item functioning on the measurement scale
itself, where practical significance is judged: the item is resolved
into one copy per group (or per cell of a factor combination), the model
is refitted, and the distance between the resolved locations is the DIF
size in logits (Andrich & Marais 2019, ch. 16: a simulated shift of 0.71
was recovered as 0.75 by exactly this method). Every pair of levels is
compared with a Wald test using the full sandwich covariance of the
resolved locations (for a between-person factor the persons behind
different levels are disjoint, but the shared calibration of the other
items still couples the estimates, so the covariance is used rather than
assumed zero), with familywise adjustment over the pairs. When person
identifiers repeat, as in a stacked repeated-measures design, the
calibration sandwich treats the rows as independent and does not carry
the covariance induced by repeated persons. The resolved point
differences and practical flags are therefore retained, but sampling
standard errors, confidence intervals, and Wald tests are withheld;
[`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md)
supplies person-level inferential tests. Differences at least
`flag_logits` in absolute size are flagged as practically significant;
half a logit is a common working criterion, to be weighed against the
test's targeting and purpose.

## Usage

``` r
dif_size(
  fit,
  item,
  by,
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

- item:

  Item name or index.

- by:

  One or more person-factor names nominated in the fit (several names
  give interaction cells), or a grouping vector/data frame with one
  entry per person.

- p_adjust:

  Familywise adjustment over the pairwise comparisons; default `"holm"`.

- alpha:

  Significance level for the adjusted probabilities.

- flag_logits:

  Absolute difference flagged as practically significant.

- min_n:

  Levels with fewer responders to the item are dropped (their resolved
  locations would be too unstable to compare), with a note.

## Value

A list of class `"rasch_dif_size"`: `levels` (resolved location and SE
per level, with its n), `pairs` (per comparison: difference in logits,
SE, z, raw and adjusted p, 95 per cent interval, `significant`,
`practical`), the settings, and any notes. Sampling-uncertainty fields
are `NA` when person IDs repeat.

## Details

For an interaction, supply several factor names: levels are then the
factor-combination cells, which is the post-hoc follow-up to a
significant factor-by-factor term in
[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md).

## References

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. Scandinavian Journal of Statistics, 6(2), 65–70.

## Examples

``` r
set.seed(1); n <- 600
d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(data.frame(X, grp = g), factors = "grp")
dif_size(fit, "I3", by = "grp")
#> DIF size for I3 by grp (resolved locations, logits)
#>  level location    se weak   n
#>      a   -0.890 0.133    0 300
#>      b    0.018 0.126    0 300
#>  level_a level_b difference    se      z p p_adj  lower  upper significant
#>        a       b     -0.907 0.204 -4.441 0     0 -1.308 -0.507           *
#>  practical
#>    >= 0.50
#> p adjusted by holm over 1 pairwise comparison(s); practical criterion 0.50 logits
```
