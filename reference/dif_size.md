# DIF differences between factor levels

Resolves an item into one parameter per factor level, refits the model,
and reports pairwise differences between the resolved locations. Wald
tests use the full sandwich covariance and are adjusted over the family
of level comparisons. Several names in `by` define factor-combination
cells for following up an interaction in
[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md).

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
`practical`, `ets`), the settings, and any notes. Sampling-uncertainty
fields are `NA` when person IDs repeat.

## Details

Let \\\delta_i\\ contain the resolved locations of item \\i\\, and let
\\\mathbf{c}\_{ab}\\ place 1 on level \\a\\, -1 on level \\b\\, and zero
elsewhere. The reported difference and its standard error are
\$\$\Delta\_{i,ab}=\mathbf{c}\_{ab}^{\mathsf T} \delta_i,\$\$
\$\$\operatorname{SE}(\Delta\_{i,ab})= \sqrt{\mathbf{c}\_{ab}^{\mathsf
T}\mathbf{V}\_i \mathbf{c}\_{ab}},\$\$ where \\\mathbf{V}\_i\\ is the
full covariance of the resolved locations. Wald probabilities are
adjusted over the pairwise family.

With repeated person identifiers, the row-level calibration covariance
does not represent within-person sampling dependence. Logit differences
and practical flags are retained, but their standard errors and Wald
tests are withheld. Use
[`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md)
for person-level inference in a repeated-measures design.

## The ETS categories

`pairs$ets` reports the A, B and C categories used at ETS, signed for
direction as they are there. They are defined from the Mantel-Haenszel
common odds ratio on the delta scale, where \\\mathrm{MH\\ D\text{-}DIF}
= -2.35\log\hat\alpha\_{MH}\\. Under the Rasch model that log odds ratio
and the difference in item location estimate the same quantity, both
being conditional on the total score, so the delta thresholds convert
exactly at 2.35 delta units to the logit: A is not significant or below
0.426 logits, C is at or above 0.638 logits and significantly beyond
0.426, and B is the remainder. The C rule tests against a non-zero null,
so it uses the standard error rather than the probability alone.

Polytomous items are classified on the same metric. Under the partial
credit model an item's difficulty decomposes as a location plus its
thresholds, so where differential functioning shifts the location and
leaves the thresholds alone – uniform functioning, which is what a
resolved location difference estimates – the signed area between the two
groups' expected score curves is the number of thresholds times that
shift, and the same cut-values apply per threshold (Golia, 2012, section
2.2, following Cohen, Kim and Baker, 1993). Where the thresholds
themselves differ across groups the shift is not a summary of the item
and the letter should not be read.

Read the letter beside the magnitude rather than instead of it, because
the categories were built to triage items for an operational bank and
not to answer whether an item is invariant. A tops out at 0.426 logits,
which is a difference of 10.6 percentage points in success at the item's
own location – plainly not invariance. What justifies calling it
negligible is its effect on the score rather than on the item: in
simulation at 3,000 persons, one item sitting at that ceiling moved the
comparison between the two groups by 0.026 logits on a ten-item test and
0.018 on a twenty-item one, while three such items moved it by 0.086 and
0.050. So the letter answers "does this item distort the total score",
and `practical` against `flag_logits` answers "is this item behaving the
same way in both groups". They are different questions.

A separate ETS convention categorises polytomous items from a
standardised mean difference at 0.17 and 0.25, but that is a statistic
in the observed-score metric rather than this one, and Zwick, Thayer and
Mazzeo (1997) record that ETS had no official polytomous policy.

## References

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. Scandinavian Journal of Statistics, 6(2), 65–70.

Zieky, M. (1993). Practical questions in the use of DIF statistics in
item development. In P. W. Holland and H. Wainer (eds), Differential
Item Functioning (pp. 337–364). Erlbaum.

Linacre, J. M. and Wright, B. D. (1989). Mantel-Haenszel DIF and PROX
are equivalent! Rasch Measurement Transactions, 3(2), 51–53.

Golia, S. (2012). Differential item functioning classification for
polytomously scored items. Electronic Journal of Applied Statistical
Analysis, 5(3), 367–373.

## See also

[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
and
[`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md).

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
#>  level_a level_b difference    se      z       p   p_adj  lower  upper
#>        a       b     -0.907 0.204 -4.441 < 0.001 < 0.001 -1.308 -0.507
#>  significant practical ets
#>            *   >= 0.50  C-
#> p adjusted by holm over 1 pairwise comparison(s); practical criterion 0.50 logits
```
