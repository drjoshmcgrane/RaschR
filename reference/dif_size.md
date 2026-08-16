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
`practical`), the settings, and any notes. Sampling-uncertainty fields
are `NA` when person IDs repeat.

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

## References

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. Scandinavian Journal of Statistics, 6(2), 65–70.

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
#>  level_a level_b difference    se      z p p_adj  lower  upper significant
#>        a       b     -0.907 0.204 -4.441 0     0 -1.308 -0.507           *
#>  practical
#>    >= 0.50
#> p adjusted by holm over 1 pairwise comparison(s); practical criterion 0.50 logits
```
