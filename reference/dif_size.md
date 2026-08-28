# DIF differences between factor levels

Resolves an item by one or more person factors and compares the
resulting locations. Several factors in `by` give pairwise comparisons
between their joint cells and can be used to quantify an interaction.

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
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) or
  [`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md).
  EFRM fits are excluded because an ordinary split refit would discard
  their frame units.

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

A list of class `"rasch_dif_size"`. `levels` contains the resolved
location, standard error and sample size for each level. `pairs`
contains logit differences, Wald statistics, confidence intervals, raw
and adjusted probabilities, and practical flags. For dichotomous items
it also contains `ets`; for polytomous items it contains the descriptive
`signed_area`. Sampling-uncertainty fields are `NA` when person
identifiers repeat.

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

## Magnitude conventions

For dichotomous items, `ets` applies the ETS A, B and C rules to the
itemwise comparison. On the logit scale the magnitude cut-points are
\\1/2.35=0.426\\ and \\1.5/2.35=0.638\\. Category A also includes an
itemwise test that is not significant. Category C requires a magnitude
of at least 0.638 and rejection of the interval null \\\|\Delta\|\leq
0.426\\; B is the remainder. The category uses the raw itemwise
probability, while `significant` uses `p_adjust` over the requested
pairwise family.

For a partial credit item with \\m_i\\ thresholds, the signed area
between the two expected-score curves has the closed form
\$\$SA\_{ab}=\int\\E_b(X\mid\theta)-E_a(X\mid\theta)\\\\d\theta
=\sum\_{k=1}^{m_i}(\delta\_{iak}-\delta\_{ibk})
=m_i(\beta\_{ia}-\beta\_{ib}).\$\$ This is returned as `signed_area`; a
positive value means that level `a` has the harder resolved item. It is
descriptive and is not given an A/B/C category: score-metric
classifications for polytomous DIF are not interchangeable with a PCM
logit difference. For pooled MFRM items, the areas use the same
precision weight for a facet cell in every group. A comparison is
withheld when the groups do not support the same observed response
categories.

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

Cohen, A. S., Kim, S.-H. and Baker, F. B. (1993). Detection of
differential item functioning in the graded response model. Applied
Psychological Measurement, 17(4), 335–350.

Raju, N. S. (1988). The area between two item characteristic curves.
Psychometrika, 53(4), 495–502.

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
#>  significant practical ets signed_area
#>            *   >= 0.50  C-            
#> p adjusted by holm over 1 pairwise comparison(s); practical criterion 0.50 logits
```
