# Equate two test calibrations through their common items

Places two calibrations on a common origin using their shared items,
then tests the shared items for drift. The reference may be a fitted
model or an item bank.

## Usage

``` r
equate_tests(fit, reference, shift = c("mean", "none"), independent = NULL)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- reference:

  A second
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  fit, or a data frame with columns `item`, `location`, and optionally
  `se`. Item names and column names must be unique. Numeric fields may
  be numeric columns, numeric text, or factors with numeric labels;
  other column classes are refused. Locations must be finite. For
  bank-based drift inference with an estimated mean shift, attach the
  bank's joint item-location covariance as a square matrix in
  `attr(reference, "cov_location")`, ordered like the bank rows (or
  named by item); marginal SEs alone do not carry the centring
  covariance. They are sufficient with `shift = "none"`. A bank treated
  as fixed may instead have zero SEs. A polytomous bank must also
  include `max`, the maximum item score. A bank covariance estimated
  from a finite number of independent sampling units may carry its
  positive residual degrees of freedom in
  `attr(reference, "df_location")`.

- shift:

  `"mean"` (default) allows a scale shift between the two analyses;
  `"none"` compares raw locations, appropriate when both analyses are
  already on a shared (anchored) scale.

- independent:

  Whether the two calibrations use independent sampling units. For two
  fitted objects this must be stated explicitly: the default `NULL`
  withholds drift tests because cross-fit covariance is otherwise
  unknown. A bank table is treated as independent unless `FALSE` is
  supplied. When `FALSE`, descriptive equating is returned but
  inferential drift columns are withheld.

## Value

A list with the comparison `table` (locations, standard errors,
difference, its standard error, t statistic, reference degrees of
freedom, raw and Holm-adjusted p, drift flag), the estimated `shift` and
its `shift_method`, the location `correlation`, the root mean square
difference after shifting (`rmsd`), the number of common items
`n_common`, the number with usable standard errors `n`, and whether
drift inference was available (`inferential`). The `note` component
records exclusions and the reason inference was withheld, where
applicable. An individual drift probability is also withheld when its
contrast has zero estimated uncertainty. Common items with unavailable
drift probabilities remain in the multiplicity family.

## Details

Let \\d_j\\ be the location difference for common item \\j\\ and \\v_j\\
its marginal variance. With `shift = "mean"` the scale shift is the
precision-weighted mean \$\$\hat s=\frac{\sum_j d_j/v_j}{\sum_j
1/v_j},\$\$ and each item is tested using \\d_j-\hat s\\ with a variance
that accounts for the estimated shift through the items' joint
covariance. If fewer than two common items have usable variances but at
least two have finite locations, the function returns their unweighted
mean difference as a descriptive fallback and records
`shift_method = "unweighted"`. An exact common anchor determines the
shift even when it is the only common item with usable uncertainty;
unavailable SEs do not override it. When the shift is estimated, drift
inference requires independent calibrations and at least three common
items with usable, positive-semidefinite joint covariance information.
With `shift = "none"`, the origin is fixed before the comparison and
each item's variance is the sum of its two marginal variances; joint
covariance information and a three-item link are then unnecessary. One
common item is sufficient for that fixed-origin comparison; estimating a
shift still requires at least two. Otherwise the function returns a
descriptive link. A fitted calibration's empirical covariance must also
pass the informative-person count, effective-support and projected-rank
checks used by
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).
Supported independent rows use the limiting normal reference. When
either covariance comes from repeated person clusters, drift
probabilities use contrast-specific Welch–Satterthwaite degrees of
freedom. The corresponding residual degrees of freedom are the number of
independent person clusters minus one. A fixed anchor contributes zero
variance and does not consume cluster degrees of freedom. Fitted
calibrations must have converged.

## Examples

``` r
set.seed(1); d <- seq(-1.5, 1.5, length.out = 8)
mk <- function() {
  X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
  colnames(X) <- paste0("I", 1:8); rasch(X)
}
eq <- equate_tests(mk(), mk(), independent = TRUE)
eq$table
#>   item location_1      se_1 location_2      se_2  difference adj_difference
#> 1   I1 -1.5760281 0.1292289 -1.4780195 0.1246379 -0.09800863    -0.10727215
#> 2   I2 -1.1398494 0.1155055 -1.1006867 0.1168127 -0.03916271    -0.04842624
#> 3   I3 -0.5314218 0.1083620 -0.8227470 0.1123970  0.29132523     0.28206171
#> 4   I4 -0.2287656 0.1063021 -0.2085726 0.1093047 -0.02019291    -0.02945643
#> 5   I5  0.3229450 0.1053953  0.1917743 0.1077741  0.13117072     0.12190720
#> 6   I6  0.5366178 0.1097245  0.7088231 0.1124051 -0.17220529    -0.18146881
#> 7   I7  1.1874751 0.1180524  1.1279739 0.1189861  0.05950128     0.05023776
#> 8   I8  1.4290269 0.1263670  1.5814546 0.1284329 -0.15242768    -0.16169120
#>     se_diff          t  df          p     p_adj drift
#> 1 0.1844083 -0.5817100 Inf 0.56076205 1.0000000 FALSE
#> 2 0.1640769 -0.2951436 Inf 0.76788421 1.0000000 FALSE
#> 3 0.1541893  1.8293214 Inf 0.06735147 0.5388118 FALSE
#> 4 0.1489597 -0.1977476 Inf 0.84324253 1.0000000 FALSE
#> 5 0.1466567  0.8312417 Inf 0.40583712 1.0000000 FALSE
#> 6 0.1551459 -1.1696652 Inf 0.24213573 1.0000000 FALSE
#> 7 0.1692070  0.2969011 Inf 0.76654200 1.0000000 FALSE
#> 8 0.1850922 -0.8735713 Inf 0.38235177 1.0000000 FALSE
```
