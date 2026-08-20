# Resolve differential item functioning by iterative item splitting

Splits DIF items one at a time, beginning with the largest estimated
effect, and refits after each split. This order addresses the artificial
DIF that a large departure can induce in otherwise invariant items
(Andrich and Hagquist 2012, 2015). Each split gives the item a separate
location and threshold structure in every factor cell. The procedure
stops when no item is flagged or the remaining anchor set reaches
`min_anchors`.

## Usage

``` r
resolve_dif(
  fit,
  factors = NULL,
  alpha = 0.05,
  p_adjust = "BH",
  min_anchors = NULL,
  max_splits = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  carrying person factors.

- factors:

  Person factors to test, as in
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md);
  defaults to every nominated factor.

- alpha:

  Significance level for the adjusted probabilities.

- p_adjust:

  Multiplicity adjustment across items each round.

- min_anchors:

  Minimum number of original items to leave unsplit; the procedure stops
  before the anchor set falls below this (pervasive DIF is not
  artificial DIF). Default `max(3, items / 4)`.

- max_splits:

  Hard cap on the number of splits. Default: the number of items.

## Value

A list of class `"rasch_resolve_dif"`: the final resolved `fit`, the
`splits` performed (order, item, factor, partial eta-squared, DIF
magnitude in logits), the `stopped` reason, and the residual `dif` table
for the final fit.

## References

Andrich, D., & Hagquist, C. (2012). Real and artificial differential
item functioning. *Journal of Educational and Behavioral Statistics*,
37(3), 387-416.

## See also

[`split_items`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md)
for a single split,
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
to remove an item instead, and
[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
for the test it resolves.

## Examples

``` r
set.seed(1); n <- 600
d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1.2      # one strong DIF item
X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(data.frame(X, grp = g), factors = "grp")
resolve_dif(fit)$splits
#>  order item factor  eta2 magnitude
#>      1   I3    grp 0.054     1.111
```
