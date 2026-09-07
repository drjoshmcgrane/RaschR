# Split items by a person factor to resolve DIF

Replaces each nominated item with one item per estimable level of a
person factor, each carrying that level's responses only (other levels
missing). A level must contain every score from zero to the fitted item
maximum; levels missing a score are omitted and recorded in the notes.
The split is refused if calibration subsequently merges a category
lacking conditional information within a retained group. Every retained
group then receives its own item location, which resolves the invariance
violation flagged by
[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md);
the distance between the split locations estimates the DIF size. The
model is refitted with the same settings.

## Usage

``` r
split_items(fit, items, by)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- items:

  Character vector naming the item(s) to split.

- by:

  The name of a person factor nominated in the fit, or a grouping vector
  with one entry per person.

## Value

A new
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) fit
in which each split item appears as `"item (level)"`, with the splits
recorded in its notes. Person and item estimates are recalculated with
the fitted grouping, keyed scoring, anchors on unchanged items, PCM
constraints and optimisation controls. An anchored item cannot be split.

## See also

[`resolve_dif`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md),
which applies this iteratively;
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md),
which removes an item rather than resolving it; and
[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md),
which identifies the items to split.

## Examples

``` r
set.seed(1); n <- 600
d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1
X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(data.frame(X, grp = g), factors = "grp")
fit2 <- split_items(fit, "I3", by = "grp")
fit2$items$item
#> [1] "I1"     "I2"     "I4"     "I5"     "I6"     "I7"     "I8"     "I3 (a)"
#> [9] "I3 (b)"
```
