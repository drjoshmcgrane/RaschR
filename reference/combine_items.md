# Combine items into subtests and re-analyse

Replaces each nominated item group by a polytomous super-item whose
score is the sum of its members, then refits the model. The function is
commonly used to examine item groups identified by
[`residual_correlations`](https://drjoshmcgrane.github.io/rasch/reference/residual_correlations.md).
Every total from zero to the sum of the component maxima must be
observed; otherwise the refit is refused rather than renumbering the
superitem score.

## Usage

``` r
combine_items(fit, groups, model = "PCM")
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- groups:

  A list of character vectors, each naming two or more items to combine;
  a single vector is also accepted.

- model:

  Model for the re-analysis; defaults to `"PCM"`, which is almost always
  required because subtests change the maximum scores.

## Value

A new
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) fit
on the combined structure, with the combinations recorded in its notes.
Person and item estimates are recalculated. Fit grouping, external
anchors on unchanged items, keyed scoring, PCM constraints and
optimisation controls are retained. Anchored items cannot be combined
because the resulting superitem has no corresponding external anchor. A
group-specific copy produced by
[`split_items()`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md)
cannot be combined: form the subtest before applying a DIF split.
Existing split-item provenance is retained for items not included in a
new subtest.

## See also

[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
to remove an item rather than combine it, and
[`residual_correlations`](https://drjoshmcgrane.github.io/rasch/reference/residual_correlations.md)
for the dependence that motivates combining.

## Examples

``` r
set.seed(1); Np <- 500; L <- 8
d <- seq(-2, 2, length.out = L)
X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
colnames(X) <- paste0("I", 1:L)
X[, 5] <- ifelse(runif(Np) < 0.9, X[, 4], X[, 5])   # dependent pair
fit <- rasch(X)
fit2 <- combine_items(fit, list(c("I4", "I5")))
fit2$items$item
#> [1] "I1"    "I2"    "I3"    "I6"    "I7"    "I8"    "I4+I5"
```
