# Plot an object characteristic curve

Plots the expected response for one object against opponent location,
with the observed mean response against each sufficiently observed
opponent. For dichotomous fits the curve is the win probability; for
ordered fits it is the expected response. With fitted position or
history effects, model values are weighted means of the fitted
expectations for each observed opponent, joined in location order.
Judge-group overlays use each group's own comparison context.

## Usage

``` r
plot_btl_icc(fit, object, group = NULL, grid = NULL, min_n = 10)
```

## Arguments

- fit:

  An object from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- object:

  Object name.

- group:

  Optional judge grouping for a DIF overlay: either one value per
  comparison row of `fit$comparisons` or a vector named by every judge
  in the fit. Observed means are then drawn separately per group, as
  [`plot_icc`](https://drjoshmcgrane.github.io/rasch/reference/plot_icc.md)
  draws person groups.

- grid:

  Opponent-location grid, in logits.

- min_n:

  An opponent's observed point is drawn only when the object (or, in the
  grouped display, that judge group) met it at least this many times;
  sparser pairs from incomplete or unbalanced designs are omitted.

## Value

Called for its plotting side effect; invisibly the names of the
opponents drawn (the ungrouped display), or `NULL` for the grouped
display.

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
p <- plogis(beta[d$a] - beta[d$b])
d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
plot_btl_icc(btl(d, "a", "b", winner = "win"), "C")
```
