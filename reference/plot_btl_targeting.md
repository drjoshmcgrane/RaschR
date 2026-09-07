# Targeting plot for a paired-comparison design

The paired-comparison counterpart of a test-information display. Every
object is a dot at its location (x) and its design information (y, the
pooled Fisher information of the comparisons it took part in), the dot
sized by how many comparisons that is. For an equal-unit fit, a
reference curve on the right axis traces the information a single *new*
comparison would carry against an opponent at each location. For a
dichotomous fit without a position effect, it peaks at gap zero. Ordered
response thresholds and position effects can change where it peaks. A
frame fit has no single reference curve because the information also
depends on the fitted panel and set units. The reference curve is also
omitted for history-dependent fits. With a position effect only, it
represents the median-location object presented first against each
opponent location.

## Usage

``` r
plot_btl_targeting(fit, grid = NULL)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- grid:

  Optional location grid for the equal-unit reference curve.

## Value

Called for its plotting side effect; invisibly `NULL`.

## See also

[`btl_information`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md),
[`btl_next_pairs`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
plot_btl_targeting(btl(d, "a", "b", "win"))
```
