# Plot Bradley-Terry-Luce object locations

Caterpillar plot of object locations with 95 per cent error bars. The
intervals use the reference degrees of freedom carried by the fitted
covariance: a t reference for judge-based uncertainty and the limiting
normal reference otherwise. An interval is omitted when its covariance
does not support inference. Objects beyond the specified fit-residual
band are marked.

## Usage

``` r
plot_btl(fit, band = 2.5)
```

## Arguments

- fit:

  An object from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- band:

  Absolute fit-residual value beyond which an object is highlighted.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pairs <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pairs[, 1], each = 30),
                b = rep(pairs[, 2], each = 30))
p <- plogis(beta[d$a] - beta[d$b])
d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
plot_btl(btl(d, "a", "b", "win"))
```
