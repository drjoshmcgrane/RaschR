# Plot the item map (location against fit residual)

Fitted columns plotted by location and fit residual, with the
conventional acceptance band at +/- 2.5. MFRM and EFRM points are
response cells; ordinary Rasch points are items.

## Usage

``` r
plot_item_map(fit, band = 2.5)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- band:

  Fit residual acceptance band.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
plot_item_map(rasch(X))
```
