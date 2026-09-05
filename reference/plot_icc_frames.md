# Plot an item's characteristic curves across frames

Plots the model expected-score curve for one item in each frame, with
observed class-interval means overlaid. Differences between the model
curves reflect the fitted frame units. A nominated non-frame person
factor separates the observed means within each frame for a DIF display.

## Usage

``` r
plot_icc_frames(fit, item, n_groups = fit$n_groups, grid = NULL, group = NULL)
```

## Arguments

- fit:

  A fitted object from
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

- item:

  Underlying item name.

- n_groups:

  Number of class intervals for the observed means.

- grid:

  Logit grid.

- group:

  Optional person grouping vector, or one or more names of non-frame
  factors nominated in the fit. Several names define their
  factor-combination cells.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
# \donttest{
# see ?rasch_efrm for a complete simulated example
# }
```
