# Drop items and refit

Removes named items and refits the analysis with the same model
specification.

## Usage

``` r
drop_items(fit, items, boot_reps = NULL)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) or
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).
  Many-facet fits are refused: remove the item's rows from the
  long-format data and refit
  [`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  instead.

- items:

  Item names to remove.

- boot_reps:

  Bootstrap replicates for an EFRM refit. The default retains the fitted
  specification; a number overrides it.

## Value

A refitted object of the same class as `fit`, carrying a note recording
which items were dropped.

## Details

The refit retains person identifiers and factors, class-interval
settings, optimisation controls, anchors, multiple-choice scoring and
PCM component constraints. An EFRM refit also retains the item-set and
crossed-frame design, linking controls and uncertainty method. The
operation is refused if it would remove every anchor, empty an item set
or leave the model unidentified.

Item removal changes both the item calibration and the person estimates.
For an EFRM it can also change the estimated frame units. Compare the
original and refitted results as a sensitivity analysis.

## See also

[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
and
[`resolve_frames`](https://drjoshmcgrane.github.io/rasch/reference/resolve_frames.md)
for frame models;
[`split_items`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md)
and
[`resolve_dif`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
for DIF; and
[`combine_items`](https://drjoshmcgrane.github.io/rasch/reference/combine_items.md)
for response dependence.

## Examples

``` r
d <- simulate_rasch(300, 8, seed = 1)
fit <- rasch(d, id = "id")
fit2 <- drop_items(fit, "I03")
nrow(fit2$items)
#> [1] 7
```
