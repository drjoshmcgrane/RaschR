# Drop items and refit

Removes named items from a fitted analysis and refits it, keeping the
original model, person identifiers, factors, and (for frame models) the
set structure and standard-error method. The result is an ordinary fit
of the same class, so every diagnostic applies to it unchanged.

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

  Bootstrap replicates for the refit. The default keeps the character of
  the fit it came from: a fit whose unit standard errors came from a
  bootstrap passes its replicate count on, a fit that has unit standard
  errors by the analytic route takes the package default, and a fit
  asked for no standard errors is refitted without them. Pass a number
  to override all three.

## Value

A refitted object of the same class as `fit`, carrying a note recording
which items were dropped.

## Details

Screening items is part of an analysis of frames rather than a step
before it. A person-group unit is estimated from the same items in every
frame, so an item with differential item functioning distorts it; an
item-set unit is estimated from the dispersion of person estimates
within each set, so an item that fits its set badly distorts that set
alone, with nothing in the other set to offset it. In simulation the
second effect is large: at eight items per set, two under-discriminating
items in one set moved a planted unit ratio of 1.40 to 1.73, and four
over-discriminating items moved it to 1.02. Dropping such an item and
comparing the units before and after is therefore a substantive
sensitivity analysis, not housekeeping.

An item that fits no set well is usually better removed than reassigned,
and it cannot be given a set of its own: a single item carries no
dispersion from which to estimate a unit.

## See also

[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md),
which identifies the items a frame model's assumption does not hold for;
[`resolve_frames`](https://drjoshmcgrane.github.io/rasch/reference/resolve_frames.md),
which gives such an item a location per frame rather than removing it;
[`split_items`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md)
and
[`resolve_dif`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md),
which resolve an item rather than remove it; and
[`combine_items`](https://drjoshmcgrane.github.io/rasch/reference/combine_items.md).

## Examples

``` r
d <- simulate_rasch(300, 8, seed = 1)
fit <- rasch(d, id = "id")
fit2 <- drop_items(fit, "I03")
nrow(fit2$items)
#> [1] 7
```
