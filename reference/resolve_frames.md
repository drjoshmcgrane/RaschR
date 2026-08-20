# Resolve items that do not hold across frames

Gives named items a separate location in each frame and refits, so they
continue to measure persons within their own frame while no longer
constraining the comparison between frames. The result is an ordinary
frame fit, so every diagnostic applies to it unchanged.

## Usage

``` r
resolve_frames(fit, items, boot_reps = NULL)
```

## Arguments

- fit:

  A fitted object from
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

- items:

  Item names to resolve.

- boot_reps:

  Bootstrap replicates for the refit, resolved as in
  [`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md):
  the refit keeps the character of the fit it came from unless a number
  is given.

## Value

A refitted object of class `"rasch_efrm"`, carrying a note for each item
resolved. The resolved versions appear in the item table as
`"item (frame)"`.

## Details

The fitted model represents an item's threshold in a frame as the frame
unit times one location shared across frames, so an item that behaves
differently in one frame is misrepresented in all of them, and the group
units absorb part of the discrepancy.
[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
tests for such items; this function and
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
are the two remedies.

Resolving is the milder one. The item is replaced by one version per
frame, named `"item (frame)"`, each answered by that frame's persons
alone. Each version keeps its own location, so the item still
contributes to the person estimates of everyone who answered it, and it
no longer contributes to the link between frames. Dropping the item
removes that contribution as well, from every frame at once.

The link is what pays for it. Person-group units are identified by the
items two frames have in common, so a resolved item leaves the units
resting on the items that remain shared, and resolving too many leaves
them unidentified. This function refuses when a set would be left with
fewer than two common items; the model's own connectivity check catches
the remaining cases.

Prefer resolving when the item measures well inside each frame and only
its comparability is in doubt, and dropping when the item is a poor
measure wherever it appears. The distinction is empirical: compare the
unit estimates and the person standard errors the two remedies produce.

## See also

[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md),
which identifies the items to resolve;
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md),
which removes an item instead; and
[`split_items`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md),
the equivalent for an ordinary fit.

## Examples

``` r
d <- simulate_efrm(n_per_group = 200, items_per_set = 6, n_sets = 2,
                   n_groups = 2, set_unit_ratio = 1.3, seed = 5)
tr <- attr(d, "truth")
fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                  id = "id", boot_reps = 0)
fit2 <- resolve_frames(fit, "S1I02", boot_reps = 0)
grep("S1I02", fit2$items$item, value = TRUE)
#> [1] "S1I02 (g1):g1" "S1I02 (g2):g2"
```
