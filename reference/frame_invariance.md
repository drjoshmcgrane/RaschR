# Test the item invariance a frame model assumes

Calibrates each frame separately and compares each item across the
frames it appears in, on two counts: whether it keeps its location once
the frame units are accounted for, and whether it discriminates alike. A
frame model assumes both, differing only by the frame's unit; this
function tests the assumption rather than imposing it.

## Usage

``` r
frame_invariance(fit, alpha = 0.05, adjust = c("holm", "none"))
```

## Arguments

- fit:

  A fitted object from
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

- alpha:

  Significance level for flagging items.

- adjust:

  Multiplicity adjustment used to flag items: `"holm"` across all
  comparisons, or `"none"` for screening. Both probabilities are
  reported regardless.

## Value

A list of class `"rasch_frame_invariance"` with `locations` (one row per
item and frame pair: locations on the common scale, their difference,
its standard error, statistic, both probabilities, and flag),
`discrimination` (the same items compared on their within-frame infit
statistics, with the Winsteps-style discrimination index for each frame
and its ratio alongside), and `summary` (per frame pair: the number of
items, the root mean squared difference, the root mean squared standard
error, their ratio, and the number of items flagged on each count).

## Details

The comparison is possible only where an item set is taken by more than
one person group, since an item must appear in at least two frames to be
compared across them. Item sets partition the items, so there is no
equivalent test across sets: screen those with the ordinary item fit
statistics within each set instead.

Each frame is refitted with
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) on
its own persons and items, giving locations in that frame's natural unit
and centred on that frame's origin. Dividing by the frame unit from the
original fit puts them on the common scale, where the model says they
should agree. The reported difference is between the two calibrations,
its standard error combines theirs (the person samples are disjoint, so
the calibrations are independent), and probabilities are Holm-adjusted
across items.

The summary compares the root mean squared difference with the root mean
squared standard error, following Humphry (2005). A root mean squared
difference materially larger than the root mean squared standard error
indicates item behaviour that the frame units do not account for,
whether or not individual items reach significance.

A location comparison cannot detect a difference in discrimination: a
steeper item still crosses one half in the same place, so its location
survives intact. The second comparison uses the within-frame fit
statistics, which carry no unit because each is computed against its own
frame's model, and treats the difference of two independent standardised
infit statistics as having variance 2. It is conservative and needs a
reasonable sample: in simulation, with 8 items and two items
discriminating half again as steeply in one frame, it detected them in
12% of replicates at 500 persons per frame and 85% at 2,000, with
false-positive rates near 1% throughout.

The two comparisons are therefore not equally sensitive, and the gap
matters because the two departures do comparable damage. At departures
that each move a unit ratio by six or seven per cent – two items shifted
a logit, against two items discriminating half again as steeply – 500
persons per frame detected the shifted items about 97% of the time and
the steeper items about 16%. Only by around 2,000 persons per frame do
the two reach comparable power. A clean result at a few hundred persons
per frame is thus much stronger evidence against differential item
functioning than against differential discrimination: it has ruled out
one departure and barely tested the other.

The discrimination table also reports a Winsteps-style index for each
frame (`disc_1`, `disc_2`) and their ratio, because a standardised
difference says only that something differs while the index says how
much and in which direction. It is fitted by maximum likelihood on each
item's own responses with the person measures and item location held at
their Rasch values, so it is relative to the frame's own model and the
unit cancels. Read it as description rather than estimate: the measures
are estimated including the item being scored, which biased the index to
about 1.19 for a true discrimination of 1.0 in simulation and attenuated
a true frame ratio of 1.5 to about 1.22. Tested on its own it is also
the weaker instrument, detecting a 1.5-fold difference in 63% of
replicates at 2,000 persons per frame against the infit comparison's
90%, which is why the test column comes from the latter.

Which item to flag is a screening decision, not a confirmatory one, and
the two call for different thresholds. `adjust` chooses: Holm across
every item and frame pair, or none. Both probabilities are reported
either way, so the choice changes only `flagged`. Screening 8 to 10
items with Holm costs between 20 and 60 points of sensitivity in
simulation, and that shows up in the repair: dropping the flagged items
from a planted unit ratio of 1.40 left the ratio at 1.479 under Holm and
1.433 unadjusted, against 1.406 for dropping the items actually planted.
The loose screen is not free – where misfit is strong it flags sound
items too, and
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
then refuses drops that would empty a set – so use `"none"` to decide
which items to examine and `"holm"` to report which ones differ.

Dropping is a complete cure where the item is found: removing the
planted items restored the ratio in every simulated departure, so what
limits the repair is detection rather than removal. Past roughly a fifth
of the items breaking invariance, no threshold rescues the ratio and the
item set itself is the problem.

## References

Humphry, S. M. (2005). *Maintaining a Common Arbitrary Unit in Social
Measurement*. PhD thesis, Murdoch University.

## See also

[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
to remove an item the test flags, and
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
for the model whose assumption is tested.

## Examples

``` r
d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 1,
                   n_groups = 2, group_unit_ratio = 1.4, seed = 2)
tr <- attr(d, "truth")
fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                  id = "id", boot_reps = 0)
frame_invariance(fit)
#> Item invariance across frames (each frame calibrated separately)
#> 
#>   set frame_1 frame_2 n_items  rmsd  rmse ratio n_location n_discrimination
#>  set1      g1      g2       8 0.199 0.200 0.994          0                0
#> 
#> rmsd/rmse above 1 indicates item behaviour the frame units do not account for
#> 
#> No item's location differs across frames at alpha = 0.05 (Holm-adjusted).
#> 
#> No item's discrimination differs across frames.
```
