# Compare fitted and generating parameters

Compares fitted parameters with the generating values from a
`simulate_*` function (carried on the data as `attr(sim, "truth")`):
item difficulties and person abilities for a Rasch fit, object locations
for a paired-comparison fit, rater severities (with item and person
measures) for a many-facet fit, set and group units for a Rasch frames
fit, and common object locations, panel and set units, and set origins
for a paired-comparison frames fit. Locations are mean-centred where the
model identifies them only up to an origin. An externally anchored Rasch
or paired-comparison fit retains its identified origin, so recovery and
bias are reported on the anchored scale. The fit must be from the
simulated model family and must have converged. For a many-facet
simulation, the planted rater facet must be identifiable uniquely by its
name or level labels. EFRM set parameters are matched by their item or
object membership, not by the spelling of the set labels; a different
fitted partition is refused.

## Usage

``` r
sim_recovery(fit, sim)
```

## Arguments

- fit:

  A fit of the simulated data
  ([`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md),
  [`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md),
  or
  [`btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)).

- sim:

  The simulated data (from a `simulate_*` function).

## Value

A list of class `"rasch_recovery"`: `summary` (per parameter type: n,
correlation, RMSE, bias) and `pieces` (the true and estimated values
behind each).

## Examples

``` r
d <- simulate_rasch(500, 12, seed = 1)
sim_recovery(rasch(d, id = "id"), d)$summary
#>         parameter   n correlation      rmse bias
#> 1 item difficulty  12   0.9980002 0.1703475   NA
#> 2  person ability 500   0.8177074 0.6756437   NA
```
