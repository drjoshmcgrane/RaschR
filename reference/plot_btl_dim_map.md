# Residual map of the leading paired-comparison bimension

Objects placed in the leading bimension plane. Reading round the swirl,
an object sits “upstream” of those it over-beats relative to the fitted
locations; a clear rotational arrangement is the second attribute, a
formless blob near the origin is noise. Point size grows with the
object's location on the primary scale.

## Usage

``` r
plot_btl_dim_map(x, ...)
```

## Arguments

- x:

  A `"rasch_btl_dim"` object.

- ...:

  Unused.

## Value

Called for its plotting side effect.

## Examples

``` r
# \donttest{
d <- simulate_btl(7, 12, reps_per_pair = 20, seed = 1)
fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
dimensions <- btl_dimensionality(fit, reps = 20)
plot_btl_dim_map(dimensions)

# }
```
