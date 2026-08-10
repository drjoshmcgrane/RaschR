# Recovery scatter of planted against recovered parameters

One true-versus-estimated panel per parameter type, with the identity
line and the correlation and RMSE.

## Usage

``` r
plot_recovery(x, ...)
```

## Arguments

- x:

  A `"rasch_recovery"` object.

- ...:

  Unused.

## Value

Called for its plotting side effect.

## Examples

``` r
# \donttest{
d <- simulate_rasch(300, 8, seed = 1)
fit <- rasch(d, id = "id")
plot_recovery(sim_recovery(fit, d))

# }
```
