# Simulate paired-comparison data

Generates dichotomous or ordered paired comparisons from the
Bradley–Terry–Luce model. Optional arguments introduce a second object
attribute, erratic judges, or within-judge dependence. Generating values
are stored in `attr(x, "truth")`.

## Usage

``` r
simulate_btl(
  n_objects = 8,
  n_judges = 12,
  reps_per_pair = 25,
  model = c("dichotomous", "polytomous", "graded"),
  n_categories = 4,
  object_sd = 1,
  second_attribute = NULL,
  erratic_judges = 0,
  dependence = NULL,
  seed = NULL
)
```

## Arguments

- n_objects, n_judges:

  Objects to scale and judges comparing them.

- reps_per_pair:

  Comparisons made of each object pair.

- model:

  `"dichotomous"` (a winner) or `"polytomous"` (a rated margin in
  `n_categories` categories; an earlier development-era value `"graded"`
  is accepted as an alias).

- n_categories:

  Categories for the polytomous model.

- object_sd:

  Spread of the object locations (evenly spaced, sum-zero).

- second_attribute:

  `NULL`, or `list(rho=)`: half the judges rank by a second object
  attribute correlated `rho` with the first. This introduces residual
  dimensionality and possible intransitivity.

- erratic_judges:

  Proportion of judges who choose at random.

- dependence:

  `NULL`, or `list(exposure=, carry_over=)`: within-judge order effects
  (a seen-before advantage and a pull from the judge's own earlier
  verdicts). Adds an `order` column. Feeds the dependence effects fitted
  by [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- seed:

  Optional RNG seed.

## Value

A data frame of class `"rasch_sim"`: `object_a`, `object_b`, `winner`
(or `response` when polytomous), `judge`, and `order` when dependence is
planted; with `attr(x, "truth")`.

## Examples

``` r
d <- simulate_btl(8, 12, erratic_judges = 0.15, seed = 1)
bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
bt$judges          # the erratic judges carry large fit residuals
#>  judge  n infit_ms outfit_ms fit_resid df_fit
#>     J1 64    1.230     1.286     1.401 63.360
#>    J10 61    0.828     0.792    -1.344 60.390
#>    J11 41    0.752     0.653    -1.860 40.590
#>    J12 62    0.911     0.877    -0.760 61.380
#>     J2 50    1.438     1.851     3.050 49.500
#>     J3 59    0.933     0.899    -0.565 58.410
#>     J4 65    1.093     1.195     1.179 64.350
#>     J5 70    1.008     1.087     0.564 69.300
#>     J6 60    0.914     0.893    -0.731 59.400
#>     J7 54    1.048     1.015     0.083 53.460
#>     J8 53    1.043     1.053     0.293 52.470
#>     J9 61    0.831     0.799    -1.422 60.390
```
