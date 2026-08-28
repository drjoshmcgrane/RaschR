# Simulate extended frame-of-reference data with differing units

Generates data whose latent unit differs across item-set by person-group
frames (Humphry 2005): a person in group g responding to an item in set
s does so at the frame unit rho = alpha_set \* phi_group scaling the
whole exponent. The planted set- and group-unit ratios are recovered by
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

## Usage

``` r
simulate_efrm(
  n_per_group = 300,
  items_per_set = 8,
  n_sets = 2,
  n_groups = 2,
  set_unit_ratio = 1.3,
  group_unit_ratio = 1,
  n_categories = 2,
  theta_sd = 1.3,
  seed = NULL
)
```

## Arguments

- n_per_group:

  Persons in each group.

- items_per_set:

  Items in each set.

- n_sets, n_groups:

  Numbers of item sets and person groups.

- set_unit_ratio, group_unit_ratio:

  Geometric span of the set and group units across their levels (1 =
  equal units, i.e. an ordinary Rasch fit).

- n_categories:

  Response categories per item: 2 (the default) gives dichotomous items;
  larger values give partial credit items whose evenly spaced thresholds
  are centred on the item locations, with the frame unit scaling the
  whole exponent as in the dichotomous case.

- theta_sd:

  Spread of person ability.

- seed:

  Optional RNG seed.

## Value

A wide data frame of class `"rasch_sim"`, containing an ID, item
columns, and group. Its truth attribute contains the item-set map
required by
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

## Examples

``` r
d <- simulate_efrm(200, 6, set_unit_ratio = 1.3, seed = 1)
tr <- attr(d, "truth")
ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                 boot_reps = 0)    # point estimates only
ef$alpha_table   # planted ratio 1.3, recovered within small-sample noise
#>   set alpha se_log_alpha
#>  set1 0.850             
#>  set2 1.176             
```
