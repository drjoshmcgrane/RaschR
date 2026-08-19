# Fit the extended frame of reference model for paired comparisons

Fits paired comparisons when judges belong to panels and objects belong
to linked sets whose units or origins can differ. It combines the
Bradley–Terry–Luce model with Humphry's extended frame of reference
structure.

## Usage

``` r
btl_efrm(
  data,
  object_a,
  object_b,
  winner,
  judge,
  panels,
  object_sets,
  response = NULL,
  ties = c("drop", "error"),
  min_link = 20,
  se_method = c("judge_bootstrap", "bootstrap", "conditional"),
  boot_reps = 200,
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- data:

  A data frame with one comparison per row.

- object_a, object_b:

  Names of the columns holding the two compared objects.

- winner:

  Name of the winner column. A value must match one of the two objects
  in that row. `"tie"` and `"draw"` mark ties; other values are treated
  as missing.

- judge:

  Name of the judge column (clusters the stage-one standard errors and
  defines the panels when `panels` is a judge attribute).

- panels:

  Either the name of a judge-attribute column in `data` or a named
  vector mapping judge to panel.

- object_sets:

  A named list mapping set names to character vectors of object names;
  every compared object must belong to exactly one set.

- response:

  Not supported: this first implementation fits dichotomous winner data
  only. Supplying it raises an informative error.

- ties:

  `"drop"` (default, removed with a note) or `"error"`.

- min_link:

  Minimum number of cross-set comparisons a set pair must supply to be
  used for linking; sets not reachable from the reference set through
  sufficient cross-set pairs raise an error.

- se_method:

  Method used for standard errors. The default, `"judge_bootstrap"`,
  resamples judges within panels and retains dependence among a judge's
  comparisons. `"bootstrap"` instead draws independent outcomes from
  fitted probabilities. Both stages are refitted. `"conditional"` uses
  analytic stage-one standard errors for `beta` and `phi`, and inverse
  observed information for `alpha` and `kappa` conditional on the
  stage-one estimates. It is faster, but does not propagate stage-one
  uncertainty into the linking parameters.

- boot_reps:

  Number of replicates for `se_method = "bootstrap"` or
  `"judge_bootstrap"`.

- maxit, tol:

  Newton iteration cap and convergence tolerance.

## Value

An object of class `"rasch_btl_efrm"`. Principal components are
`objects`, `phi_table`, `alpha_table`, `kappa_table`, `unit_omnibus`,
`frames`, `equal_unit`, `n_cross`, `notes`, and `converged`.

## Details

For object \\k\\ in set \\s\\, let \$\$v_k=\alpha_s\beta_k+\kappa_s,\$\$
where \\\beta_k\\ is its within-set location, \\\alpha_s\>0\\ is the set
unit, and \\\kappa_s\\ is the set origin. A comparison in panel \\g\\
has logit \$\$\phi_g(\beta_a-\beta_b)\$\$ for objects in the same set,
and \$\$\phi_g(v_a-v_b)\$\$ for objects in different sets. Cross-set
comparisons identify the common scale. The first set fixes \\\alpha=1\\
and \\\kappa=0\\; panel units have geometric mean one.

Estimation has two stages. Within-set comparisons estimate object
locations and panel-unit ratios. Weighted least squares reconciles the
ratios over the panel-by-set linking graph. Cross-set comparisons then
estimate the set units and origins. Unlike the person-by-item EFRM, this
linking step uses only comparison outcomes and does not require a
distribution of persons. The paired-comparison form is an extension of
Humphry's model implemented in this package.

The default judge bootstrap resamples judges within panels and refits
both stages. The parametric bootstrap draws independent outcomes from
the fitted probabilities. `se_method = "conditional"` uses analytic
stage-one errors and conditions the linking errors on stage one; it is
intended for preliminary inspection. Bootstrap failures and boundary
estimates are reported in `notes`.

With one set, the model contains panel units only. With one set and one
panel, it reduces to
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md). Omnibus
Wald tests provide inference for the unit families; individual contrasts
are Holm-adjusted follow-ups.

## References

Andrich, D. (1978). Relationships between the Thurstone and Rasch
approaches to item scaling. Applied Psychological Measurement, 2(3),
451–462.

Bradley, R. A. and Terry, M. E. (1952). Rank analysis of incomplete
block designs: I. The method of paired comparisons. Biometrika, 39,
324–345.

David, H. A. (1988). The Method of Paired Comparisons (2nd ed.).
Griffin.

Humphry, S. M. (2005). Maintaining a common arbitrary unit in social
measurement. PhD thesis, Murdoch University.

Humphry, S. M. (2012). Item set discrimination and the unit in the Rasch
model. Journal of Applied Measurement, 13(2), 165–180.

Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
Rasch model. Journal of Applied Measurement, 9(3), 249–264.

Luce, R. D. (1959). Individual Choice Behavior: A Theoretical Analysis.
Wiley.

Thurstone, L. L. (1927). A law of comparative judgment. Psychological
Review, 34, 273–286.

## See also

[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md),
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md),
[`plot_btl_units`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_units.md),
and
[`simulate_btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl_efrm.md).

## Examples

``` r
# \donttest{
d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                       set_units = c(1, 1.4), set_origins = c(0, 0.8),
                       seed = 1)
fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
                judge = "judge", panels = "panel",
                object_sets = attr(d, "truth")$object_sets,
                se_method = "conditional")
fit$alpha_table
#>   set alpha se_log_alpha     t df       p   p_adj significant
#>  set1 1.000                    11                            
#>  set2 1.607        0.077 6.140 11 < 0.001 < 0.001           *
# }
```
