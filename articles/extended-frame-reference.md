# The extended frame of reference model

``` r

library(rasch)
```

## Frames and units

Rasch’s criterion of invariant comparison holds within a specified frame
of reference (Rasch 1961). The extended frame of reference model
(Humphry 2005; Humphry and Andrich 2008) makes the frame explicit: it
allows the unit to differ across linked item-set by person-group frames.
Within each frame, the partial credit model holds in that frame’s
natural unit; dichotomous items are its one-threshold special case, and
the example below uses them.

For item \\i\\ in set \\s\\ and person \\n\\ in group \\g\\,

\\ P(X\_{ni}=x)= \frac{\exp\left\\\rho\_{sg}\left\[x\theta_n-
\sum\_{k=1}^{x}\delta\_{ik}\right\]\right\\}
{\sum\_{y=0}^{m_i}\exp\left\\\rho\_{sg}\left\[y\theta_n-
\sum\_{k=1}^{y}\delta\_{ik}\right\]\right\\}, \qquad
\rho\_{sg}=\alpha_s\phi_g. \\

Here \\\alpha_s\\ is the item-set unit and \\\phi_g\\ is the
person-group unit. A frame-defining factor sets the unit and cannot also
be tested for DIF: the factor takes a single level among the persons
responding within any frame, so no within-frame comparison remains on
which to test it. This is why
[`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
refuses such a factor, and why the fitted model cannot check its own
assumption: it holds each item at one location across frames by
construction, so the per-frame estimates it reports agree because the
model made them agree. Testing the assumption means stepping outside the
model, which is what
[`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
does – see below.

The model should be used when a frame-dependent unit is part of the
substantive measurement account, not simply because an equal-unit model
fits poorly.

``` r

d <- simulate_efrm(n_per_group = 500, items_per_set = 12,
                   set_unit_ratio = 1.30, group_unit_ratio = 1.10,
                   seed = 25)
truth <- attr(d, "truth")
```

The simulator defaults to dichotomous items;
`simulate_efrm(n_categories = 4)` generates partial credit items within
frames instead, which carry more linking information per item.

## Fit and inspect the links

``` r

fit <- rasch_efrm(
  d,
  item_sets = truth$item_sets,
  groups = "group",
  id = "id",
  boot_reps = 60
)
fit
#> rasch extended frame of reference analysis: 24 items in 2 set(s) x 2 group(s) = 4 frames, 1000 persons
#> Within-frame pairwise conditional ML: converged in 8 iterations
#> PSI 0.845, power of fit: good
#> 
#> Person group units (phi):
#>  group   phi se_log_phi
#>     g1 0.953      0.019
#>     g2 1.050      0.019
#> 
#> Item set units (alpha) and locations:
#>   set alpha se_log_alpha     mu n_items
#>  set1 0.880        0.030  0.005      12
#>  set2 1.137        0.030 -0.005      12
#> 
#> Equal-unit comparison: 2(ll_EFRM - ll_equal) = 30.368 with 1 extra unit parameter(s)
#> (composite likelihood: descriptive; informative for group units (phi))
#> Omnibus Wald tests of equal units:
#>               term df   wald       p
#>  group units (phi)  1  6.509   0.011
#>  set units (alpha)  1 18.349 < 0.001
#> Holm-adjusted exploratory unit contrasts (H0: unit = 1):
#>        parameter estimate    se      z       p   p_adj significant
#>      log phi[g1]   -0.049 0.019 -2.551   0.011   0.021           *
#>      log phi[g2]    0.049 0.019  2.551   0.011   0.021           *
#>  log alpha[set1]   -0.128 0.030 -4.284 < 0.001 < 0.001           *
#>  log alpha[set2]    0.128 0.030  4.284 < 0.001 < 0.001           *
#> 
#> Notes: person measures use the weighted score; per-group score curves replace the raw-score table (see score_curves)
```

`phi_table` reports person-group units, `alpha_table` reports item-set
units, and `set_table` reports the linked set locations. Units are
positive and are most naturally compared as ratios. The linking tables
and notes should be read before interpreting common-unit item or person
estimates.

Both recovered units sit close to the planted 1.30 and 1.10 here, but
read that as one draw rather than as the estimator’s precision. Units
are ratios of estimated dispersions and carry real sampling error: at
this design repeated draws give set-unit ratios with a standard
deviation of about 0.055 and group-unit ratios about 0.038, so analyses
routinely land several per cent from the generating value in either
direction. This example uses a seed whose draw falls near the centre of
that distribution; a different seed would land elsewhere within it.

The two units are estimated differently and behave differently. Group
units come from the within-frame conditional calibration: unbiased, with
well-calibrated intervals (no detectable bias over 60 replicates at this
design, 95% coverage). Set units come from the linking step and
additionally carry a small fixed-design offset, around one per cent at
these set lengths.

Judge a recovered unit against its standard error, then, and never
against a single alternative value alone. Whether the *estimator* is on
target is a question about repeated draws, not about one:

``` r

ratios <- vapply(seq_len(20), function(r) {
  dd <- simulate_efrm(n_per_group = 500, items_per_set = 12,
                      set_unit_ratio = 1.30, group_unit_ratio = 1.10,
                      seed = 91000 + r)
  ff <- rasch_efrm(dd, item_sets = attr(dd, "truth")$item_sets,
                   groups = "group", id = "id", boot_reps = 0)
  c(alpha = max(ff$alpha_table$alpha) / min(ff$alpha_table$alpha),
    phi = max(ff$phi_table$phi) / min(ff$phi_table$phi))
}, c(alpha = 0, phi = 0))
rowMeans(ratios)   # 1.304 and 1.107 against the planted 1.30 and 1.10
apply(ratios, 1, sd)
```

Longer item sets, not larger samples, are what tighten a set unit: the
precision of the set-linking step is limited by the number of items in
each set once a few hundred linking persons are available.

``` r

fit$phi_table
#>  group   phi se_log_phi
#>     g1 0.953      0.019
#>     g2 1.050      0.019
fit$alpha_table
#>   set alpha se_log_alpha
#>  set1 0.880        0.030
#>  set2 1.137        0.030
fit$set_table
#>   set     mu alpha n_items
#>  set1  0.005 0.880      12
#>  set2 -0.005 1.137      12
fit$linking
#> $phi_edges
#> $phi_edges[[1]]
#> [1] 2 1
#> 
#> $phi_edges[[2]]
#> [1] 2 1
#> 
#> 
#> $alpha_edges
#>  set_a set_b   n log_slope
#>   set1  set2 913     0.256
```

Group units are estimated from within-frame pairwise conditional
calibrations. Set units and locations require persons common to the
linked sets. A disconnected design cannot establish a common unit.

The set-linking step uses a method-of-moments estimator based on the
true-score variance-ratio argument in Humphry (2005), with the
true-score variance recovered by a truncated-score-moment correction:
the mean and variance of the score-to-measure map over the non-extreme
scores are exact model functions of the person location, and their
expectations are estimated through score weights whose validity does not
depend on the person distribution. This matters at any realistic length:
the classical “observed variance minus mean squared standard error”
construction under-recovers the true-score variance by more than half at
eight dichotomous items per set, and its errors do not fade with length
– in simulation it biased recovered unit ratios by about +5 per cent at
eight items and -2 to -4 per cent at sixteen to thirty-two, never
converging in the tested range, and it fails entirely below seven items.
The corrected estimator is unbiased at every tested set length, unit
ratio (1 to 2), and person distribution (skewness to 2.8, where a
normal-population marginal maximum likelihood anchor drifts), and agrees
with independent external anchors. Its remaining fixed-design offset is
small – about 0.4 per cent on the unit ratio at eight dichotomous items
per set, decaying with set length – and intervals stay calibrated to at
least five thousand linking persons; beyond that, accuracy is bounded by
the set length, not the sample. Each linking person’s pattern must span
a score range of at least four within a set (four dichotomous items; six
or more are recommended), and because the correction computes score
distributions from the fitted within-frame model, model violations
within a set bias the recovered units roughly in proportion. The
estimator is distinct from the likelihood proposed in section 5.3 of
that thesis. The crossed, multigroup, and polytomous forms are
extensions implemented in this package.

## Uncertainty and comparison with equal units

The default hybrid standard errors combine the pairwise sandwich
covariance, a person bootstrap for the set-linking stage, and
delta-method propagation. Across linking samples from 250 to 10,000
persons the reported standard errors track the empirical sampling
variability as both shrink with the square root of the sample. The
example uses 60 replicates to keep the vignette quick; a final analysis
should use enough replicates to stabilise its reported standard errors.
Set `se_method = "bootstrap"` to refit the complete pipeline on each
person resample.

``` r

fit$frames
#>   set group n_persons n_items alpha   phi   rho se_log_rho origin infit_ms
#>  set1    g1       500      12 0.880 0.953 0.838      0.035  0.005    1.012
#>  set2    g1       500      12 1.137 0.953 1.083      0.035 -0.005    1.022
#>  set1    g2       500      12 0.880 1.050 0.924      0.035  0.005    1.019
#>  set2    g2       500      12 1.137 1.050 1.193      0.035 -0.005    1.032
#>  outfit_ms fit_resid n_responses
#>      1.000     0.367        5964
#>      0.978    -0.817        5964
#>      0.995     0.156        5856
#>      0.973    -0.843        5856
fit$efrm_vs_rasch
#> $ll_efrm
#> [1] -27194
#> 
#> $ll_equal
#> [1] -27209
#> 
#> $two_delta_ll
#> [1] 30.37
#> 
#> $extra_parameters
#> [1] 1
#> 
#> $informative_for
#> [1] "group units (phi)"
#> 
#> $unit_omnibus
#>               term df   wald       p
#>  group units (phi)  1  6.509   0.011
#>  set units (alpha)  1 18.349 < 0.001
#> 
#> $unit_tests
#>        parameter estimate    se      z       p   p_adj significant
#>      log phi[g1]   -0.049 0.019 -2.551   0.011   0.021           *
#>      log phi[g2]    0.049 0.019  2.551   0.011   0.021           *
#>  log alpha[set1]   -0.128 0.030 -4.284 < 0.001 < 0.001           *
#>  log alpha[set2]    0.128 0.030  4.284 < 0.001 < 0.001           *
```

The raw equal-unit composite-likelihood difference is descriptive. The
omnibus Wald tests in `efrm_vs_rasch$unit_omnibus` test the unit
families; individual unit contrasts are Holm-adjusted follow-ups.

## Testing the invariance the model assumes

[`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
calibrates each frame separately, puts the item locations on the common
scale by dividing by that frame’s unit, and compares them item by item.
It reports two comparisons, and they are not equally sensitive: a
location difference is found about 97 per cent of the time at 500
persons per frame, a difference in discrimination about 16 per cent, and
the two reach comparable power only near 2,000. A clean result at a few
hundred persons per frame is therefore much stronger evidence against
differential item functioning than against differential discrimination.

``` r

inv <- frame_invariance(fit)
inv$summary                    # rmsd against rmse, per set and frame pair
inv$locations
```

Which items to flag is a screening decision, not a confirmatory one, and
`adjust` chooses. Holm across the ten items measured is right for
reporting that an item differs and wrong for deciding which items to
look at, where it costs between 8 and 42 points of sensitivity. Both
probabilities are reported either way, so the choice changes only which
rows are flagged.

``` r

frame_invariance(fit, adjust = "none")   # screen
```

An item the test flags is removed with
[`drop_items()`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md),
which refits and returns an ordinary fit of the same class. Rank the
items by `abs(fit_resid)` rather than applying a fixed cut-off – a cut
states detectability rather than magnitude, so it flags more items the
more persons there are – then drop the largest departure, refit, and
stop when the units no longer differ. Dropping is a complete cure where
the item is found; what limits the repair is detection, not removal.

``` r

fit2 <- drop_items(fit, "I07")
fit2$efrm_vs_rasch$unit_tests            # has the difference gone?
```

## Worked analyses on real data

Three case studies apply the frame of reference models to real datasets
and ship with the package. The two wording studies are best read
together: they find a set-unit difference of almost the same size and it
means something different in each. In the self-esteem data the
difference is carried by one ambivalent item and vanishes when that item
goes; in the height data it survives every single-item removal, and yet
leaves the person ordering essentially unchanged against a criterion
collected outside the inventory. Neither conclusion is available from
the unit test alone.

``` r

# wording effects as item-set units in a self-esteem scale
file.edit(system.file("casestudies", "wording_units_selfesteem.R",
                      package = "rasch"))
# the same question in a balanced inventory, checked against a criterion
file.edit(system.file("casestudies", "wording_units_height.R",
                      package = "rasch"))
# panel and bloc frames in paired comparisons between political parties
file.edit(system.file("casestudies", "party_blocs_crisis.R",
                      package = "rasch"))
```

## References

Humphry, S. M. (2005). *Maintaining a Common Arbitrary Unit in Social
Measurement*. PhD thesis, Murdoch University.

Humphry, S. M., and Andrich, D. (2008). Understanding the unit in the
Rasch model. *Journal of Applied Measurement*, 9(3), 249–264.

Rasch, G. (1961). On general laws and the meaning of measurement in
psychology. In *Proceedings of the Fourth Berkeley Symposium on
Mathematical Statistics and Probability* (Vol. 4, pp. 321–333).
Berkeley: University of California Press.
