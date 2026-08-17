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
natural unit.

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
which to test it.

The model should be used when a frame-dependent unit is part of the
substantive measurement account, not simply because an equal-unit model
fits poorly.

``` r

d <- simulate_efrm(n_per_group = 300, items_per_set = 8,
                   set_unit_ratio = 1.30, group_unit_ratio = 1.10,
                   seed = 12)
truth <- attr(d, "truth")
```

`simulate_efrm(n_categories = 4)` generates partial credit items within
frames instead; polytomous sets carry more linking information per item.

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
#> rasch extended frame of reference analysis: 16 items in 2 set(s) x 2 group(s) = 4 frames, 600 persons
#> Within-frame pairwise conditional ML: converged in 9 iterations
#> PSI 0.781, power of fit: reasonable
#> 
#> Person group units (phi):
#>  group   phi se_log_phi
#>     g1 0.912     0.0314
#>     g2 1.096     0.0314
#> 
#> Item set units (alpha) and locations:
#>   set alpha se_log_alpha     mu n_items
#>  set1 0.845        0.045 -0.028       8
#>  set2 1.183        0.045  0.028       8
#> 
#> Equal-unit comparison: 2(ll_EFRM - ll_equal) = 28.664 with 1 extra unit parameter(s)
#> (composite likelihood: descriptive; informative for group units (phi))
#> Omnibus Wald tests of equal units:
#>               term df   wald       p
#>  group units (phi)  1  8.516   0.004
#>  set units (alpha)  1 14.295 < 0.001
#> Holm-adjusted exploratory unit contrasts (H0: unit = 1):
#>        parameter estimate    se      z       p   p_adj significant
#>      log phi[g1]   -0.092 0.031 -2.918   0.004   0.007           *
#>      log phi[g2]    0.092 0.031  2.918   0.004   0.007           *
#>  log alpha[set1]   -0.168 0.045 -3.781 < 0.001 < 0.001           *
#>  log alpha[set2]    0.168 0.045  3.781 < 0.001 < 0.001           *
#> 
#> Notes: person measures use the weighted score; per-group score curves replace the raw-score table (see score_curves)
```

`phi_table` reports person-group units, `alpha_table` reports item-set
units, and `set_table` reports the linked set locations. Units are
positive and are most naturally compared as ratios. The linking tables
and notes should be read before interpreting common-unit item or person
estimates.

``` r

fit$phi_table
#>   group    phi se_log_phi
#> 1    g1 0.9125    0.03139
#> 2    g2 1.0959    0.03139
fit$alpha_table
#>    set  alpha se_log_alpha
#> 1 set1 0.8451      0.04451
#> 2 set2 1.1833      0.04451
fit$set_table
#>    set     mu  alpha n_items
#> 1 set1 -0.028 0.8451       8
#> 2 set2  0.028 1.1833       8
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
#>   set_a set_b   n log_slope
#> 1  set1  set2 503    0.3366
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
#>    set group n_persons n_items  alpha    phi    rho se_log_rho origin infit_ms
#> 1 set1    g1       300       8 0.8451 0.9125 0.7711    0.05446 -0.028    1.010
#> 2 set2    g1       300       8 1.1833 0.9125 1.0797    0.05446  0.028    1.041
#> 3 set1    g2       300       8 0.8451 1.0959 0.9262    0.05446 -0.028    1.043
#> 4 set2    g2       300       8 1.1833 1.0959 1.2968    0.05446  0.028    1.024
#>   outfit_ms fit_resid n_responses
#> 1    1.0036    0.6195        2344
#> 2    0.9711   -0.7375        2344
#> 3    1.0126    0.7063        2304
#> 4    0.8693   -2.0092        2304
fit$efrm_vs_rasch
#> $ll_efrm
#> [1] -6667
#> 
#> $ll_equal
#> [1] -6681
#> 
#> $two_delta_ll
#> [1] 28.66
#> 
#> $extra_parameters
#> [1] 1
#> 
#> $informative_for
#> [1] "group units (phi)"
#> 
#> $unit_omnibus
#>                term df   wald         p
#> 1 group units (phi)  1  8.516 0.0035212
#> 2 set units (alpha)  1 14.295 0.0001563
#> 
#> $unit_tests
#>         parameter estimate      se      z         p     p_adj significant
#> 1     log phi[g1] -0.09159 0.03139 -2.918 0.0035212 0.0070423        TRUE
#> 2     log phi[g2]  0.09159 0.03139  2.918 0.0035212 0.0070423        TRUE
#> 3 log alpha[set1] -0.16829 0.04451 -3.781 0.0001563 0.0006251        TRUE
#> 4 log alpha[set2]  0.16829 0.04451  3.781 0.0001563 0.0006251        TRUE
```

The raw equal-unit composite-likelihood difference is descriptive. The
omnibus Wald tests in `efrm_vs_rasch$unit_omnibus` test the unit
families; individual unit contrasts are Holm-adjusted follow-ups.

## References

Humphry, S. M. (2005). *Maintaining a Common Arbitrary Unit in Social
Measurement*. PhD thesis, Murdoch University.

Humphry, S. M., and Andrich, D. (2008). Understanding the unit in the
Rasch model. *Journal of Applied Measurement*, 9(3), 249–264.

Rasch, G. (1961). On general laws and the meaning of measurement in
psychology. In *Proceedings of the Fourth Berkeley Symposium on
Mathematical Statistics and Probability* (Vol. 4, pp. 321–333).
Berkeley: University of California Press.
