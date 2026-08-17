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
which to test it.

The model should be used when a frame-dependent unit is part of the
substantive measurement account, not simply because an equal-unit model
fits poorly.

``` r

d <- simulate_efrm(n_per_group = 500, items_per_set = 12,
                   set_unit_ratio = 1.30, group_unit_ratio = 1.10,
                   seed = 1)
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
#> PSI 0.844, power of fit: good
#> 
#> Person group units (phi):
#>  group   phi se_log_phi
#>     g1 0.959     0.0186
#>     g2 1.043     0.0186
#> 
#> Item set units (alpha) and locations:
#>   set alpha se_log_alpha     mu n_items
#>  set1 0.853        0.026 -0.008      12
#>  set2 1.172        0.026  0.008      12
#> 
#> Equal-unit comparison: 2(ll_EFRM - ll_equal) = 22.250 with 1 extra unit parameter(s)
#> (composite likelihood: descriptive; informative for group units (phi))
#> Omnibus Wald tests of equal units:
#>               term df   wald       p
#>  group units (phi)  1  5.122   0.024
#>  set units (alpha)  1 36.845 < 0.001
#> Holm-adjusted exploratory unit contrasts (H0: unit = 1):
#>        parameter estimate    se      z       p   p_adj significant
#>      log phi[g1]   -0.042 0.019 -2.263   0.024   0.047           *
#>      log phi[g2]    0.042 0.019  2.263   0.024   0.047           *
#>  log alpha[set1]   -0.159 0.026 -6.070 < 0.001 < 0.001           *
#>  log alpha[set2]    0.159 0.026  6.070 < 0.001 < 0.001           *
#> 
#> Notes: person measures use the weighted score; per-group score curves replace the raw-score table (see score_curves)
```

`phi_table` reports person-group units, `alpha_table` reports item-set
units, and `set_table` reports the linked set locations. Units are
positive and are most naturally compared as ratios. The linking tables
and notes should be read before interpreting common-unit item or person
estimates.

Compare the recovered units with the planted 1.30 and 1.10. One is
close; the other is not, and that is the point worth absorbing. Units
are ratios of estimated dispersions, and at ordinary designs they carry
real sampling error: at this one, repeated draws give set-unit ratios
with a standard deviation of about 0.055 and group-unit ratios about
0.038, so a single analysis lands several per cent from the generating
value as a matter of course. The two units are estimated differently and
behave differently: group units come from the within-frame conditional
calibration and are unbiased with well-calibrated intervals (no
detectable bias over 60 replicates at this design, 95% coverage), while
set units additionally carry a small fixed-design offset of about one
per cent at these set lengths.

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
#>   group    phi se_log_phi
#> 1    g1 0.9588    0.01858
#> 2    g2 1.0429    0.01858
fit$alpha_table
#>    set  alpha se_log_alpha
#> 1 set1 0.8533      0.02613
#> 2 set2 1.1719      0.02613
fit$set_table
#>    set        mu  alpha n_items
#> 1 set1 -0.008391 0.8533      12
#> 2 set2  0.008391 1.1719      12
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
#> 1  set1  set2 913    0.3172
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
#>    set group n_persons n_items  alpha    phi    rho se_log_rho    origin
#> 1 set1    g1       500      12 0.8533 0.9588 0.8182    0.03206 -0.008391
#> 2 set2    g1       500      12 1.1719 0.9588 1.1236    0.03206  0.008391
#> 3 set1    g2       500      12 0.8533 1.0429 0.8900    0.03206 -0.008391
#> 4 set2    g2       500      12 1.1719 1.0429 1.2222    0.03206  0.008391
#>   infit_ms outfit_ms fit_resid n_responses
#> 1    1.004    0.9871  -0.10701        5940
#> 2    1.042    1.0099   0.03523        5940
#> 3    1.022    1.0129   0.94924        5892
#> 4    1.008    0.9366  -1.76491        5892
fit$efrm_vs_rasch
#> $ll_efrm
#> [1] -27555
#> 
#> $ll_equal
#> [1] -27566
#> 
#> $two_delta_ll
#> [1] 22.25
#> 
#> $extra_parameters
#> [1] 1
#> 
#> $informative_for
#> [1] "group units (phi)"
#> 
#> $unit_omnibus
#>                term df   wald         p
#> 1 group units (phi)  1  5.122 2.363e-02
#> 2 set units (alpha)  1 36.845 1.279e-09
#> 
#> $unit_tests
#>         parameter estimate      se      z         p     p_adj significant
#> 1     log phi[g1] -0.04204 0.01858 -2.263 2.363e-02 4.726e-02        TRUE
#> 2     log phi[g2]  0.04204 0.01858  2.263 2.363e-02 4.726e-02        TRUE
#> 3 log alpha[set1] -0.15859 0.02613 -6.070 1.279e-09 5.115e-09        TRUE
#> 4 log alpha[set2]  0.15859 0.02613  6.070 1.279e-09 5.115e-09        TRUE
```

The raw equal-unit composite-likelihood difference is descriptive. The
omnibus Wald tests in `efrm_vs_rasch$unit_omnibus` test the unit
families; individual unit contrasts are Holm-adjusted follow-ups.

## Worked analyses on real data

Two case studies apply the frame of reference models to real datasets
and ship with the package:

``` r

# wording effects as item-set units in a self-esteem scale
file.edit(system.file("casestudies", "wording_units_selfesteem.R",
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
