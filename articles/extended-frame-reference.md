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
natural unit. This is a different claim from ordinary DIF: a
frame-defining factor sets the scale unit and cannot also be tested as a
separate DIF factor on the frame’s own virtual items.

The model should be used when a frame-dependent unit is part of the
substantive measurement account, not simply because an equal-unit model
fits poorly.

``` r

d <- simulate_efrm(n_per_group = 300, items_per_set = 8,
                   set_unit_ratio = 1.30, group_unit_ratio = 1.10,
                   seed = 12)
truth <- attr(d, "truth")
```

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
#> PSI 0.780, power of fit: reasonable
#> 
#> Person group units (phi):
#>  group   phi se_log_phi
#>     g1 0.912     0.0314
#>     g2 1.096     0.0314
#> 
#> Item set units (alpha) and locations:
#>   set alpha se_log_alpha     mu n_items
#>  set1 0.822        0.032 -0.027       8
#>  set2 1.217        0.032  0.027       8
#> 
#> Equal-unit comparison: 2(ll_EFRM - ll_equal) = 28.664 with 1 extra unit parameter(s)
#> (composite likelihood: descriptive; informative for group units (phi))
#> Omnibus Wald tests of equal units:
#>               term df   wald       p
#>  group units (phi)  1  8.516   0.004
#>  set units (alpha)  1 37.025 < 0.001
#> Holm-adjusted exploratory unit contrasts (H0: unit = 1):
#>        parameter estimate    se      z       p   p_adj significant
#>      log phi[g1]   -0.092 0.031 -2.918   0.004   0.007           *
#>      log phi[g2]    0.092 0.031  2.918   0.004   0.007           *
#>  log alpha[set1]   -0.196 0.032 -6.085 < 0.001 < 0.001           *
#>  log alpha[set2]    0.196 0.032  6.085 < 0.001 < 0.001           *
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
#> 1 set1 0.8217      0.03228
#> 2 set2 1.2170      0.03228
fit$set_table
#>    set       mu  alpha n_items
#> 1 set1 -0.02699 0.8217       8
#> 2 set2  0.02699 1.2170       8
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
#> 1  set1  set2 503    0.3928
```

The group units are obtained from within-frame pairwise conditional
calibrations. Set units and locations require persons common to the
linked sets. A disconnected or weakly linked design cannot establish a
common unit; increasing the number of responses within an isolated frame
does not repair the missing link.

The person-side set-linking step is an error-corrected method-of-moments
implementation based on the true-score variance-ratio argument in
Humphry (2005), rather than the separate likelihood proposal in section
5.3 of that thesis. The crossed, multigroup and polytomous
implementation is an experimental extension. Its sampling behaviour
should be checked for the intended design, preferably with the full
person bootstrap. At this vignette’s design and sample size, for
example, the recovered set-unit ratio runs a modest few per cent above
the planted value across replications – a finite-sample bias of the
moments-based link, not sampling noise alone – which is exactly the kind
of behaviour such a check reveals before the estimates are interpreted
substantively.

## Uncertainty and comparison with equal units

The default hybrid standard errors combine the pairwise sandwich
covariance, a person bootstrap for the set-linking stage, and
delta-method propagation. The example uses 60 replicates to keep the
vignette quick; a final analysis should use enough replicates to
stabilise its reported standard errors. Set `se_method = "bootstrap"` to
refit the complete pipeline on each person resample.

``` r

fit$frames
#>    set group n_persons n_items  alpha    phi    rho se_log_rho   origin
#> 1 set1    g1       300       8 0.8217 0.9125 0.7498    0.04502 -0.02699
#> 2 set2    g1       300       8 1.2170 0.9125 1.1105    0.04502  0.02699
#> 3 set1    g2       300       8 0.8217 1.0959 0.9005    0.04502 -0.02699
#> 4 set2    g2       300       8 1.2170 1.0959 1.3338    0.04502  0.02699
#>   infit_ms outfit_ms fit_resid n_responses
#> 1    1.002    0.9934    0.4098        2344
#> 2    1.049    0.9767   -0.6696        2344
#> 3    1.034    1.0026    0.5449        2304
#> 4    1.032    0.8728   -1.9346        2304
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
#> 1 group units (phi)  1  8.516 3.521e-03
#> 2 set units (alpha)  1 37.025 1.167e-09
#> 
#> $unit_tests
#>         parameter estimate      se      z         p     p_adj significant
#> 1     log phi[g1] -0.09159 0.03139 -2.918 3.521e-03 7.042e-03        TRUE
#> 2     log phi[g2]  0.09159 0.03139  2.918 3.521e-03 7.042e-03        TRUE
#> 3 log alpha[set1] -0.19642 0.03228 -6.085 1.167e-09 4.666e-09        TRUE
#> 4 log alpha[set2]  0.19642 0.03228  6.085 1.167e-09 4.666e-09        TRUE
```

The raw equal-unit composite-likelihood difference is descriptive.
Inference is carried by the omnibus Wald tests in
`efrm_vs_rasch$unit_omnibus`; individual unit contrasts are exploratory
and Holm-adjusted. None of these is a substitute for checking
frame-level fit, targeting, link strength, and the defensibility of the
frame interpretation.

## References

Andrich, D., and Marais, I. (2019). *A Course in Rasch Measurement
Theory: Measuring in the Educational, Social and Health Sciences*.
Springer.

Humphry, S. M. (2005). *Maintaining a Common Arbitrary Unit in Social
Measurement*. PhD thesis, Murdoch University.

Humphry, S. M., and Andrich, D. (2008). Understanding the unit in the
Rasch model. *Journal of Applied Measurement*, 9(3), 249–264.

Rasch, G. (1960). *Probabilistic Models for Some Intelligence and
Attainment Tests*. Copenhagen: Danish Institute for Educational
Research. (Expanded edition, 1980, Chicago: University of Chicago
Press.)

Rasch, G. (1961). On general laws and the meaning of measurement in
psychology. In *Proceedings of the Fourth Berkeley Symposium on
Mathematical Statistics and Probability* (Vol. 4, pp. 321–333).
Berkeley: University of California Press.
