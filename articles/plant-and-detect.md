# Simulation-based checks of Rasch diagnostics

``` r

library(rasch)
```

## 1. Why simulate

Simulation provides a direct way to examine whether an estimator
recovers known parameters and whether a diagnostic responds to the
departure it is intended to detect. The package therefore includes
simulators for the main model families of Rasch Measurement Theory
(Rasch 1960): `simulate_rasch`, `simulate_btl`, `simulate_mfrm`, and
`simulate_efrm`. Each can generate model-conforming data or introduce a
known departure from the model.

Every simulator returns data ready for its fitting function, with the
generating parameters and any specified departures carried on the object
as `attr(x, "truth")`. The print method reports what was planted, so a
dataset never becomes anonymous once it leaves the call that made it.

``` r

d <- simulate_rasch(n_persons = 400, n_items = 10, seed = 101)
d                                   # the print method reports the plant
#> Simulated rasch data: dichotomous, 400 persons x 10 items
#> Model-conforming (no departures planted).
names(attr(d, "truth"))            # the truth travels with the data
#>  [1] "layout"         "description"    "model"          "n_persons"     
#>  [5] "n_items"        "theta"          "theta2"         "difficulty"    
#>  [9] "thresholds"     "discrimination" "guessing"       "groups"        
#> [13] "dim_items"      "dif_items"      "careless_idx"   "style_idx"     
#> [17] "planted"
```

With no departures requested, the data conform to the model. The
remaining sections introduce one departure at a time and examine the
corresponding diagnostic. The diagnostic conventions follow Andrich and
Marais (2019).

## 2. A model-conforming baseline

It is useful to begin with parameter recovery when the model holds, so
that later departures can be interpreted against a model-conforming
reference. `sim_recovery` compares the parameters a fit recovers against
the ones the simulator planted (mean-centring the locations, which the
model identifies only up to an origin), and `plot_recovery` draws each
true-versus-estimated panel with its correlation and root-mean-square
error.

``` r

fit <- rasch(d, id = "id")
rec <- sim_recovery(fit, d)
rec
#> Parameter recovery (planted vs recovered):
#>   item difficulty  n=10   r=0.998  RMSE=0.108  bias=NA
#>   person ability   n=400  r=0.761  RMSE=0.826  bias=NA
```

``` r

plot_recovery(rec)
```

![Recovery scatter plots of planted against recovered item difficulty
and person
ability.](plant-and-detect_files/figure-html/baseline-plot-1.png)

The two panels tell different stories, and the difference is not a
defect of the estimator but the arithmetic of information. Each item
difficulty pools evidence from all 400 persons, so it is recovered
tightly: the correlation is near unity and the error is around a tenth
of a logit. Each person ability, by contrast, rests on only the ten
items that person answered, so its weighted likelihood estimate (Warm
1989) carries a standard error of roughly eight tenths of a logit; the
recovery correlation is correspondingly lower and the root-mean-square
error much larger. (The bias column reads NA by design: locations are
identified only up to the scale origin, so a mean difference from the
planted values is not an estimable quantity.) A short test measures
persons imprecisely *even when the model is exactly true* — a point
worth holding onto before reading any person-level result from a
ten-item scale.

## 3. Planting item misfit

The `discrimination` argument scales an item’s slope. A value above one
makes an item over-discriminate — its responses become more
deterministic than the model expects, a Guttman-like pattern — while a
value below one makes it under-discriminate, adding noise. The fit
residuals should separate both from the well-behaved items. Here one
central item is given a slope of 2.5 and another a slope of 0.4.

``` r

disc <- rep(1, 10); disc[5] <- 2.5; disc[6] <- 0.4
d2  <- simulate_rasch(400, 10, discrimination = disc, seed = 21)
fit2 <- rasch(d2, id = "id")
fit2$items[, c("item", "location", "infit_ms", "outfit_ms")]
#>    item location infit_ms outfit_ms
#> 1   I01  -2.5220   0.8994    0.7716
#> 2   I02  -1.8851   1.0412    0.9415
#> 3   I03  -1.5071   1.0119    0.8956
#> 4   I04  -0.7011   1.0517    0.9912
#> 5   I05  -0.3867   0.9134    0.8490
#> 6   I06   0.2056   1.2625    1.2397
#> 7   I07   0.7730   1.0445    1.0459
#> 8   I08   1.4631   1.0828    0.9437
#> 9   I09   2.1031   1.0133    0.8363
#> 10  I10   2.4572   1.0382    0.7645
```

The two planted items depart from the pack in opposite directions, as
the log-of-mean-square fit statistics (Andrich and Marais 2019, ch. 23)
require. The over-discriminating item (I05) produces standardised
residuals that are *too small* — responses more predictable than the
model allows — so its infit and outfit mean squares fall below one. The
under-discriminating item (I06) produces residuals that are too large,
so its mean squares rise above one. The remaining items sit near the
expected value of one. With a 400-by-10 draw these mean squares are
themselves estimates; what matters is the direction of departure and
that the planted items separate from the rest.

## 4. Planting DIF

Differential item functioning is planted with `dif`, naming the affected
items and the size of the group shift, with `n_groups = 2` so there is a
group to differ. Here item I06 is given a uniform shift of one logit for
the second group; the other nine items are invariant. `dif_anova`
analyses the residuals by group over automatically sized trait class
intervals (Hagquist and Andrich 2017), and should flag I06 alone.

``` r

d3   <- simulate_rasch(500, 10, dif = list(items = "I06", uniform = 1),
                       n_groups = 2, seed = 303)
fit3 <- rasch(d3, id = "id", factors = "group")
da   <- dif_anova(fit3)
da$summary[, c("item", "term", "F_uniform", "p_uniform_adj", "uniform_DIF")]
#>    item  term F_uniform p_uniform_adj uniform_DIF
#> 1   I01 group  0.164009      0.855993       FALSE
#> 2   I02 group  3.234919      0.242365       FALSE
#> 3   I03 group  0.172138      0.855993       FALSE
#> 4   I04 group  0.412059      0.855993       FALSE
#> 5   I05 group  3.694811      0.242365       FALSE
#> 6   I06 group 18.647136      0.000191        TRUE
#> 7   I07 group  0.611468      0.855993       FALSE
#> 8   I08 group  0.008686      0.925785       FALSE
#> 9   I09 group  0.138649      0.855993       FALSE
#> 10  I10 group  0.085278      0.855993       FALSE
```

Only I06 carries a significant group term after the Benjamini-Hochberg
adjustment across items; the other nine are quiet, as invariant items
should be. The analysis of variance answers *whether* there is DIF. To
ask *how much*, and whether it matters on the measurement scale,
`dif_size` resolves the item into one copy per group, refits, and
reports the distance between the resolved locations in logits (Andrich
and Marais 2019, ch. 16).

``` r

dif_size(fit3, "I06", by = "group")
#> DIF size for I06 by group (resolved locations, logits)
#>  level location    se weak   n
#>     g1    0.262 0.140    0 250
#>     g2    1.266 0.162    0 250
#>  level_a level_b difference    se      z p p_adj  lower  upper significant
#>       g1      g2     -1.004 0.234 -4.287 0     0 -1.463 -0.545           *
#>  practical
#>    >= 0.50
#> p adjusted by holm over 1 pairwise comparison(s); practical criterion 0.50 logits
```

The resolved gap recovers the planted shift of one logit and clears the
half-logit practical criterion. Flagging and magnitude are
complementary: a significant `dif_anova` term says the item is not
invariant; the `dif_size` logit difference says whether the
non-invariance is large enough to act on.

## 5. Planting local dependence

Response dependence — one item’s answer partly following another’s — is
planted with `dependence`, naming the item pairs and a strength. It
should raise the residual correlation of the planted pair above the
small negative value expected under local independence.
`residual_correlations` returns Yen’s (1984) Q3 and, following
Christensen, Makransky and Horton (2017), the adjusted Q3\* (each Q3
less the average off-diagonal value). There is no universal critical
value for adjusted Q3, so binary flagging is available only when the
analyst supplies a screening threshold and remains heuristic.

``` r

d4   <- simulate_rasch(500, 10,
                       dependence = list(pairs = list(c("I04", "I05")),
                                         strength = 1.8), seed = 41)
fit4 <- rasch(d4, id = "id")
rc   <- residual_correlations(fit4, flag = 0.20)
rc$average                          # near -1/(L-1) under independence
#> [1] -0.09841
head(rc$pairs, 3)                   # the planted pair leads the table
#>   item_a item_b       q3 q3_star flagged
#> 1    I04    I05 0.170343  0.2688    TRUE
#> 2    I02    I10 0.042300  0.1407   FALSE
#> 3    I02    I09 0.002477  0.1009   FALSE
rc$flagged
#>   item_a item_b     q3 q3_star flagged
#> 1    I04    I05 0.1703  0.2688    TRUE
```

The planted I04–I05 pair sits at the top of the sorted table and clears
the chosen heuristic screen; the average off-diagonal value is close to
the $`-1/(L-1)`$ expected when the items are locally independent. As
with DIF, screening is one thing and magnitude another.
`dependence_magnitude` puts the effect on the logit scale by the
resolution method of Andrich and Kreiner (2010): the dependent item is
split by the category of the item it follows, refitted, and the
threshold displacement read off.

``` r

dependence_magnitude(fit4, dependent = "I05", independent = "I04")
#> Response dependence of I05 on I04 (Andrich & Kreiner resolution)
#>   d = 0.760 logits (se 0.139), z = 5.48, p = < 0.001
```

The estimated dependence is around three quarters of a logit and highly
significant — a substantial displacement of I05’s threshold by the
response to I04, not a marginal correlation. (The residual principal
components of the same fit, via `plot_scree`, are the complementary tool
for a *second dimension* rather than a dependent pair.)

## 6. Paired comparisons

The paired-comparison simulator plants faults of its own.
`erratic_judges` sets a proportion of judges who choose at random; their
disorder should surface in the judge fit residuals and in the
transitivity of their choices. Here two of eight judges are erratic.

``` r

b  <- simulate_btl(n_objects = 7, n_judges = 8, reps_per_pair = 30,
                   erratic_judges = 0.25, seed = 61)
bt <- btl(b, "object_a", "object_b", winner = "winner", judge = "judge")
bt$judges[order(-bt$judges$fit_resid), ]
#>   judge  n infit_ms outfit_ms fit_resid df_fit
#> 1    J1 74   1.4569    1.5253    3.9953  73.30
#> 2    J2 75   1.3566    1.4168    3.2204  74.29
#> 6    J6 70   1.0708    1.0758    0.5876  69.33
#> 4    J4 78   0.9017    0.8810   -1.1930  77.26
#> 7    J7 77   0.8939    0.8757   -1.3396  76.27
#> 3    J3 81   0.8743    0.8520   -1.5176  80.23
#> 8    J8 92   0.8565    0.8232   -1.8897  91.12
#> 5    J5 83   0.7515    0.7151   -2.9513  82.21
```

The two erratic judges (J1 and J2) carry the largest infit and outfit
mean squares and the largest positive fit residuals; the well-behaved
judges sit below one. The judge fit residual is the paired-comparison
counterpart of person fit, and here it is the sharpest instrument.

`btl_transitivity` asks the single-dimension question in its native
paired-comparison form. A Bradley-Terry-Luce scale implies that
preferences stack into one consistent order: if A beats B and B beats C
then A should beat C. A *circular triad* (A beats B, B beats C, C beats
A) is a local contradiction, and their rate is compared with the one
quarter expected from pure guessing (Kendall and Babington Smith 1940).

``` r

tr <- btl_transitivity(bt)
tr$summary[, c("n_objects", "n_triples", "n_circular",
               "circular_rate", "consistency")]
#>   n_objects n_triples n_circular circular_rate consistency
#> 1         7        30          2       0.06667      0.7333
tr$judges[, c("judge", "n_triples", "circular_rate", "consistency")]
#>   judge n_triples circular_rate consistency
#> 1    J2        21       0.33333     -0.3333
#> 2    J1        25       0.12000      0.5200
#> 3    J6        35       0.11429      0.5429
#> 4    J3        15       0.06667      0.7333
#> 5    J8        30       0.03333      0.8667
#> 6    J4        21       0.00000      1.0000
#> 7    J5        20       0.00000      1.0000
#> 8    J7        16       0.00000      1.0000
```

Pooled across all judges the object scale is nearly transitive —
majority verdicts stack into essentially one order — so the overall
consistency is high. The per-judge index separates the individuals, with
the erratic judges falling toward the bottom of the table. With only a
handful of complete triples per judge this index is noisier than the fit
residual, but the two agree in direction. Where the concern is not judge
noise but a second attribute steering some contests,
`btl_dimensionality` provides the residual-“swirl” analogue of residual
principal components.

## 7. A small Monte Carlo power estimate

The plant-and-detect loop scales up. `sim_replicate` calls a simulator
many times with successive seeds, returning the datasets as a list, so a
diagnostic can be run across all of them and its behaviour summarised.
Asking how often a planted departure is flagged turns that into a
detection-power estimate. Here a uniform DIF of 0.8 logits is planted on
one item across ten datasets, and each is fitted and passed to
`dif_anova`.

``` r

batch <- sim_replicate(simulate_rasch, 10, n_persons = 400, n_items = 8,
                       dif = list(items = "I04", uniform = 0.8),
                       n_groups = 2, seed = 700)
flagged <- vapply(batch, function(dd) {
  s <- dif_anova(rasch(dd, id = "id", factors = "group"))$summary
  isTRUE(s$uniform_DIF[s$item == "I04"])
}, logical(1))
mean(flagged)                       # proportion of runs that flagged I04
#> [1] 0.5
```

The estimate is the proportion of the ten runs in which the planted item
was flagged after adjustment. A shift of 0.8 logits split over two
groups of 200 is detected in only about half of these runs: a real
departure of moderate size against a modest sample is caught only some
of the time. Ten replicates suffice to demonstrate the loop, not to pin
the number down — the Monte Carlo error on ten draws is large — but the
same six lines, run with a few hundred replicates and swept over sample
size and effect, are a complete power study.

## 8. Validation studies

The simulators above are also the package’s own test bench. Each release
is validated by a simulation battery that runs every model family — the
dichotomous model, the polytomous models, many-facet models, the
extended frame of reference model, and the comparative judgement models
— against its full set of diagnostics, under complete data, 25% random
missingness, and linked structural designs: some 240 checks in all.
Parameter recovery is checked against the planted values; reported
standard errors against the empirical sampling variability of the
estimates (all calibration ratios fall within 0.85–1.20, most within a
few percent of 1); significance tests against their nominal rates under
model-true data; and every diagnostic against its planted departure for
power. The identification guards are exercised in both directions:
connected linked designs fit and recover well, while genuinely
disconnected designs are refused rather than fitted.

One correction came out of the battery. `dependence_magnitude` followed
Andrich and Marais (2019, eqs. 24.9–24.11) in pooling the two resolved
items’ variances as if their estimates were independent, on the grounds
that disjoint persons answer them. In the joint refit, however, the two
estimates are negatively correlated through the shared comparator items,
and the pooled standard error is too small: with ten items, the null
rejection rate was 7.5% at a nominal 5% over 1,200 replicates. Since
release 1.14.2 the standard error is the delta-method contrast over the
fit’s sandwich covariance, which restores the nominal rate (5.4%
observed) and carries over unchanged to the polytomous case. The
estimate of *d* itself is unaffected.

## References

Andrich, D., and Kreiner, S. (2010). Quantifying response dependence
between two dichotomous items using the Rasch model. *Applied
Psychological Measurement*, 34, 181–192.

Andrich, D., and Marais, I. (2019). *A Course in Rasch Measurement
Theory*. Springer.

Christensen, K. B., Makransky, G., and Horton, M. (2017). Critical
values for Yen’s Q3. *Applied Psychological Measurement*, 41, 178–194.

Hagquist, C., and Andrich, D. (2017). Recent advances in analysis of
differential item functioning in health research using the Rasch model.
*Health and Quality of Life Outcomes*, 15, 181.

Kendall, M. G., and Babington Smith, B. (1940). On the method of paired
comparisons. *Biometrika*, 31, 324–345.

Rasch, G. (1960). *Probabilistic Models for Some Intelligence and
Attainment Tests*. Copenhagen: Danish Institute for Educational
Research. (Expanded edition, 1980, Chicago: University of Chicago
Press.)

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika*, 54, 427–450.

Yen, W. M. (1984). Effects of local item dependence on the fit and
equating performance of the three-parameter logistic model. *Applied
Psychological Measurement*, 8, 125–145.
