# Simulation-based checks of Rasch diagnostics

``` r

library(rasch)
```

## Simulating from the models

Simulation is useful when the sampling behaviour of an estimate or
diagnostic depends on the test design. The package includes simulators
for the ordinary Rasch models (Rasch 1960; Andrich and Marais 2019),
many-facet models, extended frames of reference (Humphry 2005; Humphry
and Andrich 2008), and comparative judgement of paired comparisons
(Andrich 1978; Tutz 1986). They can generate model-conforming data or
introduce a specified departure.

Each simulator stores the generating values in `attr(x, "truth")`.

``` r

d <- simulate_rasch(n_persons = 400, n_items = 10, seed = 101)
names(attr(d, "truth"))
#>  [1] "layout"         "description"    "model"          "n_persons"     
#>  [5] "n_items"        "theta"          "theta2"         "difficulty"    
#>  [9] "thresholds"     "discrimination" "guessing"       "groups"        
#> [13] "dim_items"      "dif_items"      "careless_idx"   "style_idx"     
#> [17] "planted"
```

## Parameter recovery

`sim_recovery` compares fitted parameters with the generating values.
Location parameters are centred before comparison because their origin
is arbitrary.

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

![Planted and recovered item and person
locations.](plant-and-detect_files/figure-html/baseline-plot-1.png)

Item locations pool information over persons. Each person location is
based on the items answered by that person, so person recovery is
usually less precise on a short test. This difference should be judged
against the reported standard errors rather than the raw recovery
correlations alone.

## Item misfit

The `discrimination` argument changes an item’s response slope. Values
above one produce more deterministic responses than the Rasch model
expects; values below one produce less deterministic responses.

``` r

disc <- rep(1, 10)
disc[5] <- 2.5
disc[6] <- 0.4

d2 <- simulate_rasch(400, 10, discrimination = disc, seed = 21)
fit2 <- rasch(d2, id = "id")
fit2$items[, c("item", "location", "infit_ms", "outfit_ms")]
#>  item location infit_ms outfit_ms
#>   I01   -2.522    0.899     0.772
#>   I02   -1.885    1.041     0.941
#>   I03   -1.507    1.012     0.896
#>   I04   -0.701    1.052     0.991
#>   I05   -0.387    0.913     0.849
#>   I06    0.206    1.262     1.240
#>   I07    0.773    1.044     1.046
#>   I08    1.463    1.083     0.944
#>   I09    2.103    1.013     0.836
#>   I10    2.457    1.038     0.765
```

Over-discrimination tends to give mean-square statistics below one;
under-discrimination tends to give values above one. Their sampling
variation still depends on the item location, sample, and test length.

## Differential item functioning

The `dif` argument shifts selected items for a person group. Here I06
differs by one logit in the second group.

``` r

d3 <- simulate_rasch(
  500, 10,
  dif = list(items = "I06", uniform = 1),
  n_groups = 2,
  seed = 303
)

fit3 <- rasch(d3, id = "id", factors = "group")
da <- dif_anova(fit3)
da$summary[, c("item", "term", "F_uniform", "p_uniform_adj", "uniform_DIF")]
#>  item  term F_uniform p_uniform_adj uniform_DIF
#>   I01 group     0.187         1.000            
#>   I02 group     3.197         1.000            
#>   I03 group     0.175         1.000            
#>   I04 group     0.410         1.000            
#>   I05 group     3.717         0.980            
#>   I06 group    17.499       < 0.001           *
#>   I07 group     0.573         1.000            
#>   I08 group     0.009         1.000            
#>   I09 group     0.137         1.000            
#>   I10 group     0.082         1.000
```

`dif_anova` tests invariance. `dif_size` resolves the item by group and
reports the difference between the resolved locations in logits.

``` r

dif_size(fit3, "I06", by = "group")
#> DIF size for I06 by group (resolved locations, logits)
#>  level location    se weak   n
#>     g1    0.262 0.140    0 250
#>     g2    1.266 0.162    0 250
#>  level_a level_b difference    se      z       p   p_adj  lower  upper
#>       g1      g2     -1.004 0.234 -4.287 < 0.001 < 0.001 -1.463 -0.545
#>  significant practical ets signed_area
#>            *   >= 0.50  C-            
#> p adjusted by holm over 1 pairwise comparison(s); practical criterion 0.50 logits
```

A simulation study of DIF should record false-positive rates for
invariant items as well as detection of the shifted item. Sample size,
group imbalance, targeting, test length, and shift size should be varied
separately.

## Local response dependence

The `dependence` argument makes one item’s response partly follow
another. `residual_correlations` reports Yen’s Q3 and adjusted Q3.
Because adjusted Q3 has no universal critical value, `flag` is a
screening threshold supplied by the analyst (Yen 1984; Christensen,
Makransky and Horton 2017).

``` r

d4 <- simulate_rasch(
  500, 10,
  dependence = list(
    pairs = list(c("I04", "I05")),
    strength = 1.8
  ),
  seed = 41
)

fit4 <- rasch(d4, id = "id")
rc <- residual_correlations(fit4, flag = 0.20)
head(rc$pairs, 3)
#>   item_a item_b       q3 q3_star flagged
#> 1    I04    I05 0.170343  0.2688    TRUE
#> 2    I02    I10 0.042300  0.1407   FALSE
#> 3    I02    I09 0.002477  0.1009   FALSE
```

`dependence_magnitude` resolves the dependent item by the response to
the independent item and expresses the displacement on the logit scale
(Andrich and Kreiner 2010).

``` r

dependence_magnitude(fit4, dependent = "I05", independent = "I04")
#> Response dependence of I05 on I04 (Andrich & Kreiner resolution)
#>   d = 0.760 logits (se 0.139), z = 5.48, p = < 0.001
```

## Paired comparisons

The paired-comparison simulator can introduce erratic judges, ties,
position effects, or within-judge dependence. In this example, one
quarter of the judges respond at random.

``` r

b <- simulate_btl(
  n_objects = 7,
  n_judges = 8,
  reps_per_pair = 30,
  erratic_judges = 0.25,
  seed = 61
)

bt <- btl(b, "object_a", "object_b", winner = "winner", judge = "judge")
bt$judges[order(-bt$judges$fit_resid), ]
#>  judge  n infit_ms outfit_ms fit_resid df_fit
#>     J1 74    1.457     1.525     3.995 73.295
#>     J2 75    1.357     1.417     3.220 74.286
#>     J6 70    1.071     1.076     0.588 69.333
#>     J4 78    0.902     0.881    -1.193 77.257
#>     J7 77    0.894     0.876    -1.340 76.267
#>     J3 81    0.874     0.852    -1.518 80.229
#>     J8 92    0.856     0.823    -1.890 91.124
#>     J5 83    0.752     0.715    -2.951 82.210
```

Judge fit describes agreement with the common object scale. Transitivity
is a different summary: it counts circular triads in the observed
comparisons.

``` r

tr <- btl_transitivity(bt)
tr$summary
#>  n_objects n_pairs n_triples n_circular circular_rate chance_rate consistency
#>          7      21        30          2         0.067       0.250       0.733
#>  zeta
#> 
head(tr$judges)
#>  judge n_comparisons n_triples n_circular circular_rate consistency
#>     J2            75        21          7         0.333      -0.333
#>     J1            74        25          3         0.120       0.520
#>     J6            70        35          4         0.114       0.543
#>     J3            81        15          1         0.067       0.733
#>     J8            92        30          1         0.033       0.867
#>     J4            78        21          0         0.000       1.000
```

## Repeated simulation

`sim_replicate` generates datasets with successive seeds. The same
analysis can then be applied to each dataset to estimate bias, coverage,
rejection rates, or power.

``` r

batch <- sim_replicate(
  simulate_rasch, 10,
  n_persons = 400,
  n_items = 8,
  dif = list(items = "I04", uniform = 0.8),
  n_groups = 2,
  seed = 700
)

flagged <- vapply(batch, function(dd) {
  s <- dif_anova(rasch(dd, id = "id", factors = "group"))$summary
  isTRUE(s$uniform_DIF[s$item == "I04"])
}, logical(1))

mean(flagged)
#> [1] 0.5
```

Ten replicates demonstrate the workflow but do not give a stable power
estimate. For a Monte Carlo proportion \\\hat p\\ based on \\R\\
independent replicates, the estimated Monte Carlo standard error is

\\ \operatorname{MCSE}(\hat p)= \sqrt{\frac{\hat p(1-\hat p)}{R}}. \\

The number of attempted, refused, and non-converged fits should be
reported. Bias and coverage should be calculated for each generating
condition rather than after pooling conditions with different true
values.

## Validation studies

The repository contains the simulation studies used to check parameter
recovery, standard errors, confidence-interval coverage, null rejection
rates, power, and identification guards. The scripts and result tables
are under `tools/simval/`; they are excluded from the CRAN source
package for size, and are computationally intensive to re-run. The
examples in this vignette use small runs and are intended as templates
for design-specific studies.

The principal calibration results, each carried with its script and
provenance in the result tables:

| Quantity | Design | Result |
|----|----|----|
| `lr_test` adjusted size | 500 persons, 8 items, 3 categories | 4.7% at the 0.05 level (2,000 replicates) |
| `lr_test` small-sample edge | 300 persons, 12 items, 4 categories | 6.1% among 1,927/2,000 admissible replicates |
| `dependence_magnitude` size | 800 persons, 10 items | 7.5% pooled-variance (pre-fix) to 4.8% covariance-based |
| Class-interval item fit | 8–30 dichotomous items, 600 persons | HC3 was rejected (21.9–48.3% item-wise Type I); conventional ANOVA remained approximate, whereas item-trait Holm familywise error was 4.0–7.0% from ten items onward and 12.0–17.0% with eight items (200 replicates each) |
| Repeated-measures DIF follow-up | 10:90 nuisance-cell imbalance | 5.25% size for a main effect and 5.4% for a mixed interaction (2,000 replicates each); the superseded person-frequency shortcut targeted a different contrast |
| Ordinary DIF covariance | balanced, 1:4 ability imbalance, unequal observations/person, and a three-level 1:2:3 factor | hybrid HC3-uniform/residual-ANOVA-non-uniform familywise error 4.0%, 6.4%, 4.6%, and 6.2%; full HC3 reached 22.0% and 20.2% in the two imbalanced designs and was rejected (500 replicates each) |
| Balanced homoskedastic DIF | two- and three-level factors, 10–150 observations per group-by-interval cell | hybrid familywise error 4.3–5.2%; under local alternatives the classical power advantage declined from 3.10 points at ten per cell to 0.84, 0.42, and 0.16 at 30, 75, and 150 for two levels, and from 1.72 at ten to 0.54 at 50 for three levels (5,000 paired replicates each) |
| Conditional DIF bootstrap | dichotomous data, four-category PCM and RSM data, three-level groups, and correlated person factors | preserving raw scores and refitting under the Rasch null gave acceptable global-null calibration but was usually more conservative and less powerful than the hybrid analysis; under a partial alternative it reduced, but did not remove, artificial flags on invariant items, so it was not adopted as a replacement (100–300 datasets; 99–199 bootstrap refits each) |
| DIF score purification | four-category PCM and RSM data, plus two correlated person factors | preselecting a five-item anchor scale was liberal and leave-one-out matching was rejected; the public split-and-refit procedure retained uniform-DIF power and left 4.0–5.6% familywise error among invariant items, while a strongest-item recalibration was promising for non-uniform DIF but is not yet an automatic remedy (500 replicates per refined condition) |
| BTL-DIF pairwise inference | 6 objects, 8 or 10 judges per factor level | 5.5% and 4.83% size when balanced; 5.0% with 10 raw/9.31 effective judges per level (2,000-replicate top-up); omnibus and pairwise inference are withheld below eight judges or eight effective judges per level |
| BTL-DIF omnibus under unequal precision | 8 versus 16 judges, fourfold variance ratio | classical Type I 11.7% or 1.8% depending on the allocation; HC3 Type I 5.6% and 4.0% (10,000 replicates each) |
| BTL core cluster covariance | 10 balanced judges; 20 judges with one carrying 20% | CR1 Type I 5.4% and 4.2%, coverage 94.6% and 95.8%; delete-one-judge Type I 5.6% and 4.0%, coverage 94.4% and 96.0% (500 replicates each); no default change |
| Superitem spread test | 900 persons, 8 items | 5.2% size at the binomial boundary (1,000 replicates); 100% power under the planted dependence condition (400 replicates) |
| EFRM set-unit linking | 500-600 persons, 8 items/set | hybrid Type I 4.0-5.0%, SE ratios 0.97-1.05 and coverage 0.927-0.960 under normal, bimodal and contrasting group distributions; full-bootstrap Type I 2.5% and coverage 0.975; unit probabilities require 50 persons on every contributing group or link |
| EFRM compiled-kernel parity | demonstration data, 30 seed-paired hybrid bootstrap replicates | largest absolute difference from the retained R implementation 1.30e-11 across set units, standard errors, origins, thresholds and edge likelihoods; all convergence flags agreed |
| EFRM parallel-bootstrap parity | demonstration data, 300 hybrid replicates; simulated data, 30 full-bootstrap replicates | serial, two-worker and four-worker fits used the same pre-generated samples and agreed exactly on every checked estimate and convergence flag |
| BTL-EFRM parallel-bootstrap parity | 3 sets, 2 panels, 20 judges, 200 judge-bootstrap replicates | serial and default four-worker fits used the same pre-generated judge resamples and agreed exactly on the complete reported result, apart from the recorded worker count; elapsed time fell from 17.34 to 5.47 seconds on the executing machine |
| BTL-EFRM unit tests | 12 judges, 6 objects/set | judge-bootstrap Type I 3.3-5.3% for panel units, set units and origins; independent-outcome bootstrap 3.0-6.7% (300 replicates); judge-bootstrap probabilities require six judges and 5.5 effective judges per panel and eight per set link |
| MFRM interaction omnibus | 50 or 200 persons, 6 raters, 25 df | 4.3% and 5.2% Type I error (600 fixed-truth replicates each); probabilities use the least-supported facet level and require `max(30, q + 2)` effective persons |
| MFRM multifactor DIF | 500 persons, 8 items, 6 raters | 4.7% familywise error with balanced raters and 3.8% when one group has two raters (1,000 replicates) |
| Frame-invariance bootstrap size | 500 persons/frame, 8 common items | 3.0% combined Holm familywise error; SE ratios 1.00 locations and 1.03 discrimination (300 replicates) |
| Frame-invariance bootstrap power | two affected items, 500 persons/frame | 96.3% for a one-logit location shift; 9.6% for a 1.5-fold discrimination change (120 replicates) |
| Comparative judgement contrasts | 10-50 judges, balanced | 5.0% size, 94.5% coverage (1,200 replicates) |
| Effective-judge thresholds | one judge with 15-50% of comparisons | ~9% at 4 effective, ~7% at 6-7, nominal when balanced |
| Equating familywise error | 3, 5, and 10 anchors | 4.8-5.0% under the Holm adjustment (2,000 replicates per anchor count) |
| Person-measure coverage | 10-item test, central range | 0.945-0.983; conservative in the tails |
| Tailored bootstrap | 300 persons, 8 items, 399 resamples | clean-item familywise error 0-2.5%; at least one of two hard items detected in 17.5% and 26.3% of datasets with guessing 0.15 and 0.30 (80 full-procedure replicates per effect) |
| CL-AIC model selection | PCM vs RSM; free vs PC thresholds (items and CJ) | null false selection 4.5-5.2% multi-parameter, ~17% one-parameter (the theoretical AIC rates); at the strongest tested departures, selection was 50% for PCM vs RSM and 99.5-100% for the threshold-structure comparisons |
| Paired-comparison effect tests | 8 objects, 14 judges | position/exposure nulls 5.8%/5.9%; carry-over 8.3% at 14 judges, 5.3% at 30; power 62/39/77% at 0.6 logits |
| Cross-package agreement | sirt, eRm, TAM, BradleyTerry2, VGAM, lme4 | identical-likelihood comparators at solver precision; current EFRM set-unit bias +0.0036 vs TAM +0.0008 dichotomous and +0.0035 vs +0.0020 polytomous |
| Cross-package diagnostics | eRm, TAM, psych, difR, PerFit, sirt | alpha exact; item fit r 0.97-0.99 aligned; person fit rho 0.97-0.98; DIF detection 80-88% across methods; dimensionality conservative (exact null, 67% power) vs DETECT (100%) |
| EFRM boundary conditions | 3-8 items/set; 80-1,000 persons; ratios to 3.5; targeting, missingness and non-normality | absolute bias at most 0.022 under the model; all three-item links refused; at 80 persons 11% refused and 2% did not converge; 41- and 101-point grids agreed |
| BTL-EFRM staged link | 6 objects/set, 12 judges | log set-unit bias decreased from -0.108 at 10 repetitions per pair to -0.041 at 20, -0.016 at 50 and -0.007 at 100; bootstrap coverage was 0.933-0.950 at 20 repetitions |
| Explanatory Rasch models | LLTM, LPCM, dichotomous and ordered CJ, with independent or judge-clustered comparisons | coefficient bias at most 0.005 logits, empirical SD/mean SE 0.99-1.03, coverage 0.938-0.955, and Kent-adjusted null rejection 4.2-6.0% (1,000 replicates per condition); fixed-departure Holm familywise error 4.3-4.7% and power 98.3-100% for a 0.8-logit departure (300 replicates per condition) |
| Explanatory edge cases | 300-2,000 persons; dichotomous, four-category and mixed-maximum-score items | empirical SD/mean SE 0.993-1.026, coverage 0.942-0.954, Kent-adjusted null rejection 4.3-5.8%, and no refusals or non-convergence (1,000 replicates per condition); the unscaled probability rejected 98.6-100% and is retained only as `p_naive` |
| Explanatory calibration R-squared | uninformative and true designs, 12 and 24 items | the raw coefficient averaged 0.170 and 0.085 for uninformative designs where the adjusted coefficient centred at -0.015 and -0.002; a true design gave 0.948 raw and 0.936 adjusted (300 replicates per condition) |

## References

Andrich, D. (1978). Relationships between the Thurstone and Rasch
approaches to item scaling. *Applied Psychological Measurement*, 2(3),
451–462.

Andrich, D., and Kreiner, S. (2010). Quantifying response dependence
between two dichotomous items using the Rasch model. *Applied
Psychological Measurement*, 34, 181–192.

Andrich, D., and Marais, I. (2019). *A Course in Rasch Measurement
Theory*. Springer.

Christensen, K. B., Makransky, G., and Horton, M. (2017). Critical
values for Yen’s Q3. *Applied Psychological Measurement*, 41, 178–194.

Humphry, S. M. (2005). *Maintaining a Common Arbitrary Unit in Social
Measurement*. PhD thesis, Murdoch University.

Humphry, S. M., and Andrich, D. (2008). Understanding the unit in the
Rasch model. *Journal of Applied Measurement*, 9(3), 249–264.

Rasch, G. (1960). *Probabilistic Models for Some Intelligence and
Attainment Tests*. Danish Institute for Educational Research. Expanded
edition, University of Chicago Press, 1980.

Tutz, G. (1986). Bradley-Terry-Luce models with an ordered response.
*Journal of Mathematical Psychology*, 30(3), 306–316.

Yen, W. M. (1984). Effects of local item dependence on the fit and
equating performance of the three-parameter logistic model. *Applied
Psychological Measurement*, 8, 125–145.
