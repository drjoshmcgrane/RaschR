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
#>   I01 group     0.164         0.856            
#>   I02 group     3.235         0.242            
#>   I03 group     0.172         0.856            
#>   I04 group     0.412         0.856            
#>   I05 group     3.695         0.242            
#>   I06 group    18.647       < 0.001           *
#>   I07 group     0.611         0.856            
#>   I08 group     0.009         0.926            
#>   I09 group     0.139         0.856            
#>   I10 group     0.085         0.856
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
#>  significant practical ets
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
| EFRM set-unit omnibus size | 200/group, 8 items/set, 2 sets | 9.4% (pre-fix) to 4.9% (1,200 replicates, eight designs) |
| Comparative judgement contrasts | 10-50 judges, balanced | 5.0% size, 94.5% coverage (1,200 replicates) |
| Effective-judge thresholds | one judge with 15-50% of comparisons | ~9% at 4 effective, ~7% at 6-7, nominal when balanced |
| Equating familywise error | 3, 5, and 10 anchors | 4.9-5.5% (2,000-4,000 replicates) |
| Person-measure coverage | 10-item test, central range | 0.945-0.983; conservative in the tails |
| Tailored bootstrap familywise error | 300 persons, 8 items | 0.8% over 240 full-procedure replicates |
| CL-AIC model selection | PCM vs RSM; free vs PC thresholds (items and CJ) | null false selection 4.5-5.2% multi-parameter, ~17% one-parameter (the theoretical AIC rates); detection 95-100% at strong departures |
| Paired-comparison effect tests | 8 objects, 14 judges | position/exposure nulls 5.8%/5.9%; carry-over 8.3% at 14 judges, 5.3% at 30; power 62/39/77% at 0.6 logits |
| Cross-package agreement | sirt, eRm, TAM, BradleyTerry2, VGAM, lme4 | identical-likelihood comparators at solver precision; estimator variants within 0.02-0.15 logits; unit estimators tracked against external anchors |
| Cross-package diagnostics | eRm, TAM, psych, difR, PerFit, sirt | alpha exact; item fit r 0.97-0.99 aligned; person fit rho 0.97-0.98; DIF detection 84-88% all methods; dimensionality conservative (exact null, 67% power) vs DETECT (100%) |
| EFRM unit linking (corrected) | 5-15 items/set; ratios 1-2; skew 0.5-2.8; booklet and pairwise-overlap designs | bias within ±0.015 throughout; null size 2.7-4.8%; coverage 93-99%; stays unbiased under skew where a normal-population anchor drifts |
| EFRM unit linking across sample size | 250 to 10,000 linking persons, 8 items/set (100 replicates per cell) | reported SEs track the empirical SD throughout (ratios 0.91-1.06); bias plateaus at ~0.4% on the ratio; coverage 92-99% to 5,000 persons |

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
