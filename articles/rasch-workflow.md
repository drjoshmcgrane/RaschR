# A Rasch analysis workflow

``` r

library(rasch)
```

This vignette follows a Rasch analysis from the overall summary to item
and person fit, targeting, dependence, and differential item functioning
(DIF). The order matters. A fit statistic is difficult to interpret
without knowing whether the scale separates the sample or supplies
information over the relevant part of the latent trait.

The same analyses are available in the Shiny application. Run
[`rasch::run_app()`](https://drjoshmcgrane.github.io/rasch/reference/run_app.md)
to import data, fit a model, inspect the results, and export tables or
plots without writing R analysis code.

## Fit the model required by the scoring structure

`rasch` fits the partial credit model by default. Dichotomous data are
its one-threshold special case. Set `model = "RSM"` only when the same
category threshold structure is intended to hold across items.

For person \\n\\, item \\i\\, and score \\x=0,\ldots,m_i\\, the partial
credit model is

\\ P(X\_{ni}=x)=
\frac{\exp\left\\x\theta_n-\sum\_{k=1}^{x}\delta\_{ik}\right\\}
{\sum\_{y=0}^{m_i}\exp\left\\y\theta_n-
\sum\_{k=1}^{y}\delta\_{ik}\right\\}. \\

The rating scale model constrains \\\delta\_{ik}=\beta_i+\tau_k\\.
Pairwise conditioning removes \\\theta_n\\ from the item likelihood.
Person locations are subsequently estimated by Warm’s weighted
likelihood method (Warm 1989).

The example is polytomous and contains three person groups. One item has
disordered generating thresholds, one has uniform DIF, and one item pair
has local response dependence. These departures make the diagnostic
sequence visible without changing the commands used for observed data.
Before fitting observed data, check the coding and frequency of every
category. Negative scores are read as missing; valid categories begin at
zero.

``` r

d <- simulate_rasch(
  n_persons = 600,
  n_items = 12,
  model = "PCM",
  n_categories = 4,
  difficulty = c(-1.5, 1.5),
  disordered = "I04",
  dependence = list(pairs = list(c("I10", "I11")), strength = 1.3),
  dif = list(items = "I08", uniform = 0.8),
  n_groups = 3,
  seed = 17
)

fit <- rasch(d, model = "PCM", id = "id", factors = "group")
```

## Begin with the summary and reliability

The summary establishes whether estimation converged and gives the
overall item–trait interaction, fit-residual distributions, reliability,
and the package’s qualitative assessment of the power of fit tests. It
should be read before individual item or person flags.

``` r

fit
#> rasch PCM analysis: 12 items, 600 persons
#> Pairwise conditional ML (Andrich & Luo): converged in 5 iterations
#> PSI 0.851 (no extremes 0.851), item SI 0.995, alpha 0.855, power of fit: good
#> Total item-trait chi-square 101.760 on 108 df, p = 0.651
fit_summary_table(fit)
#>                                      statistic                   value
#> 1                                        Model                     PCM
#> 2                                   Estimation pairwise conditional ML
#> 3                                    Converged                     yes
#> 4                                   Iterations                       5
#> 5                  Total item-trait chi-square                 101.760
#> 6                           Degrees of freedom                     108
#> 7                       Item-trait probability                   0.651
#> 8                              Class intervals                      10
#> 9                       Item fit residual mean                    0.15
#> 10                        Item fit residual SD                    0.87
#> 11                  Item fit residual skewness                    0.62
#> 12                  Item fit residual kurtosis                   -0.79
#> 13                    Person fit residual mean                   -0.17
#> 14                      Person fit residual SD                    0.90
#> 15                Person fit residual skewness                   -0.19
#> 16                Person fit residual kurtosis                    0.41
#> 17            Fit-location correlation (items)                   0.039
#> 18          Fit-location correlation (persons)                   0.035
#> 19 Items with Holm-adjusted chi-square p < .05                 0 of 12
#> 20                       Disordered thresholds                     I04
```

The person separation index (PSI) compares the observed variance of
person locations with their mean error variance:

\\ \mathrm{PSI}= \frac{\operatorname{Var}(\hat\theta)-
\operatorname{mean}\\\operatorname{SE}(\hat\theta)^2\\}
{\operatorname{Var}(\hat\theta)}. \\

The implementation truncates negative values at zero. `fit$psi_noext`
removes persons with extreme scores and is useful when extremes inflate
the observed spread. The separation ratio and number of strata express
the same information on more interpretable scales.

``` r

fit$psi
#> $PSI
#> [1] 0.8511
#> 
#> $separation
#> [1] 2.391
#> 
#> $strata
#> [1] 3.522
#> 
#> $var_theta
#> [1] 1.136
#> 
#> $mean_error_var
#> [1] 0.1692
#> 
#> $n
#> [1] 600
fit$psi_noext
#> $PSI
#> [1] 0.8511
#> 
#> $separation
#> [1] 2.391
#> 
#> $strata
#> [1] 3.522
#> 
#> $var_theta
#> [1] 1.136
#> 
#> $mean_error_var
#> [1] 0.1692
#> 
#> $n
#> [1] 600
fit$power_of_fit
#> [1] "good"
```

The power label is a screening judgement based on PSI, not a formal
power calculation. Low reliability weakens the ordering of persons and
the formation of distinct class intervals, so item–trait and DIF tests
may fail to detect real departures. Power also depends on sample size,
test length, targeting, category use, and trait spread. Conversely, a
large sample can make a small departure statistically significant. Fit
residuals, effect sizes, plots, and substantive importance remain
necessary.

## Examine item estimates and fit

Item locations and their standard errors should be read alongside fit
residuals, mean-square statistics, and the Holm-adjusted item–trait
probabilities. A positive residual indicates more variation than
expected; a negative residual indicates responses that are more
predictable than expected. The conventional \\\pm2.5\\ residual band is
a screening rule rather than a separate hypothesis test.

``` r

fit$items[, c(
  "item", "location", "se", "fit_resid",
  "infit_ms", "outfit_ms", "p_adj"
)]
#>  item location    se fit_resid infit_ms outfit_ms p_adj
#>   I01   -1.604 0.106    -1.103    0.949     0.899 1.000
#>   I02   -1.246 0.079     1.548    1.048     1.119 1.000
#>   I03   -0.887 0.063    -0.158    1.022     0.990 1.000
#>   I04   -0.719 0.072    -0.335    0.980     0.947 1.000
#>   I05   -0.509 0.059    -0.685    0.973     0.960 1.000
#>   I06   -0.143 0.059     1.394    1.092     1.087 1.000
#>   I07    0.077 0.057     1.501    1.088     1.108 1.000
#>   I08    0.567 0.056    -0.221    0.995     0.995 1.000
#>   I09    0.692 0.061    -0.243    1.012     0.988 1.000
#>   I10    1.023 0.077     0.068    0.993     0.997 1.000
#>   I11    1.244 0.083     0.037    1.018     0.991 1.000
#>   I12    1.505 0.094     0.048    1.021     0.985 1.000
```

``` r

plot_item_map(fit)
```

![Item locations plotted against item fit
residuals.](rasch-workflow_files/figure-html/item-fit-plot-1.png)

For a polytomous item, successive thresholds should normally increase on
the latent scale. Threshold disordering means that an intervening
category is not the most probable response over any part of the scale.
It can reflect sparse categories, indistinct category meanings, or a
scoring order that is not working as intended.

``` r

fit$thresholds
#>  id item k    tau    se anchored weak
#>   1    1 1 -3.184 0.332              
#>   2    1 2 -1.247 0.137              
#>   3    1 3 -0.380 0.100              
#>   4    2 1 -2.422 0.238              
#>   5    2 2 -1.071 0.138              
#>   6    2 3 -0.244 0.104              
#>   7    3 1 -1.626 0.189              
#>   8    3 2 -1.122 0.128              
#>   9    3 3  0.088 0.105              
#>  10    4 1 -0.990 0.245              
#>  11    4 2 -2.762 0.170              
#>  12    4 3  1.594 0.111              
#>  13    5 1 -1.666 0.163              
#>  14    5 2 -0.665 0.113              
#>  15    5 3  0.805 0.113              
#>  16    6 1 -1.549 0.142              
#>  17    6 2 -0.131 0.109              
#>  18    6 3  1.251 0.132              
#>  19    7 1 -0.837 0.129              
#>  20    7 2 -0.284 0.114              
#>  21    7 3  1.351 0.133              
#>  22    8 1 -0.309 0.112              
#>  23    8 2  0.424 0.120              
#>  24    8 3  1.585 0.152              
#>  25    9 1 -0.335 0.108              
#>  26    9 2  0.522 0.119              
#>  27    9 3  1.890 0.172              
#>  28   10 1 -0.312 0.102              
#>  29   10 2  1.200 0.124              
#>  30   10 3  2.182 0.228              
#>  31   11 1 -0.056 0.102              
#>  32   11 2  1.223 0.132              
#>  33   11 3  2.564 0.257              
#>  34   12 1  0.098 0.102              
#>  35   12 2  1.525 0.137              
#>  36   12 3  2.893 0.295
ordered <- vapply(fit$thresholds_diag, function(x) x$ordered, logical(1))
disordered_items <- names(ordered)[!ordered]
disordered_items
#> [1] "I04"
fit$thresholds_diag[disordered_items]
#> $I04
#> $I04$item
#> [1] "I04"
#> 
#> $I04$thresholds
#> [1] -0.9898 -2.7619  1.5941
#> 
#> $I04$ordered
#> [1] FALSE
#> 
#> $I04$reversed_at
#> [1] 2
#> 
#> $I04$never_modal_categories
#> [1] 1
#> 
#> $I04$category_counts
#> [1]  41  39 402 118
```

``` r

plot_ccc(fit, "I04", observed = TRUE)
```

![Category characteristic curves and observed category proportions for
item I04.](rasch-workflow_files/figure-html/category-curves-1.png)

Do not collapse categories from the ordering flag alone. Check category
frequencies, the category curves, the wording of the response options,
and the consequences for measurement. If categories are combined on
substantive grounds, refit the model and repeat the full diagnostic
sequence.

## Examine person estimates and fit

Person fit addresses the consistency of each response pattern with the
fitted scale. The following table puts the largest absolute residuals
first. Positive residuals indicate unexpectedly erratic patterns;
negative residuals indicate patterns that are unusually predictable.

``` r

person_order <- order(abs(fit$person$fit_resid), decreasing = TRUE,
                      na.last = TRUE)
head(fit$person[person_order, c(
  "id", "group", "raw", "theta", "se", "extreme",
  "fit_resid", "infit_ms", "outfit_ms"
)], 10)
#>     id group raw  theta    se extreme fit_resid infit_ms outfit_ms
#>  P0486    g3  18 -0.042 0.381            -3.454    0.202     0.212
#>  P0238    g1  23  0.691 0.395            -3.241    0.190     0.226
#>  P0001    g1  18 -0.042 0.381            -3.163    0.233     0.242
#>  P0031    g1  11 -1.062 0.397             2.999    3.209     3.955
#>  P0493    g1  21  0.391 0.387            -2.912    0.249     0.269
#>  P0327    g3  22  0.540 0.390            -2.742    0.235     0.288
#>  P0484    g1  16 -0.327 0.381            -2.610    0.307     0.311
#>  P0417    g3  15 -0.470 0.382            -2.430    0.320     0.338
#>  P0416    g2  19  0.101 0.382            -2.289    0.346     0.357
#>  P0230    g2  23  0.691 0.395            -2.288    0.343     0.349
```

``` r

plot_person_fit(fit)
```

![Person locations plotted against person fit
residuals.](rasch-workflow_files/figure-html/person-fit-plot-1.png)

An unexpected response pattern may reflect coding or data-entry errors,
careless responding, a secondary trait, or a genuine but unusual person.
It is not, by itself, a reason to remove the person. Fit residuals are
unavailable for extreme response patterns because those patterns do not
provide an interior location at which fit can be assessed.

## Examine targeting and information

Targeting concerns the match between the person distribution and the
item threshold distribution. The table reports their locations and
spread, the proportions of persons beyond the threshold range, and the
principal reliability indices.

``` r

targeting_table(fit)
#>                             statistic  value
#> 1                Person location mean -0.012
#> 2                  Person location SD  1.066
#> 3            Person location skewness  -0.03
#> 4            Person location kurtosis   0.26
#> 5  Person location mean (no extremes) -0.012
#> 6    Person location SD (no extremes)  1.066
#> 7    Item location mean (constrained) -0.000
#> 8                    Item location SD  1.017
#> 9              Item location skewness  -0.03
#> 10             Item location kurtosis  -1.24
#> 11                  Threshold minimum -3.184
#> 12                  Threshold maximum  2.893
#> 13  Persons below threshold range (%)    0.3
#> 14  Persons above threshold range (%)    0.5
#> 15                                PSI  0.851
#> 16                         Separation   2.39
#> 17                      Person strata    3.5
#> 18               PSI without extremes  0.851
#> 19                 n without extremes    600
#> 20        Item separation reliability  0.995
#> 21                  Coefficient alpha  0.855
#> 22                 n complete (alpha)    600
```

For a Rasch model, test information is the sum of the conditional
response variances:

\\ I(\theta)=\sum_i \operatorname{Var}(X_i\mid\theta), \qquad
\operatorname{SE}(\hat\theta)\approx I(\theta)^{-1/2}. \\

The person–item map below places person locations and item thresholds on
the same logit scale. The information curve shows where the test is most
precise.

``` r

plot_pimap(fit, information = TRUE)
```

![Person and item distributions with the test information curve on the
common logit
scale.](rasch-workflow_files/figure-html/targeting-map-1.png)

``` r

head(test_information(fit))
#>   theta   info   sem
#> 1  -6.0 0.1479 2.600
#> 2  -5.9 0.1630 2.477
#> 3  -5.8 0.1795 2.360
#> 4  -5.7 0.1977 2.249
#> 5  -5.6 0.2176 2.144
#> 6  -5.5 0.2395 2.044
```

The optional `WrightMap` package supplies a second form of Wright map.
The default retains one person panel and displays the estimated category
thresholds. Separate person or item panels should be requested only when
they answer a substantive question.

``` r

if (requireNamespace("WrightMap", quietly = TRUE)) {
  wright_map(fit)
}
```

![Wright map of person locations and item
thresholds.](rasch-workflow_files/figure-html/wrightmap-1.png)

If `WrightMap` is not installed, install it with
`install.packages("WrightMap")` before running this chunk.

## Examine local and trait dependence

Local response dependence occurs when two responses remain associated
after conditioning on the latent trait. Yen’s \\Q3\\ is the correlation
between two items’ standardised residuals. Because raw \\Q3\\ values
have a negative baseline in a finite test, `q3_star` subtracts the
average off-diagonal value.

``` r

q3 <- residual_correlations(fit)
q3$average
#> [1] -0.08639
head(q3$pairs, 10)
#>    item_a item_b         q3 q3_star flagged
#> 1     I10    I11  0.1266703 0.21306      NA
#> 2     I02    I05  0.0045201 0.09091      NA
#> 3     I01    I04 -0.0009814 0.08541      NA
#> 4     I06    I07 -0.0173825 0.06901      NA
#> 5     I02    I03 -0.0244845 0.06190      NA
#> 6     I05    I09 -0.0252834 0.06110      NA
#> 7     I01    I09 -0.0261819 0.06021      NA
#> 8     I02    I12 -0.0286664 0.05772      NA
#> 9     I04    I12 -0.0336436 0.05274      NA
#> 10    I04    I10 -0.0364637 0.04992      NA
```

``` r

plot_resid_cor(fit)
```

![Heatmap of adjusted residual correlations between
items.](rasch-workflow_files/figure-html/local-dependence-plot-1.png)

There is no universal critical value for adjusted \\Q3\\ (Christensen,
Makransky and Horton 2017). Its size, the response process, and the
content of the item pair matter. When theory supports treating dependent
items as one superitem,
[`combine_items()`](https://drjoshmcgrane.github.io/rasch/reference/combine_items.md)
refits the complete calibration and
[`spread_test()`](https://drjoshmcgrane.github.io/rasch/reference/spread_test.md)
compares its threshold spread with the binomial bound. The bound does
not apply to a superitem containing a polytomous component.

Trait dependence is examined through the residual components and by
comparing person estimates from opposed item subsets (Smith 2002). For
subsets \\A\\ and \\B\\, the person-level statistic is

\\ t_n=\frac{\hat\theta\_{nA}-\hat\theta\_{nB}}
{\sqrt{\operatorname{SE}(\hat\theta\_{nA})^2+
\operatorname{SE}(\hat\theta\_{nB})^2}}. \\

Under unidimensionality, about the nominated alpha level of these
comparisons should be significant. The exact binomial interval, the
residual loadings, and the score points in each subset are part of the
result.

``` r

dimensionality_test(fit)
#> $prop_significant
#> [1] 0.04498
#> 
#> $ci
#> [1] 0.02959 0.06522
#> 
#> $n
#> [1] 578
#> 
#> $n_excluded_extreme
#> [1] 22
#> 
#> $multidimensional
#> [1] FALSE
#> 
#> $split
#> [1] "residual component 1"
#> 
#> $score_points
#> positive negative 
#>       24       12 
#> 
#> $caution
#> [1] "subtests carry only 24 and 12 score points (fewer than the ~15 recommended for stable subtest estimates); read the verdict cautiously"
#> 
#> $items_positive
#> [1] "I01" "I02" "I03" "I05" "I06" "I07" "I09" "I12"
#> 
#> $items_negative
#> [1] "I04" "I08" "I10" "I11"
#> 
#> $first_eigenvalue
#> [1] 1.351
#> 
#> $paired_t
#> $paired_t$mean_difference
#> [1] -0.03322
#> 
#> $paired_t$t
#> [1] -0.9055
#> 
#> $paired_t$df
#> [1] 577
#> 
#> $paired_t$p
#> [1] 0.3656
```

``` r

plot_pca(fit)
```

![Loadings of items on the first residual
component.](rasch-workflow_files/figure-html/trait-dependence-plot-1.png)

A quiet result from short subtests is inconclusive rather than evidence
that a secondary trait is absent. The returned caution records when
either subset has fewer score points than recommended for stable person
estimates.

## Examine differential item functioning

Person factors must be nominated when the model is fitted. With a person
factor \\G\\ and trait class interval \\C\\, the residual model is

\\ z=\mu+G+C+G\mathbin{:}C+\varepsilon. \\

The factor term tests uniform DIF; the interaction with class interval
tests non-uniform DIF.
[`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
fits all nominated factors together, recognises within-person factors in
repeated designs, and applies Holm familywise correction over the
item-by-term tests.

``` r

dif <- dif_anova(fit, sizes = TRUE)
subset(dif$summary, uniform_DIF | nonuniform_DIF)
#>  item  term F_uniform p_uniform p_uniform_adj eta2_uniform uniform_DIF
#>   I08 group    20.043   < 0.001       < 0.001        0.067           *
#>  F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#>         1.578        0.109            1.000           0.026               
#>  superseded
#> 
```

``` r

plot_icc(fit, "I08", group = "group")
```

![Observed and expected item characteristic curves for item I08 by
person group.](rasch-workflow_files/figure-html/dif-plot-1.png)

For a significant factor with more than two levels, the follow-up should
estimate the relevant differences in Rasch logits rather than rely on a
generic Tukey procedure. With `sizes = TRUE`,
[`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
returns the Holm-adjusted marginal pairwise comparisons for significant
main effects and difference-in-differences for significant interactions.

``` r

dif$posthoc
#>  item  term item.1 contrast within estimate    se statistic df       p   p_adj
#>   I08 group    I08  g2 - g1          -0.039 0.139    -0.284      0.777   0.777
#>   I08 group    I08  g3 - g1           0.705 0.148     4.768    < 0.001 < 0.001
#>   I08 group    I08  g3 - g2           0.745 0.147     5.081    < 0.001 < 0.001
#>   lower upper significant practical
#>  -0.312 0.233                      
#>   0.415 0.995           *         *
#>   0.457 1.032           *         *
```

Statistical significance and practical magnitude answer different
questions. Any split must be supported by the response process and
should be applied with
[`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md),
which refits the calibration and updates the item and person estimates.
The revised fit then goes through the same summary, fit, targeting,
dependence, and DIF sequence.

## Link calibrations through common items

Two separately fitted tests can be aligned through their common items.
Drift tests are meaningful only when the calibrations use independent
sampling units; state this rather than allowing the software to assume
it.

``` r

eq <- equate_tests(current_fit, reference_fit, independent = TRUE)
eq$table
```

A saved item bank needs more than marginal standard errors for drift
inference. Its joint item-location covariance carries the dependence
created by the fitted origin and is attached as
`attr(bank, "cov_location")`, in bank row order or named by item.
Without that matrix, the package still estimates the origin shift and
reports the link descriptively. A bank intended to be fixed may instead
use zero standard errors. For polytomous items, include the maximum
score in a `max` column so that scoring-scale compatibility can be
checked.

## Report the final calibration

A report should state the model, scoring and missing-data rules, sample
and item exclusions, reliability, item and person fit, threshold
functioning, targeting, dependence, dimensionality, DIF, and the
uncertainty of person estimates. The score table supplies the
complete-response conversion from raw score to logits.

``` r

score_table(fit)
#>  score  theta    se extrapolated freq cum_pct
#>      0 -4.844 1.501                 0   0.000
#>      1 -3.637 0.874                 2   0.333
#>      2 -3.058 0.683                 1   0.500
#>      3 -2.674 0.585                 5   1.333
#>      4 -2.381 0.526                 5   2.167
#>      5 -2.137 0.487                 3   2.667
#>      6 -1.925 0.459                 6   3.667
#>      7 -1.731 0.439                 9   5.167
#>      8 -1.552 0.424                17   8.000
#>      9 -1.382 0.413                15  10.500
#>     10 -1.219 0.404                18  13.500
#>     11 -1.062 0.397                35  19.333
#>     12 -0.910 0.392                20  22.667
#>     13 -0.761 0.387                22  26.333
#>     14 -0.614 0.384                22  30.000
#>     15 -0.470 0.382                29  34.833
#>     16 -0.327 0.381                41  41.667
#>     17 -0.185 0.380                28  46.333
#>     18 -0.042 0.381                33  51.833
#>     19  0.101 0.382                26  56.167
#>     20  0.245 0.384                29  61.000
#>     21  0.391 0.387                32  66.333
#>     22  0.540 0.390                40  73.000
#>     23  0.691 0.395                29  77.833
#>     24  0.848 0.402                29  82.667
#>     25  1.009 0.409                15  85.167
#>     26  1.177 0.418                21  88.667
#>     27  1.353 0.429                20  92.000
#>     28  1.539 0.443                11  93.833
#>     29  1.738 0.461                16  96.500
#>     30  1.955 0.483                 5  97.333
#>     31  2.194 0.512                 5  98.167
#>     32  2.467 0.552                 7  99.333
#>     33  2.790 0.609                 1  99.500
#>     34  3.198 0.702                 3 100.000
#>     35  3.780 0.882                 0 100.000
#>     36  4.951 1.485                 0 100.000
```

## References

Andrich, D. (1978). Relationships between the Thurstone and Rasch
approaches to item scaling. *Applied Psychological Measurement*, 2(3),
451–462.

Andrich, D., and Marais, I. (2019). *A Course in Rasch Measurement
Theory: Measuring in the Educational, Social and Health Sciences*.
Springer.

Christensen, K. B., Makransky, G., and Horton, M. (2017). Critical
values for Yen’s Q3: Identification of local dependence in the Rasch
model using residual correlations. *Applied Psychological Measurement*,
41(3), 178–194.

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics*, 6(2), 65–70.

Humphry, S. M. (2005). *Maintaining a Common Arbitrary Unit in Social
Measurement*. PhD thesis, Murdoch University.

Humphry, S. M., and Andrich, D. (2008). Understanding the unit in the
Rasch model. *Journal of Applied Measurement*, 9(3), 249–264.

Rasch, G. (1960). *Probabilistic Models for Some Intelligence and
Attainment Tests*. Copenhagen: Danish Institute for Educational
Research. (Expanded edition, 1980, Chicago: University of Chicago
Press.)

Smith, E. V. Jr. (2002). Detecting and evaluating the impact of
multidimensionality using item fit statistics and principal component
analysis of residuals. *Journal of Applied Measurement*, 3(2), 205–231.

Tutz, G. (1986). Bradley-Terry-Luce models with an ordered response.
*Journal of Mathematical Psychology*, 30(3), 306–316.

Torres Irribarra, D., and Freund, R. (2025). *WrightMap: IRT item-person
map with ConQuest integration*. R package version 1.5.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika*, 54(3), 427–450.

Yen, W. M. (1984). Effects of local item dependence on the fit and
equating performance of the three-parameter logistic model. *Applied
Psychological Measurement*, 8(2), 125–145.
