# Many-facet Rasch measurement

``` r

library(rasch)
```

## When to use the MFRM

The many-facet Rasch model (Linacre 1989) extends the Rasch model (Rasch
1960) to responses jointly indexed by a person, an item, and one or more
measurement facets such as rater, task, or occasion. Use `rasch_mfrm`
for such designs. A positive facet parameter denotes greater severity.
Person-group variables such as sex or treatment are not facets: carry
them as person factors and assess them with `dif_anova`.

For facet levels \\f_1,\ldots,f_Q\\, the model is

\\ P(X\_{ni\mathbf{f}}=x)= \frac{\exp\left\\x\theta_n-
\sum\_{k=1}^{x}\left(\delta\_{ik}+\sum\_{q=1}^{Q}\rho\_{qf_q}\right)\right\\}
{\sum\_{y=0}^{m_i}\exp\left\\y\theta_n-
\sum\_{k=1}^{y}\left(\delta\_{ik}+\sum\_{q=1}^{Q}\rho\_{qf_q}\right)\right\\}.
\\

Item thresholds have a common sum-zero origin, and the levels of each
facet sum to zero.

``` r

d <- simulate_mfrm(n_persons = 60, n_items = 4, n_raters = 5,
                   rater_severity_sd = 0.7, seed = 8)
fit <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                  facets = "rater")
fit
#> rasch multiple ratings analysis: 4 items x 5 rater level(s) = 20 response cells, 60 persons
#> Pairwise conditional ML: converged in 7 iterations
#> PSI 0.923, power of fit: excellent
#> 
#> Facet 'rater' severities (logits):
#>  level severity    se   n fit_resid
#>     R1    0.756 0.076 240     0.080
#>     R2    0.585 0.076 240    -0.154
#>     R3   -0.981 0.102 240     0.195
#>     R4   -0.438 0.086 240     0.011
#>     R5    0.077 0.080 240    -0.598
#> (pooled fit residuals and their df on fit$facet_effects)
#> 
#> Notes: item I1:R3: only 0 response(s) in category 0; threshold(s) 1 and the item location are weakly determined (SE reported as NA) -- consider pc_components or collapsing categories; item I1:R3: category 0/1 has only 0/4 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I1:R4: only 2 response(s) in category 0; threshold(s) 1 and the item location are weakly determined (SE reported as NA) -- consider pc_components or collapsing categories; item I1:R4: category 0 has only 2 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I1:R5: category 0 has only 4 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I2:R1: category 3 has only 4 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I2:R2: category 3 has only 7 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I2:R3: category 0 has only 4 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I3:R1: category 3 has only 4 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I3:R2: category 3 has only 6 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I3:R3: category 0 has only 5 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I4:R1: only 0 response(s) in category 3; threshold(s) 3 and the item location are weakly determined (SE reported as NA) -- consider pc_components or collapsing categories; item I4:R1: category 3 has only 0 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I4:R2: only 2 response(s) in category 3; threshold(s) 3 and the item location are weakly determined (SE reported as NA) -- consider pc_components or collapsing categories; item I4:R2: category 2/3 has only 6/2 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; item I4:R5: category 3 has only 4 response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories; a universal raw-score conversion is not defined across the expanded facet response cells; use the design-specific information curves
```

## Read the structural parameters

The model estimates item thresholds and facet severities together.
Internally, each observed item-by-facet combination is a virtual item
whose thresholds are the item thresholds shifted by the relevant facet
effects. Pairwise conditioning removes the person parameter before
calibration.

``` r

fit$item_effects
#>  item location    se   n infit_ms outfit_ms fit_resid fit_resid_pooled  df_fit
#>    I1   -0.988 0.074 300    1.045     1.059     0.164            0.657 281.250
#>    I2   -0.327 0.061 300    1.054     1.021     0.111            0.352 281.250
#>    I3    0.386 0.066 300    0.990     0.989    -0.095            0.000 281.250
#>    I4    0.930 0.090 300    0.923     0.878    -0.553           -1.142 281.250
fit$facet_effects$rater
#>  level severity    se   n infit_ms outfit_ms fit_resid fit_resid_pooled df_fit
#>     R1    0.756 0.076 240    1.075     1.013     0.080            0.243    225
#>     R2    0.585 0.076 240    0.997     0.956    -0.154           -0.299    225
#>     R3   -0.981 0.102 240    0.945     1.075     0.195            0.736    225
#>     R4   -0.438 0.086 240    1.099     1.010     0.011            0.210    225
#>     R5    0.077 0.080 240    0.905     0.877    -0.598           -1.172    225
head(fit$item_thresholds)
#>  item k    tau    se
#>    I1 1 -2.038 0.259
#>    I1 2 -0.993 0.175
#>    I1 3  0.068 0.145
#>    I2 1 -1.448 0.222
#>    I2 2 -0.463 0.171
#>    I2 3  0.928 0.172
```

``` r

plot_facets(fit, facet = "rater")
```

![Rater severity estimates with confidence
intervals.](many-facet_files/figure-html/facets-1.png)

The design must connect facet levels through common items and persons. A
facet nested within an item or a person-disjoint block can be confounded
with item location. `rasch_mfrm` checks the structural rank and
informative co-observation graph and stops when the decomposition is not
identified.

## Item-by-facet interaction

The additive model assumes that severity differences are invariant
across items. When the design and substantive question require it,
`interaction` adds item-by-level terms with double sum-to-zero
constraints.

``` r

fit_interaction <- rasch_mfrm(
  d, person = "person", item = "item", score = "score",
  facets = "rater", interaction = "rater"
)
head(fit_interaction$interaction_effects)
#>  item level  gamma    se      z     p p_adj significant
#>    I1    R1 -0.157 0.179 -0.878 0.383 1.000            
#>    I2    R1  0.041 0.132  0.311 0.757 1.000            
#>    I3    R1  0.059 0.156  0.379 0.706 1.000            
#>    I4    R1  0.057 0.144  0.394 0.695 1.000            
#>    I1    R2 -0.071 0.145 -0.490 0.626 1.000            
#>    I2    R2 -0.086 0.154 -0.555 0.581 1.000
fit_interaction$interaction_test
#>  facet df  wald     f df2     p min_effective_persons minimum_required
#>  rater 12 5.993 0.406  48 0.954                    60               30
#>  inference_available
#>                    *
```

The interaction model retains equal discrimination, but comparisons
among raters become item-dependent. A material interaction therefore
qualifies the claim of invariant rater severity rather than merely
improving fit. The joint Wald test is the primary test of the
interaction family; individual cells are exploratory and adjusted by
Holm’s method. The test uses the least-supported facet level and
withholds probabilities when that level has fewer than `max(30, q + 2)`
persons or effective persons, where `q` is the omnibus degrees of
freedom.

## Diagnostics

The returned object is also a `rasch` fit in which each item-by-facet
cell enters as its own column of the response matrix (a *virtual item*),
so the fit, targeting, dependence, dimensionality, and plotting
functions remain available. MFRM margin tables additionally report an
equal-cell fit residual and a response-weighted pooled residual. Their
weighting differs, so both the design and the location of any misfit
should guide interpretation.

## References

Linacre, J. M. (1989). *Many-Facet Rasch Measurement*. Chicago: MESA
Press.

Rasch, G. (1960). *Probabilistic Models for Some Intelligence and
Attainment Tests*. Copenhagen: Danish Institute for Educational
Research. (Expanded edition, 1980, Chicago: University of Chicago
Press.)
