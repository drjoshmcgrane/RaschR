# Changelog

## rasch 1.12.0

This update extends the model suite, adds mixed-design DIF analysis, and
strengthens identification and uncertainty checks.

### Differential item functioning

- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  fits several person factors jointly. It supports additive or factorial
  effects and uses Type II sums of squares.
- Repeated measurements are analysed with the person as the sampling
  unit. Between-person and within-person terms use their respective
  error strata, with a Greenhouse–Geisser correction for multilevel
  within-person factors.
- [`dif_contrasts()`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md)
  provides planned logit contrasts, including between-by-within
  interactions.
  [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  separates descriptive logit magnitude from repeated-measures
  inference.
- MFRM residuals can be pooled to their underlying items for DIF
  analysis. Person factors not used to define an EFRM frame can also be
  tested.
- [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  splits items iteratively, beginning with the largest confirmed effect,
  while retaining a minimum anchor set.

### Facet, frame, and paired-comparison models

- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  accepts several facets and optional item-by-facet interaction. The
  interaction family has an omnibus test with adjusted cell follow-ups.
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  accepts crossed person-group factors and returns their GLS factorial
  decomposition. Hybrid and full-bootstrap standard errors propagate
  uncertainty in the estimated units.
- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  supports ordered paired comparisons, judge-clustered covariance,
  position effects, and within-judge exposure and carry-over effects.
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  tests object invariance over one or more judge factors using a
  judge-level split-plot analysis.
- [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  fits linked object sets and judge panels with set and panel units.
  Judge bootstrap inference refits both estimation stages.
- Paired-comparison functions now cover common-object equating,
  transitivity, residual bimensions, design information, and adaptive
  pair selection.

### Statistical corrections and identification

- [`dependence_magnitude()`](https://drjoshmcgrane.github.io/rasch/reference/dependence_magnitude.md)
  uses the full covariance of the resolved thresholds when calculating
  its standard error.
- The EFRM item-set unit linking recovers the true-score variance by a
  truncated-score-moment correction. The previous construction (observed
  variance minus mean squared standard error) under-recovers the
  true-score variance on short tests and biased recovered unit ratios
  upward by about five per cent at eight dichotomous items per set; the
  corrected estimator is unbiased there, confirmed against an external
  TAM slope-group anchor.
- The hybrid EFRM set-unit covariance includes uncertainty from the
  within-frame calibration.
- Judge-clustered BTL inference requires enough nominal and effective
  judges, residual cluster degrees of freedom, and a sufficiently
  balanced workload.
- BTL-EFRM unit tests use judge-limited F and t reference distributions.
- Standard errors are withheld for an item’s thresholds when any
  response category is critically sparse.
- Item, facet, frame, and paired-comparison models check connectedness,
  rank, and separation before reporting estimates. Structurally
  unidentified parameters are refused; identified but imprecise
  parameters are marked.
- Equating tests require independent calibrations and the joint
  covariance of banked locations. Links without the required covariance
  remain descriptive.
- Information curves are calculated for item sets and facet designs that
  can actually be administered together.
- Classical statistics use complete responders by default.
  Available-case results are labelled exploratory.

### Follow-ups, comparisons, projects, and reports

- [`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  tests the item invariance a frame model assumes rather than imposing
  it. The fitted model holds one location per item, shared across frames
  and scaled by the frame unit, so the assumption cannot be checked from
  the fit; the function calibrates each frame separately, puts the
  locations on the common scale, and compares them item by item,
  reporting the root mean squared difference against the root mean
  squared standard error. Its `adjust` argument chooses between Holm
  across all comparisons and no adjustment, because screening for items
  to examine and reporting a difference are different jobs: on eight to
  ten items the adjustment costs between 20 and 60 points of
  sensitivity, and simulation shows that carrying through to the
  repaired unit ratio. Both probabilities are reported either way, and
  the printed output names the rule it applied. The application exposes
  the choice as a switch.
- Case study: `inst/casestudies/wording_units_height.R` applies the
  item-set units to a balanced inventory of 26 items, 13 worded in each
  direction, and reads it against a criterion collected outside the
  inventory. It is the counterpart to the self-esteem study rather than
  a repeat of it: there a significant set-level difference turns out to
  be one ambivalent item, while here the difference survives every
  single-item removal and still leaves the person ordering unchanged. A
  unit difference can be decisive and inconsequential at the same time,
  and only the criterion says which.
- The item summary reports `disc`, the slope that maximises an item’s
  own likelihood with the person measures and its thresholds held at the
  values the model gave them, following the index Winsteps reports as
  DISCRIM and generalised here to polytomous items. It is a description
  of how steeply an item sorts the people the model has already located,
  not an estimate of a discrimination parameter: it runs high, because
  the measures are estimated from a set that includes the item being
  scored, and its ordering is more dependable than its level.
  [`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  reports it per frame alongside the comparison it tests, which stays on
  the fit statistics.
- Results tables print in decimals rather than exponents. Tables read
  directly off a fitted object carried no print method, so base R
  formatted them and a probability below 1e-4 arrived as `4.00e-83`
  where the package’s own print methods read `< 0.001`. They now share
  one formatting vocabulary, with the values untouched. Saved tables and
  the application’s downloads keep full precision without the exponent.
- [`drop_items()`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
  removes items from a fitted analysis and refits it, keeping the model,
  person identifiers, factors, and – for frame models – the set
  structure and standard-error method. The application offers the same
  action on the selected item, recorded as an undoable step. For frame
  models this is a sensitivity analysis rather than housekeeping: a set
  unit is estimated from the dispersion its own items produce, so an
  item that fits its set badly moves the unit that decides whether the
  sets differ.
- [`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
  provides post-hoc pairwise DIF comparisons after a significant omnibus
  term: resolved item-location differences in logits with the joint
  contrast covariance, Holm-adjusted, with trend and adjacent-level
  contrasts for ordered factors and simple-effect and
  difference-of-differences contrasts for interactions.
- The application presents automatic model comparisons beside each
  analysis (partial credit against rating scale, free against
  principal-component thresholds, additive against interactive
  many-facet, frame models against their equal-unit restrictions, and
  paired-comparison effect terms), separating formal tests,
  composite-likelihood information criteria, and descriptive fit
  changes. Structure-altering procedures show before-and-after summaries
  instead.
- Analyses can be saved as `.rasch` project files and reopened exactly;
  [`report_document()`](https://drjoshmcgrane.github.io/rasch/reference/report_document.md)
  renders self-contained HTML and editable Word reports for every model
  family.

### Simulation, documentation, and interface

- New simulation functions cover each model family and retain the
  generating parameters for recovery studies.
  [`sim_replicate()`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md),
  [`sim_apply()`](https://drjoshmcgrane.github.io/rasch/reference/sim_apply.md),
  and
  [`sim_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/sim_recovery.md)
  support repeated simulation.
- Six vignettes cover the main Rasch workflow, mixed-design DIF, MFRM,
  EFRM, paired comparisons, and simulation studies.
- The Shiny application includes the extended model suite, model
  comparison, downloadable tables and plots, and the R call
  corresponding to each fit.
- The package and function documentation now state the fitted models,
  identification constraints, uncertainty methods, and principal
  references more directly.

## rasch 1.11.7

CRAN release: 2026-07-30

- [`print()`](https://rdrr.io/r/base/print.html) preserves the reference
  distribution used by saved BTL fits: current and transitional results
  are labelled `t`, while older results without cluster degrees of
  freedom retain their original `z` label.

## rasch 1.11.6

- BTL print methods read both current and earlier dependence-statistic
  columns. Current clustered statistics are labelled `t`, in accordance
  with their t reference distribution.

## rasch 1.11.5

- Ordered paired-comparison margins must be numeric or an ordered
  factor; unordered categorical margins are rejected.
- Invalid graded responses and malformed secondary-dimension
  specifications now produce errors rather than being coerced or
  dropped.
- Clustered dependence and position statistics are labelled `t`.

## rasch 1.11.4

- Score validation reads factor labels and rejects non-integer,
  non-numeric, or non-finite values.
- Ordered BTL responses require an ordered factor or integer scores.
- Clustered dependence and position tests use a t reference with `G - 1`
  degrees of freedom.
- [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md)
  validates secondary-trait correlations and item sets.

## rasch 1.11.3

- All estimators share the same integer-score validation.
- [`item_moments()`](https://drjoshmcgrane.github.io/rasch/reference/item_moments.md)
  uses a log-sum-exp calculation for wide category ranges.
- BTL object-separation reliability is withheld when the clustered
  covariance is rank deficient.
- [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  requires each judge to belong to one panel.
- Undefined reliability and omnibus statistics are reported as `NA`.
- Secondary-trait simulation retains its requested mean and standard
  deviation, and MLE score-table calculations use the correct
  common-unit score equation.

## rasch 1.11.2

- The BTL dimensionality reference simulates count-weighted data at the
  unordered-pair level and returns the leading strength from each
  replicate.
- Judge-bootstrap EFRM inference requires more than one judge per panel
  and notes panels with fewer than five judges.
- BTL fits report when the number of judge clusters cannot support a
  full-rank clustered covariance.

## rasch 1.11.1

- BTL dimensionality reference draws now reproduce count-weighted
  binomial or multinomial sampling.
- Pairwise chi-square degrees of freedom include position and dependence
  parameters; untestable designs return `NA`.
- Judge-clustered covariance uses the CR1 small-sample factor.
- `btl_efrm(se_method = "judge_bootstrap")` resamples judges within
  panels and refits both stages.
- The interpretation of
  [`btl_next_pairs()`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)
  as a ranking rule is stated in its documentation.

## rasch 1.11.0

- Threshold and item-location covariance is transformed consistently
  after recentering items with different maximum scores.
- Warm WLE uses the common-discrimination score equation in which the
  common discrimination cancels.
- Untestable item-trait statistics return `NA`, and the total test
  includes testable items only.
- Equating drift tests account for the estimated origin shift and joint
  covariance of common-item locations.
- Judge-clustered inference requires more than one judge and reports a
  caution for fewer than ten clusters.
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  uses the judge as the sampling unit and count-weighted opponent bands.
- [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  and [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  warn when estimation has not converged.

## rasch 1.10.3

- [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  reports unknown `id`, `factors`, and `items` columns as errors.
- Fractional scores are rejected rather than truncated.
- MFRM rows with missing design identifiers are omitted with a note.
- [`equate_tests()`](https://drjoshmcgrane.github.io/rasch/reference/equate_tests.md)
  excludes common items without usable locations or standard errors from
  weighted linking and drift inference.
- [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  escapes data-derived labels and notes.

## rasch 1.10.2

- The Shiny application adds consistent hover labels to scalograms,
  residual heatmaps, equating plots, tailored analysis, and
  paired-comparison displays.

## rasch 1.10.1

- Person- and item-fit plots in the Shiny application show identifiers,
  locations, and fit residuals on hover.

## rasch 1.10.0

- [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  adds composite-likelihood AIC and BIC based on the Godambe effective
  parameter count.
- Model comparison now covers BTL position, threshold, and dependence
  specifications. The Shiny comparison page supports these fits.
- Optional cross-package tests compare results with `eRm`, `sirt`, and
  `psychotools`.

## rasch 1.9.3

- Pairwise estimation now checks that the observed item graph, together
  with any anchors, identifies a common scale.
- Thresholds adjoining critically sparse categories are marked weak.
  Their threshold and item-location standard errors are withheld.

## rasch 1.9.2

- Case study: `inst/casestudies/party_blocs_crisis.R` applies
  [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  to the Tuebingen 2009 party-preference data.
- Sets without stable panel-ratio information are excluded from the unit
  reconciliation and refitted at the reconciled panel units.
- Boundary-unstable bootstrap parameters receive `NA` standard errors
  with the number of boundary replicates reported.
- BTL-EFRM convergence is assessed by the gradient per comparison.

## rasch 1.9.1

- The Shiny Frames page supports BTL-EFRM panel and object-set
  definitions, unit tables, frame fit, and unit plots.

## rasch 1.9.0

- [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  fits the paired-comparison extension of the extended frame of
  reference model, with panel units, object-set units, and set origins.
- Bootstrap standard errors refit both stages. Conditional standard
  errors remain available for descriptive work.
- [`plot_btl_units()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_units.md)
  and
  [`simulate_btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl_efrm.md)
  support display and simulation.

## rasch 1.8.0

- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) adds
  anchored estimation and a first-position effect.
- [`btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)
  and
  [`plot_btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_equate.md)
  provide common-object linking and drift tests for paired-comparison
  calibrations.
- [`btl_information()`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md),
  [`plot_btl_targeting()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_targeting.md),
  and
  [`btl_next_pairs()`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)
  provide design information and greedy next-pair selection.
- Count-weighted BTL sandwich covariance now reproduces expanded-data
  covariance. Equating uses each calibration’s stored covariance.
- The Shiny Targeting and Equating pages support paired comparisons.

## rasch 1.7.1

- Rasch simulation layers now retain all previously specified DIF,
  dimensionality, response-style, dependence, and guessing terms.
- PCM simulation uses item-specific threshold structures and applies
  stricter validation to DIF, guessing, and disordered-threshold
  specifications.
- BTL dimensionality references include fitted within-judge dependence.
- MFRM simulation and recovery use the recorded item and rater
  parameters consistently.
- The Simulate page links recovery output to the current generated
  dataset.
- A simulation vignette and package logo were added.

## rasch 1.7.0

- Simulation functions add population-distribution controls, response
  styles, speededness, and MFRM halo effects.
- [`sim_replicate()`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md),
  [`sim_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/sim_recovery.md),
  and
  [`plot_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/plot_recovery.md)
  support repeated simulation and parameter-recovery summaries.
- The Shiny Simulate page exposes the population controls and recovery
  output.

## rasch 1.6.1

- The Shiny application adds a Simulate page for the four data layouts
  and their model-departure controls.

## rasch 1.6.0

- [`simulate_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md)
  gains `n_categories` for partial credit items within frames, with
  planted thresholds recorded in the truth attribute.
- [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md),
  [`simulate_btl()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl.md),
  [`simulate_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_mfrm.md),
  and
  [`simulate_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md)
  generate data from the package’s model families and can introduce
  nominated departures. Generating parameters are stored in the returned
  data.

## rasch 1.5.1

- [`plot_btl_judge_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_judge_map.md)
  now displays individual matchups.
  [`judge_pair_surprise()`](https://drjoshmcgrane.github.io/rasch/reference/judge_pair_surprise.md)
  returns the corresponding residuals.

## rasch 1.5.0

- [`judge_surprise()`](https://drjoshmcgrane.github.io/rasch/reference/judge_surprise.md)
  and
  [`plot_btl_judge_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_judge_map.md)
  compare a judge’s object-level preferences with the consensus object
  scale. The display is available from the Shiny Judge fit page.

## rasch 1.4.0

- [`btl_transitivity()`](https://drjoshmcgrane.github.io/rasch/reference/btl_transitivity.md)
  reports circular triads and Kendall’s consistency coefficient for
  suitable paired-comparison designs.
- [`btl_dimensionality()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
  decomposes the skew-symmetric residual preference matrix and compares
  its leading component with a model-based reference.
- New plots display BTL transitivity, scree, and residual maps.

## rasch 1.3.1

- The package was renamed from its development name, `rmt`, to `rasch`.
  Result classes use the `rasch_` prefix.

## rasch 1.3.0

- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) adds
  count-weighted exposure and carry-over effects, separation handling,
  and
  [`plot_btl_dependence()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dependence.md).
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  carries fitted dependence effects into its residual analysis and
  handles aggregated comparison counts.
- Mixed-design DIF terms are evaluated in their corresponding error
  strata. Repeated-person logit contrasts are provided by
  [`dif_contrasts()`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md).
- Residual parallel analysis uses data simulated from the fitted model.
- The Shiny application retains the settings used by each DIF analysis
  and applies consistent fit flags.

## rasch 1.2.0

- [`plot_pca_biplot()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pca_biplot.md)
  draws the item loadings on the first two residual principal components
  on equal axes.
- [`residual_correlations()`](https://drjoshmcgrane.github.io/rasch/reference/residual_correlations.md)
  now also returns the adjusted-Q3 `star_matrix` and
  [`plot_resid_cor()`](https://drjoshmcgrane.github.io/rasch/reference/plot_resid_cor.md)
  can draw raw Q3 or adjusted Q3\*.
- The Shiny trait and local-dependence pages pair tables with their
  plots and allow the original data to be restored after restructuring.
- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  is now the single DIF analysis-of-variance function. One factor is
  analysed one-way; several factors are fitted jointly. It supports
  repeated-measures and mixed designs.
- [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  resolves DIF iteratively by item splitting.

## rasch 1.0.0

First stable release.

### Models

- [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  fits dichotomous, partial credit, and rating scale models by pairwise
  conditional maximum likelihood, with Warm WLE person estimates.
- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  fits additive and item-by-facet many-facet models.
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  fits the extended frame of reference model.
- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) fits
  dichotomous and ordered paired-comparison models.

### Diagnostics

- Item and person fit, category functioning, targeting, reliability,
  test information, residual dimensionality, and local dependence.
- DIF analysis, item splitting, tailored analysis, common-item equating,
  and anchored calibration.
- BTL judge fit, invariance, and paired-comparison diagnostics.

### Display and reporting

- Base-graphics functions cover item, person, threshold, targeting,
  information, residual, and paired-comparison displays.
- [`fit_summary_table()`](https://drjoshmcgrane.github.io/rasch/reference/fit_summary_table.md)
  and
  [`targeting_table()`](https://drjoshmcgrane.github.io/rasch/reference/targeting_table.md)
  return the headline statistics;
  [`save_outputs()`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md)
  and
  [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  export results.
- [`run_app()`](https://drjoshmcgrane.github.io/rasch/reference/run_app.md)
  launches the Shiny interface and shows the R call corresponding to
  each analysis.
