# rasch 1.12.0

This update extends the model suite, adds mixed-design DIF analysis, and
strengthens identification and uncertainty checks.

## Differential item functioning

* `dif_anova()` fits several person factors jointly. It supports additive or
  factorial effects and uses Type II sums of squares.
* Repeated measurements are analysed with the person as the sampling unit.
  Between-person and within-person terms use their respective error strata,
  with a Greenhouse--Geisser correction for multilevel within-person factors.
* `dif_contrasts()` provides planned logit contrasts, including
  between-by-within interactions. `dif_size()` separates descriptive logit
  magnitude from repeated-measures inference.
* MFRM residuals can be pooled to their underlying items for DIF analysis.
  Person factors not used to define an EFRM frame can also be tested.
* `resolve_dif()` splits items iteratively, beginning with the largest
  confirmed effect, while retaining a minimum anchor set.

## Facet, frame, and paired-comparison models

* `rasch_mfrm()` accepts several facets and optional item-by-facet
  interaction. The interaction family has an omnibus test with adjusted cell
  follow-ups.
* `rasch_efrm()` accepts crossed person-group factors and returns their GLS
  factorial decomposition. Hybrid and full-bootstrap standard errors propagate
  uncertainty in the estimated units.
* `btl()` supports ordered paired comparisons, judge-clustered covariance,
  position effects, and within-judge exposure and carry-over effects.
* `btl_dif()` tests object invariance over one or more judge factors using a
  judge-level split-plot analysis.
* `btl_efrm()` fits linked object sets and judge panels with set and panel
  units. Judge bootstrap inference refits both estimation stages.
* Paired-comparison functions now cover common-object equating, transitivity,
  residual bimensions, design information, and adaptive pair selection.

## Statistical corrections and identification

* `dependence_magnitude()` uses the full covariance of the resolved
  thresholds when calculating its standard error.
* The EFRM item-set unit linking recovers the true-score variance by a
  truncated-score-moment correction. The previous construction (observed
  variance minus mean squared standard error) under-recovers the
  true-score variance on short tests and biased recovered unit ratios
  upward by about five per cent at eight dichotomous items per set; the
  corrected estimator is unbiased there, confirmed against an external
  TAM slope-group anchor.
* The hybrid EFRM set-unit covariance includes uncertainty from the
  within-frame calibration.
* Judge-clustered BTL inference requires enough nominal and effective judges,
  residual cluster degrees of freedom, and a sufficiently balanced workload.
* BTL-EFRM unit tests use judge-limited F and t reference distributions.
* Standard errors are withheld for an item's thresholds when any response
  category is critically sparse.
* Item, facet, frame, and paired-comparison models check connectedness, rank,
  and separation before reporting estimates. Structurally unidentified
  parameters are refused; identified but imprecise parameters are marked.
* Equating tests require independent calibrations and the joint covariance of
  banked locations. Links without the required covariance remain descriptive.
* Information curves are calculated for item sets and facet designs that can
  actually be administered together.
* Classical statistics use complete responders by default. Available-case
  results are labelled exploratory.

## Follow-ups, comparisons, projects, and reports

* `frame_invariance()` tests the item invariance a frame model assumes
  rather than imposing it. The fitted model holds one location per item,
  shared across frames and scaled by the frame unit, so the assumption
  cannot be checked from the fit; the function calibrates each frame
  separately, puts the locations on the common scale, and compares them
  item by item, reporting the root mean squared difference against the root
  mean squared standard error. Its `adjust` argument chooses between Holm
  across all comparisons and no adjustment, because screening for items to
  examine and reporting a difference are different jobs: on eight to ten
  items the adjustment costs between 20 and 60 points of sensitivity, and
  simulation shows that carrying through to the repaired unit ratio. Both
  probabilities are reported either way, and the printed output names the
  rule it applied. The application exposes the choice as a switch.
* `drop_items()` removes items from a fitted analysis and refits it,
  keeping the model, person identifiers, factors, and -- for frame models
  -- the set structure and standard-error method. The application offers
  the same action on the selected item, recorded as an undoable step. For
  frame models this is a sensitivity analysis rather than housekeeping: a
  set unit is estimated from the dispersion its own items produce, so an
  item that fits its set badly moves the unit that decides whether the
  sets differ.
* `dif_posthoc()` provides post-hoc pairwise DIF comparisons after a
  significant omnibus term: resolved item-location differences in logits
  with the joint contrast covariance, Holm-adjusted, with trend and
  adjacent-level contrasts for ordered factors and simple-effect and
  difference-of-differences contrasts for interactions.
* The application presents automatic model comparisons beside each
  analysis (partial credit against rating scale, free against
  principal-component thresholds, additive against interactive
  many-facet, frame models against their equal-unit restrictions, and
  paired-comparison effect terms), separating formal tests,
  composite-likelihood information criteria, and descriptive fit changes.
  Structure-altering procedures show before-and-after summaries instead.
* Analyses can be saved as `.rasch` project files and reopened exactly;
  `report_document()` renders self-contained HTML and editable Word
  reports for every model family.

## Simulation, documentation, and interface

* New simulation functions cover each model family and retain the generating
  parameters for recovery studies. `sim_replicate()`, `sim_apply()`, and
  `sim_recovery()` support repeated simulation.
* Six vignettes cover the main Rasch workflow, mixed-design DIF, MFRM, EFRM,
  paired comparisons, and simulation studies.
* The Shiny application includes the extended model suite, model comparison,
  downloadable tables and plots, and the R call corresponding to each fit.
* The package and function documentation now state the fitted models,
  identification constraints, uncertainty methods, and principal references
  more directly.

# rasch 1.11.7

* `print()` preserves the reference distribution used by saved BTL fits:
  current and transitional results are labelled `t`, while older results
  without cluster degrees of freedom retain their original `z` label.

# rasch 1.11.6

* BTL print methods read both current and earlier dependence-statistic
  columns. Current clustered statistics are labelled `t`, in accordance
  with their t reference distribution.

# rasch 1.11.5

* Ordered paired-comparison margins must be numeric or an ordered factor;
  unordered categorical margins are rejected.
* Invalid graded responses and malformed secondary-dimension specifications
  now produce errors rather than being coerced or dropped.
* Clustered dependence and position statistics are labelled `t`.

# rasch 1.11.4

* Score validation reads factor labels and rejects non-integer, non-numeric,
  or non-finite values.
* Ordered BTL responses require an ordered factor or integer scores.
* Clustered dependence and position tests use a t reference with
  `G - 1` degrees of freedom.
* `simulate_rasch()` validates secondary-trait correlations and item sets.

# rasch 1.11.3

* All estimators share the same integer-score validation.
* `item_moments()` uses a log-sum-exp calculation for wide category ranges.
* BTL object-separation reliability is withheld when the clustered covariance
  is rank deficient.
* `btl_efrm()` requires each judge to belong to one panel.
* Undefined reliability and omnibus statistics are reported as `NA`.
* Secondary-trait simulation retains its requested mean and standard
  deviation, and MLE score-table calculations use the correct common-unit
  score equation.

# rasch 1.11.2

* The BTL dimensionality reference simulates count-weighted data at the
  unordered-pair level and returns the leading strength from each replicate.
* Judge-bootstrap EFRM inference requires more than one judge per panel and
  notes panels with fewer than five judges.
* BTL fits report when the number of judge clusters cannot support a
  full-rank clustered covariance.

# rasch 1.11.1

* BTL dimensionality reference draws now reproduce count-weighted binomial
  or multinomial sampling.
* Pairwise chi-square degrees of freedom include position and dependence
  parameters; untestable designs return `NA`.
* Judge-clustered covariance uses the CR1 small-sample factor.
* `btl_efrm(se_method = "judge_bootstrap")` resamples judges within panels
  and refits both stages.
* The interpretation of `btl_next_pairs()` as a ranking rule is stated in
  its documentation.

# rasch 1.11.0

* Threshold and item-location covariance is transformed consistently after
  recentering items with different maximum scores.
* Warm WLE uses the common-discrimination score equation in which the common
  discrimination cancels.
* Untestable item-trait statistics return `NA`, and the total test includes
  testable items only.
* Equating drift tests account for the estimated origin shift and joint
  covariance of common-item locations.
* Judge-clustered inference requires more than one judge and reports a
  caution for fewer than ten clusters.
* `btl_dif()` uses the judge as the sampling unit and count-weighted opponent
  bands.
* `rasch()` and `btl()` warn when estimation has not converged.

# rasch 1.10.3

* `rasch()` reports unknown `id`, `factors`, and `items` columns as errors.
* Fractional scores are rejected rather than truncated.
* MFRM rows with missing design identifiers are omitted with a note.
* `equate_tests()` excludes common items without usable locations or standard
  errors from weighted linking and drift inference.
* `report_html()` escapes data-derived labels and notes.

# rasch 1.10.2

* The Shiny application adds consistent hover labels to scalograms, residual
  heatmaps, equating plots, tailored analysis, and paired-comparison displays.

# rasch 1.10.1

* Person- and item-fit plots in the Shiny application show identifiers,
  locations, and fit residuals on hover.

# rasch 1.10.0

* `compare_fits()` adds composite-likelihood AIC and BIC based on the
  Godambe effective parameter count.
* Model comparison now covers BTL position, threshold, and dependence
  specifications. The Shiny comparison page supports these fits.
* Optional cross-package tests compare results with `eRm`, `sirt`, and
  `psychotools`.

# rasch 1.9.3

* Pairwise estimation now checks that the observed item graph, together with
  any anchors, identifies a common scale.
* Thresholds adjoining critically sparse categories are marked weak. Their
  threshold and item-location standard errors are withheld.

# rasch 1.9.2

* Case study: `inst/casestudies/party_blocs_crisis.R` applies `btl_efrm()`
  to the Tuebingen 2009 party-preference data.
* Sets without stable panel-ratio information are excluded from the unit
  reconciliation and refitted at the reconciled panel units.
* Boundary-unstable bootstrap parameters receive `NA` standard errors with
  the number of boundary replicates reported.
* BTL-EFRM convergence is assessed by the gradient per comparison.

# rasch 1.9.1

* The Shiny Frames page supports BTL-EFRM panel and object-set definitions,
  unit tables, frame fit, and unit plots.

# rasch 1.9.0

* `btl_efrm()` fits the paired-comparison extension of the extended frame of
  reference model, with panel units, object-set units, and set origins.
* Bootstrap standard errors refit both stages. Conditional standard errors
  remain available for descriptive work.
* `plot_btl_units()` and `simulate_btl_efrm()` support display and simulation.

# rasch 1.8.0

* `btl()` adds anchored estimation and a first-position effect.
* `btl_equate()` and `plot_btl_equate()` provide common-object linking and
  drift tests for paired-comparison calibrations.
* `btl_information()`, `plot_btl_targeting()`, and `btl_next_pairs()` provide
  design information and greedy next-pair selection.
* Count-weighted BTL sandwich covariance now reproduces expanded-data
  covariance. Equating uses each calibration's stored covariance.
* The Shiny Targeting and Equating pages support paired comparisons.

# rasch 1.7.1

* Rasch simulation layers now retain all previously specified DIF,
  dimensionality, response-style, dependence, and guessing terms.
* PCM simulation uses item-specific threshold structures and applies stricter
  validation to DIF, guessing, and disordered-threshold specifications.
* BTL dimensionality references include fitted within-judge dependence.
* MFRM simulation and recovery use the recorded item and rater parameters
  consistently.
* The Simulate page links recovery output to the current generated dataset.
* A simulation vignette and package logo were added.

# rasch 1.7.0

* Simulation functions add population-distribution controls, response styles,
  speededness, and MFRM halo effects.
* `sim_replicate()`, `sim_recovery()`, and `plot_recovery()` support repeated
  simulation and parameter-recovery summaries.
* The Shiny Simulate page exposes the population controls and recovery output.

# rasch 1.6.1

* The Shiny application adds a Simulate page for the four data layouts and
  their model-departure controls.

# rasch 1.6.0

* `simulate_efrm()` gains `n_categories` for partial credit items within
  frames, with planted thresholds recorded in the truth attribute.
* `simulate_rasch()`, `simulate_btl()`, `simulate_mfrm()`, and
  `simulate_efrm()` generate data from the package's model families and can
  introduce nominated departures. Generating parameters are stored in the
  returned data.

# rasch 1.5.1

* `plot_btl_judge_map()` now displays individual matchups.
  `judge_pair_surprise()` returns the corresponding residuals.

# rasch 1.5.0

* `judge_surprise()` and `plot_btl_judge_map()` compare a judge's object-level
  preferences with the consensus object scale. The display is available from
  the Shiny Judge fit page.

# rasch 1.4.0

* `btl_transitivity()` reports circular triads and Kendall's consistency
  coefficient for suitable paired-comparison designs.
* `btl_dimensionality()` decomposes the skew-symmetric residual preference
  matrix and compares its leading component with a model-based reference.
* New plots display BTL transitivity, scree, and residual maps.

# rasch 1.3.1

* The package was renamed from its development name, `rmt`, to `rasch`.
  Result classes use the `rasch_` prefix.

# rasch 1.3.0

* `btl()` adds count-weighted exposure and carry-over effects, separation
  handling, and `plot_btl_dependence()`.
* `btl_dif()` carries fitted dependence effects into its residual analysis
  and handles aggregated comparison counts.
* Mixed-design DIF terms are evaluated in their corresponding error strata.
  Repeated-person logit contrasts are provided by `dif_contrasts()`.
* Residual parallel analysis uses data simulated from the fitted model.
* The Shiny application retains the settings used by each DIF analysis and
  applies consistent fit flags.

# rasch 1.2.0

* `plot_pca_biplot()` draws the item loadings on the first two residual
  principal components on equal axes.
* `residual_correlations()` now also returns the adjusted-Q3 `star_matrix`
  and `plot_resid_cor()` can draw raw Q3 or adjusted Q3*.
* The Shiny trait and local-dependence pages pair tables with their plots and
  allow the original data to be restored after restructuring.
* `dif_anova()` is now the single DIF analysis-of-variance function. One
  factor is analysed one-way; several factors are fitted jointly. It supports
  repeated-measures and mixed designs.
* `resolve_dif()` resolves DIF iteratively by item splitting.

# rasch 1.0.0

First stable release.

## Models

* `rasch()` fits dichotomous, partial credit, and rating scale models by
  pairwise conditional maximum likelihood, with Warm WLE person estimates.
* `rasch_mfrm()` fits additive and item-by-facet many-facet models.
* `rasch_efrm()` fits the extended frame of reference model.
* `btl()` fits dichotomous and ordered paired-comparison models.

## Diagnostics

* Item and person fit, category functioning, targeting, reliability, test
  information, residual dimensionality, and local dependence.
* DIF analysis, item splitting, tailored analysis, common-item equating, and
  anchored calibration.
* BTL judge fit, invariance, and paired-comparison diagnostics.

## Display and reporting

* Base-graphics functions cover item, person, threshold, targeting,
  information, residual, and paired-comparison displays.
* `fit_summary_table()` and `targeting_table()` return the headline
  statistics; `save_outputs()` and `report_html()` export results.
* `run_app()` launches the Shiny interface and shows the R call corresponding
  to each analysis.
