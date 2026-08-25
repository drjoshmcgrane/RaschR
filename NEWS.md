# rasch 1.12.1

* Every location axis is labelled at whole logits when the span allows,
  through one shared tick rule; the previous defaults could leave an axis
  extreme between labels. The person-item map begins and ends its
  proportion axis on labelled ticks, closes the axis corner, and drops its
  unused bottom margin.
* The Wright map and the person-item map mark the person and threshold
  means with dashed lines in their distributions' colours.
* `plot_person_fit()` and `plot_item_map()` can display the standardised
  infit or outfit in place of the fit residual, under the same +/-2.5
  band, and annotate the flagged count with its percentage. The person
  table reports `infit_z` alongside `outfit_z`.
* The observed points of a paired-comparison ICC no longer abort the
  display when every comparator falls below the informativeness
  threshold; the model curve and the omission note draw on their own.
* A package-wide sweep for reporting-table misalignments found and
  corrected three defects. A numeric anchor index now resolves against
  the data as supplied, where previously a dropped constant item shifted
  the anchor to the wrong item. The class-interval detail refuses an
  item whose only responders are extreme, instead of failing obscurely
  and aborting exports. `dependence_magnitude()` withholds its standard
  error and probability when the resolved thresholds are weakly
  identified, reporting the magnitude descriptively, as the thresholds
  themselves already were.
* An undefeated or winless comparison object is now reported in the
  object table at an extrapolated location, its score moved half a point
  inside the boundary against the calibrated scale, with `extreme = TRUE`
  and its standard error withheld -- the reporting practice already used
  for extreme person measures. The row takes no part in estimation,
  inference, or equating. Validated against `sirt::btm`: identical
  likelihood to machine precision on clean replicates, with the boundary
  policies agreeing in direction.
* The heaviest examples are smaller, and CRAN runs fewer scenario test
  blocks, keeping the check well inside the incoming pretest budget.

# rasch 1.12.0

## Models and inference

* `wright_map()` sends fitted person and item estimates to `WrightMap`. It
  supports several person distributions and the item-panel layout introduced
  in WrightMap 1.5, including the person-group and item-set structure of EFRM
  fits.
* `rasch_explanatory()` fits the linear logistic test model and linear partial
  credit model from continuous, categorical or ordinal item- or
  threshold-level predictors. Formulae may include selected interactions.
  `explanatory_test()` compares the restrictions with
  a free calibration using the Kent adjustment; `explanatory_diagnostics()`
  and `relax_explanatory()` support fixed item and threshold departures.
  Refitted departures propagate to item and person estimates and are retained
  through item deletion, DIF splitting, superitem construction and
  response-dependence resolution. Keyed option responses remain available
  after item deletion, splitting and fixed-departure refits.
* `btl_explanatory()` applies a fixed explanatory design to object locations
  in dichotomous or ordered comparative judgements. It supports the same
  model comparison, Holm-adjusted diagnostics and fixed departures while
  retaining the nominated ordered-response threshold structure.
* A worked case study on the documentation site uses the verbal aggression
  data to develop and check an explanatory partial credit model.
* `explanatory_test()` now places the Kent-calibrated probability in both
  `p` and `p_kent`. The unscaled composite-likelihood probability is named
  `p_naive` so it cannot be mistaken for the inferential result. The table
  also reports calibration R-squared, with an adjusted counterpart whose
  null expectation is near zero, against the free threshold or object
  calibration.
* Pairwise conditional calibrations now use the remaining Newton move as a
  second convergence check. This prevents numerical false refusals at large
  sample sizes without changing the estimates.
* Holm adjustment for item-fit statistics now excludes items whose tests are
  unavailable. Their probabilities remain `NA`.

* Sparse-unit safeguards now use the sampling units that inform each test.
  MFRM interaction tests use the least-supported facet level; EFRM unit tests
  require adequate persons on every group or set link; BTL-EFRM judge
  bootstraps require adequate effective judges in every panel or link; and
  frame-invariance tests exclude weak frame calibrations.
* BTL-EFRM judge bootstraps now distinguish refit errors from non-convergence
  and report the underlying worker error when parallel refits fail.

* `rasch_mfrm()` supports several facets and an optional item-by-facet
  interaction. Omnibus and cell follow-up tests use the fitted joint
  covariance.
* `rasch_efrm()` supports crossed person-group factors and reports their GLS
  factorial decomposition. Set-unit linking uses a finite-grid
  semiparametric likelihood, with a separate nuisance distribution for each
  observed person group. Hybrid standard errors retain the joint uncertainty
  of the within-frame calibration and set link; full person-bootstrap
  inference remains available. The convergence flag covers both estimation
  stages, and non-converged links are excluded from bootstrap covariance
  calculations. EFRM data require one response row per person.
* The repeated semiparametric linking calculations in the EFRM bootstrap now
  use a compiled numerical kernel. Bootstrap replicates can also be distributed
  over a reproducible, cross-platform worker cluster. The Shiny application
  runs EFRM fits in a background process, defaults to four workers where the
  system permits, records the bootstrap seed, reports progress and permits the
  fit to be cancelled without retaining a partial result.
* BTL-EFRM judge bootstraps likewise default to four workers where available.
  A fixed seed gives the same result for any worker count. The application runs
  these fits in the background and supports progress reporting and cancellation.
* `frame_invariance()` compares item locations and discrimination across
  separately calibrated frames. The conditional method tests locations and
  reports discrimination descriptively. The person-within-frame bootstrap
  provides inference for both, with one combined Holm family.
* MFRM and EFRM summaries report item estimates separately from the
  item-by-facet or item-by-frame response cells used in estimation.
  Coefficient alpha is not reported for the expanded response-cell matrix.
  EFRM DIF tests pool residual evidence by item and exclude the person factors
  that define the frames.
* `btl()`, `btl_dif()` and `btl_efrm()` add ordered paired comparisons,
  judge-clustered inference, judge-factor DIF, linked object sets and judge
  panels. Paired-comparison diagnostics now include equating, transitivity,
  residual dimensions, design information and adaptive pair selection.
* Carry-over probabilities are withheld below 30 judges. `btl_equate()` uses
  Welch--Satterthwaite degrees of freedom when fitted calibrations have a
  finite number of judge clusters. Conditional BTL-EFRM unit probabilities
  are withheld; the application defaults to the judge bootstrap. BTL-EFRM
  judge bootstraps use finite-judge references, whereas its independent-outcome
  parametric bootstrap uses normal and chi-square references.

## Differential item functioning

* Confirmatory multiplicity defaults are now consistently Holm familywise
  adjustments across item fit, DIF, equating and the application. BH remains
  available where false-discovery-rate screening is explicitly requested.
  BTL DIF uses HC3 covariance for unequal judge workloads and withholds
  omnibus probabilities below eight judges or eight effective judges in a
  factor cell.

* Item-fit documentation now distinguishes the principal item-trait test from
  the supplementary class-interval ANOVA and notes the limits of both in short
  administrations. HC3 was evaluated for item fit and was not adopted.

* `dif_anova()` fits several person factors jointly using Type II sums of
  squares. Repeated measurements use the person as the sampling unit and
  separate between- and within-person error strata. Multiplicity adjustment
  covers the complete family of uniform and non-uniform DIF tests rather than
  treating each term as a separate family; `btl_dif()` follows the same rule.
  Uniform between-person terms now use HC3 covariance. Class-interval
  interactions retain the residual-ANOVA reference used for non-uniform DIF.
* `dif_contrasts()` and `dif_posthoc()` provide planned and post-hoc logit
  contrasts, including simple effects and difference-in-differences for
  interactions. MFRM follow-ups pool the fitted facet cells of an underlying
  item; resolved EFRM follow-ups are withheld because an ordinary split would
  discard the frame units.
  The residual-mean Tukey table has been removed from `dif_anova()`;
  `dif_posthoc()` is the supported follow-up for multilevel terms.
* Repeated-measures DIF follow-ups use the full design-cell weights in their
  person-level tests. Reported resolved estimates and probabilities therefore
  address the same marginal contrast when nuisance factors are imbalanced.
* `dif_size()` reports resolved pairwise logit differences. Dichotomous
  items receive the itemwise ETS A/B/C classification. Polytomous items
  report the PCM signed expected-score area descriptively, without importing
  an incompatible score-metric classification.
* `resolve_dif()` splits confirmed DIF items iteratively while retaining a
  minimum anchor set. Automatic splitting is restricted to uniform DIF;
  non-uniform DIF remains visible for item review. MFRM residuals can be
  pooled to their source items, and EFRM factors that do not define frames can
  be tested.
* `btl_dif()` retains anchors and fitted dependence terms in its resolution
  refit. Resolved pairwise inference is withheld unless each factor cell has
  at least eight effective judges; pairwise degrees of freedom use the two
  cells' effective counts, and the pairwise table reports the raw and effective
  support for both cells. BTL-EFRM fits require a frame-specific analysis rather
  than the equal-unit resolution model.

## Diagnostics and model changes

* Identification checks now cover item, facet, frame and paired-comparison
  graphs, rank, separation and sparse categories. Unidentified estimates are
  refused; identified but weak estimates are marked or have inference
  withheld.
* `dependence_magnitude()` uses the joint covariance of resolved thresholds.
  Equating tests require independent calibrations and the covariance of
  banked locations.
* `spread_test()` applies the binomial least-upper-bound only to superitems
  formed entirely from dichotomous components. It now distinguishes a point
  estimate below the bound from adjusted one-sided evidence of dependence. Its
  significance level and multiplicity adjustment are available in the
  application. The component structure is retained through subsequent item
  splits and removals.
* The tailored-analysis bootstrap resamples complete persons, including all
  rows of a repeated-measures record.
* `drop_items()`, `resolve_frames()`, DIF splitting and superitem
  construction refit the active model and update downstream item and person
  estimates. Refit specifications retain anchors, keyed scoring, threshold
  constraints, factors and frame-linking controls; a non-converged downstream
  calibration is not returned as a completed analysis.
* Classical whole-test statistics and the Guttman scalogram are withheld when
  an item is represented by several facet or frame response cells. They remain
  available for a one-cell-per-item reduction.
* MFRM characteristic and information curves now combine facet conditions
  administered to the same person. Distinct rating designs receive separate
  curves rather than being added into a test no person received.
* Automatic model comparisons are available for the main model families.
  Structural changes are accompanied by before-and-after item and person
  summaries.

## Application and documentation

* The Shiny application uses responsive control and result columns, compact
  explainers for outputs and options, scalable plots and downloadable tables.
  Plot controls sit below the plot, and related item curves may be overlaid.
* Analyses can be saved as `.rasch` projects and reopened. Reports can be
  produced as self-contained HTML, Word or PDF documents; the R code for each
  displayed result is available in the application.
* The application covers the extended model suite, including model comparison,
  DIF follow-ups, frame-invariance checks and refitted structural changes.
* The manuals and vignettes have been revised to state the fitted models,
  estimands, identification requirements and uncertainty methods directly.
* The shipped EFRM and BTL-EFRM case studies now use the current linking and
  uncertainty methods.
* `plot_scree()` and `plot_btl_scree()` label their component axes at whole
  components only, instead of overprinting the default axis.

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
