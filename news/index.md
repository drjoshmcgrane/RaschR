# Changelog

## rasch 1.12.1

- [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  now labels unnamed response matrices consistently, and direct
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  calls reject an empty anchor table. Available-case item-rest
  correlations exclude respondents with no observed rest score.

- Frame-invariance bootstraps keep singleton person-group strata in
  their original group. Paired-comparison order effects are refused when
  the comparison design confounds them with the object locations.

- Automatic DIF follow-ups now use the adjustment method and
  significance level supplied to
  [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md).
  [`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
  also checks item selectors before fitting the contrast family.

- A saved app analysis now restores its data roles, estimation controls
  and embedded supporting data, so it can be re-estimated after it is
  reopened. Simulation recovery covers paired-comparison Extended
  Frames, including object locations, panel and set units, and set
  origins.

- Shiny background fits are tied to the data and analysis that launched
  them. A completed EFRM or paired-comparison frame fit is discarded if
  that context has changed, and opening a saved analysis cancels work
  still in progress. A new Comparative Judgement fit clears earlier
  requested DIF and frame results, while reopening a saved fit retains
  the results stored with it.

- The Shiny application can simulate ordinary and explanatory Rasch,
  Comparative Judgement, Multiple Ratings, Extended Frames and paired-
  comparison Extended Frames data. Model parameters can be varied and
  model-specific departures planted.
  [`simulate_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md)
  now supports item drift, careless response and missingness, while
  [`simulate_btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl_efrm.md)
  can plant erratic judges. A positive planted proportion affects at
  least one observation whenever the requested departures can coexist;
  incompatible mixtures are diagnosed. Rasch simulation truth records
  speeded persons and the cells selected for missingness.

- [`weighted_person_estimates()`](https://drjoshmcgrane.github.io/rasch/reference/weighted_person_estimates.md)
  provides supplementary person measures from externally imposed item or
  item-set weights. It leaves the fitted calibration and the ordinary
  measures used for fit, reliability, targeting and DIF unchanged. Its
  Warm correction and sandwich standard error now both use the
  variability of the weighted score; the earlier correction treated
  weights as replicated observations. A new first-listed vignette gives
  the data structure required by every model family.

- Simulation truth now retains person identifiers, so recovery matches
  person estimates by ID after rows are reordered. A planted second
  trait has its requested realised correlation. Paired-comparison
  simulators balance comparisons across the declared judges and refuse
  designs with too few comparisons. MFRM interaction probabilities
  require adequate support in every item-by-facet cell rather than only
  at the pooled facet level.

- EFRM convergence now covers the fitted nonparametric masses as well as
  the set transformation. Linking uncertainty requires at least 30
  usable resamples and more than half of those requested; the requested,
  usable and failed counts are retained in the fit. The full-bootstrap
  attempt remains recorded when too few refits succeed and hybrid
  standard errors are returned. BTL-EFRM applies the same
  usable-resample rule and accounting. Fixed-iteration EFRM linking no
  longer recalculates an unused likelihood at every intermediate mass
  update.

- Automatic DIF resolution now acts on the adjusted omnibus result.
  Pairwise contrasts describe the location of a multifactor effect but
  do not impose a second significance test. Thin or incompatible factor
  cells are still refused. A prior manual split now counts as one
  resolved source item rather than several anchors, and its provenance
  survives later item dropping or subtest formation. The reported
  residual count is the number of source items, not the number of
  item-by-factor terms.

- BTL-DIF retains the validated Welch reference for two-cell contrasts
  and uses the least-supported effective-judge cell as a conservative
  reference for contrasts spanning more cells. In 500 balanced four-cell
  null fits, the new reference rejected 3.6%, against 6.2% for the
  superseded pooled-count rule. BTL-EFRM no longer counts a set unit
  fixed after an identification failure as an estimated parameter.

- The application uses Holm-adjusted probabilities for its DIF, frame
  invariance, threshold-spread, explanatory and supplementary item-fit
  decisions. BTL judge-group factors that vary within judge are refused
  rather than reduced to the first comparison row.

- Subtests must be formed before DIF splitting. A group-specific split
  copy cannot be combined because it does not provide a common item
  across groups.

- BTL-EFRM now refits every object set after its panel units are
  reconciled. The reported object locations, expected probabilities,
  composite likelihood and equal-unit comparison therefore come from the
  same fitted parameters. Previously, stable sets retained their
  independently optimised locations while their probabilities were
  evaluated at the reconciled units.

- Multifactor DIF magnitudes now match the adjusted omnibus estimand.
  Ordinary and paired-comparison main effects average complete factor
  cells equally over nuisance factors; interactions use differences
  between differences. DIF contrasts are withheld when the compared
  groups do not share an observed score structure. ETS classifications
  use the adjusted probabilities for both significance and departure
  beyond category A.

- Structural refits preserve score categories. Subtests, DIF splits and
  EFRM frame resolutions are refused when a required score category is
  absent, rather than allowing ordinary data preparation to renumber the
  scores. Resampling refits treat category loss as an unsuccessful
  replicate. Supplied app anchors and scoring keys now fail closed when
  malformed, inapplicable or unmatched.

- Crossed EFRM factorial tests and BTL order-effect tests now retain raw
  probabilities but use Holm-adjusted probabilities for decisions. In
  fresh null simulations, crossed-EFRM familywise rejection was 5.55%
  over 2,000 fits (5.93% pooled over 3,000); the three marginal rates
  were 5.25–5.75%. BTL familywise rejection was 5.9% over 1,000 fits
  with position, exposure and carry-over fitted together.

- [`simulate_btl()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl.md)
  now constructs a second object attribute with the requested realised
  correlation, rather than obtaining that correlation only in
  expectation. The largest error over 1,200 short and long object sets
  was 3.9e-16; dimensionality power was 87% in the re-run of the strong
  design. Structured simulator options are checked before their
  components are used.

- DIF splitting omits factor levels for which an item has fewer than two
  observed categories, and refuses a split unless at least two levels
  remain. EFRM refits and wide MFRM conversion use collision-free
  internal names; named BTL-EFRM panel maps take precedence over
  data-column names. Numeric `NaN` responses are refused rather than
  treated as ordinary missing data.

- EFRM and BTL-EFRM retain raw unit-test probabilities but use
  Holm-adjusted probabilities for decisions. Omnibus tests are adjusted
  across the reported unit families. BTL-EFRM follow-up contrasts form
  one family across panel units, set units and set origins, rather than
  three separately adjusted tables. The application uses the same
  adjusted probabilities. In null simulations, EFRM omnibus-family
  rejection was 3.5% and follow-up-family rejection was 1.6% among 489
  analysed fits. After the reconciled-panel BTL-EFRM refit, the
  corresponding rates were 3.9% and 3.0% over 1,000 null fits in the
  six-judge-per-panel caution design. A 500-fit supported-design top-up
  gave 4.6% raw set-unit rejection and 0.934 interval coverage.

- Character role arguments are resolved from matching column names
  rather than their length. Named judge maps are matched by judge even
  when their length happens to equal the number of comparisons. Missing
  DIF identifiers remain as separate analysis units, conflicting
  external factor columns are refused, and scores outside R’s integer
  range no longer become missing on coercion. Paired-comparison equating
  now permits its documented descriptive two-object link while
  continuing to require three objects for drift tests.

- Common-item and common-object equating now label an unweighted
  descriptive shift when fewer than two common items or objects have
  usable variances. The functions no longer describe that fallback as
  precision-weighted or say that observations used in it were excluded
  from the shift.

- Every location axis is labelled at whole logits when the span allows,
  through one shared tick rule; the previous defaults could leave an
  axis extreme between labels. The person-item map begins and ends its
  proportion axis on labelled ticks, closes the axis corner, and drops
  its unused bottom margin.

- The Wright map and the person-item map mark the person and threshold
  means with dashed lines in their distributions’ colours.

- [`plot_person_fit()`](https://drjoshmcgrane.github.io/rasch/reference/plot_person_fit.md)
  and
  [`plot_item_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_item_map.md)
  can display the standardised infit or outfit in place of the fit
  residual, under the same +/-2.5 band, and annotate the flagged count
  with its percentage. The person table reports `infit_z` alongside
  `outfit_z`.

- Input and selection boundaries are hardened package-wide, closing two
  further review rounds. Every estimator and the DIF family validate
  their arguments through shared checks: iteration caps, tolerances,
  class-interval counts, reference sample sizes, significance levels,
  adjustment methods, practical thresholds, linking minima, component
  and replication counts, and plot bins all reject fractional,
  non-finite, or out-of-range values instead of silently truncating.
  Selection can no longer alter the analysis silently: EFRM refuses
  misspelled or missized id, factor, and item inputs exactly as ordinary
  Rasch does, matrix input honours the items argument in both,
  multiple-choice keys refuse fractional option scores and duplicate
  item entries, replication counts must be whole and finite, explanatory
  threshold labels and predictor rows are validated, duplicate column
  names are refused in direct estimation, and unknown names in anchor,
  distractor, and object-set requests are errors rather than silent
  drops. Pooled MFRM items flow through the DIF follow-ups under the
  stricter item resolver.

- Nine defect families from an adversarial review are corrected.
  Wide-format many-facet data now scores factor columns by their labels,
  where reordered factor levels previously shifted item locations. Item,
  threshold, and anchor indices are validated: a fractional or unknown
  index errors instead of silently truncating to a different item, and
  the numeric fitting controls reject fractional interval counts and
  non-finite references. Duplicate named mappings – comparison anchors,
  item-set maps, judge-panel maps, and named judge factors – are refused
  instead of silently taking one of the conflicting values. A boundary
  comparison object keeps its extrapolated location when a dependence
  effect was dropped after its removal, reports count-weighted
  comparisons, stays inside the plotted range, and is refused by with
  its calibrated companions unaffected. Explanatory departure
  probabilities are withheld for items whose thresholds the calibration
  marks as weak, and the judge diagnostics report count-weighted
  comparison totals.

- A row that is dropped can no longer define the analysis it is dropped
  from. The paired-comparison response scale is derived after zero-count
  and unusable rows are removed, where a single zero-count row could
  take an otherwise identical model from three categories to six; margin
  levels and the presence of ties follow the kept rows for the same
  reason. An ordered factor’s declared levels remain the stated scale
  whatever the weights are, so an empty declared extreme is still
  refused as the identifiability question it is. A tie carries no
  margin, so its margin value no longer opens win and loss categories
  nothing was judged in, which had left both extremes empty and stopped
  the fit; and a response column that is entirely missing reports no
  usable comparisons rather than failing on an impossible vector length.

- A missing person identifier is unknown, not shared. The DIF family no
  longer reads missing identifiers as repeats – which declared a
  repeated-measures design and changed every test, and in
  [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  withheld every standard error and Wald test in the table – and the
  tailored bootstrap resamples each unidentified row as its own person
  rather than clustering them into one. A genuinely repeated design
  still withholds its Wald inference.

- Item and object banks are read through their labels. A factor column
  of locations, standard errors or maximum scores was passed to
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html), which returns
  level codes, so a bank of -2.5, 0.25 and 4.0 equated as 1, 2 and 3.
  Numeric text remains accepted; invalid text and non-numeric classes
  are refused. An attached joint covariance now completes individual
  missing standard errors rather than doing so only when the bank
  omitted the whole column. Bank, predictor, key and anchor tables also
  require unique column names.

- A judging sequence must order the comparisons it describes. Values
  repeating within a judge left the order of those comparisons to the
  row order of the data, so the same data read in a different order
  carried different exposure and carry-over covariates; repeated and
  non-finite sequence values are now refused. A retained non-tie margin
  must be an ordered factor or a finite positive magnitude; zero denotes
  a tie. Margins on ties and excluded rows do not define the response
  scale. Logical, complex and Date columns are refused. A replication
  count must be real, since coercing a complex one discards its
  imaginary part in silence.

- A paired-comparison call must state one outcome. Supplying both
  `winner` and `response` fitted the response and ignored the winner, so
  changing every winner left the fit identical; the combination is now
  refused, in
  [`btl_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/btl_explanatory.md)
  too. Every column-role argument names exactly one existing column
  before it is dereferenced, in the comparison, frame-adjusted
  comparison, explanatory comparison and many-facet entry points alike,
  and in the many-facet case through the wide entry as well as the long
  one, which needs at least one item column. Many-facet person, item,
  score, facet and person-factor roles must also be distinct.

- Paired-comparison equating uses the documented precision-weighted
  origin shift when two common objects have usable standard errors.
  Three remain necessary for object-level drift inference, but that
  inferential threshold no longer changes the descriptive link
  estimator.

- Named item-set, judge-panel, judge-factor and Wright-map assignments
  must cover their fitted units exactly. Missing entries no longer
  discard units, and extra entries no longer pass as unnoticed spelling
  errors. Empty or blank panels are refused. Repeated-measure person and
  time roles must be distinct.

- Simulation counts, parameters and seeds are read only as plain numeric
  values. Factors and classed or complex vectors can no longer be
  interpreted through their storage codes, and replicate seeds cannot
  overflow the integer range. A zero-effect BTL dependence specification
  is treated as no planted dependence.

- The low-level
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  and
  [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  estimators now refuse negative, non-consecutive, constant and
  all-missing item-score columns. Their pair tables require at least two
  observed categories numbered from zero; negative values were
  previously omitted by
  [`tabulate()`](https://rdrr.io/r/base/tabulate.html) rather than
  diagnosed. The main
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  entry point also refuses non-finite scores rather than treating them
  as ordinary non-numeric missing entries.

- MFRM and EFRM exports and app displays draw their observed residual
  scree without requesting the parallel reference that is unavailable
  for their virtual response-cell designs.

- A multiple-choice scoring table refuses a missing or blank item name
  in every key form – the option/score table, the item/key table and the
  named vector – as well as a missing or blank option, and a key refuses
  a blank value; both had scored their item zero throughout and then
  dropped it as constant under a misleading message. A saved analysis
  file’s schema is read as stored rather than coerced, so “1”, TRUE, 1.5
  and a factor’s level code are no longer accepted as schema 1. An
  analysis file that declares a model type that is not one character –
  several, missing, or a factor whose integer code would select a branch
  by level order – is reported as unsupported instead of failing inside
  a length-one condition.

- External person factors supplied as a data frame are checked for
  constancy within person exactly as named columns are, instead of the
  first row’s value being kept. A frame group given by value no longer
  removes a person factor whose name matches one of the group labels,
  and a named panel map must give one stated panel for every judge.

- Set names are trimmed and validated in both frame families: a
  whitespace-only set name, or a blank set in the item-to-set map, is
  refused rather than fitted.

- An empty frame definition is refused rather than dropped: an item set
  or object set naming no identifier would have been fitted away,
  answering a different question from the one asked. A paired-comparison
  DIF call needs at least one judge factor.

- Administration patterns are built with explicit person-by-set
  dimensions. A single respondent – one person in a group, or a
  one-person fit – simplified to a vector that was then read transposed,
  so the design could name sets or facet cells the person never saw.

- A person factor may not take a name the fitted person table generates
  for itself: a factor called `class_interval` silently replaced the
  intervals every fit statistic is computed over, and one called `theta`
  stopped the fit. The generated names are reserved centrally. A frame
  design refuses repeated group columns, and the generated name for
  items no set lists is refused when a nominated set already uses it.

- The standalone HTML report names the estimator the fit used, as the
  fit summary table does, instead of always reporting pairwise
  conditional estimation.

- Names are carried as values rather than parsed out of labels: an
  item-set name containing the label separator keeps its items in the
  score curves; a predictor level no item or object carries is dropped
  before the design is built, where it added an all-zero column and made
  an identified model look rank-deficient; and a resolved comparison
  copy whose generated name already belongs to another object is refused
  with its magnitude withheld, rather than the two silently merging. A
  report table longer than its display limit carries the omission note
  as a caption, so it still renders as a table.

- EFRM score curves are keyed by the administration as well as the
  group: people in one group who sat different item sets have different
  maximum scores and different expected totals, and previously shared
  one curve. The table gains `design` and `n_persons` columns.

- A response style redraws from the probabilities the response was drawn
  under, so planted local dependence survives it, and a style of zero
  strength or zero prevalence is no longer recorded as an active effect.
  In a chain of dependence pairs the middle item’s expectation now
  includes its own carry-over, where the residual it passed on otherwise
  carried the first pair’s shift as a systematic mean into the second.

- An object set aside at a response boundary keeps its predictor row, so
  an explanatory comparison model no longer fails when one object always
  wins or always loses; a predictor row for an object the comparisons
  never mention is still an error.

- [`plot_pcc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pcc.md)
  and
  [`plot_kidmap()`](https://drjoshmcgrane.github.io/rasch/reference/plot_kidmap.md)
  refuse a person identifier that appears in more than one row, naming
  the rows, rather than silently drawing the first. Repeated identifiers
  are the ordinary case in stacked and racked longitudinal data.
  [`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
  validates the repeated-measures identifier where it is described
  rather than failing later on length.

- [`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  computes the covariance of the centred location differences as C{V1 +
  V2}C’, not C{V1 + V2}C. The centring matrix is not symmetric when the
  compared items differ in maximum score, so the standard errors,
  statistics, p values, rmse and ratio were wrong in that case: a
  polytomous item’s standard error was inflated (22% in a five-category
  example) and every dichotomous item’s deflated, flagging short items
  and hiding long ones. Equal maximum scores were unaffected.

- Derived fits keep the controls they were built from.
  [`lr_test()`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)‘s
  rating-scale refit carries the reference sample size and the person
  factors, so both parameterisations’ item-trait statistics are on one
  scale and the refit can be used for follow-up analyses; the fourth
  step of
  [`tailored_analysis()`](https://drjoshmcgrane.github.io/rasch/reference/tailored_analysis.md)
  carries the reference sample size its three siblings use; and a
  subtest total is no longer read as missing when it happens to equal a
  missing-data code, which had deleted complete responses and renumbered
  the categories that remained.

- [`chisq_detail()`](https://drjoshmcgrane.github.io/rasch/reference/chisq_detail.md)
  reports the rescaled per-interval components beside the standardised
  ones, so the detail and the item table reconcile under `adjust_N`, and
  gives the adjustment factor and the unadjusted total.

- Item and object drift are refused for explanatory calibrations, where
  a location is a function of its predictors: a drifted item is smeared
  over every item sharing its design cell and the standard errors belong
  to the design coefficients. The dimensionality magnitude is refused
  for the same reason, since its subtest refit frees every superitem. An
  explanatory comparison design now centres on the objects actually
  calibrated, so its locations sit on the model’s own origin.

- An inestimable DIF analysis is refused rather than returned malformed,
  [`plot_frames()`](https://drjoshmcgrane.github.io/rasch/reference/plot_frames.md)
  draws without intervals when the units carry no standard error instead
  of failing on an empty range, the frame-invariance summary counts item
  comparisons and items separately, and the observed points of the item
  displays follow the fit’s own per-item class intervals, so a graphical
  fit check shows the intervals of the test it illustrates.

- Simulation plants what it records, in three further cases: speededness
  needs a not-reached tail, a dependence source may not be regenerated
  as a later target, and a bias planted on a rater who answers at random
  is refused. Frame group units follow their group labels rather than
  the sorted level order, which attached the wrong unit to each group
  from ten groups upward.
  [`sim_apply()`](https://drjoshmcgrane.github.io/rasch/reference/sim_apply.md)
  requires one atomic scalar per replicate.

- Reports and exports check their arguments before they write: the
  output path, title, plot dimensions and resolution are validated
  first, so a bad size can no longer leave a populated folder that reads
  as a complete export. A plot archive is written fresh rather than
  appended to, and the report closes only the devices it opened. Report
  tables state how many rows were omitted and print small probabilities
  in the package’s own vocabulary rather than as an impossible zero.

- [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  treats presentation as data and the position covariate as a model
  term: paired-comparison fits differing only in which object was
  presented first, or in the judging sequence, are no longer reported as
  the same data, while a plain fit and a position-effect fit of the same
  comparisons still are.

- Four procedures no longer relax an explanatory restriction in silence.
  [`lr_test()`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)
  refuses an explanatory fit, whose rating re-parameterisation would
  drop the design and leave the two models not nested;
  [`tailored_analysis()`](https://drjoshmcgrane.github.io/rasch/reference/tailored_analysis.md)
  refuses one, because the tailored recalibration would differ from the
  original by the design as well as the tailoring; and
  [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  refuses an explanatory comparison fit, whose resolved copies would
  either be forced equal by the design or estimated without it. The
  parallel scree reference now analyses each simulated draw under the
  model that was fitted, where an explanatory calibration was previously
  compared against an unrestricted refit.

- Simulation plants what it records. An item named as the second element
  of several dependence pairs now carries every one of them, where a
  later pair regenerated the item and erased the earlier dependence; a
  dichotomous item cannot be disordered, so the request is refused with
  a warning instead of recorded as planted; and a set or group unit
  ratio must be 1 when there is only one set or group, since a ratio
  between frames cannot be planted in a single frame.

- A failed plot export no longer closes the caller’s graphics device: it
  closes only a device it opened itself. Item names that sanitise to the
  same file stem now keep separate files, where the later plot silently
  overwrote the earlier and the export still looked complete.

- [`plot_pcc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pcc.md)
  draws the person characteristic curve from the fitted model’s own
  expectations. For a polytomous, many-facet, or explanatory fit it
  previously drew a dichotomous logistic curve, which ignores the
  thresholds and the frame’s units; the curve is now the expected
  proportion of the maximum score across the fitted item locations, and
  a dichotomous fit with a common discrimination keeps its exact
  logistic form.

- Selection can no longer alter a multiplicity family in silence.
  Duplicate items, objects, contrast names, and contrast cells are
  refused, because a repeated hypothesis quietly changes the Holm
  adjustment; a grouping given to the DIF family must name fitted
  factors or supply one value per person, rather than being recycled
  into a grouping the fit never contained; dimensionality subsets must
  be free of duplicates and disjoint; and a person-factor frame must
  carry unique, non-empty names whichever input branch assembles it.

- Sustained adversarial review closed the remaining input and display
  boundaries. Role selection is unambiguous in
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  and
  [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  alike: a vector as long as the person count is read by value even when
  its labels collide with column names, and a vector whose values match
  a data column exactly is refused, with `items` named as the
  resolution, where an item whose responses happened to agree with a
  role vector was previously dropped from the fit without notice. A
  repeated-measures identifier must carry one entry per fitted row: a
  short vector was recycled, understating the sampling units and
  changing every test in the DIF table. A person with a missing or blank
  frame group is refused rather than expanded into rows missing from
  every set. Reshaping refuses missing or blank person identifiers and
  occasions and requires one existing column name for the person and
  occasion; paired comparisons, frame-adjusted comparisons, and
  many-facet data refuse blank objects, judges, panels, persons, items,
  and facet levels, since a whitespace label is not a level but was
  calibrated as one; and every stated flag must be TRUE or FALSE.
  Display controls are checked before anything is drawn: an evaluation
  grid needs at least two finite locations, class intervals and bins
  must be whole numbers, limits must be two finite ascending values,
  label sizes and colour caps must be positive and finite,
  single-person, single-item, and single-facet displays take exactly one
  name, and limits admitting no thresholds report an empty range instead
  of drawing an empty panel. Batch plot exports and the HTML report name
  the plots they could not draw, where a file was previously written
  with the failures omitted; a batch in which nothing could be drawn is
  now an error rather than a returned path to an archive that was never
  created, and the archive is confirmed on disk before the path is
  returned. Export device dimensions are validated before any device is
  opened.

- The observed points of a paired-comparison ICC no longer abort the
  display when every comparator falls below the informativeness
  threshold; the model curve and the omission note draw on their own.

- A package-wide sweep for reporting-table misalignments found and
  corrected three defects. A numeric anchor index now resolves against
  the data as supplied, where previously a dropped constant item shifted
  the anchor to the wrong item. The class-interval detail refuses an
  item whose only responders are extreme, instead of failing obscurely
  and aborting exports.
  [`dependence_magnitude()`](https://drjoshmcgrane.github.io/rasch/reference/dependence_magnitude.md)
  withholds its standard error and probability when the resolved
  thresholds are weakly identified, reporting the magnitude
  descriptively, as the thresholds themselves already were.

- An undefeated or winless comparison object is now reported in the
  object table at an extrapolated location, its score moved half a point
  inside the boundary against the calibrated scale, with
  `extreme = TRUE` and its standard error withheld – the reporting
  practice already used for extreme person measures. The row takes no
  part in estimation, inference, or equating. Validated against
  [`sirt::btm`](https://rdrr.io/pkg/sirt/man/btm.html): identical
  likelihood to machine precision on clean replicates, with the boundary
  policies agreeing in direction.

- The heaviest examples are smaller, and CRAN runs fewer scenario test
  blocks, keeping the check well inside the incoming pretest budget.

## rasch 1.12.0

CRAN release: 2026-08-24

### Models and inference

- [`wright_map()`](https://drjoshmcgrane.github.io/rasch/reference/wright_map.md)
  sends fitted person and item estimates to `WrightMap`. It supports
  several person distributions and the item-panel layout introduced in
  WrightMap 1.5, including the person-group and item-set structure of
  EFRM fits.

- [`rasch_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_explanatory.md)
  fits the linear logistic test model and linear partial credit model
  from continuous, categorical or ordinal item- or threshold-level
  predictors. Formulae may include selected interactions.
  [`explanatory_test()`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_test.md)
  compares the restrictions with a free calibration using the Kent
  adjustment;
  [`explanatory_diagnostics()`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_diagnostics.md)
  and
  [`relax_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/relax_explanatory.md)
  support fixed item and threshold departures. Refitted departures
  propagate to item and person estimates and are retained through item
  deletion, DIF splitting, superitem construction and
  response-dependence resolution. Keyed option responses remain
  available after item deletion, splitting and fixed-departure refits.

- [`btl_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/btl_explanatory.md)
  applies a fixed explanatory design to object locations in dichotomous
  or ordered comparative judgements. It supports the same model
  comparison, Holm-adjusted diagnostics and fixed departures while
  retaining the nominated ordered-response threshold structure.

- A worked case study on the documentation site uses the verbal
  aggression data to develop and check an explanatory partial credit
  model.

- [`explanatory_test()`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_test.md)
  now places the Kent-calibrated probability in both `p` and `p_kent`.
  The unscaled composite-likelihood probability is named `p_naive` so it
  cannot be mistaken for the inferential result. The table also reports
  calibration R-squared, with an adjusted counterpart whose null
  expectation is near zero, against the free threshold or object
  calibration.

- Pairwise conditional calibrations now use the remaining Newton move as
  a second convergence check. This prevents numerical false refusals at
  large sample sizes without changing the estimates.

- Holm adjustment for item-fit statistics now excludes items whose tests
  are unavailable. Their probabilities remain `NA`.

- Sparse-unit safeguards now use the sampling units that inform each
  test. MFRM interaction tests use the least-supported item-by-level
  cell; EFRM unit tests require adequate persons on every group or set
  link; BTL-EFRM judge bootstraps require adequate effective judges in
  every panel or link; and frame-invariance tests exclude weak frame
  calibrations.

- BTL-EFRM judge bootstraps now distinguish refit errors from
  non-convergence and report the underlying worker error when parallel
  refits fail.

- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  supports several facets and an optional item-by-facet interaction.
  Omnibus and cell follow-up tests use the fitted joint covariance.

- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  supports crossed person-group factors and reports their GLS factorial
  decomposition. Set-unit linking uses a finite-grid semiparametric
  likelihood, with a separate nuisance distribution for each observed
  person group. Hybrid standard errors retain the joint uncertainty of
  the within-frame calibration and set link; full person-bootstrap
  inference remains available. The convergence flag covers both
  estimation stages, and non-converged links are excluded from bootstrap
  covariance calculations. EFRM data require one response row per
  person.

- The repeated semiparametric linking calculations in the EFRM bootstrap
  now use a compiled numerical kernel. Bootstrap replicates can also be
  distributed over a reproducible, cross-platform worker cluster. The
  Shiny application runs EFRM fits in a background process, defaults to
  four workers where the system permits, records the bootstrap seed,
  reports progress and permits the fit to be cancelled without retaining
  a partial result.

- BTL-EFRM judge bootstraps likewise default to four workers where
  available. A fixed seed gives the same result for any worker count.
  The application runs these fits in the background and supports
  progress reporting and cancellation.

- [`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  compares item locations and discrimination across separately
  calibrated frames. The conditional method tests locations and reports
  discrimination descriptively. The person-within-frame bootstrap
  provides inference for both, with one combined Holm family.

- MFRM and EFRM summaries report item estimates separately from the
  item-by-facet or item-by-frame response cells used in estimation.
  Coefficient alpha is not reported for the expanded response-cell
  matrix. EFRM DIF tests pool residual evidence by item and exclude the
  person factors that define the frames.

- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md),
  [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  and
  [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  add ordered paired comparisons, judge-clustered inference,
  judge-factor DIF, linked object sets and judge panels.
  Paired-comparison diagnostics now include equating, transitivity,
  residual dimensions, design information and adaptive pair selection.

- Carry-over probabilities are withheld below 30 judges.
  [`btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)
  uses Welch–Satterthwaite degrees of freedom when fitted calibrations
  have a finite number of judge clusters. Conditional BTL-EFRM unit
  probabilities are withheld; the application defaults to the judge
  bootstrap. BTL-EFRM judge bootstraps use finite-judge references,
  whereas its independent-outcome parametric bootstrap uses normal and
  chi-square references.

### Differential item functioning

- Confirmatory multiplicity defaults are now consistently Holm
  familywise adjustments across item fit, DIF, equating and the
  application. BH remains available where false-discovery-rate screening
  is explicitly requested. BTL DIF uses HC3 covariance for unequal judge
  workloads and withholds omnibus probabilities below eight judges or
  eight effective judges in a factor cell.

- Item-fit documentation now distinguishes the principal item-trait test
  from the supplementary class-interval ANOVA and notes the limits of
  both in short administrations. HC3 was evaluated for item fit and was
  not adopted.

- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  fits several person factors jointly using Type II sums of squares.
  Repeated measurements use the person as the sampling unit and separate
  between- and within-person error strata. Multiplicity adjustment
  covers the complete family of uniform and non-uniform DIF tests rather
  than treating each term as a separate family;
  [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  follows the same rule. Uniform between-person terms now use HC3
  covariance. Class-interval interactions retain the residual-ANOVA
  reference used for non-uniform DIF.

- [`dif_contrasts()`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md)
  and
  [`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
  provide planned and post-hoc logit contrasts, including simple effects
  and difference-in-differences for interactions. MFRM follow-ups pool
  the fitted facet cells of an underlying item; resolved EFRM follow-ups
  are withheld because an ordinary split would discard the frame units.
  The residual-mean Tukey table has been removed from
  [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md);
  [`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
  is the supported follow-up for multilevel terms.

- Repeated-measures DIF follow-ups use the full design-cell weights in
  their person-level tests. Reported resolved estimates and
  probabilities therefore address the same marginal contrast when
  nuisance factors are imbalanced.

- [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  reports resolved pairwise logit differences. Dichotomous items receive
  the itemwise ETS A/B/C classification. Polytomous items report the PCM
  signed expected-score area descriptively, without importing an
  incompatible score-metric classification.

- [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  splits confirmed DIF items iteratively while retaining a minimum
  anchor set. Automatic splitting is restricted to uniform DIF;
  non-uniform DIF remains visible for item review. MFRM residuals can be
  pooled to their source items, and EFRM factors that do not define
  frames can be tested.

- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  retains anchors and fitted dependence terms in its resolution refit.
  Resolved pairwise inference is withheld unless each factor cell has at
  least eight effective judges; pairwise degrees of freedom use the two
  cells’ effective counts, and the pairwise table reports the raw and
  effective support for both cells. BTL-EFRM fits require a
  frame-specific analysis rather than the equal-unit resolution model.

### Diagnostics and model changes

- Identification checks now cover item, facet, frame and
  paired-comparison graphs, rank, separation and sparse categories.
  Unidentified estimates are refused; identified but weak estimates are
  marked or have inference withheld.
- [`dependence_magnitude()`](https://drjoshmcgrane.github.io/rasch/reference/dependence_magnitude.md)
  uses the joint covariance of resolved thresholds. Equating tests
  require independent calibrations and the covariance of banked
  locations.
- [`spread_test()`](https://drjoshmcgrane.github.io/rasch/reference/spread_test.md)
  applies the binomial least-upper-bound only to superitems formed
  entirely from dichotomous components. It now distinguishes a point
  estimate below the bound from adjusted one-sided evidence of
  dependence. Its significance level and multiplicity adjustment are
  available in the application. The component structure is retained
  through subsequent item splits and removals.
- The tailored-analysis bootstrap resamples complete persons, including
  all rows of a repeated-measures record.
- [`drop_items()`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md),
  [`resolve_frames()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_frames.md),
  DIF splitting and superitem construction refit the active model and
  update downstream item and person estimates. Refit specifications
  retain anchors, keyed scoring, threshold constraints, factors and
  frame-linking controls; a non-converged downstream calibration is not
  returned as a completed analysis.
- Classical whole-test statistics and the Guttman scalogram are withheld
  when an item is represented by several facet or frame response cells.
  They remain available for a one-cell-per-item reduction.
- MFRM characteristic and information curves now combine facet
  conditions administered to the same person. Distinct rating designs
  receive separate curves rather than being added into a test no person
  received.
- Automatic model comparisons are available for the main model families.
  Structural changes are accompanied by before-and-after item and person
  summaries.

### Application and documentation

- The Shiny application uses responsive control and result columns,
  compact explainers for outputs and options, scalable plots and
  downloadable tables. Plot controls sit below the plot, and related
  item curves may be overlaid.
- Analyses can be saved as `.rasch` projects and reopened. Reports can
  be produced as self-contained HTML, Word or PDF documents; the R code
  for each displayed result is available in the application.
- The application covers the extended model suite, including model
  comparison, DIF follow-ups, frame-invariance checks and refitted
  structural changes.
- The manuals and vignettes have been revised to state the fitted
  models, estimands, identification requirements and uncertainty methods
  directly.
- The shipped EFRM and BTL-EFRM case studies now use the current linking
  and uncertainty methods.
- [`plot_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md)
  and
  [`plot_btl_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_scree.md)
  label their component axes at whole components only, instead of
  overprinting the default axis.

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
