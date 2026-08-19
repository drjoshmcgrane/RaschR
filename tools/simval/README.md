# Simulation validation battery

Version note: development builds 1.12.1 through 1.14.2 were never
released; their content ships publicly as the 1.12.0 CRAN update. The
version numbers and commit hashes below refer to those internal
development states and remain the accurate provenance record.

Scripts and results for the release validation batteries summarised in the
plant-and-detect vignette ("Validation studies"). Nothing here ships in the
CRAN package (see `.Rbuildignore`); everything runs from the package root
with `Rscript`, loading the in-tree package via `pkgload::load_all(".")`.

## Layout

- `round1/` — the August 2026 battery (rasch 1.14.2): one directory per
  model/diagnostic family (`dichotomous`, `pcm_rsm`, `mfrm`, `efrm`, `btl`,
  `btl-efrm`, `equating`, `dif`, `dimensionality`, `mc`, `tailored`), each
  containing the study scripts with their seeds inline. `dm_*.R` are the
  follow-up studies that diagnosed and verified the `dependence_magnitude`
  standard-error correction (null calibration, covariance diagnosis at
  1,200 replicates, post-fix confirmation).
- `results/round1-checks.csv` — the structured result of every round-1
  check: area, condition, metric, observed, expected, pass/fail, notes
  (including Monte Carlo error where computed), and provenance columns.
  231 of 234 checks passed. The three failures are retained deliberately:
  two are the `dependence_magnitude` null-calibration failure that led to
  the 1.14.2 standard-error correction (those rows describe PRE-fix
  behaviour; the post-fix confirmations are `round1/dm_post_fix.R` and the
  1,200-replicate diagnosis `round1/dm_worker.R`), and one is a
  mis-designed power scenario in the battery itself (non-uniform DIF
  planted on an item whose discrimination crossover falls outside the
  observed trait range -- undetectable by construction; the redesigned
  check appears as a passing row). Round 1 was executed 2026-08-10/11
  against the 1.14.1 release-candidate R code (commits 433e2a5..326f0b9),
  before the 1.14.2 fix; the `area` column maps to `round1/`
  subdirectories as follows: dichotomous-core -> `dichotomous/`,
  PCM/RSM -> `pcm_rsm/`, MFRM -> `mfrm/`, EFRM -> `efrm/`,
  BTL/Recovery -> `btl/`, BTL-EFRM -> `btl-efrm/`,
  equate_tests -> `equating/`, DIF -> `dif/`,
  residual_correlations/dependence/dimensionality -> `dimensionality/`.
- `results/dependence-magnitude-fix.csv` — the full lifecycle of the one
  round-1 defect as structured rows: pre-fix Type I 7.5% (1,200
  replicates), 5.4% with the corrected formula on the same draws, 4.8%
  post-fix on fresh seeds, 3.5% partial credit; each row names its script
  and the package state it ran against.
- `parse_check.R` — parses every script in the battery and fails loudly on
  any syntax error; run it (and a smoke subset) before trusting the
  reproducibility claim.

## Known limitations surfaced by the battery

- Small-sample undercoverage for extreme items: with 20-100 persons per
  form, item-threshold coverage runs 0.88-0.94 (SE ratios up to ~1.5) for
  items placed 3 logits from the person mean. The fully crossed reference
  arm shows the same undercoverage as the linked-booklet arm at every
  sample size, and both are nominal by 600 persons -- this is
  small-sample behaviour of the pairwise-conditional standard errors for
  poorly-informed items, not a cost of the booklet structure, and it is
  compounded at the smallest sizes by conditioning on the 4-20% of fits
  the identification guard refuses (`structural-missingness.csv`,
  weak_links rows).
- Pooled threshold rows in the sparse-category scenarios mix items of
  very different precision; their emp_sd/mean_se exceeds 1 by Jensen's
  inequality even when every item calibrates. The per-item rows carry the
  calibration-relevant figures.

Exact provenance state of the result files: the three studies regenerated
in the final release round (`btl-clustered.csv`, `btl-share-sweep.csv`,
`efrm-fix-sweep.csv`) carry R-tree hash `4a0e2bf8b357`, matching commit
`8de917b` exactly; `lr-smalln-topup.csv` carries the dirty-tree hash of
its execution state (`42945dd80663`): the `lr_test` estimator code it
exercised was unchanged relative to the release tree, and the hash
difference comes from documentation comments and unrelated BTL source
edits in flight at execution time. Files produced before the tree-hash
mechanism existed carry none. Documentation-only edits under `R/` move the tree hash of
later commits without touching any estimator, so a result file's hash
identifies the sources it ran against, not necessarily HEAD.

Second-round rows carry exact provenance automatically: `sv_row()` stamps
each row with the study script (`options(simval.script = ...)`), the
package git SHA (`+dirty` when the R code differs from HEAD), and the run
date. `round1-checks.csv` predates this mechanism; its rows carry the
executed date, the commit range, and a per-row `script_dir` mapping to the
`round1/` directory that produced them (exact per-script attribution was
not recorded by the round-1 orchestration and is documented at directory
granularity deliberately).
- `studies/` — the second-round studies (custom Wald/contrast tests,
  clustered comparative-judgement inference, tailored bootstrap
  calibration, `lr_test` size and power, equating multiplicity, structural
  missingness, WLE coverage), each writing its result table to `results/`.
  The `custom-wald-tests` and `btl-clustered` tables reflect re-runs after
  the 1.14.2 statistical corrections (their pre-fix rows motivated those
  corrections and are superseded); `tailored-bootstrap` documents the
  bootstrap's strong conservatism at feasible replicate counts, which is
  the procedure's documented design.
- `results/efrm-fix-sweep.csv` and `results/efrm-fix/` — the extended
  frame set-unit correction: pre-fix 9.4% null rejection (1,000
  replicates), the corrected hybrid at 4.9% across eight design cells, the
  seed-paired marginal-vs-joint draw comparison (identical decisions in
  all 400 replicates), and the full-bootstrap benchmark arm (ratio 0.963)
  showing the corrected hybrid prices the same uncertainty.
- `results/round2-followups.csv` — the adjudication trail for every
  round-2 suspect: the MFRM q=25 exoneration (600 fixed-truth replicates
  per cell), the btl_efrm origin-test correction (8.5% pre-fix to 5.5%
  post-fix at 400 replicates each), the PCM item-level guard validation,
  and the concentration guard's field behaviour.
- `results/comparison-validation.csv` — the release gate for the app's
  automatic model-comparison cards: every comparison surface validated
  under its null model and at least two departure magnitudes. Citation
  rows point surfaces already validated elsewhere (`lr_test`, the EFRM,
  BTL-EFRM, and MFRM omnibus tests) at their provenance CSVs. New cells
  (400 replicates per selection condition): CL-AIC null false selection
  5.2% for PCM-vs-RSM and 4.5% for free-vs-two-component thresholds
  (matching the theoretical multi-parameter AIC rates), rising to ~17%
  for the one-parameter comparative-judgement threshold comparison
  (the familiar P(chi-sq_1 > 2) = 15.7% AIC property); detection 95-100%
  at the stronger departures; CL-BIC selects the smaller model almost
  always under nulls with little power against mild departures. Effect
  tests: position/exposure nulls 5.8%/5.9% (800 replicates); carry-over
  8.3% at 14 judges falling to 5.3% at 30 judges (the few-cluster
  elevation round 1 flagged as a soft note, now characterised); power at
  0.6 logits 62/39/77%.
- `results/cross-package-validation.csv` also carries (post-correction
  rerun): the corrected EFRM item-set ratio anchored at −0.005 vs TAM's
  +0.001; the person-group unit phi at +0.003 vs a per-group
  `lme4::glmer` coefficient-slope anchor at −0.003 (phi's first external
  check); the crossed 2x2 design clean on both units (alpha +0.009, phi
  +0.0004 against per-cell GLMM anchors); the polytomous item side at
  +0.008 vs TAM `GPCM.design` set-constrained slopes at +0.002 (note:
  TAM's `2PL.groups` frees the category loadings on polytomous data and
  is unsound as a GPCM anchor); and `rasch_mfrm` vs `tam.mml.mfr` with
  item, rater, and threshold correlations at 0.9999+ and agreement to
  ~0.01 logits.
- `harness.R` — shared reporting helpers for the second-round studies: one
  row per scenario with bias, empirical SD, mean reported SE, SE ratio,
  95% coverage, Type I / familywise error or power, refusal and
  convergence rates, and the Monte Carlo standard error of each rate.
- `results/cross-package-validation.csv` — estimates checked against
  independent implementations on shared datasets (25 replicates per
  cell). Identical-likelihood comparators agree to solver precision:
  `btl` dichotomous vs `BradleyTerry2::BTm` to 2.6e-9 logits, polytomous
  comparative judgement vs constraint-matched `VGAM::vglm(acat)` to
  ~1e-7 with log-likelihoods equal to 5e-12. Estimator-variant
  comparators agree to the expected order: `sirt::rasch.pairwise` mean
  max-difference 0.02 logits, `eRm` full conditional ML 0.07
  (dichotomous) and 0.15 (PCM thresholds, worst 0.34 at sparse
  extremes), with equal truth RMSE. The EFRM item-set unit ratio shows
  its documented few-per-cent upward bias (+0.046 in the log at 8
  dichotomous items per set) against an essentially unbiased
  `TAM::tam.mml.2pl` slope-group anchor (+0.001); the `btl_efrm`
  panel-unit ratio is unbiased within Monte Carlo error and tracks a
  per-panel intercept-free adjacent-category anchor.
- `results/alpha-correction.csv` — the truncated-score-moment correction
  of the EFRM item-set unit linking. Diagnosis: the naive construction
  (observed variance minus mean squared standard error) under-recovers
  the true-score variance by more than half per set at 8 dichotomous
  items — the reported SE^2 overstates the actual error variance and WLE
  shrinkage makes the errors covary negatively with the true values —
  and its imperfect cancellation between sets biased the log unit ratio
  by +0.046 at 8 items and −0.036 at 15 (a sign flip: the distortions
  rebalance with length, so the old estimator was uncontrolled, not
  merely short-test-biased). Corrected estimator on the default hybrid
  path: null bias −0.002 with omnibus size 4.8% and coverage 95.3% (400
  replicates), unbiased at 5, 8, and 15 items per set, two-group frame
  grid size 3.5%, and the TAM-anchored same-seed rerun at −0.005 against
  the anchor's +0.001 (pre-correction: +0.046 on identical seeds).
- `results/alpha-correction-designs.csv` — design robustness of the
  corrected estimator: unit-ratio sweep 1.0–2.0 (bias within ±0.015
  everywhere, a mild +0.01 at ratios ≥1.5), person-skewness
  dose-response 0.5–2.8 (corrected flat at +0.008 to +0.014 while the
  TAM normal-population anchor drifts monotonically to −0.026),
  polytomous frames via `simulate_efrm(n_categories = 4)` (null size
  3.3%, power 100% at ratio 1.4), three sets linked by pairwise-only
  person overlap (size 2.7%), booklet-style within-set missingness
  exercising the per-pattern correction moments (size 3.3%), and 15
  items per set (coverage 98%).
- `results/alpha-correction-limits.csv` — deliberate stress tests mapping
  where the corrected set-unit estimator breaks, and how. Loud failures
  (informative refusals): sets whose patterns span a score range under 4
  (three dichotomous items previously collapsed to a silent −0.27 with
  near-zero variance — now guarded and refused). Silent-but-moderate:
  mistargeting by 1.5–2.5 logits ≈ +0.04 (common offsets do not cancel),
  ratio 3.5 ≈ +0.05, t3 tails ≈ +0.03; silent-and-large: widely
  bimodal persons (modes ±2.2) ≈ +0.18 — flagged in the documentation.
  Inherited model violations degrade in proportion (discrimination
  jitter −0.017/−0.035 at sd 0.25/0.5; guessing −0.017); ability-
  dependent missingness measured benign (+0.013) after fixing a crash it
  exposed. Citation rows carry the fixed-design floor curve
  (`alpha-floor-curve.R`, N=20,000): corrected +0.012 → +0.006 from 8 to
  32 items (~1/sqrt(I)), while the raw construction returns NaN at 4–6
  items and still sits at −0.02 to −0.04 at 16–32 — it never converges
  in the tested range.
- `results/cross-package-diagnostics.csv` — the diagnostics checked
  against independent implementations, at the level each comparison
  supports. Formula parity: alpha vs `psych::alpha` to 4e-15. Value
  level: infit/outfit vs `eRm::itemfit` r 0.97/0.99 once the mean-square
  convention is aligned (our E[z^2] divisor sits 8-10% above eRm's /n),
  with 25/25 same-direction flags on a planted over-discriminating item;
  PSI vs `eRm::SepRel` to 0.004 under matched conventions; q3/q3_star vs
  TAM Q3/aQ3 r 0.92 native, 0.97 at shared person estimates (top-1
  localisation of a planted dependent pair is weak for both: 2/25 vs
  8/25 — the flag rule, not top rank, is the operative criterion).
  Rank level: person fit vs PerFit lzstar/U3 |Spearman| 0.97-0.98, 87%
  worst-decile overlap, equal planted-careless detection (0.76 vs 0.77).
  Decision level: uniform DIF detection 84/88/80% (rasch/Waldtest/MH,
  BH-aligned) with comparable false-positive rates; dimensionality —
  ours holds an exact null (0 false flags) with 67% power at a balanced
  planted second dimension where DETECT flags 100% (also with clean
  nulls) and the quasi-exact Tmd/T11 tests pair full power with 13-20%
  null false-flag rates. Citation rows name the diagnostics with no
  external parallel (dependence_magnitude, spread_test, the tailored
  bootstrap, the comparative judgement family, equate_tests), which
  remain simulation-validated only.
- `results/alpha-n-sweep.csv` — the corrected set-unit estimator across
  linking samples N = 250 to 10,000 at the shipped configuration,
  uniform 100 replicates per cell (coverage resolved to ±0.02; the
  original 60-replicate cells were topped up under their fixed seeds and
  combined exactly): reported SEs track the empirical SD at every size
  (ratios 0.91–1.06) while both shrink with the square root of N; the
  bias sits flat in the +0.001 to +0.008 band with the large-N plateau
  at ~+0.004 in log alpha[set2] -- which is +0.008 on the log RATIO, so
  about 0.8% on the ratio scale; the per-set log is half the log ratio,
  an easy factor-of-two trap when quoting these figures -- instead of
  vanishing; coverage runs
  92–99% through N = 5,000 with a mild slip to 91% at N = 10,000,
  consistent with the offset-to-SE arithmetic (0.32 SE units predicts
  ~94%). The apparent sharp erosion in the initial 60-replicate pass
  (0.867 at N = 10,000) did not survive uniform replication. Operating
  consequence in the rasch_efrm help: calibrated to at least five
  thousand linking persons; beyond that, accuracy is bounded by the
  design, not the sample.
- `studies/alpha-bootstrap-pointest.R` — the answer to "could a bootstrap
  make the unit estimate better?": no. Bagging and bootstrap bias
  correction both fail to beat the shipped point estimate on RMSE (0.090
  and 0.092 against 0.0895 over 20 datasets x 30 resamples), bias
  correction being worse on bias as well because it removes the
  small-sample Jensen term that partly offsets the fixed-design floor.
  Resampling redistributes information already present; the constraint is
  the linking channel's width, so set length is the lever.
- `results/humphry-item-side.csv` — is the item-side variance-ratio
  argument (Humphry 2005, eq. 2.27) biased? Not detectably. Over 15 design
  cells (8 to 40 items, 200 to 5,000 persons per frame, planted unit ratio
  1.30, 30 replicates each) the raw SD ratio's pooled log bias is +0.00008
  with a standard error of 0.0025, and no cell reaches two standard
  errors. The attenuation mechanism is real but arithmetic puts it out of
  reach: item standard errors contribute one to two per cent of the
  observed item variance, and the two frames' error variances differ in
  the offsetting direction. Correcting is pointless here — subtracting
  `mean(se^2)`, or the covariance-correct `tr(V)/(K-1)` that accounts for
  the sum-zero constraint, moves the estimate by about 0.004. Contrast the
  person side, where the error share exceeds half and the naive
  construction biased the ratio by five per cent.
- `results/common-item-channel.csv` — when two frames share items, three
  channels estimate the unit ratio from the same data: the bilinear
  conditional fit `rasch_efrm` uses for phi, the SD ratio of two
  within-frame calibrations, and their regression slope. Conditional ML
  and the SD ratio are indistinguishable — both unbiased, with sampling
  standard deviations agreeing to three decimals across all six cells. The
  slope channel is attenuated as errors-in-variables predicts (−0.013 to
  −0.034, worsening as the item grid densifies and the true item variance
  falls), which is the control confirming the comparison measures what it
  claims. Implication for design: a bridge design putting some items in
  both set contexts would estimate alpha on an unbiased channel needing no
  correction and no new estimator, at a log-ratio SD of 0.087 against
  about 0.105 for the person-side link at equal budget, improving to 0.024
  at 40 items. Untested caveat: the study gave each frame independent
  persons; a bridge design shares persons across contexts, so the two
  calibrations' errors correlate. That caveat is settled in
  `bridge-item-design.csv`, and the bridge idea does not survive it.
- `results/bridge-item-design.csv` — a bridge design for item-set units
  requires the same person to answer an item in both set contexts, which
  means re-administration. Sharing persons is harmless on its own: with
  responses conditionally independent given theta the SD ratio stays
  unbiased and its sampling standard deviation is no worse than under
  independent persons (0.067 against 0.087 at 8 items and 250 persons),
  since correlated errors partly cancel in a ratio. Carry-over between the
  two administrations is what breaks it: half a logit inflates the unit
  ratio by about 12 per cent and a full logit by 17, in every cell, with no
  decay in sample size or item count — a person repeating their first
  answer makes the second administration look more consistent, which the
  calibration reads as a larger unit. For scale, the person-side link's
  corrected defect was 5 per cent. Since item sets are defined by item
  properties, an item belongs to one set, so a bridge means literal
  re-administration and conditional independence is not credible. The
  person-side link is therefore the only sound route to item-set units,
  and the truncated-score-moment correction was necessary rather than
  avoidable. The dependence would at least be detectable: paired
  administrations show a large Q3 in `residual_correlations`, though
  detection only tells you to abandon the bridge.
- `results/humphry-pgd-replication.csv` — a simulation of the
  person-group-discrimination design in Humphry (2005, ch. 4): 12 common
  items linking Year 5 and Year 7, calibrated separately, the unit ratio
  read off the ratio of the common items' location standard deviations.
  Planted 1.306 (his phi_5 = 0.875, phi_7 = 1.143); recovered 1.303 to
  1.311 in every cell, bias within 0.004 of zero at 200 replicates. The
  ability gap between year levels, swept from 0 to 1 logit, changes
  nothing, so differential targeting is not a threat to the estimator.
  Sampling SD is 0.04 at his N of 980 and 0.02 at 5,000, which reconciles
  his own two figures: the full-population 1.22 and the sample 1.30 differ
  by 0.064 in the log, about 1.6 SD of their difference — sampling
  variation, not a discrepancy needing explanation.
- `results/humphry-pgd-misfit.csv` — the exposure that design does carry.
  Item-level departures on the common items enter the dispersion as if
  they were unit differences: two of twelve items with a one-logit uniform
  DIF shift move the recovered ratio from 1.306 to 1.403, and four items
  with discrimination multiplied by 1.5 take it to 1.442. The direction is
  not fixed — four items with moderate alternating-sign DIF pull it down
  to 1.289, their shifts partly cancelling within the spread. Robust
  dispersion measures do not rescue it: the median absolute deviation is
  worse than the standard deviation in five of six cells (1.712 against
  1.442 under differential discrimination), because a MAD over 12 items is
  itself a poor dispersion estimator and that noise outweighs its
  contamination resistance. The discipline the method needs is screening
  the common items for DIF and misfit before computing the ratio, which is
  what his own RMSD 0.24 against RMSE 0.12 diagnostic was detecting.
- `results/channel-head-to-head.csv` — Humphry's item-side estimator and
  ours on identical simulated data (same persons, same items, two frames
  differing only in unit, conditional independence given theta). Both are
  unbiased everywhere, so the comparison is of efficiency: the item-side
  channel is 1.96 times more precise at 8 items and 980 persons, 1.46
  times at 12 items, and 1.15 to 1.33 times by 20 to 40 items. Adding
  persons sharpens item locations and helps his channel; adding items
  sharpens person estimates and helps ours. The person-side route is
  therefore costly on short sets and close to equivalent on long ones.
- `results/misfit-both-channels.csv` — ordinary item misfit, over- and
  under-discriminating items equally in both frames, as distinct from the
  frame-specific departures in `humphry-pgd-misfit.csv`. It largely
  cancels in a ratio, on both channels: four of twelve items with
  discrimination doubled or halved attenuate the recovered ratio by 2.3
  per cent (item side) or 1.0 (person side), and scattering every item's
  discrimination log-normally costs about 1 per cent. Compare +10.4 per
  cent for the same number of items when the discrimination differs
  ACROSS frames. The screening a unit ratio requires is therefore for
  frame-specific departures, not for poor fit as such.
- `results/humphry-isd-replication.csv` — his ITEM-SET discrimination
  study replicated on its own design (4 sets of 10 items spanning -4 to 4,
  N = 1000, planted ISDs 0.604/0.906/1.209/1.511). ISD is estimated
  person-side in the thesis -- "a matrix of log ratios of standard
  deviations for common persons across the sets", corrected by equation
  2.29, var(WLE) minus the mean squared standard error -- which is the
  construction this package replaced. On his design the uncorrected ratio
  is attenuated to 2.089 against a planted 2.502 end-to-end; equation 2.29
  overshoots to 2.698 (+7.8 per cent); the truncated-score-moment
  correction lands at 2.483 (-0.8 per cent). Note that the product
  constraint fixes the mean of log alpha, so mean bias is zero by
  construction for every estimator and only the SPREAD can be wrong --
  which is why the end-to-end ratio is the discriminating statistic.
  Caveat: the per-set SDs here (0.74/1.15/1.57/1.99) sit about 12 per cent
  below his Table 3.10 (0.84/1.30/1.74/2.14). The thesis is internally
  inconsistent about the person spread (his expected SDs imply 1.51, the
  text reports a generated 1.76) and he used RUMM2020's WLEs, so a design
  detail differs; the ordering of the three estimators is unaffected since
  all three run on identical data.
- `results/pgd-ours-vs-his.csv` — the like-for-like comparison on a
  common-item linking design. Given common items and disjoint person
  groups this package does not use its person-side link at all: it
  estimates the person-group unit by conditional ML on the bilinear
  threshold structure. Run against Humphry's SD ratio on his own design
  (12 WALNA common items, N = 980 per year, ability gap 0.5, planted phi
  ratio 1.306), the two are indistinguishable on clean data: 1.306
  against 1.308, sampling SD 0.040 against 0.041. Under item-level
  departures conditional ML is consistently the worse of the two, by 3 to
  5 percentage points (1.458 against 1.406 with two DIF items; 1.480
  against 1.460 with four differentially discriminating items). Weighting
  by information is the reason: an item whose discrimination is inflated
  in one frame carries more information there, so conditional ML leans on
  the items that mislead it, where a dispersion weights items by squared
  distance from the mean. The conclusion is symmetrical. On a common-item
  design his estimator is at least as good as ours, and our advantage is
  confined to designs where item sets partition the items and his channel
  does not exist. It also settles a tempting change: switching phi to a
  dispersion ratio would buy a few points under contamination that
  screening should remove anyway, at the cost of the conditional standard
  errors an SD ratio cannot provide.
- `results/alpha-set-misfit.csv` — what actually threatens an item-set
  unit, and the largest effect measured anywhere in this battery. Sets
  partition the items, so no item appears in two sets and DIF across sets
  is undefined; the hazard is misfit concentrated in ONE set, which the
  other carries nothing to offset. With 8 items per set and a planted
  ratio of 1.40: two over-discriminating items in set 1 recover 1.17, two
  under-discriminating items recover 1.73, and four over-discriminating
  items recover 1.02 — a real 40 per cent unit difference read as none.
  The same misfit spread evenly across both sets cancels almost exactly
  (1.41 against 1.42 clean). The sign follows the mechanism: an
  over-discriminating item disperses its own set's person estimates, which
  reads as a larger unit for that set. The under-discriminating case is
  the practical one — it is the wording case study's Q8, an ambivalent
  item filed among the negatively worded ones, whose removal moved that
  analysis from 1.24 to 1.03. Set membership is thus the most
  consequential decision a user makes here, it is a falsifiable hypothesis
  rather than a given, and the diagnostic is ordinary item fit within each
  set.
- `results/misfit-repair.csv` — does the diagnose-and-drop workflow put a
  planted unit ratio back? Plant misfit, read the diagnostic, drop what it
  flags with `drop_items()`, refit, and compare against dropping the
  planted items regardless of what was flagged. Planted ratio 1.40,
  N = 500 per group, 200 replicates.

  Dropping is always a sufficient cure: the oracle recovers 1.406 to 1.430
  in every cell, against clean references of 1.408 and 1.423. Nothing is
  lost by removing an item, so the binding constraint is never the repair
  — it is whether the diagnostic finds the item, and recovery tracks
  sensitivity almost one for one.

  Sensitivity ranges from 91 to 22 per cent across the departures. Items
  with DIF across person frames are found 91 per cent of the time and the
  loop closes completely (1.448 damaged, 1.413 repaired, 1.406 oracle).
  Items that merely discriminate differently across frames are found 40
  per cent of the time and the loop half closes (1.530, 1.479, 1.406).
  Under-discriminating items concentrated in one item set are found 22 per
  cent of the time, the item fit test flags nothing at all in 114 of 200
  replicates, and the repair is nearly worthless (1.694, 1.638, 1.430).

  The multiplicity adjustment is the wrong instrument for this job.
  Screening 8 to 10 items with Holm costs 20 to 60 points of sensitivity,
  and the loose screens recover more: `|infit z| > 2` lifts the
  under-discriminating cell from 1.638 to 1.486 (sensitivity 22 to 95 per
  cent), and unadjusted probabilities lift the differential-discrimination
  cell from 1.479 to 1.433 (40 to 82 per cent). The loose screen is not
  free. Where misfit is strong it over-flags: `|infit z| > 2` flags 17 per
  cent of sound items in the over-discriminating cell and `drop_items()`
  refuses 48 of 200 drops for emptying a set, so the surviving mean rests
  on 152 replicates and is not comparable with the rest.
  `fit-residual-screens.csv` later replaced this screen with the
  standardised fit residual, which detects as much without the
  over-flagging; read that entry before adopting `|infit z| > 2`.

  Screening should therefore be separated from confirming, which
  `frame_invariance()` did not allow at the time of this run: it flagged on
  Holm-adjusted probabilities only. Its `adjust` argument now offers both.

  One cell resists every screen. Four of ten items breaking invariance
  leaves 1.663 damaged against 1.410 oracle, and the best screen reaches
  only 1.602 while flagging a quarter of the sound items. Past roughly a
  fifth of the items, no amount of screening substitutes for a set of
  items that holds together.
- `results/fit-residual-screens.csv` — which standardised fit statistic
  should drive a screen. Both the frame comparison and the item-set screen
  used `infit_z`, the cube-root standardisation of a mean square; the
  package also computes the log-transformed fit residual
  `f (log y2 - log f) / sqrt(v)`. The answer differs by which comparison is
  meant, so the two were tested separately.

  Across frames, `infit_z` should stay. It detects two items discriminating
  1.8 times as steeply in 67, 95 and 100 per cent of replicates at 500,
  1,000 and 2,000 persons per frame against the fit residual's 50, 85 and
  99, at type I rates of 3.0, 4.3 and 5.5 per cent against 2.0, 2.8 and
  4.0. The extra power is largely the difference between a test at nominal
  size and a conservative one, and on ranking — whether the planted items
  are the largest departures, which is what a screen depends on — they are
  indistinguishable (80.5 against 78.5 per cent at 500 persons). An earlier
  single dataset suggested the fit residual ranked better; it does not
  generalise.

  Within a set the fit residual is the statistic to use, though the cut
  this study used does not generalise -- see
  `fit-residual-threshold-n.csv`, which supersedes the threshold advice
  below while leaving the choice between statistics standing. At the 500
  persons measured here, the question is not which test is more powerful
  but which is better calibrated at the conventional cut of 2. Against a planted 1.40 distorted to 1.69 by two
  under-discriminating items, `|fit_resid| > 2` recovered 1.439 where
  `|infit z| > 2` reached 1.486 and the chi-square test reached 1.638, with
  an oracle of 1.430. In the over-discriminating cell it recovered 1.376
  against 1.321 and 1.337, oracle 1.430. It detects as much as `infit_z`
  (88 and 93 per cent against 86 and 95) while flagging 6 and 2 per cent of
  sound items against 17 and 4, and that gap decides the repair:
  `drop_items()` refused 48 of 200 drops under `infit_z` for emptying a set
  against 4 under the fit residual.
- `results/fit-residual-threshold-n.csv` — why a fixed cut on a fit
  statistic does not survive a real sample size, and what to do instead.
  Run against the Rosenberg Self-Esteem data behind the wording case study
  (6,000 respondents, ten items), the `|fit_resid| > 2` screen recommended
  by `fit-residual-screens.csv` selects seven of the ten items, `|infit z|
  > 2` selects eight, and the chi-square test selects all ten, so
  `drop_items()` refuses every one for emptying a set. The recommendation
  did not survive contact with the data it was written for.

  The tempting explanation — real items never fit exactly while simulated
  ones do — is wrong. The cut degrades with N even when the sound items
  are generated from the model exactly: they clear it 0.8 times out of 14
  at 500 persons, 3.1 at 2,000 and 7.0 at 6,000. Two of eight items
  discriminating twice as steeply forces the fitted model to a compromise
  under which the rest genuinely depart, and that departure is fixed in
  size, so only its detectability grows. A fixed threshold on any fit
  statistic is a statement about power, not magnitude. Carried through to
  the repair it is worse than useless: the drop it implies is refused in
  61 per cent of replicates at 2,000 persons and 100 per cent at 6,000.

  Ranking survives what thresholding does not. The two planted items are
  the two largest departures in 100 per cent of replicates at 2,000
  persons and above when the others fit exactly, and 75 to 78 per cent
  when every item carries its own small slope departure. On the
  self-esteem data the ranking puts Q8 first at 22.5 and Q4 second at
  12.4, and a free-slope model fitted to the same respondents ranks the
  same two lowest — the two orderings agree on the extreme three items.
  Dropping Q8 alone moves the unit ratio from 1.24 to 1.03 at p = 0.63.

  The stopping rule has to be the unit test rather than the ratio, and this
  dataset shows why: dropping the second-ranked item as well does not
  improve on 1.03 but returns a significant difference of 1.16, still
  favouring the positive set. Once the sets no longer differ in unit there is nothing left
  for a further drop to explain, so the advice is to order by
  `abs(fit_resid)`, drop the largest, refit, and stop when the unit test
  goes quiet. `inst/casestudies/wording_units_selfesteem.R` runs exactly
  this sequence, so the figures above are reproducible rather than prose.
- `sha-map-2026-08-16.txt` — commit-ID map (old, new) from the 2026-08-16
  message-only history rewrite. `package_sha` values stamped in result
  tables before that date are pre-rewrite IDs; look them up in the first
  column to find the corresponding commit in the current history. File
  contents were untouched by the rewrite, and the `r_tree_md5` content
  hashes remain directly verifiable.

## Conventions

- True generating parameters are held fixed across replicates within a
  scenario; only responses (and persons, where stated) are redrawn.
  Redrawing the truth each replicate folds between-replicate truth
  variation into the empirical SD and fakes a calibration failure.
- Principal null-calibration claims target >= 1,000 replicates (Monte
  Carlo SE ~0.7 percentage points at a true 5% rate); power and secondary
  conditions use fewer, with the Monte Carlo error reported alongside.
- Every study records refusals and non-convergences rather than silently
  dropping them.
