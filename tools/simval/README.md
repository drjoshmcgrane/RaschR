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
