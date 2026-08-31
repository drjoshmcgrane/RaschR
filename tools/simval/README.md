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

The current-estimator studies run on 2026-08-21 carry R-tree hash
`360e3609691b`: `alpha-correction-limits.csv`,
`alpha-npml-coverage.csv`, `mfrm-pooled-dif.csv`,
`tailored-bootstrap-topup.csv`, `btl-equating-clustered.csv`, and
`cross-package-validation.csv`. Their script hashes identify the exact study
files used. `btl-efrm-current.csv` and `btl-efrm-bias-sweep.csv` were rerun on
2026-08-27 after the reconciled-panel refit; both carry R-tree hash
`dc534f567649`. Their current results are described below.

The explanatory-model study run on 2026-08-23 is in
`studies/explanatory-models.R`, with results in
`results/explanatory-models.csv`. It covers LLTM and LPCM coefficients,
dichotomous and ordered comparative judgement, judge-clustered covariance,
Kent-adjusted comparisons with free calibrations, and Holm-adjusted fixed
departure diagnostics. Each principal condition used 1,000 replicates; the
diagnostic conditions used 300. The result rows carry script hash
`e72c33409a15c6ecc04bc5fb91f413ca` and R-tree hash `dd5f1154d99f`.

The edge-case extension is in `studies/explanatory-edge-cases.R`, with results
in `results/explanatory-edge-cases.csv`. Four LLTM/LPCM conditions cover 300
to 2,000 persons, mixed maximum scores, and four-category items. Across 1,000
replicates per condition, coefficient bias was at most 0.0049 logits,
empirical SD/mean SE was 0.993--1.026, coverage was 0.942--0.954, and
Kent-adjusted null rejection was 4.3--5.8%. The unscaled probability rejected
98.6--100% and is retained only as `p_naive`. No fit was refused or failed to
converge. Mean calibration R-squared was 0.872--0.994 under the correctly
specified generating models. The rows carry script hash
`42edd1a6a6fda75a19006c01953c2dae` and R-tree hash `954bcba447fe`.

The compiled EFRM linking kernel was checked against the retained R
implementation with the same seed and 30 hybrid bootstrap replicates. Across
set units, their standard errors, origins, edge likelihoods and thresholds,
the largest absolute difference was 1.30e-11; every convergence flag agreed.
The study is `studies/efrm-cpp-parity.R` and its six result rows are in
`results/efrm-cpp-parity.csv`. They carry script hash
`1a6468c78f37e505d3c01870f8a1c693` and R-tree hash `2ff015352a0d`.

Parallel EFRM bootstrap execution was then checked from a fresh isolated
installation. Serial, two-worker and default four-worker hybrid fits used the
same 300 pre-generated resamples; the full-bootstrap comparison used 30. Every
checked estimate and convergence flag was identical. On the executing machine,
two and four workers reduced the hybrid elapsed time from 86.0 seconds to 52.4
and 33.6 seconds; two workers reduced the full-bootstrap time from 8.6 to 6.4
seconds. These timings are contextual, not general performance guarantees.
The study is `studies/efrm-parallel-parity.R`, with results in
`results/efrm-parallel-parity.csv`; the rows carry script hash
`64574bd4400ad20f64b2746b0b1453ed` and R-tree hash `f25f578582a2`.

The BTL-EFRM judge bootstrap was checked separately from a fresh isolated
installation. Serial and default four-worker fits used the same 200 judge
resamples and agreed exactly after excluding the recorded worker count. On the
executing machine, elapsed time fell from 17.34 to 5.47 seconds (3.17 times
faster). The timing is machine-specific. The study is
`studies/btl-efrm-parallel-parity.R`, with results in
`results/btl-efrm-parallel-parity.csv`; the row carries script hash
`db8d37e3a52d75b7a342fa6e43d7f8f8` and R-tree hash `eba5f653300c`.

`results/item-fit-bootstrap.csv` records the 2026-08-30 execution of the
item-fit study. Its performance rows are conditional on analysed replicates,
but that execution combined refusals and other failures in the
`n_nonconv` field and did not retain inner-bootstrap counts. The current
`studies/item-fit-bootstrap.R` distinguishes refusals, non-convergence and
other errors, carries `B`, `B_used`, `B_nonconverged` and `B_errors`, and
treats an entirely unavailable adjusted family as unavailable rather than as
no rejection. The existing CSV retains its original script hash and is not
presented as a rerun of that accounting revision.

The score-conditional person-fit and fitted-design comparative-judgement
bootstraps are checked in `studies/person-cj-fit-bootstrap.R`. The supported
four-category follow-ups are `studies/person-fit-bootstrap-pcm-topup.R` and
`studies/cj-fit-bootstrap-polytomous-topup.R`; their result files have the
same stems under `results/`. Dichotomous person-fit familywise error was 5.0%.
The 500-person PCM top-up gave 4.25% marginal error and 7.07% familywise error
(MCSE 1.9 points) among 184/200 analysed fits. Fourteen datasets had fewer than
half of their refits estimable, one observed calibration was singular, and one
complete adjusted family was unavailable. None of 39,601 inner refits failed
to converge; 5,952 otherwise failed. The 200-dataset polytomous CJ top-up gave 3.0% total-test,
2.0% pair-family, 7.0% object-family and 4.5% judge-family error, with no
refusals; all 39,800 inner refits were usable. All three studies carry R-tree
hash `ed5c283cbe90`; their script hashes are recorded in the result rows and
match the study files.

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
- `results/efrm-fix-sweep.csv` and `results/efrm-fix/` — covariance checks
  for the superseded moment-based set link. They established the need to
  propagate within-frame calibration uncertainty, which the current
  semiparametric link retains, but they do not validate its point estimator.
- `results/round2-followups.csv` — the adjudication trail for every
  round-2 suspect: the MFRM q=25 exoneration (600 fixed-truth replicates
  per cell), the btl_efrm origin-test correction (8.5% pre-fix to 5.5%
  post-fix at 400 replicates each), the PCM item-level guard validation,
  and the concentration guard's field behaviour.
- `results/audit3-fixes.csv` — checks the corrections from the third audit.
  Under a 10:90 nuisance-cell imbalance, the repeated-measures DIF follow-up
  held its 5% size (5.25%, 2,000 replicates) when it retained the equal-cell
  estimand; the superseded person-frequency shortcut rejected every dataset
  because it targeted a different contrast. The corresponding mixed
  time-by-group interaction had 5.4% size. BTL-DIF pairwise size was 5.5%
  with eight judges per level and 4.83% with ten. In the five-level design
  with four judges per level, the superseded global degrees of freedom gave
  9.45% pairwise rejection; current public inference is withheld in that
  region. At the binomial spread boundary, the raw estimate fell below the
  bound in 47.1% of datasets, while the current one-sided test rejected 5.2%
  (1,000 replicates). The planted-dependence condition had 100% power.
- `results/audit3-btl-imbalance-topup.csv` — fresh-seed adjudication of the
  mildly imbalanced BTL-DIF cell. With ten raw and 9.31 effective judges per
  level, the current pair-specific reference rejected exactly 5.0% over 2,000
  datasets. A diagnostic minimum-cell reference rejected 3.55% and was not
  adopted for the two-cell comparison. There were no refused or non-converged
  fits.
- `results/btl-dif-multicell-df.csv` — the separate four-cell case. In 500
  balanced 2 by 2 null datasets with 12 judges per cell, the weakest-cell
  reference rejected 3.6% and gave 0.964 interval coverage. The superseded
  pooled-count extension rejected 6.2% and covered 0.938. Multi-cell
  contrasts therefore use the conservative weakest-cell rule while the
  validated two-cell Welch rule is unchanged.
- `results/btl-dif-hc3.csv` — null calibration of the between-judge residual
  test under a fourfold variance ratio. With 8 versus 16 judges, the classical
  equal-variance test rejected 11.72% or 1.78%, depending on which group was
  less precise; HC3 gave 5.57% and 4.04%. With eight judges per group, HC3
  gave 4.13% (10,000 replicates per condition).
- `results/dif-hc3.csv` — ordinary between-person DIF under balanced groups,
  a 1:4 ability imbalance, and unequal observations per person. The adopted
  hybrid (HC3 for uniform terms, residual ANOVA for class-interval
  interactions) gave 4.0%, 6.4%, and 4.6% Holm familywise rejection. Applying
  HC3 to every term gave 22.0% under the ability imbalance and was rejected.
  At a planted 0.6-logit shift, hybrid power was 29.6% against 23.6% for the
  classical analysis (500 replicates per condition).
- `results/dif-hc3-multilevel.csv` — a three-level factor with group sizes in
  the ratio 1:2:3 and different group locations. Hybrid Holm familywise
  rejection was 6.2%, compared with 6.4% for the classical analysis and 20.2%
  when HC3 was also applied to the class-interval interaction. At a planted
  0.6-logit shift, hybrid power was 36.6%, compared with 32.8% classically
  (500 replicates per condition).
- `results/dif-hc3-homoskedastic.csv` — balanced group-by-interval cells with
  independent normal errors and a common variance. Hybrid Holm familywise
  rejection remained between 4.3% and 5.2%. With ten observations per cell,
  HC3 reduced planted-item power from 21.7% to 18.8% for a two-level factor
  and from 23.5% to 21.9% for a three-level factor. The difference was 0.9
  percentage points with 30 observations per cell; fixed-effect power then
  approached its ceiling (5,000 replicates per condition).
- `results/dif-hc3-homoskedastic-local-power.csv` — the same comparison with
  effects reduced as cell size increased, preventing the larger designs from
  reaching the power ceiling. For two levels, the classical power advantage
  declined from 3.10 percentage points at ten observations per cell to 0.84,
  0.42, and 0.16 points at 30, 75, and 150. For three levels it declined from
  1.72 points at ten per cell to 0.54 at 50. This confirms a real but diminishing
  efficiency cost for HC3 when the classical assumptions hold exactly (5,000
  paired replicates per condition).
- `results/dif-conditional-bootstrap.csv` — a conditional Rasch null bootstrap
  that preserves each raw score and refits the model 199 times per dataset. In
  the balanced null, item-wise rejection was 4.79% versus 4.96% for uniform DIF
  and 4.04% versus 4.00% for non-uniform DIF (current hybrid versus bootstrap);
  familywise rejection was 6.33% for both. With 1:4 group sizes and a
  0.8-logit ability difference, the corresponding rates were 4.75% versus
  4.42%, 4.46% versus 4.08%, and 7.0% versus 6.0%. The bootstrap reduced power
  from 29.7% to 24.3% for a 0.6-logit uniform shift. Adjusted power against a
  centred 0.7 slope departure was low for both procedures (5.7% versus 4.3%).
  These results do not support replacing the current hybrid (300 datasets per
  condition).
- `results/dif-conditional-bootstrap-extended.csv` — the same comparison for
  four-category PCM and RSM data, a three-level group, and two correlated
  person factors. For a response vector \(x\) with raw score \(r\), the
  polytomous sampler draws from
  \(P(X=x\mid r) \propto \exp\{-\sum_i\sum_{k=1}^{x_i}\tau_{ik}\}\), and every
  draw is checked against the conditioned score before the model is refitted.
  Global-null rejection was broadly consistent with 5% in all designs. The
  bootstrap was usually a little more conservative and a little less
  powerful than the hybrid analysis. It did not remove artificial flags on
  invariant items when another item truly had DIF (100 datasets and 99
  bootstrap refits per condition).
- `results/dif-conditional-bootstrap-confirm.csv` — fresh-seed confirmation
  with 200 datasets and 199 bootstrap refits. Under the imbalanced global
  null, Holm familywise rejection was 3.0% versus 2.5% for the PCM and 4.0%
  versus 2.5% for the RSM (hybrid versus bootstrap). With a 1.4 slope
  departure on one item, false flags among the other five items occurred in
  19.0% versus 13.5% of PCM datasets and 14.5% versus 11.0% of RSM datasets.
  The bootstrap attenuates score contamination but does not solve it.
- `results/dif-score-purification.csv` — an initial comparison of leave-one-out,
  fixed-anchor, and iterative matching scores. Leave-one-out testing was
  rejected because null familywise error reached 28--34%. Re-estimating the
  full analysis from a five-item anchor scale also lost too much uniform-DIF
  power. This screen motivated the staged comparison below (100 datasets per
  condition).
- `results/dif-score-purification-refined.csv` — 500-replicate comparison of
  anchor-based class intervals, full anchor recalibration, strongest-item
  iteration, and the public split-and-refit workflow. The unmodified hybrid
  gave global-null Holm familywise error of 4.0% for the PCM, 3.8% for the RSM,
  and 6.0% with two correlated person factors. Preselecting a five-item anchor
  scale was liberal (7.2--10.6%) and is not a valid default.

  For a uniform 0.6-logit shift, `resolve_dif()` split the planted item in
  84.4% of PCM and 86.6% of RSM datasets, against initial detection of 85.0%
  and 87.4%. It split an invariant item in 0.6% and 1.8%, and the final
  invariant-item familywise rates were 4.0% and 5.6%. The existing
  split-and-refit procedure therefore supplies an effective purification step
  for uniform DIF.

  For a centred 1.4 slope departure, correct-term power was 89.0% for the PCM
  and 86.8% for the RSM, while familywise flags among invariant items rose to
  14.8% and 16.8%. A strongest-item, one-at-a-time procedure selected the
  planted item first in 97.6% and 96.0% of datasets; it selected an invariant
  item first in 0.2% and 1.0%, and ever excluded one in 3.4% and 5.2%. After
  anchor recalibration, remaining false flags occurred in 1.6% and 2.0%.
  Retesting the selected item on the short anchor scale needlessly reduced
  power. With two correlated person factors, non-target factor error remained
  controlled (1.6%) but correct non-uniform power was only 12.0%. This is a
  power limit, not a calibration defect. No method from this study has been
  installed as an automatic non-uniform-DIF remedy.
- `results/item-fit-hc3.csv` — sensitivity study for class-interval item fit.
  HC3 was rejected: item-wise null rejection ranged from 21.9% to 48.3% over
  8--30 items, against 5.6--17.0% for the conventional ANOVA. The
  conventional ANOVA remained approximate (Holm familywise rejection 7.5%
  at 30 items and 11.0--31.5% at 8--15 items). The item-trait test was
  calibrated from ten items onward (4.0--7.0% familywise) but not with eight
  items (12.0--17.0%). These results support the short-test qualification in
  `?rasch`, not an HC3 replacement (200 replicates per condition).
- `results/item-fit-interval-count.csv` — reducing the requested number of
  class intervals did not repair short-test calibration and generally reduced
  power. The interval-count change was therefore rejected (100 replicates per
  condition).
- `results/btl-cluster-jackknife.csv` — CR1 versus delete-one-judge covariance
  for the core BTL fit. CR1 Type I was 5.4% with ten balanced judges and 4.2%
  with one of twenty judges carrying 20% of the work; jackknife rates were
  5.6% and 4.0%. The concentrated design below the public effective-judge
  guard gave 6.4% for CR1 and 5.4% for the jackknife, both within Monte Carlo
  uncertainty of 5%. CR1 therefore remains the default (500 replicates per
  condition).
- `results/btl-equating-clustered.csv` — common-object drift under two
  independent 12-judge panels. Welch--Satterthwaite probabilities with Holm
  adjustment gave 4.4% familywise rejection over 1,000 null replicates; the
  superseded normal reference gave 7.4% on the same fitted samples.
- `results/btl-btm-agreement.csv` — `btl()` against `sirt::btm` on shared
  dichotomous comparisons with eps = 0 and fixed home advantage: identical
  likelihood, mean max-difference 1.1e-15 and worst 3.1e-15 over 23 clean
  replicates. Replicates with a fully extreme object are characterised
  separately: btl sets the unidentified object aside and reports an
  extrapolated boundary location (score half a point inside the boundary,
  SE withheld), whereas btm keeps it, diverging at eps = 0 and shrinking
  it under its default eps = 0.3. The two policies agreed in direction in
  every case, differing by 1.6 logits on average in the flat region of
  the likelihood — a policy difference, not an estimation difference.
- `results/equating-holm-refresh.csv` — the null familywise cells of the
  equating study re-run after `equate_tests()` moved from BH to Holm
  adjustment, mirroring the original design: 4.8-5.0% at 3, 5, and 10
  anchors (2,000 replicates each, no refusals). The null rows of the
  round-2 equating-multiplicity table are BH-era and superseded by this
  table; its drift-power rows also predate the switch and read at most
  slightly high for the current function.
- `results/explanatory-r2-adjusted.csv` — sampling behaviour of the
  calibration R-squared reported by `explanatory_test()`. Under
  uninformative designs the raw coefficient averaged 0.170 (12 items) and
  0.085 (24 items) while the rank-based adjusted coefficient centred at
  -0.015 and -0.002; a true 12-item design gave 0.948 raw and 0.936
adjusted (300 replicates per condition).

- `results/person-external-weights.csv` — externally weighted person
  estimation with fixed generating calibrations. The 18 conditions cover
  equal, moderate, strong and zero weights; dichotomous and partial credit
  items; three person locations; and differing model units. Across 5,000
  persons per condition, absolute bias was at most 0.016 logits,
  empirical SD/mean reported SE was 0.940--1.004, and 95% coverage was
  0.941--0.978. No estimate was refused. The study is
  `studies/person-external-weights.R`.
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
  8.3% at 14 judges falling to 5.3% at 30 judges, which is why its
  probability is now withheld below 30 judges; power at 0.6 logits was
  62/39/77%.
  The EFRM log set-unit bias is +0.0036 against TAM's +0.0008 for
  dichotomous data and +0.0035 against +0.0020 for polytomous data. In a
  crossed two-set by two-group design, EFRM log-alpha bias is +0.0141 and
  log-phi bias +0.0004. The person-group unit is +0.003 against a per-group
  `lme4::glmer` coefficient-slope anchor at −0.003. `rasch_mfrm` and
  `tam.mml.mfr` item and rater estimates correlate above 0.9999.
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
  extremes), with equal truth RMSE. Current EFRM comparisons are summarised
  above. The `btl_efrm` panel-unit ratio is unbiased within Monte Carlo error
  and tracks a
  per-panel intercept-free adjacent-category anchor.
- `results/alpha-correction.csv`, `alpha-correction-designs.csv`,
  `alpha-n-sweep.csv` and `studies/alpha-bootstrap-pointest.R` evaluate the
  superseded moment-based set link. They are retained to document why it was
  replaced and must not be cited as evidence for the current estimator. The
  scratch script that produced `alpha-n-sweep.csv` was not retained, so that
  file is a historical, non-reproducible record; its recorded path and hash
  identify the missing script rather than an executable repository study.
- `results/alpha-correction-limits.csv` — point-estimator stress tests for the
  current finite-grid semiparametric link. The study covers targeting, unit
  ratios, short sets, small samples, heavy-tailed and bimodal populations,
  missingness, guessing and within-set discrimination departures. Sets with
  fewer than four score steps are refused by design. Absolute bias is at most
  0.022 under the model. At 80 persons, 11% of datasets are refused and 2%
  do not converge; the 41- and 101-point grid results are effectively equal.
- `results/alpha-npml-coverage.csv` — sampling calibration of the current
  estimator under normal, wide bimodal and deliberately different group
  distributions. Raw marginal hybrid set-unit Type I is 4.0--5%, SE ratios
  are 0.97--1.05, and coverage is 0.927--0.960. Common-scale item SE ratios are
  0.97--1.04. The complete bootstrap is mildly conservative under the null
  (2.5% rejection, 0.975 coverage) and calibrated under the planted ratio
  (SE ratio 0.99, coverage 0.938).
- `results/efrm-unit-multiplicity-supported.csv` — complete-family null
  calibration for the EFRM decision policy, with 200 persons per group and
  six items per set. Of 500 attempted fits, 489 were analysed and 11 were
  refused. At least one raw omnibus probability was below 0.05 in 10.4% of
  analysed fits. Holm familywise rejection was 3.5% for the omnibus family
  and 1.6% for the separate individual-contrast family. The script is
  `studies/efrm-unit-multiplicity-supported.R`; the rows carry script hash
  `6341acc22941fb7041fb08ec8f8bef0c` and R-tree hash `7356d398664f`.
- `results/frame-unit-multiplicity.csv` — the same complete-family check for
  BTL-EFRM, and an EFRM boundary design. The BTL-EFRM rows predate the
  reconciled-panel refit and are superseded; they are retained as provenance,
  not as evidence for the current estimator. The four-item,
  100-person-per-group EFRM boundary attempted 1,000 fits: 493
  were refused, 24 did not converge, and 483 were analysed. Conditional on
  analysis, its Holm rates were 5.2% and 2.5%. The script is
  `studies/frame-unit-multiplicity.R`; the rows carry script hash
  `45abbc7d7938f9fc2416d2dbbf356eeb` and the same R-tree hash.
- `results/audit-adjusted-dependence.csv` — null familywise error and power
  for the crossed-EFRM and BTL dependence decisions, the finite-object
  correlation check, and the affected dimensionality-power cell. Crossed
  EFRM rejected 6.7% in the first 1,000 null fits; the independent top-up in
  `results/crossed-efrm-factorial-topup.csv` gave 5.55% over 2,000 fresh fits,
  with marginal rates 5.25--5.75% and no refusals or non-convergence. The
  pooled familywise rate is 5.93% over 3,000 fits. BTL dependence familywise
  rejection was 5.9% over 1,000 fits. The simulator reproduced requested
  finite-object correlations to 3.9e-16 and dimensionality power was 87%.
  The scripts are `studies/audit-adjusted-dependence.R` and
  `studies/crossed-efrm-factorial-topup.R`; their hashes are
  `3f1af2fcb6098b5996bda3e0987763f6` and
  `70e9130564c43819760d330e925f53b1`, and both result sets carry R-tree hash
  `8a6cc825b06a`. The final tree differs only in the subsequent compatibility
  change that permits an empty optional simulator list; every generating and
  fitting path used by these studies is unchanged.
- `results/btl-efrm-current.csv` — current judge- and independent-outcome
  bootstrap calibration for BTL-EFRM, rerun after the reconciled-panel refit.
  Over 300 null fits, raw marginal judge-bootstrap Type I was 3.7% for panel
  units, 7.3% for set units and 5.0% for origins. The corresponding
  independent-outcome rates were 4.3%, 6.7% and 3.7%. Set-unit coverage was
  0.890 with the judge bootstrap and 0.933 with the independent-outcome
  bootstrap; coverage for the other units was 0.923--0.963. This design has
  six judges per panel and lies in the documented caution band.
- `results/btl-efrm-multiplicity-current.csv` — 1,000 current-estimator null
  fits in that caution-band design. Raw marginal Type I was 3.3% for panel
  units, 6.7% for set units and 6.0% for origins. Holm familywise error was
  3.9% across the three omnibus decisions and 3.0% across the individual
  follow-ups. Set-unit coverage was 0.900; no fit was refused or failed to
  converge. The script is `studies/btl-efrm-multiplicity-current.R`; its rows
  carry script hash `a3635aff9d2c9ed79f437686651f3dc0` and R-tree hash
  `8b2530afb990`.
- `results/btl-efrm-supported-topup.csv` — the supported-design top-up with
  12 judges per panel and the public default of 200 resamples. Across 500
  null fits, raw set-unit Type I was 4.6%, the empirical-SD/mean-SE ratio was
  0.992 and coverage was 0.934. The Holm-adjusted set-unit rate was 2.0% in
  the omnibus family and 1.6% in the follow-up family. There were no refusals
  or non-convergences. The script is
  `studies/btl-efrm-supported-topup.R`; its rows carry script hash
  `1f58bf65d9a4d192773c296596a78135` and the same R-tree hash.
- `results/btl-efrm-bias-sweep.csv` — finite-sample attenuation from the
  staged BTL-EFRM set link, rerun after the reconciled-panel refit. Log
  set-unit bias decreases from −0.106 at 10 repetitions per pair to −0.041
  at 20, −0.016 at 50 and −0.006 at 100 (500 datasets per cell). The
  caution-band bootstrap study above shows that this finite-sample
  attenuation can also reduce Wald coverage.
- `results/coherence-fixes.csv` — direct checks of the repaired BTL-EFRM
  fit and the multifactor DIF estimands. Across three information levels,
  reported BTL-EFRM log likelihoods agreed with likelihoods reconstructed
  from every stored fitted probability to numerical precision. Mean object
  RMSE fell from 0.401 to 0.232 and 0.123 as repetitions increased. In a
  correlated 3:1:1:3 two-factor design, ordinary DIF marginal magnitudes had
  biases −0.021 and +0.017 logits, and paired-comparison DIF biases +0.026
  and +0.032; 95% coverage was 0.92--0.96. The study is
  `studies/coherence-fixes.R`; its 29 rows carry script hash
  `25e3f34eddee8bb9d71f77d3c79e8a31` and R-tree hash `dc534f567649`.
- `results/mfrm-pooled-dif.csv` — null calibration after putting uniform and
  non-uniform item tests in the one multiplicity family used for decisions.
  Familywise rejection is 4.7% with balanced raters and 3.8% when the second
  group has only two raters (1,000 attempted datasets per cell).
- `results/tailored-bootstrap-topup.csv` — full automatic anchor selection
  repeated inside each of 399 person-bootstrap draws. Clean-item familywise
  error is 2.5% at guessing 0.15 and 0% at 0.30. At least one of two planted
  hard items is detected in 17.5% and 26.3% of datasets respectively, showing
  that calibration is conservative and power is limited for this design.
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
  calibration reads as a larger unit. Since item sets are defined by item
  properties, an item belongs to one set, so a bridge means literal
  re-administration and conditional independence is not credible. The
  person-side link is therefore the practical route to item-set units. The
  dependence would at least be detectable: paired
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
- `results/channel-head-to-head.csv` — a historical comparison of Humphry's
  item-side estimator and the superseded score-moment person-side link on
  identical simulated data (same persons, same items, two frames differing
  only in unit, conditional independence given theta). Both are unbiased
  everywhere, so the comparison is of efficiency: the item-side
  channel is 1.96 times more precise at 8 items and 980 persons, 1.46
  times at 12 items, and 1.15 to 1.33 times by 20 to 40 items. Adding
  persons sharpens item locations and helps the item-side channel; adding
  items sharpens person estimates and helps the score-moment channel. This
  comparison does not validate the current semiparametric link.
- `results/misfit-both-channels.csv` — a historical comparison of ordinary
  item misfit in the item-side and superseded score-moment channels. Four of
  twelve items with discrimination doubled or halved attenuate the recovered
  ratio by 2.3 per cent (item side) or 1.0 per cent (score-moment side), and
  scattering every item's discrimination log-normally costs about 1 per cent.
  It does not validate the current semiparametric link; the corresponding
  current-estimator departures are in `alpha-correction-limits.csv`.
- `results/humphry-isd-replication.csv` — historical comparison of the
  earlier item-set estimators on Humphry's item-set discrimination
  study replicated on its own design (4 sets of 10 items spanning -4 to 4,
  N = 1000, planted ISDs 0.604/0.906/1.209/1.511). ISD is estimated
  person-side in the thesis -- "a matrix of log ratios of standard
  deviations for common persons across the sets", corrected by equation
  2.29, var(WLE) minus the mean squared standard error -- which is the
  construction this package replaced. On his design the uncorrected ratio
  is attenuated to 2.089 against a planted 2.502 end-to-end; equation 2.29
  overshoots to 2.698 (+7.8 per cent); the superseded score-moment
  correction lands at 2.483 (-0.8 per cent). Note that the product
  constraint fixes the mean of log alpha, so mean bias is zero by
  construction for every estimator and only the SPREAD can be wrong --
  which is why the end-to-end ratio is the discriminating statistic.
  Caveat: the per-set SDs here (0.74/1.15/1.57/1.99) sit about 12 per cent
  below his Table 3.10 (0.84/1.30/1.74/2.14). The thesis is internally
  inconsistent about the person spread (his expected SDs imply 1.51, the
  text reports a generated 1.76) and he used RUMM2020's WLEs, so a design
  detail differs; the ordering of the three estimators is unaffected since
  all three run on identical data. These results do not validate the current
  semiparametric estimator.
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
  lost by removing an item, so the binding constraint is never the repair —
  it is whether the diagnostic finds the item. Recovery is not a simple
  function of sensitivity across the whole table, because the two damage
  directions push the ratio opposite ways and cancel when pooled; read the
  cells within a direction, where the cell with the lowest sensitivity is
  also the cell whose repair leaves the most damage behind.

  Sensitivity ranges from 91 to 22 per cent across the departures. Items
  with DIF across person frames are found 91 per cent of the time and the
  loop closes completely (1.448 damaged, 1.413 repaired, 1.406 oracle).
  Items that merely discriminate differently across frames are found 40
  per cent of the time and the loop half closes (1.530, 1.479, 1.406).
  Under-discriminating items concentrated in one item set are found 22 per
  cent of the time, the item fit test flags nothing at all in 114 of 200
  replicates, and the repair is nearly worthless (1.694, 1.638, 1.430).

  The multiplicity adjustment is the wrong instrument for this job.
  Screening ten items with Holm costs 8 to 42 points of sensitivity,
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
  only 1.602 while flagging a quarter of the sound items. Two of ten and
  four of ten are the only contamination levels this study runs, so where
  between them screening stops substituting for a coherent item set is not
  established here -- only that at two of ten it substitutes partially and
  at four of ten it does not.
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
- `results/resolve-versus-drop.csv` — when an item breaks frame invariance,
  is it better removed or resolved? `misfit-repair.csv` showed that dropping
  a flagged item restores a planted unit ratio; it did not ask what the
  repair costs. Dropping takes the item out of every frame, so it stops
  contributing to any person's measure, including in the frames it behaved
  well in. Resolving gives it a location per frame: it stops linking the
  frames, which is what the diagnosis found wrong with it, and goes on
  measuring the person inside their own frame.

  One item shifted 1.2 logits in frame 2, twelve items, 500 persons per
  frame, planted group-unit ratio 1.40, 200 replicates.

  On the unit the two remedies are indistinguishable: 1.400 dropped against
  1.401 resolved, from 1.386 damaged and 1.399 clean. Both are unbiased, and
  the damage itself is modest — one differentially functioning item in twelve
  moves the ratio by about one per cent, so the repair is not what
  distinguishes them.

  The cost does. Dropping raises the mean person standard error from 0.7268
  to 0.7592, a 4.5 per cent loss of precision paid by every respondent, and
  lowers the correlation between the person estimates and the locations that
  generated them from 0.8501 to 0.8374. Resolving leaves both where the
  clean analysis had them, 0.7267 and 0.8498 — indistinguishable from never
  having had the problem. It buys this with one parameter per extra frame
  and by removing the item from the link, so the group units then rest on
  the items that remain common; `resolve_frames()` refuses when the remaining
  frame graph or information matrix no longer identifies those units.

  So the two remedies are not a trade-off on this evidence: where an item
  measures well within each frame and only its comparability fails,
  resolving dominates. Dropping earns its place where the item is a poor
  measure wherever it appears, which is a different diagnosis than the one
  `frame_invariance()` makes.
- `results/chained-linking.csv` — does a unit ratio recover when the two
  frames share no items at all? A vertical design rarely gives every pair of
  year levels a common block: year 3 and year 5 share one anchor, year 5 and
  year 7 share a different one, and years 3 and 7 share nothing. The model
  accepts this, because identification of the person-group units runs over
  the CONNECTED COMPONENTS of the graph whose edges are group pairs sharing
  at least two items within one set. A chain is enough. Breaking it is
  refused, naming the components: "relative units (phi) are unidentified
  between: year3+year5 | year7".

  Accepting a design is not recovering from it, and the year 3 to year 7
  comparison is made entirely through year 5. It recovers anyway. Against
  planted ratios of 1.250 directly linked and 1.562 chained, at 300, 700 and
  2,000 persons per year the direct ratio came back at 1.238, 1.262 and
  1.249, and the chained one at 1.536, 1.573 and 1.567. No bias worth the
  name in either, at any of the three sizes.

  What the chain costs is precision, not accuracy. The chained ratio's
  empirical standard deviation runs about 1.4 times the direct one at every
  sample size — 0.171 against 0.119, 0.115 against 0.080, and 0.068 against
  0.051. Contrast standard errors use \eqn{c^T V c}. They track the empirical
  spread closely from 700 persons per year; at 300, the indirect-link mean SE
  is conservative (0.374 against 0.171).
- `results/frame-invariance-conditional-topup.csv` — the 2,000-replicate
  null check that led to conditional discrimination probabilities being
  withdrawn. The standardised-infit comparison produced 7.1% combined Holm
  familywise error, with rejection strongly dependent on item position. The
  descriptive infit and fitted-slope columns remain useful, but they are not
  treated as tests.
- `results/frame-invariance-bootstrap.csv` — calibration of the replacement
  bootstrap inference at 500 persons per frame. Across 300 null replicates,
  empirical-SD/mean-SE ratios were 1.002 for locations and 1.034 for log
  discrimination ratios, coverage was 0.948 and 0.954, and the combined Holm
  familywise error was 3.0%. At two planted items, power was 96.3% for a
  one-logit location shift and 9.6% for a 1.5-fold discrimination change over
  120 replicates. The latter comparison is valid but weak at this design.
- `results/frame-invariance-power.csv` — the conditional location study,
  with 1,000 null and 400 departure replicates at 500, 1,000 and 2,000
  persons per frame. Empirical-SD/mean-SE ratios were 0.908--0.916, coverage
  was 0.966--0.968, and Holm familywise error was 1.9--2.6%. Holm-adjusted
  power for each of two one-logit shifts was 95.1% at 500 persons and 100% at
  1,000 and 2,000. Because separate frame origins centre the common items,
  unshifted items can carry non-zero relative contrasts when a few items
  move; the study records that rate rather than treating it as an ordinary
  false-positive rate.
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
