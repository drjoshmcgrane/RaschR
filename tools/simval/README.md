# Simulation validation battery

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
  frame set-unit correction: pre-fix 9.4% null rejection (1,200
  replicates), the corrected hybrid at 4.9% across nine design cells, the
  seed-paired marginal-vs-joint draw comparison (identical decisions in
  all 400 replicates), and the full-bootstrap benchmark arm (ratio 0.963)
  showing the corrected hybrid prices the same uncertainty.
- `results/round2-followups.csv` — the adjudication trail for every
  round-2 suspect: the MFRM q=25 exoneration (600 fixed-truth replicates
  per cell), the btl_efrm origin-test correction (8.5% pre-fix to 5.5%
  post-fix at 400 replicates each), the PCM item-level guard validation,
  and the concentration guard's field behaviour.
- `harness.R` — shared reporting helpers for the second-round studies: one
  row per scenario with bias, empirical SD, mean reported SE, SE ratio,
  95% coverage, Type I / familywise error or power, refusal and
  convergence rates, and the Monte Carlo standard error of each rate.

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
