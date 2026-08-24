# CRAN comments for rasch 1.12.0

## Summary of this update

This is the first update since rasch 1.11.7 was accepted on 2026-07-30.

The release extends differential item functioning analysis to multiple
between-person and within-person factors, with factorial terms and the
corresponding repeated-measures error structure. It adds explanatory
(linear logistic) variants of the Rasch and comparative judgement models,
item-invariance testing and repair for the frame-of-reference models, and
an interface to the WrightMap package (in Suggests, used conditionally).

This is the first release containing compiled code: two inner loops of the
extended-frame estimator are implemented in C++ via Rcpp (LinkingTo: Rcpp),
with no system requirements beyond a C++ compiler.

Corrections concern standard errors, reference distributions, and the
conditions under which inferential results are reported. Familywise
multiplicity control was standardised on Holm across the package. The
estimators now apply explicit checks for connectedness, separation, rank,
sparse categories, and weak identification, and withhold probabilities
rather than report unstable ones.

The package title, description, help pages, eight vignettes, and graphical
interface have also been revised for this release.

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.6.1
* win-builder: R-release and R-devel
* GitHub Actions: macOS, Windows and Ubuntu; R devel, release and oldrel

## R CMD check results

0 errors | 0 warnings | 0 notes locally; win-builder adds the
incoming-feasibility NOTE only.

This is a resubmission. The previous submission was archived by the
pretest for overall check time on Windows (32 minutes). The test suite
now runs a core subset on CRAN that exercises every estimator and
inference path once; the complete suite runs whenever NOT_CRAN is true,
locally and in continuous integration on three operating systems. The
heaviest vignette was also reduced. The full check now completes well
inside the pretest budget.

The two words flagged by the incoming spell check, 'Ponocny' and 'Tutz',
are author surnames from references in the Description.

Some check services report the BugReports URL,
https://github.com/drjoshmcgrane/rasch/issues, as possibly invalid. The
issue tracker is enabled and the repository is public: the GitHub API
reports `has_issues: true` for it, and the repository root and
documentation site both resolve normally. GitHub returns 404 to
unauthenticated non-browser requests for the /issues and /pulls paths
alike, and pull requests cannot be disabled on a public repository, so
the status reflects that behaviour rather than a broken link.

The package was built from source with its eight vignettes before checking.
Simulation-intensive and bootstrap-calibration tests use `skip_on_cran()`;
they run locally and in continuous integration.

## Current CRAN status

Version 1.11.7 is currently OK on all CRAN check flavours.

## Reverse dependencies

There are no known reverse dependencies.
