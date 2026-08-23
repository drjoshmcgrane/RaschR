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

* local: macOS (aarch64-apple-darwin20), R 4.5.1
* GitHub Actions: macOS, Windows and Ubuntu; R devel, release and oldrel

## R CMD check results

0 errors | 0 warnings | 2 notes

The first note reports the BugReports URL,
https://github.com/drjoshmcgrane/rasch/issues, as possibly invalid. The
issue tracker is enabled and the repository is public: the GitHub API
reports `has_issues: true` for it, and the repository root and
documentation site both resolve normally. GitHub returns 404 to
unauthenticated non-browser requests for the /issues and /pulls paths
alike, and pull requests cannot be disabled on a public repository, so
the status reflects that behaviour rather than a broken link.

The second note states that the locally installed HTML Tidy is not recent
enough to validate the HTML manual; it is specific to the local check
environment and unrelated to the package. The PDF manual was built and
checked successfully.

The package was built from source with its eight vignettes before checking.
Simulation-intensive and bootstrap-calibration tests use `skip_on_cran()`;
they run locally and in continuous integration.

## Current CRAN status

Version 1.11.7 is currently OK on all CRAN check flavours.

## Reverse dependencies

There are no known reverse dependencies.
