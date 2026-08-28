# CRAN comments for rasch 1.12.1

## Summary of this update

This is a maintenance update to rasch 1.12.0, and chiefly a correctness
release. A sustained internal review of the
package found and corrected a number of defects that could report a wrong
number rather than an error. The most consequential is the covariance of
the centred location differences in `frame_invariance()`, which was
computed without the transpose of the centring matrix: when compared items
differed in maximum score, every standard error, test statistic and
probability in that table was wrong, over-flagging the dichotomous items
and hiding the polytomous one.

Procedures that compare a fit with a refit now refuse an explanatory fit
rather than silently releasing its design restriction, and the parallel
scree reference refits under the model that was fitted. Derived fits carry
the estimation controls they were built from. Simulation now plants what
its recorded generating values claim, and refuses requests it cannot plant.
Exports no longer report success they did not achieve: an archive in which
nothing could be drawn is an error rather than a path to a file that was
never written.

Input and selection boundaries are hardened throughout, so that a
mis-specified role, identifier, key or display control is refused where it
is written instead of changing the analysis in silence.

The application gains simulation of explanatory models and supplementary
weighted person measures, and a saved analysis now restores its data roles
and estimation controls. A data-structures vignette has been added, making
eight in total.

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.6.1
* win-builder: R-release and R-devel
* GitHub Actions: macOS, Windows and Ubuntu; R devel, release and oldrel

## R CMD check results

0 errors | 0 warnings | 0 notes locally; win-builder adds the
incoming-feasibility NOTE only.

The check-time measures introduced at 1.12.0 are retained: the test suite
runs a small core on CRAN that exercises every estimator once, and the
complete suite runs whenever NOT_CRAN is true, locally and in continuous
integration on three operating systems. Simulation-intensive and
bootstrap-calibration tests use `skip_on_cran()`.

Any words flagged by the incoming spell check are author surnames from the
references cited in the help pages.

Some check services report the BugReports URL,
https://github.com/drjoshmcgrane/rasch/issues, as possibly invalid. The
issue tracker is enabled and the repository is public: the GitHub API
reports `has_issues: true` for it, and the repository root and
documentation site both resolve normally. GitHub returns 404 to
unauthenticated non-browser requests for the /issues and /pulls paths
alike, and pull requests cannot be disabled on a public repository, so
the status reflects that behaviour rather than a broken link.

The package was built from source with its eight vignettes before checking.

## Current CRAN status

Version 1.12.0 is currently OK on all CRAN check flavours.

## Reverse dependencies

There are no known reverse dependencies.
