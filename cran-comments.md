# CRAN comments for rasch 1.14.2

## Summary of this update

This is the first update since rasch 1.11.7 was accepted on 2026-07-30.

It arrives promptly because post-release auditing identified statistical
correctness issues that should not remain in the released version: some
standard errors and significance flags could be reported for parameters the
data did not identify (for example differential-item-functioning magnitudes
resting on near-empty categories, and information curves pooled over item
sets no respondent took together). These are corrected, alongside checks
for connectedness, separation, rank, and weakly identified parameters
across all model families. A release-wide simulation validation battery
(recovery, standard error calibration, null rates, and power for every
model family and diagnostic, under complete and missing data) was run
before submission; the one defect it found, an understated standard error
in the response-dependence magnitude estimate, is corrected in this
version.

The update also extends differential item functioning analysis to multiple
between-person and within-person factors (factorial terms with the
appropriate repeated-measures error structure), expands the many-facet,
extended frame of reference, and paired-comparison model families, adds
simulation and recovery tools, and adds five vignettes covering the main
analysis workflows.

The package title and short description have been revised to state the
purpose more directly. Technical details remain in the function
documentation and vignettes.

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.5.1
* GitHub Actions: macOS, Windows and Ubuntu; R devel, release and oldrel

## R CMD check results

0 errors | 0 warnings | 1 note

The note states that the locally installed HTML Tidy is not recent enough
to validate the HTML manual; it is specific to the local check environment
and unrelated to the package. The PDF manual was built and checked
successfully.

The package was built from source with its six vignettes before checking.
Simulation-intensive and bootstrap-calibration tests use `skip_on_cran()`;
they run locally and in continuous integration.

## Current CRAN status

Version 1.11.7 is currently OK on all CRAN check flavours.

## Reverse dependencies

There are no known reverse dependencies.
