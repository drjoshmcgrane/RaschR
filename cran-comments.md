# CRAN comments for rasch 1.12.0

## Summary of this update

This is the first update since rasch 1.11.7 was accepted on 2026-07-30.

The release extends differential item functioning analysis to multiple
between-person and within-person factors, with factorial terms and the
corresponding repeated-measures error structure. It also expands the
many-facet, extended frame of reference, and paired-comparison models, and
adds simulation and recovery functions.

Corrections concern standard errors, reference distributions, and the
conditions under which inferential results are reported. The affected areas
include response-dependence magnitude, extended-frame linking, and
judge-clustered paired-comparison inference. The estimators now apply more
explicit checks for connectedness, separation, rank, sparse categories, and
weak identification.

The package title, description, help pages, six vignettes, and graphical
interface have also been revised for this release.

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
