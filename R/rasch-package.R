#' rasch: Models and Diagnostics for Rasch Measurement Theory
#'
#' Implements Rasch Measurement Theory for the construction and evaluation of
#' measurement scales. The package is organised around sufficiency and
#' invariance: the total score is sufficient for the person parameter, which
#' permits person and item parameters to be separated in conditional
#' estimation; and, within a specified frame of reference, comparisons
#' between persons do not depend on which items are used, while comparisons
#' between items do not depend on which persons respond -- Rasch's criterion
#' of invariant comparison (specific objectivity; Rasch 1961, 1977). These
#' are requirements the data must meet for measurement, and the diagnostic
#' suite examines whether they do.
#'
#' The model suite (Rasch 1960; Andrich and Marais 2019) includes the
#' dichotomous Rasch model, partial credit and rating scale models,
#' many-facet Rasch models, the extended frame of reference model (Humphry
#' and Andrich 2008), and comparative judgement models for paired
#' comparisons -- the conditional form of the dichotomous Rasch model
#' (Andrich 1978). Functions are
#' provided for item and person estimation, fit, targeting, dimensionality,
#' local response dependence, differential item functioning, threshold
#' functioning, reliability, equating, repeated measurements, simulation, and
#' graphical reporting. See \code{\link{rasch}}, \code{\link{rasch_mfrm}},
#' \code{\link{rasch_efrm}}, and \code{\link{btl}} for the principal fitting
#' functions. Results can be exported with \code{\link{save_outputs}}.
#'
#' @section Point-and-click graphical interface:
#' The package includes a guided Shiny application, launched with
#' \code{\link{run_app}}. It supports data import, assignment of item, person,
#' group, rater, and comparison roles, model fitting, diagnostics, plots, and
#' export through a graphical interface. R is needed to install the package
#' and launch the application, but users can then complete the principal
#' analyses without writing R code. The app displays the corresponding R call
#' for each analysis so that work conducted through the interface remains
#' transparent and reproducible.
#'
#' @section Choosing a model:
#' Use \code{\link{rasch}} for persons responding to a common set of items.
#' Dichotomous data are a special case of its partial credit model; set
#' \code{model = "RSM"} when the same category-threshold structure is to be
#' shared across items. Use \code{\link{rasch_mfrm}} when observations also
#' contain facets such as raters, tasks, or occasions whose additive severity
#' is part of the measurement model. Use \code{\link{rasch_efrm}} only when
#' the substantive model allows the unit itself to differ across linked
#' item-set by person-group frames. For paired-comparison data, use
#' \code{\link{btl}}; \code{\link{btl_efrm}} supplies the corresponding
#' linked-frame extension.
#'
#' These models address different designs and invariance claims; they are not
#' a ladder of progressively better-fitting alternatives. In particular,
#' person factors ordinarily belong in an invariance analysis with
#' \code{\link{dif_anova}}, whereas rating facets belong in the MFRM and
#' frame-defining factors belong in the EFRM.
#'
#' @section Typical workflow:
#' Begin by checking coding, category use, missingness, and design
#' connectedness. Missing responses can be accommodated when the observed
#' response graph identifies a common scale and the missingness process is
#' ignorable for the model; informative missingness requires an explicit
#' sensitivity analysis. Fit the model required by the design, then examine targeting,
#' item and person fit, threshold ordering, local response dependence, and
#' dimensionality. Test invariance across all relevant person factors jointly;
#' repeated-measures factors should be identified as within-person in
#' \code{\link{dif_anova}}. Investigate the size and substantive meaning of a
#' departure before modifying the scale. Refit after any defensible change and
#' report the final model, exclusions or resolutions, uncertainty, targeting,
#' and remaining limitations. The package vignettes give worked workflows for
#' each main model family.
#'
#' @keywords internal
#' @import stats
#' @import graphics
#' @import grDevices
#' @importFrom utils write.csv combn
#'
#' @references
#' Rasch, G. (1960). Probabilistic Models for Some Intelligence and
#' Attainment Tests. Copenhagen: Danish Institute for Educational Research.
#' (Expanded edition, 1980, Chicago: University of Chicago Press.)
#'
#' Rasch, G. (1961). On general laws and the meaning of measurement in
#' psychology. In Proceedings of the Fourth Berkeley Symposium on
#' Mathematical Statistics and Probability (Vol. 4, pp. 321--333).
#' Berkeley: University of California Press.
#'
#' Rasch, G. (1977). On specific objectivity: An attempt at formalizing the
#' request for generality and validity of scientific statements. Danish
#' Yearbook of Philosophy, 14, 58--94.
#'
#' Andrich, D. (1978). Relationships between the Thurstone and Rasch
#' approaches to item scaling. Applied Psychological Measurement, 2(3),
#' 449--460.
#'
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#'
#' Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
#' Rasch model. Journal of Applied Measurement, 9(3), 249--264.
"_PACKAGE"
