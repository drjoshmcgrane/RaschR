#' rasch: An R Implementation of Rasch Measurement Theory
#'
#' Implements Rasch Measurement Theory for the construction and evaluation of
#' measurement scales. The package is organised around sufficiency and
#' invariance: the relevant score contains the information required for person
#' measurement, and comparisons between persons and items should remain stable
#' across relevant samples and conditions.
#'
#' The model suite includes the dichotomous Rasch model, partial credit and
#' rating scale models, many-facet Rasch models, the extended frame of reference
#' model, and Bradley-Terry-Luce models for paired comparisons. Functions are
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
"_PACKAGE"
