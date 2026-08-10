# rasch :: Shiny launcher

#' Launch the rasch point-and-click graphical interface
#'
#' Opens the guided Shiny application for users who prefer a graphical
#' workflow. The app supports data import; assignment of item, person, group,
#' rater, and comparison roles; model selection and fitting; interactive
#' diagnostics, tables, and plots; and one-click export. The corresponding R
#' call is shown for each analysis, but users do not need to write R analysis
#' code to use the principal workflows.
#'
#' @param ... Passed to \code{shiny::runApp}.
#' @return Called for its side effect of launching the app.
#' @examples
#' if (interactive()) run_app()
#' @export
run_app <- function(...) {
  for (pkg in c("shiny", "bslib", "DT", "bsicons")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("the rasch app needs the '", pkg, "' package: install.packages(\"", pkg, "\")")
  }
  dir <- system.file("shiny", package = "rasch")
  if (dir == "") stop("app not found: reinstall rasch")
  shiny::runApp(dir, ...)
}
