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
#' The app's display packages ('shiny', 'bslib', 'DT', 'bsicons') are
#' suggested rather than required by the package. If any are missing,
#' \code{run_app} lists them all and, in an interactive session, offers to
#' install them before launching.
#'
#' @param ... Passed to \code{shiny::runApp}.
#' @return Called for its side effect of launching the app.
#' @examples
#' if (interactive()) run_app()
#' @export
run_app <- function(...) {
  .app_require(c("shiny", "bslib", "DT", "bsicons"))
  dir <- system.file("shiny", package = "rasch")
  if (dir == "") stop("app not found: reinstall rasch")
  shiny::runApp(dir, ...)
}

# Check the app's display packages in one pass. Reports everything missing
# at once, and in an interactive session offers to install the whole set
# (with explicit consent) rather than sending the user through one
# stop-install-retry loop per package.
.app_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1),
                          quietly = TRUE)]
  if (!length(missing)) return(invisible(TRUE))
  cmd <- sprintf("install.packages(c(%s))",
                 paste(sprintf('\"%s\"', missing), collapse = ", "))
  if (interactive()) {
    ans <- utils::askYesNo(sprintf(
      "The rasch app needs %d more package(s): %s. Install now?",
      length(missing), paste(missing, collapse = ", ")), default = TRUE)
    if (isTRUE(ans)) {
      utils::install.packages(missing)
      still <- missing[!vapply(missing, requireNamespace, logical(1),
                               quietly = TRUE)]
      if (!length(still)) return(invisible(TRUE))
      stop("installation did not complete for: ",
           paste(still, collapse = ", "), "; try ", cmd, call. = FALSE)
    }
  }
  stop("the rasch app needs ", paste(missing, collapse = ", "),
       "; run ", cmd, call. = FALSE)
}
