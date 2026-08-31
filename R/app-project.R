# rasch :: Shiny project files

# The graphical interface stores complete analysis sessions as ordinary RDS
# files with a small identifying header. Validation lives here so it can be
# exercised without a running Shiny session. Schema 2 adds an integrity seal:
# schema 1 stored data and fits side by side but had no way to establish that
# they came from the same saved analysis.

.app_fit_family <- function(fit) {
  if (inherits(fit, "rasch_btl")) return("btl")
  if (inherits(fit, "rasch_efrm")) return("efrm")
  if (inherits(fit, "rasch_mfrm")) return("mfrm")
  if (inherits(fit, "rasch")) return("rasch")
  NA_character_
}

.app_scalar_text <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.validate_app_fit <- function(fit, what = "fitted model") {
  fail <- function(message) stop(what, " ", message, call. = FALSE)
  family <- .app_fit_family(fit)
  if (is.na(family)) fail("is not a fitted rasch model")

  if (identical(family, "btl")) {
    if (!is.data.frame(fit$objects) || nrow(fit$objects) < 2L ||
        !all(c("object", "location") %in% names(fit$objects)) ||
        !is.character(fit$objects$object) || anyNA(fit$objects$object) ||
        any(!nzchar(fit$objects$object)) || anyDuplicated(fit$objects$object))
      fail("has an invalid object calibration")
    if (!is.data.frame(fit$comparisons) || !nrow(fit$comparisons) ||
        !all(c("object_a", "object_b", "response", "weight", "judge") %in%
             names(fit$comparisons)))
      fail("has an invalid comparison design")
    if (!is.data.frame(fit$pairs) ||
        !all(c("object_a", "object_b") %in% names(fit$pairs)))
      fail("has an invalid pair-fit table")
    if (length(fit$m) != 1L || !is.numeric(fit$m) || !is.finite(fit$m) ||
        fit$m < 1 || fit$m != floor(fit$m))
      fail("has an invalid response scale")
    if (length(fit$n_comparisons) != 1L ||
        !is.numeric(fit$n_comparisons) || !is.finite(fit$n_comparisons) ||
        fit$n_comparisons < 1)
      fail("has an invalid comparison count")
    if (length(fit$converged) != 1L || !is.logical(fit$converged) ||
        is.na(fit$converged))
      fail("has no valid convergence record")

    # The extended-frame result has its own fitted-probability machinery.
    # Ordinary and explanatory BTL fits must retain the row-wise probabilities
    # and public refit specification used by the fitted-model bootstrap.
    if (inherits(fit, "rasch_btl_efrm")) {
      if (!is.data.frame(fit$phi_table) || !is.data.frame(fit$alpha_table) ||
          !is.data.frame(fit$kappa_table) || !is.data.frame(fit$frames))
        fail("has an incomplete frame calibration")
    } else {
      if (!is.matrix(fit$fitted_prob) ||
          !identical(dim(fit$fitted_prob),
                     c(nrow(fit$comparisons), as.integer(fit$m) + 1L)) ||
          any(!is.finite(fit$fitted_prob)))
        fail("has invalid fitted comparison probabilities")
      if (!is.list(fit$refit_spec))
        fail("has no valid refit specification")
    }
    return(invisible(family))
  }

  X <- fit$X
  if (!is.matrix(X) || length(dim(X)) != 2L || any(dim(X) < 1L))
    fail("has no valid response matrix")
  N <- nrow(X); L <- ncol(X)
  if (!.app_scalar_text(fit$model)) fail("has no valid model name")
  if (!is.numeric(fit$m) || length(fit$m) != L ||
      any(!is.finite(fit$m)) || any(fit$m < 1) || any(fit$m != floor(fit$m)))
    fail("has an invalid item response scale")
  if (!is.data.frame(fit$items) || nrow(fit$items) != L ||
      !all(c("item", "location") %in% names(fit$items)) ||
      !is.character(fit$items$item) || anyNA(fit$items$item) ||
      any(!nzchar(fit$items$item)) || anyDuplicated(fit$items$item))
    fail("has an invalid item calibration")
  if (!is.data.frame(fit$thresholds) ||
      !all(c("item", "k", "tau") %in% names(fit$thresholds)) ||
      nrow(fit$thresholds) != sum(fit$m))
    fail("has an invalid threshold calibration")
  if (!is.list(fit$tau_list) || length(fit$tau_list) != L ||
      !identical(as.integer(lengths(fit$tau_list)), as.integer(fit$m)) ||
      any(!vapply(fit$tau_list, function(z)
        is.numeric(z) && all(is.finite(z)), logical(1))))
    fail("has invalid item threshold vectors")
  if (!is.data.frame(fit$person) || nrow(fit$person) != N ||
      !all(c("id", "raw", "theta") %in% names(fit$person)))
    fail("has an invalid person calibration")
  if (!is.matrix(fit$residuals) || !identical(dim(fit$residuals), dim(X)))
    fail("has an invalid residual matrix")
  if (!is.list(fit$moments) || !is.matrix(fit$moments$E) ||
      !is.matrix(fit$moments$V) ||
      !identical(dim(fit$moments$E), dim(X)) ||
      !identical(dim(fit$moments$V), dim(X)))
    fail("has invalid fitted response moments")
  if (!is.list(fit$est) || length(fit$est$converged) != 1L ||
      !is.logical(fit$est$converged) || is.na(fit$est$converged))
    fail("has no valid estimation record")
  if (!is.null(fit$factors) &&
      (!is.data.frame(fit$factors) || nrow(fit$factors) != N))
    fail("has an invalid person-factor design")
  invisible(family)
}

.app_project_binding <- function(project) {
  x <- project
  attr(x, "rasch_project_legacy") <- NULL
  x$binding <- NULL
  .fit_boot_md5(x)
}

.seal_app_project <- function(project) {
  if (!is.list(project)) stop("`project` must be a list", call. = FALSE)
  project$schema <- 2L
  project$binding <- .app_project_binding(project)
  project
}

.validate_app_project <- function(project) {
  fail <- function(message) stop(message, call. = FALSE)
  if (!is.list(project) || !identical(project$format, "rasch-shiny-project"))
    fail("not a rasch analysis file")
  # Read the schema as stored: coercion would accept "2", TRUE, 2.5 or a
  # factor's level code as schema 2.
  schema <- project$schema
  if (length(schema) != 1L || !is.numeric(schema) || is.na(schema) ||
      schema != 2L)
    fail(paste("unsupported rasch analysis-file schema; this version needs",
               "a schema-2 file with data-to-fit integrity information"))
  if (!(is.data.frame(project$data) || is.matrix(project$data)) ||
      nrow(project$data) < 1L || ncol(project$data) < 1L)
    fail("the analysis file does not contain a valid source dataset")
  data_names <- colnames(project$data)
  if (is.null(data_names) || anyNA(data_names) ||
      any(!nzchar(data_names)) || anyDuplicated(data_names))
    fail("the analysis file has invalid source-data column names")

  base_problem <- tryCatch({
    .validate_app_fit(project$base_fit, "the saved base fit")
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(base_problem)) fail(base_problem)
  base_family <- .app_fit_family(project$base_fit)
  if (!is.null(project$model_type)) {
    # A factor would pass a membership test and then reach switch() as its
    # integer code, selecting a branch by level order.
    if (!.app_scalar_text(project$model_type) ||
        !(project$model_type %in% c("rasch", "mfrm", "efrm", "btl")))
      fail("the analysis file names an unsupported model type")
    if (!identical(base_family, project$model_type))
      fail(sprintf(paste("the analysis file declares model type '%s' but",
                         "contains a fit of class '%s'"),
                   project$model_type, class(project$base_fit)[1]))
  }
  for (field in c("rasch_steps", "btl_steps", "kept_fits", "kept_fit_code",
                  "simulation", "results", "settings", "resources"))
    if (!is.null(project[[field]]) && !is.list(project[[field]]))
      fail(sprintf("the analysis file has an invalid %s field", field))

  validate_history <- function(history, family, field) {
    if (!length(history)) return(invisible(NULL))
    for (i in seq_along(history)) {
      entry <- history[[i]]
      metadata_ok <- is.list(entry) && .app_scalar_text(entry$type) &&
        .app_scalar_text(entry$label) && is.list(entry$details) &&
        (is.null(entry$code) || .app_scalar_text(entry$code)) &&
        .app_scalar_text(entry$created)
      if (!metadata_ok)
        fail(sprintf("the analysis file has invalid %s history metadata at entry %d",
                     field, i))
      problem <- tryCatch({
        .validate_app_fit(entry$fit,
                          sprintf("the %s history fit at entry %d", field, i))
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(problem) ||
          !identical(.app_fit_family(entry$fit), family))
        fail(sprintf("the analysis file has an invalid fitted-model history in %s at entry %d%s",
                     field, i,
                     if (is.null(problem)) "" else paste0(": ", problem)))
    }
    invisible(NULL)
  }

  is_btl <- identical(base_family, "btl")
  rasch_history <- project$rasch_steps %||% list()
  btl_history <- project$btl_steps %||% list()
  if ((is_btl && length(rasch_history)) ||
      (!is_btl && length(btl_history)))
    fail("the analysis file has fitted-model history for the inactive model family")
  validate_history(if (is_btl) btl_history else rasch_history,
                   base_family, if (is_btl) "btl_steps" else "rasch_steps")
  history <- if (is_btl) btl_history else rasch_history
  active_fit <- if (length(history)) history[[length(history)]]$fit
                else project$base_fit

  # Results belong to the active fit, which is the final structural change.
  # In particular, a bootstrap null from an earlier fit must not be restored
  # beside later DIF splits, superitems or paired-comparison frame changes.
  bootstrap <- project$results$bootstrap
  if (!is.null(bootstrap)) {
    if (!is.list(bootstrap) || is.null(bootstrap$bs))
      fail("the analysis file has an invalid saved bootstrap result")
    problem <- tryCatch({
      .validate_fit_bootstrap(bootstrap$bs, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved bootstrap result does not belong to the active fit:",
                 problem))
  }

  if (!.app_scalar_text(project$binding))
    fail("the analysis file has no valid data-to-fit integrity information")
  expected <- .app_project_binding(project)
  if (!identical(project$binding, expected))
    fail(paste("the analysis file's source data, fitted models or results have",
               "changed since they were saved"))
  invisible(project)
}

.save_app_project <- function(project, file) {
  .validate_app_project(project)
  saveRDS(project, file, version = 3, compress = "xz")
  invisible(file)
}

.read_app_project <- function(file) {
  project <- readRDS(file)
  legacy <- is.list(project) &&
    identical(project$format, "rasch-shiny-project") &&
    length(project$schema) == 1L && is.numeric(project$schema) &&
    !is.na(project$schema) && project$schema == 1L
  if (legacy) {
    # Schema 1 did not record an integrity binding. It can be checked
    # structurally and upgraded, but its original data-to-fit relationship
    # cannot be authenticated retrospectively.
    project <- .seal_app_project(project)
  }
  .validate_app_project(project)
  if (legacy) {
    attr(project, "rasch_project_legacy") <- TRUE
    warning(paste("this schema-1 analysis predates data-to-fit integrity",
                  "checks; it passed structural validation and has been",
                  "upgraded in memory; save it again to retain schema 2"),
            call. = FALSE)
  }
  project
}
