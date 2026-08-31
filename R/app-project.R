# rasch :: Shiny project files

# The graphical interface stores complete analysis sessions as ordinary RDS
# files with a small identifying header. Keeping validation here makes the
# format testable without a running Shiny session and leaves room for explicit
# schema migration in later releases.
.validate_app_project <- function(project) {
  fail <- function(message) stop(message, call. = FALSE)
  if (!is.list(project) || !identical(project$format, "rasch-shiny-project"))
    fail("not a rasch analysis file")
  # the schema is read as stored: coercing it would accept "1", TRUE, 1.5
  # and a factor's level code as schema 1
  schema <- project$schema
  if (length(schema) != 1L || !is.numeric(schema) || is.na(schema) ||
      schema != 1L)
    fail("unsupported rasch analysis-file schema")
  if (!(is.data.frame(project$data) || is.matrix(project$data)))
    fail("the analysis file does not contain a valid source dataset")
  if (!(inherits(project$base_fit, "rasch") ||
        inherits(project$base_fit, "rasch_btl")))
    fail("the analysis file does not contain a fitted rasch model")
  if (!is.null(project$model_type)) {
    # a factor would pass the membership test and then reach switch() as its
    # integer code, selecting a branch by level order
    if (!is.character(project$model_type) ||
        length(project$model_type) != 1L || is.na(project$model_type) ||
        !(project$model_type %in% c("rasch", "mfrm", "efrm", "btl")))
      fail("the analysis file names an unsupported model type")
    expected_class <- switch(project$model_type,
      rasch = "rasch", mfrm = "rasch_mfrm", efrm = "rasch_efrm",
      btl = "rasch_btl")
    if (!inherits(project$base_fit, expected_class))
      fail(sprintf(paste("the analysis file declares model type '%s' but",
                         "contains a fit of class '%s'"),
                   project$model_type, class(project$base_fit)[1]))
  }
  for (field in c("rasch_steps", "btl_steps", "kept_fits", "kept_fit_code",
                  "simulation", "results", "settings", "resources"))
    if (!is.null(project[[field]]) && !is.list(project[[field]]))
      fail(sprintf("the analysis file has an invalid %s field", field))

  # Results belong to the active fit, which is the final structural change
  # when a saved history is present. In particular, a bootstrap null from an
  # earlier fit must not be restored beside later DIF splits, superitems or
  # paired-comparison frame changes.
  is_btl <- inherits(project$base_fit, "rasch_btl")
  history <- if (is_btl) project$btl_steps else project$rasch_steps
  active_fit <- project$base_fit
  if (length(history)) {
    last <- history[[length(history)]]
    expected <- if (is_btl) "rasch_btl" else "rasch"
    if (!is.list(last) || !inherits(last$fit, expected))
      fail("the analysis file has an invalid fitted-model history")
    active_fit <- last$fit
  }
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
  invisible(project)
}

.save_app_project <- function(project, file) {
  .validate_app_project(project)
  saveRDS(project, file, version = 3, compress = "xz")
  invisible(file)
}

.read_app_project <- function(file) {
  project <- readRDS(file)
  .validate_app_project(project)
  project
}
