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
  is.character(x) && is.null(dim(x)) && is.null(oldClass(x)) &&
    length(x) == 1L && !is.na(x) && nzchar(trimws(x))
}

.validate_app_fit <- function(fit, what = "fitted model") {
  fail <- function(message) stop(what, " ", message, call. = FALSE)
  family <- .app_fit_family(fit)
  if (is.na(family)) fail("is not a fitted rasch model")

  if (identical(family, "btl")) {
    if (!is.data.frame(fit$objects) || nrow(fit$objects) < 2L ||
        !all(c("object", "location") %in% names(fit$objects)) ||
        !is.character(fit$objects$object) || anyNA(fit$objects$object) ||
        any(!nzchar(trimws(fit$objects$object))) ||
        anyDuplicated(fit$objects$object))
      fail("has an invalid object calibration")
    if (!is.data.frame(fit$comparisons) || !nrow(fit$comparisons) ||
        !all(c("object_a", "object_b", "response", "weight", "judge") %in%
             names(fit$comparisons)))
      fail("has an invalid comparison design")
    if (!is.data.frame(fit$pairs) ||
        !all(c("object_a", "object_b") %in% names(fit$pairs)))
      fail("has an invalid pair-fit table")
    if (length(fit$m) != 1L || !is.numeric(fit$m) || is.complex(fit$m) ||
        !is.null(dim(fit$m)) || !is.null(oldClass(fit$m)) ||
        !is.finite(fit$m) ||
        fit$m < 1 || fit$m != floor(fit$m))
      fail("has an invalid response scale")
    if (length(fit$n_comparisons) != 1L ||
        !is.numeric(fit$n_comparisons) || is.complex(fit$n_comparisons) ||
        !is.null(dim(fit$n_comparisons)) ||
        !is.null(oldClass(fit$n_comparisons)) ||
        !is.finite(fit$n_comparisons) ||
        fit$n_comparisons < 1)
      fail("has an invalid comparison count")
    if (length(fit$converged) != 1L || !is.logical(fit$converged) ||
        !is.null(dim(fit$converged)) || !is.null(oldClass(fit$converged)) ||
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
      any(!nzchar(trimws(fit$items$item))) || anyDuplicated(fit$items$item))
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
  attr(x, "rasch_project_legacy_dropped") <- NULL
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
  if (length(schema) != 1L || !is.numeric(schema) || is.complex(schema) ||
      !is.null(dim(schema)) || !is.null(oldClass(schema)) || is.na(schema) ||
      schema != 2L)
    fail(paste("unsupported rasch analysis-file schema; this version needs",
               "a schema-2 file with data-to-fit integrity information"))
  if (!(is.data.frame(project$data) || is.matrix(project$data)) ||
      nrow(project$data) < 1L || ncol(project$data) < 1L)
    fail("the analysis file does not contain a valid source dataset")
  data_names <- colnames(project$data)
  if (is.null(data_names) || anyNA(data_names) ||
      any(!nzchar(trimws(data_names))) || anyDuplicated(data_names))
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

  # Exact indexing matters here: `$dif` partially matches `dif_bootstrap`
  # when the primary result is absent. A primary result is tied to the active
  # fit even when no bootstrap was requested; the project seal proves only
  # that the bundle has not changed since it was written, not that two
  # separately supplied fitted objects belong together.
  primary_dif <- if (is_btl) project$results[["btl_dif"]] else
    project$results[["dif"]]
  if (!is.null(primary_dif)) {
    problem <- tryCatch({
      if (is_btl) .validate_btl_dif_result(primary_dif, active_fit)
      else .validate_dif_result(primary_dif, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved primary DIF analysis does not belong to the active fit:",
                 problem))
    # Ordinary Rasch DIF is recomputed reactively after the controls are
    # restored. Pin the stored analysis to those controls even when it has no
    # bootstrap, so the displayed run is the one that was validated.
    if (!is_btl) {
      settings <- project$settings %||% list()
      effects <- settings$dif_effects %||% "main"
      alpha <- settings$dif_alpha %||% 0.05
      if (!.app_scalar_text(effects) ||
          !identical(primary_dif$effects, effects) ||
          length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
          !isTRUE(all.equal(primary_dif$alpha, alpha, tolerance = 0)))
        fail(paste("the saved primary DIF analysis does not match the",
                   "restored DIF settings"))
    }
  }

  btl_dif_meta <- project$results[["btl_dif_meta"]]
  if (!is.null(btl_dif_meta)) {
    if (!is_btl || is.null(primary_dif) || !is.list(btl_dif_meta) ||
        !.app_scalar_text(btl_dif_meta$judge_col) ||
        !btl_dif_meta$judge_col %in% data_names)
      fail("the analysis file has invalid saved Comparative Judgement DIF display metadata")
  }

  dif_bootstrap <- project$results$dif_bootstrap
  if (!is.null(dif_bootstrap)) {
    if (!is.list(dif_bootstrap) || is.null(dif_bootstrap$db))
      fail("the analysis file has an invalid saved DIF bootstrap result")
    if (is.null(primary_dif))
      fail(paste("the saved DIF bootstrap has no accompanying primary DIF",
                 "analysis"))
    problem <- tryCatch({
      .validate_dif_bootstrap(dif_bootstrap$db, active_fit, primary_dif)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved DIF bootstrap result does not belong to the active fit and DIF analysis:",
                 problem))
  }

  dimensionality <- project$results$dimensionality
  if (!is.null(dimensionality)) {
    problem <- tryCatch({
      if (is_btl) .validate_btl_dimensionality(dimensionality, active_fit)
      else .validate_scree_result(dimensionality, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved dimensionality analysis does not belong to the active fit:",
                 problem))
  }
  subtest <- project$results$subtest
  if (!is.null(subtest)) {
    if (is_btl)
      fail("the saved person-subset dimensionality test accompanies a paired-comparison fit")
    problem <- tryCatch({
      .validate_dimensionality_test(subtest, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved person-subset dimensionality test does not belong to the active fit:",
                 problem))
  }
  invariance <- project$results$frame_invariance
  if (!is.null(invariance)) {
    problem <- tryCatch({
      .validate_frame_invariance(invariance, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved frame-invariance analysis does not belong to the active fit:",
                 problem))
  }

  if (!.app_scalar_text(project$binding))
    fail("the analysis file has no valid data-to-fit integrity information")
  unsigned_project <- project
  attr(unsigned_project, "rasch_project_legacy") <- NULL
  attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
  unsigned_project$binding <- NULL
  if (!.fit_boot_hash_matches(project$binding, unsigned_project))
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
  dropped <- character(0)
  # Early schema-2 projects contain valid signed bootstrap arrays but predate
  # external (leave-one-out) maxT standardisation. Their adjusted
  # probabilities must not be displayed under the corrected algorithm. The
  # original project seal is checked before the derived result is removed.
  old_maxt <- !legacy && is.list(project) &&
    length(project$schema) == 1L && is.numeric(project$schema) &&
    !is.na(project$schema) && project$schema == 2L &&
    is.list(project$results) &&
    is.list(project$results$bootstrap) &&
    is.list(project$results$bootstrap$bs) &&
    is.null(project$results$bootstrap$bs$algorithm)
  old_frame_invariance <- !legacy && is.list(project) &&
    length(project$schema) == 1L && is.numeric(project$schema) &&
    !is.na(project$schema) && project$schema == 2L &&
    is.list(project$results) &&
    is.list(project$results$frame_invariance) &&
    !all(c("family_n", "boot_reps", "boot_reps_used",
           "boot_reps_nonconverged", "boot_reps_errors",
           "boot_minimum_usable") %in%
         names(project$results$frame_invariance))
  if (isTRUE(old_maxt)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    project$results$bootstrap <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped, "fit bootstrap (superseded maxT adjustment)")
  }
  # Earlier schema-2 projects can carry frame-invariance results produced
  # before complete multiplicity/bootstrap accounting and the strict
  # all-frame comparison rule. Authenticate the bundle before omitting that
  # derived result.
  if (isTRUE(old_frame_invariance)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    project$results$frame_invariance <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped, "frame-invariance analysis (superseded inference)")
  }
  if (legacy) {
    # Schema 1 did not record an integrity binding. It can be checked
    # structurally and upgraded, but its original data-to-fit relationship
    # cannot be authenticated retrospectively. Results written before the
    # current result fingerprints cannot be validated against the retained
    # fit, so omit those results while preserving the data, fit and history.
    results <- project$results
    has_signature <- function(x)
      is.list(x) && is.character(x$result_signature) &&
        length(x$result_signature) == 1L && !is.na(x$result_signature)
    if (is.list(results)) {
      if (!is.null(results$bootstrap) &&
          (!is.list(results$bootstrap) ||
           !has_signature(results$bootstrap$bs))) {
        results$bootstrap <- NULL
        dropped <- c(dropped, "fit bootstrap")
      }
      primary_dropped <- FALSE
      for (nm in c("dif", "btl_dif")) {
        if (!is.null(results[[nm]]) && !has_signature(results[[nm]])) {
          results[[nm]] <- NULL
          dropped <- c(dropped, if (nm == "btl_dif")
            "Comparative Judgement DIF" else "DIF")
          primary_dropped <- TRUE
        }
      }
      if (isTRUE(primary_dropped)) results$btl_dif_meta <- NULL
      if (!is.null(results$dif_bootstrap) &&
          (isTRUE(primary_dropped) || !is.list(results$dif_bootstrap) ||
           !has_signature(results$dif_bootstrap$db))) {
        results$dif_bootstrap <- NULL
        dropped <- c(dropped, "DIF bootstrap")
      }
      project$results <- results
    }
    project <- .seal_app_project(project)
  }
  .validate_app_project(project)
  if (legacy) {
    attr(project, "rasch_project_legacy") <- TRUE
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("this schema-1 analysis predates data-to-fit integrity",
                  "checks; it passed structural validation and has been",
                  "upgraded in memory; save it again to retain schema 2",
                  if (length(dropped)) paste0("; unverifiable saved results ",
                    "were omitted: ", paste(unique(dropped), collapse = ", "))
                  else ""),
            call. = FALSE)
  }
  if (isTRUE(old_maxt)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved fit bootstrap used the earlier maxT",
                  "standardisation and was omitted; recompute it before",
                  "reporting adjusted bootstrap probabilities"),
            call. = FALSE)
  }
  if (isTRUE(old_frame_invariance)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved frame-invariance analysis used earlier",
                  "comparison or bootstrap-accounting rules and was omitted;",
                  "recompute it before reporting frame-invariance inference"),
            call. = FALSE)
  }
  project
}
