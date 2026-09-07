test_that("tailored results are authenticated, restored and exported", {
  set.seed(921)
  n <- 220L
  difficulty <- seq(-1.4, 1.6, length.out = 6L)
  theta <- rnorm(n)
  probability <- 0.25 + 0.75 * plogis(outer(theta, difficulty, "-"))
  X <- matrix(rbinom(n * 6L, 1L, probability), n, 6L,
              dimnames = list(NULL, paste0("I", 1:6)))
  fit <- rasch(X)
  tailored <- tailored_analysis(fit, chance = 0.25)

  expect_s3_class(tailored, "rasch_tailored")
  expect_identical(tailored$algorithm, "tailored-four-stage-2")
  expect_null(tailored$anchor_items_requested)
  expect_no_error(.validate_tailored_result(tailored, fit))
  carried <- fit
  attr(carried, "report_tailored") <- tailored
  expect_true(.fit_boot_signature_matches(tailored$fit_signature, carried))

  changed_fit <- rasch(X[rev(seq_len(nrow(X))), , drop = FALSE])
  expect_error(.validate_tailored_result(tailored, changed_fit),
               "different fitted model")
  changed_result <- tailored
  changed_result$table$shift[1L] <- changed_result$table$shift[1L] + 0.1
  expect_error(.validate_tailored_result(changed_result, fit),
               "internally inconsistent")

  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), settings = list(),
    results = list(guessing = tailored)))
  expect_no_error(.validate_app_project(project))

  # Results from releases before tailored-result fingerprints cannot be tied
  # to the active calibration. A valid enclosing project is reopened with
  # only that derived result omitted.
  old <- unclass(tailored)
  old$result_signature <- NULL
  old$fit_signature <- NULL
  old_project <- project
  old_project$results$guessing <- old
  old_project <- .seal_app_project(old_project)
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(old_project, path)
  expect_warning(restored <- .read_app_project(path),
                 "predates current result provenance")
  expect_null(restored$results$guessing)
  expect_match(attr(restored, "rasch_project_legacy_dropped"),
               "tailored analysis")
  expect_no_error(.validate_app_project(restored))

  old_selection <- tailored
  old_selection$anchor_items_requested <- NULL
  old_selection$result_signature <- NULL
  old_selection$result_signature <- .fit_boot_md5(old_selection)
  old_selection_project <- project
  old_selection_project$results$guessing <- old_selection
  old_selection_project <- .seal_app_project(old_selection_project)
  saveRDS(old_selection_project, path)
  expect_warning(restored_selection <- .read_app_project(path),
                 "predates current result provenance")
  expect_null(restored_selection$results$guessing)

  old_bootstrap <- project
  old_bootstrap$results$guessing$algorithm <- "tailored-four-stage-1"
  old_bootstrap$results$guessing$se_method <- "bootstrap"
  old_bootstrap <- .seal_app_project(old_bootstrap)
  saveRDS(old_bootstrap, path)
  expect_warning(restored_bootstrap <- .read_app_project(path),
                 "predates current result")
  expect_null(restored_bootstrap$results$guessing)
  expect_match(attr(restored_bootstrap, "rasch_project_legacy_dropped"),
               "superseded bootstrap")

  out <- tempfile("tailored-output-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  files <- suppressWarnings(save_outputs(
    fit, out, formats = "png", dpi = 72, item_plots = FALSE,
    tailored = tailored))
  expect_true(file.path(out, "tables", "tailored_item_shifts.csv") %in%
                files)
  expect_true(file.path(out, "tables", "tailored_analysis.csv") %in% files)
  saved <- utils::read.csv(file.path(out, "tables", "tailored_analysis.csv"),
                           stringsAsFactors = FALSE)
  expect_identical(saved$algorithm, tailored$algorithm)
  expect_identical(saved$fit_fingerprint,
                   tailored$fit_signature$fingerprint)
  expect_identical(saved$result_signature, tailored$result_signature)
  expect_identical(saved$anchor_selection, "automatic")
  expect_true(is.na(saved$anchor_items_requested))

  supplied <- tailored_analysis(fit, chance = 0.25,
                                anchor_items = tailored$anchor_items)
  expect_identical(supplied$anchor_items_requested, supplied$anchor_items)
  expect_no_error(.validate_tailored_result(supplied, fit))
  falsified_selection <- tailored
  falsified_selection$anchor_items_requested <- tailored$anchor_items
  falsified_selection$result_signature <- NULL
  falsified_selection$result_signature <- .fit_boot_md5(falsified_selection)
  expect_error(.validate_tailored_result(falsified_selection, fit),
               "internally inconsistent")

  html <- tempfile(fileext = ".html")
  on.exit(unlink(html), add = TRUE)
  expect_identical(suppressWarnings(report_html(
    fit, html, dpi = 72, tailored = tailored)), html)
  rendered <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(rendered, "Tailored analysis for guessing", fixed = TRUE)
  expect_match(rendered, "Item shifts are descriptive", fixed = TRUE)

  wrong_out <- tempfile("wrong-tailored-output-")
  on.exit(unlink(wrong_out, recursive = TRUE), add = TRUE)
  expect_error(save_outputs(
    changed_fit, wrong_out, formats = "png", item_plots = FALSE,
    tailored = tailored), "different fitted model")
  expect_false(dir.exists(wrong_out))
})

test_that("tailored bootstrap keeps zero changes and classifies failures", {
  fit <- rasch(simulate_rasch(100, 4, seed = 929)[, sprintf("I%02d", 1:4)])
  items <- fit$items$item
  maxima <- fit$m
  # A positive chance floor below every fitted observed-cell probability
  # makes the no-tailoring case deterministic without altering the fit.
  chance <- min(fit$moments$E[!is.na(fit$X)], na.rm = TRUE) / 2
  expect_equal(.tailored_boot_refit(fit, chance, NULL, items, maxima),
               rep(0, length(items)))
  expect_identical(.fit_boot_status(.tailored_boot_refit(
    fit, chance, NULL, rev(items), maxima)), "error")
  expect_identical(.fit_boot_status(.tailored_boot_refit(
    fit, chance, NULL, items, maxima + 1L)), "error")

  for (stage in c("tailored", "origin-equated")) {
    ans <- with_mocked_bindings(
      .tailored_boot_refit(fit, 0.99, NULL, items, maxima),
      tailored_analysis = function(...)
        .tailored_nonconvergence(paste(stage, "did not converge")),
      .package = "rasch")
    expect_identical(.fit_boot_status(ans), "nonconverged")
  }
})
