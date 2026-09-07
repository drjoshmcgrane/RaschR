test_that("tailored scoring fits retain anchors and refuse unsupported refits", {
  set.seed(726)
  n <- 400L
  X <- matrix(rbinom(n * 8L, 1,
    plogis(outer(rnorm(n), seq(-2, 2, length.out = 8), "-"))), n, 8,
    dimnames = list(NULL, paste0("I", 1:8)))
  fit <- rasch(X)
  tailored <- tailored_analysis(fit)
  fixed <- tailored$anchored
  expect_identical(fixed$est$n_parameters, 0L)
  expect_true(all(fixed$thresholds$anchored))
  expect_true(fixed$refit_spec$fixed_calibration)
  expect_identical(fixed$refit_spec$calibration_source, "tailored")
  expect_identical(fixed$refit_spec$anchors$item, colnames(X))
  expect_identical(fixed$refit_spec$anchors$k, fixed$thresholds$k)
  expect_identical(fixed$refit_spec$anchors$tau, fixed$thresholds$tau)
  expect_no_error(.validate_tailored_result(tailored, fit))
  expect_true(.require_refittable_calibration(fit))
  expect_no_error(score_table(fixed))
  expect_no_error(residual_pca(fixed))
  expect_no_error(.scree_analysis(fixed, parallel = FALSE))
  observed_scree <- .scree_analysis(fixed, parallel = FALSE)
  expect_no_error(.validate_scree_result(observed_scree, fixed))
  observed_dim <- dimensionality_test(fixed,
    items_positive = paste0("I", 1:4), items_negative = paste0("I", 5:8))
  expect_no_error(.validate_dimensionality_test(observed_dim, fixed))
  weighted <- weighted_person_estimates(fixed, setNames(rep(1, 8), colnames(X)))
  expect_equal(weighted$theta, fixed$person$theta, tolerance = 1e-7)

  # No generator or calibration should run before refusing. This also protects
  # callers that select many bootstrap workers or request long simulations.
  local_mocked_bindings(
    .fit_gen_conditional = function(...) stop("unexpected null generation"),
    .pcml_fit = function(...) stop("unexpected item recalibration"))
  legacy <- fixed
  legacy$refit_spec <- NULL
  legacy <- unserialize(serialize(legacy, NULL))
  before <- fixed$thresholds
  for (f in list(fixed, legacy)) {
    expect_identical(.estimation_label(f),
                      "fixed item calibration; person scoring")
    summary <- fit_summary_table(f)
    expect_identical(summary$value[summary$statistic == "Estimation"],
                      "fixed item calibration; person scoring")
    expect_error(drop_items(f, "I8"), "fully anchored calibration",
                  class = "rasch_refusal")
    expect_error(combine_items(f, list(c("I1", "I2"))),
                  "fully anchored calibration")
    expect_error(split_items(f, "I1", rep(c("A", "B"), length.out = n)),
                  "fully anchored calibration")
    expect_error(.rasch_refit(f, f$X), "fully anchored calibration")
    expect_error(tailored_analysis(f), "fully anchored calibration")
    expect_error(fit_bootstrap(f, B = 20, workers = 1),
                  "fully anchored calibration")
    expect_error(dif_bootstrap(f, B = 20, workers = 1),
                  "fully anchored calibration")
    expect_error(.validate_fit_bootstrap(list(), f),
                  "fully anchored calibration")
    expect_error(.validate_dif_bootstrap(list(), f),
                  "fully anchored calibration")
    expect_error(.dif_boot_refit_ordinary(f$X, f, NULL),
                  "fully anchored calibration")
    expect_error(.scree_reference(f, k = 2, reps = 20),
                  "fully anchored calibration")
    expect_error(dimensionality_test(f, items_positive = paste0("I", 1:4),
      items_negative = paste0("I", 5:8), B = 20, workers = 1),
      "fully anchored calibration")
    expect_error(dimensionality_magnitude(f,
      list(paste0("I", 1:4), paste0("I", 5:8))),
      "fully anchored calibration")
  }
  expect_identical(fixed$thresholds, before)
  # Even an authenticated saved result must not supply a refit reference
  # unsupported for its calibration. Use signed wrappers to isolate this
  # model-eligibility check from the separate provenance checks.
  old_scree <- observed_scree
  attr(old_scree, "parallel") <- TRUE
  attr(old_scree, "result_signature") <- NULL
  attr(old_scree, "result_signature") <- .fit_boot_md5(old_scree)
  expect_error(.validate_scree_result(old_scree, fixed),
                "fully anchored calibration")
  old_dim <- observed_dim
  old_dim$bootstrap <- list(B = 20L)
  attr(old_dim, "result_signature") <- NULL
  attr(old_dim, "result_signature") <- .fit_boot_md5(old_dim)
  expect_error(.validate_dimensionality_test(old_dim, fixed),
                "fully anchored calibration")
})

test_that("partial anchoring still permits supported calibration refits", {
  set.seed(729)
  X <- matrix(rbinom(1200, 1, .5), 200, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  anchors <- data.frame(item = "I1", k = 1L, tau = .25)
  fit <- rasch(X, anchors = anchors)
  expect_gt(fit$est$n_parameters, 0)
  expect_true(.require_refittable_calibration(fit))
  reduced <- drop_items(fit, "I6")
  expect_identical(reduced$refit_spec$anchors, fit$refit_spec$anchors)
  expect_equal(reduced$items$location[reduced$items$item == "I1"], .25)
  refit <- .fit_refit_residuals(X, fit$model, fit$refit_spec$anchors,
                               fit$m, 60, 1e-8)
  expect_false(inherits(refit, "rasch_fit_boot_failure"))
  expect_gt(refit$est$n_parameters, 0)
  expect_equal(refit$est$thr$tau[1L], .25)
})
