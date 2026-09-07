test_that("failed EFRM stage one withholds every exposed covariance", {
  set.seed(1)
  n <- 240L
  X <- matrix(rbinom(n * 4L, 1L, plogis(rnorm(n))), n, 4L)
  colnames(X) <- paste0("I", 1:4)
  d <- data.frame(X, group = rep(c("A", "B"), each = n / 2L))
  expect_warning(f <- rasch_efrm(
    d, list(all = colnames(X)), "group", boot_reps = 0, maxit = 1L),
    "did NOT converge")
  expect_false(f$est$stage1_converged)
  expect_false(f$est$converged)
  expect_false(f$unit_support$phi_inference)
  expect_false(f$unit_support$alpha_inference)
  expect_true(all(is.na(f$est$cov_tau)))
  expect_true(all(is.na(f$est$thr$se)))
  for (nm in setdiff(names(f$unit_cov), "method"))
    if (!is.null(f$unit_cov[[nm]]))
      expect_true(all(is.na(f$unit_cov[[nm]])), info = nm)
  expect_true(all(is.na(f$score_curves$sem)))
  expect_true(all(is.finite(f$score_curves$expected_score)))
  expect_match(paste(f$notes, collapse = " "),
               "probabilities are withheld because the conditional calibration did not converge")
})

test_that("converged weak EFRM group units remain descriptive", {
  # With almost no item spread, the relative unit can be estimated but is
  # much too uncertain for its normal/Wald reference. This ordinary fit
  # previously warned about weak identification while still reporting p.
  set.seed(16)
  n <- 2000L
  X <- matrix(rbinom(n * 2L, 1L, 0.5), n, 2L)
  colnames(X) <- c("I1", "I2")
  d <- data.frame(X, group = rep(c("A", "B"), each = n / 2L))
  expect_warning(f <- rasch_efrm(
    d, list(all = colnames(X)), "group", boot_reps = 0),
    "only weakly identified")
  expect_true(f$est$converged)
  expect_true(all(f$phi_table$se_log_phi > 5))
  expect_true(all(f$unit_support$group$n_persons >= 50))
  expect_false(f$unit_support$phi_inference)
  expect_true(all(is.finite(f$unit_cov$cov_log_phi)))
  # The one-set alpha is fixed at one, not estimated through a group link.
  expect_true(f$unit_support$alpha_inference)
  expect_true(all(is.finite(f$phi_table$phi)))
  for (part in c("unit_tests", "unit_omnibus")) {
    expect_true(all(is.na(f$efrm_vs_rasch[[part]]$p)))
    expect_true(all(is.na(f$efrm_vs_rasch[[part]]$p_adj)))
    expect_true(all(is.na(f$efrm_vs_rasch[[part]]$significant)))
  }
  interval_drawn <- FALSE
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  testthat::with_mocked_bindings(
    plot_frames(f),
    segments = function(...) interval_drawn <<- TRUE,
    .package = "rasch")
  expect_false(interval_drawn)
})

test_that("EFRM weak-unit decisions use the returned bootstrap SE", {
  d <- simulate_efrm(n_per_group = 150, items_per_set = 5,
                     n_sets = 1, n_groups = 2, seed = 4)
  original_solve <- rasch:::.efrm_solve
  # Only the intermediate reported SE changes; all estimates and covariance
  # matrices, including the actual full-bootstrap draws, remain unmodified.
  testthat::local_mocked_bindings(
    .efrm_solve = function(...) {
      z <- original_solve(...)
      z$se_log_phi[] <- 6
      z
    }, .package = "rasch")
  expect_no_warning(f <- rasch_efrm(
    d, attr(d, "truth")$item_sets, "group", se_method = "bootstrap",
    boot_reps = 30, workers = 1, seed = 41))
  expect_identical(f$se_method, "bootstrap")
  expect_true(f$est$converged)
  expect_true(all(f$phi_table$se_log_phi < 5))
  expect_true(f$unit_support$phi_inference)
  expect_true(all(is.finite(f$efrm_vs_rasch$unit_tests$p)))
  expect_false(any(grepl("weakly identified unit", f$notes)))
})

test_that("EFRM set inference follows the group units used in its link", {
  d <- simulate_efrm(n_per_group = 200, items_per_set = 4,
                     n_sets = 2, n_groups = 2, seed = 93)
  original_solve <- rasch:::.efrm_solve
  testthat::local_mocked_bindings(
    .efrm_solve = function(...) {
      z <- original_solve(...)
      z$se_log_phi[] <- 6
      z
    }, .package = "rasch")
  expect_warning(f <- rasch_efrm(
    d, attr(d, "truth")$item_sets, "group", boot_reps = 40,
    workers = 1, seed = 49), "only weakly identified")
  expect_true(f$est$converged)
  expect_true(all(f$unit_support$set$n_common_persons >= 50))
  expect_true(all(is.finite(f$alpha_table$se_log_alpha)))
  expect_true(all(is.finite(f$unit_cov$cov_log_alpha)))
  expect_false(f$unit_support$phi_inference)
  expect_false(f$unit_support$alpha_inference)
  expect_true(all(is.na(f$efrm_vs_rasch$unit_omnibus$p_adj)))
  expect_true(all(is.na(f$efrm_vs_rasch$unit_tests$p_adj)))
  expect_match(paste(f$notes, collapse = " "), "their links depend on")
})

test_that("EFRM failed set links retain only stage-one uncertainty", {
  d <- simulate_efrm(n_per_group = 120, items_per_set = 6,
                     n_sets = 2, n_groups = 1, seed = 1901)
  original_mass <- rasch:::efrm_fit_weights_cpp
  testthat::local_mocked_bindings(
    efrm_fit_weights_cpp = function(...) {
      z <- original_mass(...)
      z$converged <- FALSE
      z
    }, .package = "rasch")
  expect_warning(f <- rasch_efrm(
    d, attr(d, "truth")$item_sets, "group", boot_reps = 0, workers = 1),
    "nuisance masses")
  expect_true(f$est$stage1_converged)
  expect_false(f$est$converged)
  expect_true(all(is.finite(f$unit_cov$cov_dtilde)))
  expect_true(all(is.finite(f$unit_cov$cov_log_phi)))
  expect_true(all(is.na(f$unit_cov$cov_delta)))
  expect_true(all(is.na(f$score_curves$sem)))
  expect_true(all(is.finite(f$score_curves$expected_score)))
})
