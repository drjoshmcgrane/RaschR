.dm_sample_fixture <- function() {
  set.seed(917)
  n <- 1200L
  theta <- c(rnorm(n / 2, sd = 0.6), rnorm(n / 2, sd = 3))
  delta <- rep(seq(-1, 1, length.out = 4), 2)
  X <- sapply(delta, function(d) rbinom(n, 1, plogis(theta - d)))
  colnames(X) <- paste0("I", 1:8)
  X[601:1200, c(1, 5)] <- NA
  X
}

test_that("dimensionality PSI compares the same complete response rows", {
  fit <- rasch(.dm_sample_fixture())
  original <- fit
  z <- dimensionality_magnitude(fit, list(paste0("I", 1:4), paste0("I", 5:8)))
  keep <- complete.cases(fit$X) & complete.cases(z$refit$X)
  p <- fit$person; q <- z$refit$person
  keep <- keep & is.finite(p$theta) & is.finite(q$theta) &
    is.finite(p$se) & p$se >= 0 & is.finite(q$se) & q$se >= 0
  psi <- function(theta, se) max(1 - mean(se^2) / var(theta), 0)
  r1 <- psi(p$theta[keep], p$se[keep])
  r2 <- psi(q$theta[keep], q$se[keep])
  old_c2 <- 2 * (fit$psi$PSI / z$refit$psi$PSI - 1) / (6 / 7)
  expect_gt(old_c2, 4)
  tab <- as.data.frame(z$table)
  expect_equal(tab$run1[1], r1)
  expect_equal(tab$subtest[1], r2)
  expect_equal(tab$c2[1], max(2 * (r1 / r2 - 1) / (6 / 7), 0))
  expect_lt(tab$c2[1], 0.5)
  expect_equal(tab$n, c(600L, 600L))
  expect_equal(tab$n_excluded, c(600L, 600L))
  expect_equal(tab$run1[2], fit$alpha$alpha)
  expect_equal(tab$subtest[2], z$refit$alpha$alpha)
  expect_identical(fit, original)
})

test_that("partial subscale coverage is not a matched complete-test sample", {
  X <- .dm_sample_fixture()
  # These people retain one complete subscale and can be scored after the
  # refit, but that is still different coverage from their seven original items.
  X[601:1200, 5] <- X[601:1200, 6]
  fit <- rasch(X)
  z <- dimensionality_magnitude(fit, list(paste0("I", 1:4), paste0("I", 5:8)))
  expect_gt(sum(is.finite(z$refit$person$theta)), 600L)
  expect_equal(z$table$n[1], 600L)
  expect_equal(z$table$n_excluded[1], 600L)
})

test_that("complete data retain their original reliability comparison", {
  set.seed(918)
  X <- matrix(rbinom(400 * 8, 1, .5), 400, 8,
              dimnames = list(NULL, paste0("I", 1:8)))
  fit <- rasch(X)
  z <- dimensionality_magnitude(fit, list(paste0("I", 1:4), paste0("I", 5:8)))
  expect_equal(z$table$run1, c(fit$psi$PSI, fit$alpha$alpha))
  expect_equal(z$table$subtest, c(z$refit$psi$PSI, z$refit$alpha$alpha))
  expect_equal(z$table$n, c(400L, 400L))
  expect_equal(z$table$n_excluded, c(0L, 0L))
  expect_identical(z$algorithm, "complete-panel-1")

  refit <- z$refit
  refit$person$se[1] <- Inf
  matched <- with_mocked_bindings(
    dimensionality_magnitude(fit, list(paste0("I", 1:4), paste0("I", 5:8))),
    combine_items = function(...) refit, .package = "rasch")
  expect_equal(matched$table$n, c(399L, 400L))
  expect_equal(matched$table$run1[1],
    max(1 - mean(fit$person$se[-1]^2) / var(fit$person$theta[-1]), 0))
  refit$person$se <- c(z$refit$person$se[1:2], rep(Inf, 398))
  matched <- with_mocked_bindings(
    dimensionality_magnitude(fit, list(paste0("I", 1:4), paste0("I", 5:8))),
    combine_items = function(...) refit, .package = "rasch")
  expect_equal(matched$table$n[1], 2L)
  expect_true(is.na(matched$table$c2[1]))

  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), settings = list(),
    results = list(dimension_magnitude = z)))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  .save_app_project(project, path)
  expect_identical(.read_app_project(path)$results$dimension_magnitude, z)

  old <- project
  old$results$dimension_magnitude$algorithm <- NULL
  old <- .seal_app_project(old)
  saveRDS(old, path)
  expect_warning(restored <- .read_app_project(path), "earlier reliability")
  expect_null(restored$results$dimension_magnitude)
  expect_identical(restored$base_fit, fit)
  expect_match(attr(restored, "rasch_project_legacy_dropped"),
               "unmatched reliability samples")
  expect_no_error(.validate_app_project(restored))

  old$results$dimension_magnitude$table$run1[1] <- .99
  saveRDS(old, path)
  expect_error(.read_app_project(path), "changed since they were saved")

  old$schema <- 1L
  old$binding <- NULL
  saveRDS(old, path)
  expect_warning(restored <- .read_app_project(path), "schema-1")
  expect_null(restored$results$dimension_magnitude)
})

test_that("no complete response panel gives unavailable magnitudes", {
  set.seed(919)
  # Three overlapping subscales keep the calibration connected without
  # giving any row a complete response panel.
  X <- matrix(rbinom(400 * 12, 1, .5), 400, 12,
              dimnames = list(NULL, paste0("I", 1:12)))
  X[seq(1, 400, 3), 1] <- NA
  X[seq(2, 400, 3), 5] <- NA
  X[seq(3, 400, 3), 9] <- NA
  fit <- rasch(X)
  z <- dimensionality_magnitude(fit,
    list(paste0("I", 1:4), paste0("I", 5:8), paste0("I", 9:12)))
  expect_equal(z$table$n, c(0L, 0L))
  expect_equal(z$table$n_excluded, c(400L, 400L))
  expect_true(all(is.na(z$table$c2)))
  expect_true(all(is.na(z$table$rho)))
})
