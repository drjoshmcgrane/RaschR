.explanatory_keyed_fixture <- function() {
  set.seed(223)
  X <- matrix(rbinom(1200, 1, .5), 200, 6)
  colnames(X) <- paste0("I", 1:6)
  X[7, ] <- NA
  raw <- matrix(ifelse(is.na(X), NA_character_,
                       ifelse(X == 1, "A", "B")), nrow(X),
                 dimnames = dimnames(X))
  predictors <- data.frame(item = colnames(X), w = rep(0:2, each = 2))
  fit <- rasch_explanatory(raw, predictors, ~ w,
    factors = data.frame(group = rep(c("a", "b"), 100)),
    key = setNames(rep("A", 6), colnames(X)))
  list(fit = fit, X = X, predictors = predictors)
}

test_that("observed explanatory refits align keyed responses with person rows", {
  f <- .explanatory_keyed_fixture()$fit
  rows <- rev(which(is.finite(f$person$theta)))
  source <- f$X[rows, , drop = FALSE]
  source[1, 1] <- NA
  actual <- .explanatory_refit_modified(f, source, person_rows = rows)
  expected <- f$mc$raw[rows, , drop = FALSE]
  expected[is.na(source)] <- NA_character_
  expect_equal(actual$mc$raw, expected, ignore_attr = TRUE)
  expect_identical(actual$mc$map, f$mc$map)
  expect_identical(actual$person$id, f$person$id[rows])
  expect_equal(actual$factors, f$factors[rows, , drop = FALSE],
               ignore_attr = TRUE)
  expect_identical(actual$X, source)
  for (bad in list(NA_real_, 0, .5, 201, Inf, "1"))
    expect_error(.explanatory_refit_modified(f, source[1, , drop = FALSE],
                                            person_rows = bad),
                 "valid fitted row indices")
})

test_that("simulated explanatory refits do not inherit observed answer options", {
  f <- .explanatory_keyed_fixture()$fit
  X <- 1L - f$X
  actual <- .explanatory_refit_modified(f, X, inherit_mc = FALSE)
  expect_true(actual$est$converged)
  expect_null(actual$mc)
  expect_identical(actual$X, X)
  expect_equal(actual$est$B, f$est$B)
  expect_error(.explanatory_refit_modified(f, X),
               "cannot inherit observed multiple-choice answers for changed scores")
  # The DIF worker must opt out too, even when it retains every person row.
  real_refit <- .explanatory_refit_modified
  result <- testthat::with_mocked_bindings(
    .dif_boot_refit_explanatory(X, f, list()),
    .explanatory_refit_modified = function(..., inherit_mc = TRUE) {
      expect_false(inherit_mc)
      real_refit(..., inherit_mc = inherit_mc)
    },
    dif_anova = function(fit, ...) {
      expect_null(fit$mc)
      expect_identical(fit$X, X)
      data.frame(p = .5)
    }, .package = "rasch")
  expect_equal(result$dif$p, .5)
})

test_that("keyed explanatory scree references survive blank person rows", {
  d <- .explanatory_keyed_fixture()
  numeric_fit <- rasch_explanatory(d$X, d$predictors, ~ w,
                                   factors = d$fit$factors)
  expect_equal(sum(is.finite(d$fit$person$theta)), 199L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  keyed <- plot_scree(d$fit, n_components = 2, reps = 20, seed = 224)
  numeric <- plot_scree(numeric_fit, n_components = 2, reps = 20, seed = 224)
  expect_true(all(keyed$n_reference == 20L))
  expect_true(all(keyed$n_reference_errors == 0L))
  expect_true(all(keyed$n_reference_nonconverged == 0L))
  for (nm in c("eigenvalue", "reference_mean", "reference_critical",
                "parallel_p", "parallel_p_adj", "parallel_significant"))
    expect_equal(keyed[[nm]], numeric[[nm]], tolerance = 1e-12)
})
