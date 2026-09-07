test_that("explanatory matching preserves exact item names before trimming", {
  items <- c(" I1", "I1", "I2 ")
  expect_identical(.explanatory_match_items(items, items), items)
  expect_identical(.explanatory_match_items(factor(items), items), items)
  expect_identical(.explanatory_match_items("I2", items), "I2 ")
  expect_error(.explanatory_match_items(" I1 ", items), "ambiguous item name")
  expect_true(is.na(.explanatory_match_items(NA_character_, items)))
  expect_true(is.na(.explanatory_match_items(NaN, "NaN")))
})

test_that("spaced item names work through explanatory estimation and relaxation", {
  set.seed(222)
  X <- matrix(rbinom(3200, 1, .5), 400, 8)
  colnames(X) <- paste0("I", 1:8)
  p <- data.frame(item = colnames(X), x = 0:7)
  base <- rasch_explanatory(X, p, ~ x)
  # The first two names have the same trimmed form but identify distinct items.
  colnames(X)[1:3] <- c(" I1", "I1", "I3 ")
  p$item <- colnames(X)
  actual <- rasch_explanatory(X, p[8:1, ], ~ x)
  expect_identical(actual$items$item, colnames(X))
  expect_identical(actual$explanatory$metadata$item, colnames(X))
  expect_equal(actual$est$B, base$est$B)
  expect_equal(actual$est$coefficients, base$est$coefficients)
  expect_equal(actual$est$cov_beta, base$est$cov_beta)
  expect_equal(actual$est$loglik, base$est$loglik)
  expect_equal(actual$person$theta, base$person$theta)
  relaxed <- relax_explanatory(actual, " I1")
  base_relaxed <- relax_explanatory(base, "I1")
  expect_identical(relaxed$explanatory$relaxations$item, " I1")
  expect_equal(relaxed$thresholds$tau, base_relaxed$thresholds$tau)
  expect_equal(relaxed$est$loglik, base_relaxed$est$loglik)
  expect_error(relax_explanatory(actual, " I1 "), "ambiguous item name")
  p$item[3] <- "I3" # An unambiguous whitespace-normalised match is supported.
  fallback <- rasch_explanatory(X, p, ~ x)
  expect_equal(fallback$est$B, actual$est$B)
  expect_identical(fallback$explanatory$source_predictors$item[3], "I3 ")
  p$item[1] <- "I1"
  expect_error(rasch_explanatory(X, p, ~ x), "exactly one row per item")
})

test_that("threshold-level metadata retains distinct spaced item names", {
  items <- c(" I1", "I1", "I3 ", "I4")
  X <- matrix(rep(0:2, 4), 3, 4, dimnames = list(NULL, items))
  p <- data.frame(item = rep(items, each = 2), threshold = rep(1:2, 4),
                   x = rep(0:3, each = 2))
  actual <- .explanatory_metadata(p[8:1, ], ~ x + threshold, X, "threshold")
  expect_identical(actual$metadata$item, p$item)
  ordinary <- X
  colnames(ordinary) <- paste0("I", 1:4)
  p$item <- rep(colnames(ordinary), each = 2)
  expected <- .explanatory_metadata(p, ~ x + threshold, ordinary, "threshold")
  expect_equal(actual$B, expected$B)
})
