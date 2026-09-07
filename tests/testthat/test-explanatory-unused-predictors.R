test_that("unused ordinal metadata leaves Rasch explanatory fits unchanged", {
  set.seed(220)
  X <- matrix(rbinom(3200, 1, .5), 400, 8)
  colnames(X) <- paste0("I", 1:8)
  p <- data.frame(item = colnames(X), x = 0:7)
  base <- rasch_explanatory(X, p, ~ x)
  p$unused <- ordered(rep("A", 8), levels = c("A", "B"))
  p$absent <- ordered(rep(NA_character_, 8), levels = c("A", "B"))
  actual <- rasch_explanatory(X, p, ~ x)
  expect_true(actual$est$converged)
  expect_identical(actual$X, base$X)
  expect_identical(actual$est$B, base$est$B)
  expect_equal(actual$est$coefficients, base$est$coefficients)
  expect_equal(actual$est$cov_beta, base$est$cov_beta)
  expect_equal(actual$est$loglik, base$est$loglik)
  expect_equal(actual$person$theta, base$person$theta)
  expect_equal(explanatory_test(actual), explanatory_test(base))
  expect_equal(explanatory_diagnostics(actual), explanatory_diagnostics(base))
  expect_true(all(c("unused", "absent") %in%
                    names(actual$explanatory$metadata)))
  expect_error(rasch_explanatory(X, p, ~ x + unused),
               "ordinal predictor needs at least two observed levels")
  expect_error(rasch_explanatory(X, p, ~ x + absent), "missing values")
})

test_that("unused ordinal metadata leaves CJ explanatory fits unchanged", {
  set.seed(221)
  p <- data.frame(object = LETTERS[1:8], x = 0:7)
  pairs <- t(combn(p$object, 2))
  d <- data.frame(a = rep(pairs[, 1], each = 30),
                  b = rep(pairs[, 2], each = 30))
  d$winner <- ifelse(runif(nrow(d)) < .5, d$a, d$b)
  base <- btl_explanatory(d, p, ~ x, "a", "b", winner = "winner")
  p$unused <- ordered(rep("A", 8), levels = c("A", "B"))
  p$absent <- ordered(rep(NA_character_, 8), levels = c("A", "B"))
  actual <- btl_explanatory(d, p, ~ x, "a", "b", winner = "winner")
  expect_true(actual$converged)
  expect_identical(actual$comparisons, base$comparisons)
  expect_identical(actual$location_design, base$location_design)
  expect_equal(actual$object_coefficients, base$object_coefficients)
  expect_equal(actual$cov_beta, base$cov_beta)
  expect_equal(actual$loglik, base$loglik)
  expect_equal(actual$objects$location, base$objects$location)
  expect_equal(explanatory_test(actual), explanatory_test(base))
  expect_equal(explanatory_diagnostics(actual), explanatory_diagnostics(base))
  expect_true(all(c("unused", "absent") %in%
                    names(actual$explanatory$metadata)))
  expect_error(btl_explanatory(d, p, ~ x + unused, "a", "b",
                               winner = "winner"),
               "ordinal predictor needs at least two observed levels")
  expect_error(btl_explanatory(d, p, ~ x + absent, "a", "b",
                               winner = "winner"), "missing values")
})

test_that("selected ordinal interactions retain adjacent contrasts", {
  p <- expand.grid(x = 0:3, grade = c("low", "middle", "high"))
  p$grade <- ordered(p$grade, levels = c("low", "middle", "high"))
  p$item <- paste0("I", seq_len(nrow(p)))
  X <- matrix(rep(0:1, nrow(p)), 2, nrow(p),
               dimnames = list(NULL, p$item))
  q <- p
  contrasts(q$grade) <- .explanatory_ordinal_contrasts(q$grade)
  expected <- model.matrix(~ x * grade, q)
  expected <- .explanatory_centre(expected)[, -1, drop = FALSE]
  p$unused <- ordered(rep("A", nrow(p)))
  r <- .explanatory_metadata(p, ~ x * grade, X)
  expect_equal(r$B, expected, ignore_attr = TRUE)
  expect_true(any(grepl("adjacent_", colnames(r$B))))
  names(p)[names(p) == "item"] <- "object"
  b <- .btl_explanatory_design(p, ~ x * grade, p$object)
  expect_equal(b$B, expected, ignore_attr = TRUE)
  expect_identical(colnames(b$B), colnames(r$B))
  expect_error(.btl_explanatory_design(p, ~ . - object, p$object),
               "ordinal predictor needs at least two observed levels")
})
