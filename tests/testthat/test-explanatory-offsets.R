test_that("explanatory centring retains the weighted threshold origin", {
  m <- c(1L, 3L, 2L, 1L, 4L, 2L)
  thr <- threshold_index(m)
  weights <- 1 / (length(m) * m[thr$item])
  mm <- cbind(`(Intercept)` = 1, x = (0:5)[thr$item], k = thr$k)
  expected <- .explanatory_projector(m, thr) %*% mm
  actual <- .explanatory_centre(mm, weights)
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_equal(drop(crossprod(weights, actual)), c(`(Intercept)` = 0, x = 0, k = 0),
               tolerance = 1e-12)
  shifted <- sweep(mm, 2L, c(1e12, -1e12, 1e12), `+`)
  expect_identical(.explanatory_centre(shifted, weights), actual)
  expect_true(all(actual[, 1L] == 0))
})

test_that("Rasch explanatory fits are invariant to a predictor offset", {
  set.seed(917)
  x <- as.numeric(0:7)
  q <- data.frame(item = paste0("I", 1:8), x = x)
  theta <- rnorm(300)
  X <- vapply(x, function(v) rbinom(length(theta), 1,
    plogis(theta - .4 * (v - 3.5))), integer(length(theta)))
  colnames(X) <- q$item
  base <- rasch_explanatory(X, q, ~ x)
  for (offset in c(-1e12, 1e12)) {
    qs <- q; qs$x <- qs$x + offset
    expect_identical(outer(x, x, `-`), outer(qs$x, qs$x, `-`))
    f <- rasch_explanatory(X, qs, ~ x)
    expect_true(f$est$converged)
    expect_identical(f$est$B, base$est$B)
    expect_equal(f$est$beta, base$est$beta, tolerance = 1e-10)
    expect_equal(f$est$cov_beta, base$est$cov_beta, tolerance = 1e-10)
    expect_equal(f$thresholds$tau, base$thresholds$tau, tolerance = 1e-10)
    expect_equal(f$person$theta, base$person$theta, tolerance = 1e-10)
    expect_equal(f$est$loglik, base$est$loglik, tolerance = 1e-10)
    expect_equal(explanatory_test(f), explanatory_test(base), tolerance = 1e-10)
    expect_equal(explanatory_diagnostics(f)$p_adj,
                 explanatory_diagnostics(base)$p_adj, tolerance = 1e-10)
  }
  q$x <- 1e12
  expect_error(rasch_explanatory(X, q, ~ x), "no estimable variation")
})

test_that("CJ explanatory fits are invariant to a predictor offset", {
  set.seed(918)
  q <- data.frame(object = LETTERS[1:8], x = as.numeric(0:7))
  beta <- setNames(.4 * (q$x - 3.5), q$object)
  pr <- t(combn(q$object, 2))
  d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  base <- btl_explanatory(d, q, ~ x, "a", "b", winner = "win")
  for (offset in c(-1e12, 1e12)) {
    qs <- q; qs$x <- qs$x + offset
    f <- btl_explanatory(d, qs, ~ x, "a", "b", winner = "win")
    expect_true(f$converged)
    expect_identical(f$location_design, base$location_design)
    expect_equal(f$objects$location, base$objects$location, tolerance = 1e-10)
    expect_equal(f$cov_beta, base$cov_beta, tolerance = 1e-10)
    expect_equal(f$object_coefficients, base$object_coefficients, tolerance = 1e-10)
    expect_equal(f$loglik, base$loglik, tolerance = 1e-10)
    expect_equal(explanatory_test(f), explanatory_test(base), tolerance = 1e-10)
    expect_equal(explanatory_diagnostics(f)$p_adj,
                 explanatory_diagnostics(base)$p_adj, tolerance = 1e-10)
  }
  q$x <- 1e12
  expect_error(btl_explanatory(d, q, ~ x, "a", "b", winner = "win"),
                "no estimable variation")
})

test_that("centring does not change an interaction-only restriction", {
  X <- matrix(rep(c(0L, 1L), 40), 10, 8)
  colnames(X) <- LETTERS[1:8]
  q <- data.frame(item = colnames(X), x = as.numeric(0:7), z = rep(c(1, 2), 4))
  design <- .explanatory_metadata(q, ~ x:z, X)
  mm <- model.matrix(~ x:z, q)
  expected <- sweep(mm, 2L, colMeans(mm))[, "x:z", drop = FALSE]
  expect_equal(design$B, expected, ignore_attr = TRUE)
  qs <- q; qs$x <- qs$x + 10
  shifted <- .explanatory_metadata(qs, ~ x:z, X)
  # Adding an offset to x in x:z alone changes the model by adding an
  # unmodelled z main effect. It must not be silently centred away beforehand.
  expect_gt(max(abs(shifted$B - design$B)), 1)
})
