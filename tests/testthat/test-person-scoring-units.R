test_that("common-unit WLE locations and SEs are scale equivariant", {
  tau <- list(c(-1.2, -.3), c(-.5, .9), c(.2, 1.1))
  reference <- person_wle(tau)
  for (d in c(1e-200, 1e-9, 1e-3, 1e3, 1e9, 1e12, 1e200)) {
    actual <- person_wle(lapply(tau, function(t) t / d), disc = d)
    expect_true(all(is.finite(actual$theta)))
    expect_true(all(is.finite(actual$se) & actual$se > 0))
    expect_equal(actual$theta * d, reference$theta, tolerance = 1e-8)
    expect_equal(actual$se * d, reference$se, tolerance = 1e-8)
    # Check the weighted score equation on the unit-one representation,
    # independently of the root search's reported convergence.
    for (r in 0:6) {
      mo <- lapply(tau, item_moments, theta = actual$theta[r + 1L] * d)
      E <- sum(vapply(mo, `[[`, 0, "E"))
      V <- sum(vapply(mo, `[[`, 0, "V"))
      m3 <- sum(vapply(mo, `[[`, 0, "mu3"))
      expect_lt(abs(r - E + m3 / (2 * V)), 1e-7)
    }
  }
})

test_that("ML conversion and extrapolated scores retain their measurement unit", {
  set.seed(9581)
  X <- matrix(rbinom(300 * 6, 1, .5), 300, 6,
               dimnames = list(NULL, paste0("I", 1:6)))
  X[1:2, ] <- rep(c(0L, 1L), 6)
  X[3:4, 1:4] <- rep(c(0L, 1L), 4)
  X[3:4, 5:6] <- NA_integer_
  fit <- rasch(X)
  base_ml <- score_table(fit, method = "mle")
  base_ext <- lapply(c("wle", "mle"), function(method)
    score_table(fit, method = method, extremes = "extrapolated"))
  base_persons <- person_extrapolated(fit)

  for (d in c(1e-200, 1e-9, 1e9, 1e12, 1e200)) {
    # Express the scoring fields of the same calibration in another unit:
    # d * (theta/d - tau/d) gives exactly the original response model.
    scaled <- fit
    scaled$disc <- rep(d, ncol(X))
    scaled$tau_list <- lapply(fit$tau_list, function(t) t / d)
    scaled$score_table$theta <- fit$score_table$theta / d
    scaled$score_table$se <- fit$score_table$se / d
    scaled$person$theta <- fit$person$theta / d
    scaled$person$se <- fit$person$se / d
    ml <- score_table(scaled, method = "mle")
    expect_identical(is.na(ml$theta), is.na(base_ml$theta))
    expect_equal(ml$theta * d, base_ml$theta, tolerance = 1e-8)
    expect_equal(ml$se * d, base_ml$se, tolerance = 1e-8)
    for (j in 1:2) {
      actual <- score_table(scaled, method = c("wle", "mle")[j],
                             extremes = "extrapolated")
      expect_equal(actual$theta * d, base_ext[[j]]$theta, tolerance = 1e-8)
      expect_equal(actual$se * d, base_ext[[j]]$se, tolerance = 1e-8)
      expect_identical(actual$extrapolated, base_ext[[j]]$extrapolated)
    }
    persons <- person_extrapolated(scaled)
    expect_equal(persons$theta_extrapolated * d,
                 base_persons$theta_extrapolated, tolerance = 1e-8)
    expect_equal(persons$se_extrapolated * d,
                 base_persons$se_extrapolated, tolerance = 1e-8)
  }
})
