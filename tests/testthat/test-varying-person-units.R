.varying_unit_fixture <- function() {
  X <- as.matrix(expand.grid(A = 0:1, B = 0:2, C = 0:3))
  X <- rbind(X, c(NA, 0, 0), c(NA, 1, 2), c(1, NA, 3),
             c(0, NA, NA), c(NA, NA, NA))
  tau <- list(A = -.3, B = c(-.7, .5), C = c(-1, .2, 1))
  list(X = X, tau_list = tau, disc = c(.6, 1.1, 1.5),
       person = data.frame(id = seq_len(nrow(X))), factors = NULL,
       est = list(converged = TRUE))
}

test_that("unequal-unit WLE preserves locations, SEs and the score equation", {
  f <- .varying_unit_fixture()
  reference <- .efrm_person_estimates(f$X, f$tau_list, f$disc)
  for (scale in c(1e-200, 1e-9, 1, 1e9, 1e12, 1e200)) {
    actual <- .efrm_person_estimates(
      f$X, lapply(f$tau_list, function(t) t / scale), f$disc * scale)
    expect_identical(is.finite(actual$theta), rowSums(!is.na(f$X)) > 0)
    expect_equal(actual$theta * scale, reference$theta, tolerance = 1e-8)
    expect_equal(actual$se * scale, reference$se, tolerance = 1e-8)
    expect_equal(actual$weighted_score / scale, reference$weighted_score,
                 tolerance = 1e-12)
    expect_identical(actual$extreme, reference$extreme)
    for (i in which(is.finite(actual$theta))) {
      cols <- which(!is.na(f$X[i, ]))
      r <- f$disc[cols]
      mo <- lapply(cols, function(j) item_moments(
        actual$theta[i] * scale, f$tau_list[[j]], f$disc[j]))
      E <- vapply(mo, `[[`, 0, "E")
      V <- vapply(mo, `[[`, 0, "V")
      m3 <- vapply(mo, `[[`, 0, "mu3")
      expect_lt(abs(sum(r * (f$X[i, cols] - E)) +
        sum(r^3 * m3) / (2 * sum(r^2 * V))), 1e-7)
      expect_equal(actual$se[i] * scale, 1 / sqrt(sum(r^2 * V)),
                   tolerance = 1e-10)
    }
  }
})

test_that("external weights retain their equation under changes of frame unit", {
  f <- .varying_unit_fixture()
  class(f) <- "rasch"
  for (w in list(c(A = 1, B = 1, C = 1), c(A = .5, B = 2, C = 0),
                 c(A = 1, B = 3, C = .25))) {
    reference <- weighted_person_estimates(f, w)
    if (all(w == 1)) {
      ordinary <- .efrm_person_estimates(f$X, f$tau_list, f$disc)
      expect_equal(reference$theta, ordinary$theta, tolerance = 1e-8)
      expect_equal(reference$se, ordinary$se, tolerance = 1e-8)
    }
    for (scale in c(1e-200, 1e-9, 1, 1e9, 1e12, 1e200)) {
      scaled <- f
      scaled$tau_list <- lapply(f$tau_list, function(t) t / scale)
      scaled$disc <- f$disc * scale
      actual <- weighted_person_estimates(scaled, w)
      expect_identical(is.finite(actual$theta),
                       rowSums(!is.na(f$X[, w > 0, drop = FALSE])) > 0)
      expect_equal(actual$theta * scale, reference$theta, tolerance = 1e-8)
      expect_equal(actual$se * scale, reference$se, tolerance = 1e-8)
      expect_equal(actual$weighted_score / scale, reference$weighted_score,
                   tolerance = 1e-12)
      expect_identical(actual$extreme, reference$extreme)
      for (i in which(is.finite(actual$theta))) {
        cols <- which(!is.na(f$X[i, ]) & w > 0)
        r <- f$disc[cols]; q <- w[cols]
        mo <- lapply(cols, function(j) item_moments(
          actual$theta[i] * scale, f$tau_list[[j]], f$disc[j]))
        E <- vapply(mo, `[[`, 0, "E")
        V <- vapply(mo, `[[`, 0, "V")
        m3 <- vapply(mo, `[[`, 0, "mu3")
        H <- sum(q * r^2 * V); J <- sum(q^2 * r^2 * V)
        expect_lt(abs(sum(q * r * (f$X[i, cols] - E)) +
          J * sum(q * r^3 * m3) / (2 * H^2)), 1e-7)
        expect_equal(actual$se[i] * scale, sqrt(J) / H,
                     tolerance = 1e-10)
      }
    }
  }
})
