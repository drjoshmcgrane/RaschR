check_rasch_units <- function(base, scaled, factors) {
  expect_true(scaled$est$converged)
  expect_identical(scaled$est$cluster_inference, base$est$cluster_inference)
  expect_equal(scaled$est$thr$tau, base$est$thr$tau, tolerance = 1e-7)
  expect_equal(scaled$est$cov_tau, base$est$cov_tau, tolerance = 1e-7)
  expect_equal(scaled$est$beta * factors, base$est$beta, tolerance = 1e-7)
  expect_equal(scaled$est$cov_beta * outer(factors, factors),
               base$est$cov_beta, tolerance = 1e-7)
  expect_equal(scaled$est$H_beta / outer(factors, factors),
               base$est$H_beta, tolerance = 1e-7)
  expect_equal(scaled$est$loglik, base$est$loglik, tolerance = 1e-9)
  expect_equal(scaled$est$coefficients$p_adj, base$est$coefficients$p_adj,
               tolerance = 1e-7)
  expect_equal(explanatory_test(scaled)$p, explanatory_test(base)$p,
               tolerance = 1e-7)
  expect_equal(explanatory_test(scaled)$df, explanatory_test(base)$df)
  expect_equal(explanatory_test(scaled)$r_squared_adj,
               explanatory_test(base)$r_squared_adj, tolerance = 1e-7)
  expect_equal(scaled$person$theta, base$person$theta, tolerance = 1e-7)
  expect_true(all(is.finite(scaled$est$coefficients$p_adj)))
  expect_equal(.cl_ic(scaled), .cl_ic(base), tolerance = 1e-7)
}

test_that("LLTM convergence and inference are independent of predictor units", {
  set.seed(622)
  q <- data.frame(item = paste0("I", 1:12),
                   x = seq(-1.5, 1.5, length.out = 12))
  delta <- 2 * q$x + 3 * sin(3 * q$x)
  th <- rnorm(800)
  X <- vapply(delta, function(d) rbinom(length(th), 1, plogis(th - d)),
               integer(length(th)))
  colnames(X) <- q$item
  base <- rasch_explanatory(X, q, ~ x)
  for (s in c(1e-12, 1e-6, 1e9)) {
    qs <- q; qs$x <- qs$x * s
    f <- rasch_explanatory(X, qs, ~ x)
    check_rasch_units(base, f, setNames(s, "x"))
  }
  qs <- q; qs$x <- qs$x * 1e9
  relaxed <- relax_explanatory(base, "I1")
  scaled_relaxed <- relax_explanatory(rasch_explanatory(X, qs, ~ x), "I1")
  expect_equal(explanatory_test(scaled_relaxed)$p, explanatory_test(relaxed)$p,
               tolerance = 1e-7)
  expect_equal(explanatory_diagnostics(rasch_explanatory(X, qs, ~ x))$p_adj,
               explanatory_diagnostics(base)$p_adj, tolerance = 1e-7)
  design <- .explanatory_metadata(qs, ~ x, base$X)
  short <- .pcml_design(base$X, design$B, maxit = 1)
  expect_false(short$converged)
  expect_true(all(is.na(short$coefficients$p_adj)))
  for (s in c(1e-200, 1e200)) {
    qs <- q; qs$x <- qs$x * s
    expect_error(.explanatory_metadata(qs, ~ x, base$X),
                  "outside the representable covariance range")
  }
})

test_that("LPCM mixed-unit designs retain the person-cluster covariance", {
  set.seed(624)
  q <- data.frame(item = paste0("I", 1:8), x = rep(c(-1, -.3, .3, 1), 2))
  theta <- rep(rnorm(160), each = 2)
  X <- vapply(q$x, function(x) vapply(theta, function(th)
    sample.int(3, 1, prob = item_moments(th, c(-.8, .8) + .4 * x)$P) - 1L,
    integer(1)), integer(length(theta)))
  colnames(X) <- q$item
  id <- rep(seq_len(160), each = 2)
  base <- rasch_explanatory(X, q, ~ x + threshold, id = id)
  expect_identical(base$est$cluster_support$correction,
                   "linearised delete-one-person jackknife")
  for (s in c(1e-9, 1e9)) {
    qs <- q; qs$x <- qs$x * s
    f <- rasch_explanatory(X, qs, ~ x + threshold, id = id)
    factors <- setNames(rep(1, length(base$est$beta)), names(base$est$beta))
    factors["x"] <- s
    check_rasch_units(base, f, factors)
  }
})

test_that("explanatory CJ preserves locations and uncertainty across units", {
  set.seed(623)
  q <- data.frame(object = LETTERS[1:8], x = seq(-1.5, 1.5, length.out = 8))
  beta <- setNames(1.2 * q$x, q$object)
  pr <- t(combn(q$object, 2))
  d <- data.frame(a = rep(pr[, 1], each = 60), b = rep(pr[, 2], each = 60))
  eta <- beta[d$a] - beta[d$b]
  d$winner <- ifelse(runif(nrow(d)) < plogis(eta), d$a, d$b)
  d$grade <- vapply(eta, function(v)
    sample.int(3, 1, prob = item_moments(v + .3, c(-.6, .6))$P) - 1L,
    integer(1))
  d$judge <- rep(seq_len(20), length.out = nrow(d))
  for (graded in c(FALSE, TRUE)) {
    args <- if (graded) list(response = "grade", judge = "judge", position = TRUE)
            else list(winner = "winner")
    fit <- function(qs) do.call(btl_explanatory,
      c(list(data = d, predictors = qs, formula = ~ x,
             object_a = "a", object_b = "b"), args))
    base <- fit(q)
    for (s in c(1e-12, 1e-6, 1e9)) {
      qs <- q; qs$x <- qs$x * s
      f <- fit(qs)
      expect_true(f$converged)
      expect_true(all(is.finite(f$object_coefficients$p_adj)))
      expect_equal(f$fitted_prob, base$fitted_prob, tolerance = 1e-7)
      expect_equal(f$objects$location, base$objects$location, tolerance = 1e-7)
      expect_equal(f$cov_beta, base$cov_beta, tolerance = 1e-7)
      expect_equal(f$loglik, base$loglik, tolerance = 1e-9)
      expect_equal(f$object_coefficients$estimate * s,
                   base$object_coefficients$estimate, tolerance = 1e-7)
      expect_equal(f$object_coefficients$se * s,
                   base$object_coefficients$se, tolerance = 1e-7)
      expect_equal(f$object_coefficients$p_adj, base$object_coefficients$p_adj,
                   tolerance = 1e-7)
      factors <- c(s, rep(1, nrow(base$sensitivity) - 1L))
      expect_equal(f$sensitivity / outer(factors, factors), base$sensitivity,
                   tolerance = 1e-7)
      expect_equal(f$cov_parameters * outer(factors, factors), base$cov_parameters,
                   tolerance = 1e-7)
      expect_equal(explanatory_test(f)$p, explanatory_test(base)$p,
                   tolerance = 1e-7)
      expect_equal(explanatory_test(f)$df, explanatory_test(base)$df)
      expect_equal(explanatory_test(f)$r_squared_adj,
                   explanatory_test(base)$r_squared_adj, tolerance = 1e-7)
      expect_equal(.cl_ic(f), .cl_ic(base), tolerance = 1e-7)
      if (graded) {
        expect_equal(f$thresholds$tau, base$thresholds$tau, tolerance = 1e-7)
        expect_equal(f$dependence$estimate, base$dependence$estimate,
                     tolerance = 1e-7)
      }
    }
  }
  qs <- q; qs$x <- qs$x * 1e-6
  base <- btl_explanatory(d, q, ~ x, "a", "b", winner = "winner")
  scaled <- btl_explanatory(d, qs, ~ x, "a", "b", winner = "winner")
  expect_equal(explanatory_diagnostics(scaled)$p_adj,
               explanatory_diagnostics(base)$p_adj, tolerance = 1e-7)
  expect_equal(explanatory_test(relax_btl_explanatory(scaled, "A"))$p,
               explanatory_test(relax_btl_explanatory(base, "A"))$p,
               tolerance = 1e-7)
  design <- .btl_explanatory_design(qs, ~ x, q$object)
  short <- suppressWarnings(btl(d, "a", "b", winner = "winner",
                                  .object_design = design, maxit = 1))
  expect_false(short$converged)
  expect_true(all(is.na(short$object_coefficients$p_adj)))
})
