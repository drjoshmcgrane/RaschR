test_that("Rasch explanatory fits retain distinct coincident coefficient labels", {
  d <- simulate_rasch(180, 8, seed = 671)
  q <- data.frame(item = sprintf("I%02d", 1:8),
    x = rep(c("A", "B"), 4), xB = rep(0:3, each = 2))
  expect_identical(colnames(model.matrix(~ x + xB, q)),
    c("(Intercept)", "xB", "xB"))
  renamed <- q
  names(renamed)[3L] <- "amount"
  fit <- rasch_explanatory(d, q, ~ x + xB)
  reference <- rasch_explanatory(d, renamed, ~ x + amount)
  expect_identical(fit$est$coefficients$term, c("xB", "xB.1"))
  expect_identical(colnames(fit$est$B), fit$est$coefficients$term)
  expect_identical(colnames(fit$explanatory$model_matrix),
    fit$est$coefficients$term)
  expect_identical(rownames(fit$est$cov_beta), fit$est$coefficients$term)
  expect_identical(fit$explanatory$formula, ~ x + xB)
  expect_identical(fit$explanatory$source_predictors, q)
  expect_equal(fit$thresholds$tau, reference$thresholds$tau)
  expect_equal(fit$person$theta, reference$person$theta)
  expect_equal(fit$est$coefficients$p_adj, reference$est$coefficients$p_adj)
  expect_equal(explanatory_test(fit), explanatory_test(reference))
  expect_equal(explanatory_diagnostics(fit), explanatory_diagnostics(reference))
  relaxed <- relax_explanatory(fit, "I01")
  relaxed_reference <- relax_explanatory(reference, "I01")
  expect_equal(relaxed$thresholds$tau, relaxed_reference$thresholds$tau)
  expect_equal(explanatory_test(relaxed), explanatory_test(relaxed_reference))
  expect_false(anyDuplicated(relaxed$est$coefficients$term) > 0L)
})

test_that("CJ explanatory fits retain distinct coincident coefficient labels", {
  d <- simulate_btl(8, 12, 15, seed = 672)
  q <- data.frame(object = paste0("O", 1:8),
    x = rep(c("A", "B"), 4), xB = rep(0:3, each = 2))
  renamed <- q
  names(renamed)[3L] <- "amount"
  fit <- btl_explanatory(d, q, ~ x + xB, "object_a", "object_b",
    winner = "winner", judge = "judge")
  reference <- btl_explanatory(d, renamed, ~ x + amount,
    "object_a", "object_b", winner = "winner", judge = "judge")
  expect_identical(fit$object_coefficients$term, c("xB", "xB.1"))
  expect_identical(colnames(fit$location_design), fit$object_coefficients$term)
  expect_identical(colnames(fit$explanatory$active_B),
    fit$object_coefficients$term)
  expect_equal(fit$objects$location, reference$objects$location)
  expect_equal(fit$object_coefficients$p_adj, reference$object_coefficients$p_adj)
  expect_equal(explanatory_test(fit), explanatory_test(reference))
  expect_equal(explanatory_diagnostics(fit), explanatory_diagnostics(reference))
  relaxed <- relax_btl_explanatory(fit, "O1")
  relaxed_reference <- relax_btl_explanatory(reference, "O1")
  expect_equal(relaxed$objects$location, relaxed_reference$objects$location)
  expect_equal(explanatory_test(relaxed), explanatory_test(relaxed_reference))
  expect_false(anyDuplicated(relaxed$object_coefficients$term) > 0L)
})

test_that("new departures cannot reuse an existing predictor coefficient label", {
  d <- simulate_rasch(180, 8, seed = 673)
  q <- data.frame(item = sprintf("I%02d", 1:8),
    departure_location = factor(rep(c("other", "[I01]"), 4),
      levels = c("other", "[I01]")))
  fit <- rasch_explanatory(d, q, ~ departure_location)
  expect_identical(fit$est$coefficients$term, "departure_location[I01]")
  relaxed <- relax_explanatory(fit, "I01")
  expect_identical(relaxed$est$coefficients$term,
    c("departure_location[I01]", "departure_location[I01].1"))
  expect_identical(colnames(relaxed$est$B), relaxed$est$coefficients$term)
  expect_no_error(explanatory_test(relaxed))
  expect_no_error(explanatory_diagnostics(fit))

  dc <- simulate_btl(8, 12, 15, seed = 674)
  qc <- data.frame(object = paste0("O", 1:8),
    departure = factor(rep(c("other", "[O1]"), 4),
      levels = c("other", "[O1]")))
  fc <- btl_explanatory(dc, qc, ~ departure, "object_a", "object_b",
    winner = "winner", judge = "judge")
  expect_identical(fc$object_coefficients$term, "departure[O1]")
  rc <- relax_btl_explanatory(fc, "O1")
  expect_identical(rc$object_coefficients$term,
    c("departure[O1]", "departure[O1].1"))
  expect_identical(colnames(rc$explanatory$active_B), rc$object_coefficients$term)
  expect_no_error(explanatory_test(rc))
  expect_no_error(explanatory_diagnostics(fc))
})
