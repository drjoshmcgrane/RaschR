test_that("combining explanatory items cannot silently ignore the model", {
  set.seed(9560)
  n <- 500L
  predictors <- data.frame(item = paste0("I", 1:8),
                            feature = rep(0:1, each = 4))
  theta <- rnorm(n)
  X <- vapply(.6 * predictors$feature, function(d)
    rbinom(n, 1, plogis(theta - d)), integer(n))
  colnames(X) <- predictors$item
  explanatory <- rasch_explanatory(X, predictors, ~ feature)
  ordinary <- explanatory$reference_fit
  groups <- lapply(seq(1, 7, 2), function(i) paste0("I", i:(i + 1)))

  for (fit in list(explanatory, ordinary)) {
    for (bad in list(NULL, NA_character_, c("PCM", "RSM"),
                     matrix("PCM"), factor("PCM"), 1))
      expect_error(combine_items(fit, groups, model = bad),
                   "`model` must name one model", fixed = TRUE)
    expect_error(combine_items(fit, groups, model = "invalid"),
                 "arg.*should be one of")
  }
  expect_error(combine_items(explanatory, groups, model = "RSM"),
               'model = "RSM" is not supported for explanatory fits',
               fixed = TRUE)

  # Four two-item sums make RSM a legitimate request for an ordinary fit.
  # Explanatory superitems are instead fully relaxed by the default workflow.
  combined <- combine_items(explanatory, groups)
  rsm <- combine_items(ordinary, groups, model = "RSM")
  direct_rsm <- rasch(combined$X, model = "RSM")
  expect_s3_class(combined, "rasch_explanatory")
  expect_equal(unname(combined$m), rep(2L, 4))
  expect_equal(combined$est$n_parameters, 7L)
  expect_identical(rsm$model, "RSM")
  expect_equal(rsm$est$n_parameters, 4L)
  expect_equal(rsm$est$loglik, direct_rsm$est$loglik, tolerance = 1e-10)
  expect_equal(rsm$person$theta, direct_rsm$person$theta, tolerance = 1e-10)
  expect_equal(rsm$X, combined$X)

  # The ordinary default still fits PCM; it must not inherit an input RSM's
  # common threshold restriction when unequal superitem maxima are created.
  base_rsm <- rasch(X, model = "RSM")
  mixed <- combine_items(base_rsm, c("I1", "I2"))
  expect_identical(mixed$model, "PCM")
  expect_error(combine_items(ordinary, c("I1", "I2"), model = "RSM"),
               "RSM requires equal max score")
})
