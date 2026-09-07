test_that("explanatory relaxation retains split and superitem provenance", {
  set.seed(4105)
  n <- 500
  predictors <- data.frame(item = paste0("I", 1:10),
                           feature = rep(0:1, each = 5))
  theta <- rnorm(n)
  X <- sapply(.65 * predictors$feature, function(d)
    rbinom(n, 1, plogis(theta - d)))
  colnames(X) <- predictors$item
  fit <- rasch_explanatory(X, predictors, ~ feature,
    factors = data.frame(group = rep(c("A", "B"), length.out = n)))

  split <- split_items(fit, "I1", by = "group")
  relaxed <- relax_explanatory(split, "I2")
  expect_identical(relaxed$X, split$X)
  expect_identical(relaxed$split_map, split$split_map)
  expect_length(unique(.split_source_map(relaxed)), 10L)
  expect_equal(.n_unsplit_sources(.split_source_map(relaxed)), 9L)

  combined <- combine_items(fit, c("I1", "I2"))
  relaxed_combined <- relax_explanatory(combined, "I3")
  expect_identical(relaxed_combined$X, combined$X)
  expect_identical(relaxed_combined$subtest_map, combined$subtest_map)
  expect_identical(relaxed_combined$subtest_binary, combined$subtest_binary)
  expect_length(relaxed_combined$subtest_map, 1L)
  expect_equal(spread_test(relaxed_combined), spread_test(combined))
})
