test_that("the rating-scale comparison retains keyed scoring and split history", {
  set.seed(9240)
  n <- 500
  theta <- rnorm(n)
  X <- sapply(seq(-1, 1, length.out = 6), function(b)
    vapply(theta, function(t)
      sample(0:2, 1, prob = item_moments(t, b + c(-.6, .6))$P), 0L))
  colnames(X) <- paste0("I", 1:6)
  raw <- matrix(c("C", "B", "A")[X + 1L], nrow(X), dimnames = dimnames(X))
  key <- do.call(rbind, lapply(colnames(X), function(it)
    data.frame(item = it, option = c("A", "B", "C"), score = c(2L, 1L, 0L))))
  fit <- rasch(raw, key = key,
    factors = data.frame(group = rep(c("a", "b"), length.out = n)))
  result <- lr_test(fit)
  rsm <- result$fit_rsm
  direct <- rasch(raw, key = key, model = "RSM", factors = fit$factors)
  expect_identical(rsm$X, fit$X)
  expect_equal(rsm$mc, fit$mc)
  # lr_test tags nested data frames for printing; scoring values are unchanged.
  expect_identical(as.data.frame(rsm$refit_spec$key),
                   as.data.frame(fit$refit_spec$key))
  expect_equal(rsm$est$loglik, direct$est$loglik, tolerance = 1e-10)
  expect_equal(distractor_analysis(rsm), distractor_analysis(direct))
  # The scoring recipe must survive a further refit, not just the first table.
  dropped <- drop_items(rsm, "I6")
  expect_equal(colnames(dropped$mc$raw), colnames(dropped$X))
  expect_no_error(distractor_analysis(dropped))

  split <- split_items(fit, "I1", by = "group")
  split_rsm <- lr_test(split)$fit_rsm
  expect_identical(split_rsm$X, split$X)
  expect_identical(split_rsm$split_map, split$split_map)
  expect_equal(.n_unsplit_sources(.split_source_map(split_rsm)), 5L)
  expect_equal(split_rsm$mc, split$mc)
})

test_that("the rating-scale comparison retains superitem definitions", {
  fit <- rasch(simulate_rasch(500, 6, seed = 9241))
  groups <- split(colnames(fit$X), rep(1:3, each = 2))
  combined <- combine_items(fit, groups)
  rsm <- lr_test(combined)$fit_rsm
  expect_identical(rsm$X, combined$X)
  expect_identical(rsm$subtest_map, combined$subtest_map)
  expect_identical(rsm$subtest_binary, combined$subtest_binary)
  expect_equal(spread_test(rsm), spread_test(combined))
})
