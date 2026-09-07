test_that("legacy restrictions cannot disappear during downstream refits", {
  set.seed(773)
  X <- matrix(rbinom(2100, 1, .5), 350, 6,
    dimnames = list(NULL, paste0("I", 1:6)))
  fixed <- rasch(X, anchors = data.frame(item = "I1", k = 1, tau = 2))
  average <- rasch(X, anchors = data.frame(
    item = c("I1", "I2"), k = NA_integer_, tau = c(1, 2), average = TRUE))
  d <- simulate_rasch(300, 6, model = "PCM", n_categories = 4, seed = 773)
  pc <- rasch(d, pc_components = 2)
  expect_true(fixed$est$converged)
  expect_true(average$est$converged)
  expect_true(pc$est$converged)
  expect_false(any(average$est$thr$anchored))
  expect_gt(NROW(average$est$anchors), 0L)

  # A refit with intact settings remains supported and retains restrictions.
  reduced <- drop_items(fixed, "I6")
  expect_equal(reduced$items$location[reduced$items$item == "I1"], 2)
  reduced_pc <- drop_items(pc, "I06")
  expect_equal(reduced_pc$est$n_components, 2)
  reduced_average <- drop_items(average, "I6")
  expect_equal(mean(reduced_average$items$location[
    reduced_average$items$item %in% c("I1", "I2")]), 1.5)

  for (metadata in list(NULL, list())) {
    for (f in list(fixed, average, pc)) {
      f$refit_spec <- metadata
      f <- unserialize(serialize(f, NULL))
      message <- if (.has_pc_thresholds(f)) "saved component settings" else
        "saved anchor settings"
      nm <- colnames(f$X)
      expect_error(drop_items(f, nm[6]), message, class = "rasch_refusal")
      expect_error(combine_items(f, nm[1:2]), message)
      expect_error(split_items(f, nm[1], rep(c("A", "B"), length.out = nrow(f$X))),
        message)
      expect_error(.rasch_refit(f, f$X), message)
      expect_error(tailored_analysis(f), message)
      expect_error(fit_bootstrap(f, B = 20, workers = 1), message)
      expect_error(dif_bootstrap(f, B = 20, workers = 1), message)
      expect_error(.scree_reference(f, 2, 20), message)
      expect_error(.dim_bootstrap_check(f), message)
      expect_error(.validate_fit_bootstrap(list(), f), message)
      expect_error(.validate_dif_bootstrap(list(), f), message)
      # Scoring and descriptive diagnostics do not recalibrate the model.
      expect_no_error(score_table(f))
      expect_no_error(residual_pca(f))
    }
  }
  # Older fixed-threshold objects may retain only the threshold flags.
  fixed$refit_spec <- NULL
  fixed$est$anchors <- NULL
  expect_error(drop_items(fixed, "I6"), "saved anchor settings")

  # Absence of metadata alone is not evidence of a restriction.
  ordinary <- rasch(X)
  expected <- drop_items(ordinary, "I6")
  ordinary$refit_spec <- NULL
  actual <- drop_items(ordinary, "I6")
  expect_equal(actual$items$location, expected$items$location)
  expect_equal(actual$person$theta, expected$person$theta)
})

test_that("legacy anchored recovery retains its identified origin", {
  d <- simulate_rasch(300, 6, difficulty = c(-1, 1), seed = 775)
  truth <- attr(d, "truth")
  f <- rasch(d, anchors = data.frame(item = "I01", k = 1L,
    tau = unname(truth$difficulty["I01"]) + 1))
  expected <- sim_recovery(f, d)
  expect_true(all(is.finite(expected$summary$bias)))
  expect_gt(abs(expected$summary$bias[1L]), .5)
  for (metadata in list(NULL, list())) {
    f$refit_spec <- metadata
    actual <- sim_recovery(f, d)
    expect_equal(actual$summary, expected$summary)
    expect_equal(actual$pieces, expected$pieces)
  }
})
