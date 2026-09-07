test_that("tailored fits retain keys, censor option data and preserve split records", {
  set.seed(9350)
  n <- 600
  theta <- rnorm(n)
  raw <- sapply(seq(-2, 2, length.out = 8), function(d) {
    x <- rbinom(n, 1, .15 + .85 * plogis(theta - d))
    ifelse(x == 1, "A", sample(c("B", "C", "D"), n, replace = TRUE))
  })
  colnames(raw) <- paste0("I", 1:8)
  fit <- rasch(raw, key = setNames(rep("A", 8), colnames(raw)),
    factors = data.frame(group = rep(c("a", "b"), length.out = n)))
  split <- split_items(fit, "I1", by = "group")
  for (source in list(fit, split)) {
    result <- tailored_analysis(source)
    for (nm in c("tailored", "origin_equated", "anchored")) {
      out <- result[[nm]]
      expect_equal(out$mc$map, source$mc$map)
      # Returned tables receive a printing class; the scoring recipe is exact.
      expect_identical(as.data.frame(out$refit_spec$key),
                       as.data.frame(source$refit_spec$key))
      expect_identical(out$split_map, source$split_map)
      expect_equal(.n_unsplit_sources(.split_source_map(out)),
                    .n_unsplit_sources(.split_source_map(source)))
      expected <- source$mc$raw
      expected[is.na(out$X[, colnames(expected), drop = FALSE])] <- NA_character_
      expect_identical(out$mc$raw, expected)
      expect_no_error(distractor_analysis(out))
    }
    keyed <- colnames(result$tailored$mc$raw)
    removed <- !is.na(source$mc$raw) &
      is.na(result$tailored$X[, keyed, drop = FALSE])
    expect_true(any(removed))
    expect_true(all(is.na(result$tailored$mc$raw[removed])))
    expect_identical(result$origin_equated$mc$raw, source$mc$raw)
    expect_identical(result$anchored$mc$raw, source$mc$raw)
    expect_true(result$anchored$refit_spec$fixed_calibration)
  }
  # A later keyed refit must not restore censored observations.
  keyed_result <- tailored_analysis(fit)
  numeric_result <- tailored_analysis(rasch(fit$X, id = fit$person$id,
                                             factors = fit$factors))
  expect_equal(keyed_result$table, numeric_result$table, tolerance = 1e-10)
  expect_equal(keyed_result$anchored$person$theta,
                numeric_result$anchored$person$theta, tolerance = 1e-10)
  ta <- keyed_result$tailored
  dropped <- drop_items(ta, "I8")
  expect_identical(dropped$X, ta$X[, colnames(dropped$X), drop = FALSE])
  expect_identical(is.na(dropped$mc$raw), is.na(dropped$X))
  expect_no_error(distractor_analysis(dropped))
})
