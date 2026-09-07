test_that("dropping an item cannot rescore retained Rasch or explanatory items", {
  d <- simulate_rasch(300, 4, model = "PCM", n_categories = 3,
                      n_groups = 2, seed = 9791)
  d[d$I01 == 2, c("I02", "I03")] <- 2
  fit <- rasch(d, id = "id", factors = "group")
  expect_identical(unname(fit$m["I01"]), 2L)
  keep <- c("I01", "I02", "I03")
  # All categories are present before fitting. Removing I04 nevertheless
  # makes category 2 of I01 conditionally uninformative.
  expect_no_error(.require_score_structure(fit$X[, keep], fit$m[keep]))
  raw <- rasch(fit$X[, keep])
  expect_lt(raw$m["I01"], fit$m["I01"])
  expect_error(drop_items(fit, "I04"), "cannot preserve the fitted score structure")

  predictors <- data.frame(item = sprintf("I%02d", 1:4), feature = 1:4)
  ex <- rasch_explanatory(d, predictors, ~ feature, id = "id", factors = "group")
  expect_error(drop_items(ex, "I04"), "cannot preserve the fitted score structure")
})

test_that("EFRM structural refits preserve frame-specific observed scores", {
  d <- simulate_rasch(300, 4, model = "PCM", n_categories = 3,
                      n_groups = 2, seed = 9791)
  d[d$I01 == 2, c("I02", "I03")] <- 2
  fit <- rasch_efrm(d, item_sets = list(s = sprintf("I%02d", 1:4)),
                    groups = "group", id = "id", boot_reps = 0)
  source <- .efrm_source_matrix(fit, c("I01", "I02", "I03"))
  direct <- rasch_efrm(data.frame(source, group = d$group),
    item_sets = list(s = colnames(source)), groups = "group", boot_reps = 0)
  expect_lt(.efrm_source_maxima(direct)["I01"],
            .efrm_source_maxima(fit)["I01"])
  expect_false(isTRUE(all.equal(.efrm_source_matrix(direct), source,
                                check.attributes = FALSE)))
  expect_error(drop_items(fit, "I04", boot_reps = 0),
               "EFRM structural refit cannot preserve the fitted score structure")

  # A single-group boundary that retains its original scores must remain
  # usable; the guard tests actual rescoring rather than the input pattern.
  one <- simulate_rasch(300, 4, model = "PCM", n_categories = 3,
                        n_groups = 2, seed = 9791)
  one[one$group == "g1" & one$I01 == 2, c("I02", "I03")] <- 2
  supported <- rasch_efrm(one, item_sets = list(s = sprintf("I%02d", 1:4)),
                          groups = "group", id = "id", boot_reps = 0)
  expect_no_error(drop_items(supported, "I04", boot_reps = 0))
})

test_that("dependence resolution cannot compare recoded threshold copies", {
  set.seed(9792)
  X <- matrix(sample(0:2, 5000, TRUE), 1000, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  X[X[, 1] == 0 & X[, 2] == 2, 3:5] <- 2
  fit <- rasch(X)
  expect_true(all(fit$m == 2L))
  expect_true(all(vapply(0:2, function(k)
    setequal(X[X[, 1] == k, 2], 0:2), logical(1))))
  expect_error(dependence_magnitude(fit, "I2", "I1"),
    "cannot preserve the fitted score structure for I2|I1=0", fixed = TRUE)
})

test_that("superitem refits retain every intended summed score", {
  set.seed(9793)
  X <- matrix(sample(0:2, 5000, TRUE), 1000, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  X[X[, 1] == 2 & X[, 2] == 2, 3:5] <- 2
  fit <- rasch(X)
  expect_true(all(fit$m == 2L))
  expect_setequal(rowSums(X[, 1:2]), 0:4)
  expect_error(combine_items(fit, list(c("I1", "I2"))),
    "cannot preserve the fitted score structure for I1+I2", fixed = TRUE)
})
