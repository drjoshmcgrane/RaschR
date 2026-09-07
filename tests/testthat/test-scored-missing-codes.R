test_that("raw missing codes do not erase scores assigned by a key", {
  set.seed(9561)
  n <- 400L
  raw <- matrix(sample(c("A", "B", "C"), n * 6, TRUE), n, 6,
                 dimnames = list(NULL, paste0("M", 1:6)))
  key <- do.call(rbind, lapply(colnames(raw), function(it)
    data.frame(item = it, option = c("A", "B", "C"), score = 2:0)))
  for (code in c(0L, 1L, 2L)) {
    x <- raw
    x[1, 1] <- sprintf("%02d", code)
    x[2, 2] <- "-2"
    fit <- rasch(x, key = key, na_codes = code)
    expected <- matrix(unname(c(A = 2L, B = 1L, C = 0L)[x]), n, 6,
                        dimnames = dimnames(x))
    direct <- rasch(expected, na_codes = integer(0))
    expect_identical(fit$X, direct$X)
    expect_equal(sum(is.na(fit$X)), 2L)
    expect_true(all(fit$m == 2L))
    expect_true(all(is.na(fit$mc$raw[cbind(1:2, 1:2)])))
    expect_equal(fit$est$loglik, direct$est$loglik, tolerance = 1e-10)
    expect_equal(fit$est$cov_tau, direct$est$cov_tau, tolerance = 1e-10)
    expect_equal(fit$person$theta, direct$person$theta, tolerance = 1e-10)
  }

  # With binary scoring, a raw missing code of zero must not remove every
  # incorrect answer and leave a constant item.
  x <- raw; x[1, 1] <- "0"
  binary <- rasch(x, key = setNames(rep("A", 6), colnames(x)), na_codes = 0)
  expected <- 1L * (x == "A"); expected[1, 1] <- NA_integer_
  expect_identical(binary$X, expected)
  expect_equal(sum(is.na(binary$X)), 1L)
})

test_that("mixed keyed and numeric refits use the prepared score scale", {
  set.seed(9562)
  n <- 500L
  raw <- matrix(sample(c("A", "B", "C"), n * 6, TRUE), n, 6,
                 dimnames = list(NULL, paste0("M", 1:6)))
  raw[1, 1] <- "01"
  key <- do.call(rbind, lapply(colnames(raw), function(it)
    data.frame(item = it, option = c("A", "B", "C"), score = 2:0)))
  data <- data.frame(raw, N = sample(0:3, n, TRUE), check.names = FALSE)
  expected <- data.frame(matrix(unname(c(A = 2L, B = 1L, C = 0L)[raw]),
                                 n, 6, dimnames = dimnames(raw)),
                          N = data$N, check.names = FALSE)
  expected$N[expected$N == 1L] <- NA_integer_
  factors <- data.frame(group = rep(c("A", "B"), length.out = n))
  fit <- rasch(data, key = key, na_codes = "01", factors = factors)
  direct <- rasch(expected, na_codes = integer(0), factors = factors)
  expect_identical(fit$X, direct$X)
  expect_identical(is.na(fit$X[, "N"]), data$N == 1L)
  expect_true(any(fit$X[, "N"] == 1L, na.rm = TRUE))
  expect_equal(fit$est$cov_tau, direct$est$cov_tau, tolerance = 1e-10)

  compare_refit <- function(actual, reference) {
    expect_identical(actual$X, reference$X)
    expect_equal(actual$est$loglik, reference$est$loglik, tolerance = 1e-9)
    expect_equal(actual$person$theta, reference$person$theta, tolerance = 1e-9)
  }
  dropped <- drop_items(fit, "M6")
  compare_refit(dropped, drop_items(direct, "M6"))
  expect_identical(dropped$mc$raw, fit$mc$raw[, -6, drop = FALSE])
  compare_refit(drop_items(dropped, "M5"), drop_items(drop_items(direct, "M6"), "M5"))
  compare_refit(split_items(fit, "M1", "group"),
                split_items(direct, "M1", "group"))
  compare_refit(combine_items(fit, c("M1", "M2")),
                combine_items(direct, c("M1", "M2")))

  anchors <- data.frame(item = "M1", k = 1:2, tau = direct$tau_list[[1]])
  anchored <- rasch(data, key = key, na_codes = "01", factors = factors,
                    anchors = anchors)
  anchored_direct <- rasch(expected, na_codes = integer(0), factors = factors,
                           anchors = anchors)
  compare_refit(anchored, anchored_direct)
  compare_refit(drop_items(anchored, "M6"), drop_items(anchored_direct, "M6"))

  predictors <- data.frame(item = names(data), feature = seq_along(data))
  explanatory <- rasch_explanatory(data, predictors, ~ feature + threshold,
    key = key, na_codes = "01", factors = factors)
  explanatory_direct <- rasch_explanatory(expected, predictors,
    ~ feature + threshold, na_codes = integer(0), factors = factors)
  compare_refit(explanatory, explanatory_direct)
  compare_refit(drop_items(explanatory, "M6"),
                drop_items(explanatory_direct, "M6"))

  # The issue also occurs without a key: removing raw category 1 renumbers
  # 0, 2, 3 to 0, 1, 2. Refitting must retain the newly assigned score 1.
  numeric_data <- matrix(sample(0:3, n * 6, TRUE), n, 6,
                          dimnames = list(NULL, paste0("I", 1:6)))
  numeric_fit <- rasch(numeric_data, na_codes = 1)
  numeric_drop <- drop_items(numeric_fit, "I6")
  expect_identical(numeric_drop$X, numeric_fit$X[, -6, drop = FALSE])
  expect_true(all(numeric_drop$m == 2L))
  compare_refit(numeric_drop,
    rasch(numeric_fit$X[, -6, drop = FALSE], na_codes = integer(0)))
})

test_that("EFRM item removal does not reapply raw missing codes", {
  d <- simulate_efrm(n_per_group = 250, items_per_set = 5, n_sets = 1,
                     n_groups = 1, n_categories = 3, seed = 9563)
  sets <- attr(d, "truth")$item_sets
  items <- unlist(sets, use.names = FALSE)
  # Raw 0, 2, 4 become prepared 0, 1, 2; the declared raw code 1 is absent.
  d[items] <- lapply(d[items], function(x) 2L * x)
  fit <- rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                    na_codes = 1, se_method = "hybrid", boot_reps = 0)
  source <- .efrm_source_matrix(fit)
  removed <- items[5]
  keep <- setdiff(items, removed)
  refit <- drop_items(fit, removed, boot_reps = 0)
  expect_identical(.efrm_source_matrix(refit), source[, keep, drop = FALSE])
  expect_true(all(refit$m == 2L))
  manual_data <- data.frame(id = fit$person$id, source[, keep, drop = FALSE],
                            group = fit$factors$group, check.names = FALSE)
  manual <- rasch_efrm(manual_data,
    item_sets = split(keep, fit$set_of[keep]), groups = "group", id = "id",
    na_codes = integer(0), se_method = "hybrid", boot_reps = 0)
  expect_equal(refit$est$loglik, manual$est$loglik, tolerance = 1e-9)
  expect_equal(refit$person$theta, manual$person$theta, tolerance = 1e-9)
})
