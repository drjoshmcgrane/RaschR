test_that("fixed-origin equating uses marginal variances without a bank covariance", {
  set.seed(4601)
  N <- 450; L <- 6
  theta <- rnorm(N)
  delta <- seq(-1.1, 1.1, length.out = L)
  X <- matrix(rbinom(N * L, 1, plogis(outer(theta, delta, "-"))), N, L)
  colnames(X) <- paste0("I", seq_len(L))
  fit <- rasch(X)
  bank <- fit$items[, c("item", "location", "se", "max")]
  bank$location <- bank$location + seq(0.05, 0.30, length.out = L)

  # No joint covariance is attached. It is required when a shift is
  # estimated, but not when an anchored origin is fixed before comparison.
  mean_link <- equate_tests(fit, bank, shift = "mean")
  fixed_link <- equate_tests(fit, bank, shift = "none")
  expect_false(mean_link$inferential)
  expect_true(fixed_link$inferential)
  expected_se <- sqrt(fit$items$se^2 + bank$se^2)
  expect_equal(fixed_link$table$t,
               (fit$items$location - bank$location) / expected_se,
               tolerance = 1e-12)
  expect_equal(fixed_link$table$p_adj,
               p.adjust(fixed_link$table$p, "holm", n = L))

  # Unavailable common items remain in the declared family.
  sparse_bank <- bank
  sparse_bank$se[-1L] <- NA_real_
  sparse_link <- equate_tests(fit, sparse_bank, shift = "none")
  expect_true(sparse_link$inferential)
  expect_identical(sum(is.finite(sparse_link$table$p)), 1L)
  expect_equal(sparse_link$table$p_adj[1],
               p.adjust(sparse_link$table$p[1], "holm", n = L))

  one <- bank[1L, , drop = FALSE]
  one_link <- equate_tests(fit, one, shift = "none")
  expect_true(one_link$inferential)
  expect_identical(one_link$n_common, 1L)
  expect_true(is.finite(one_link$table$p_adj))
  expect_error(equate_tests(fit, one, shift = "mean"),
               "at least two common items")
})

test_that("estimated-shift equating does not invent joint fit covariance", {
  set.seed(4607)
  N <- 350; L <- 5
  theta <- rnorm(N)
  delta <- seq(-1, 1, length.out = L)
  X <- matrix(rbinom(N * L, 1,
    plogis(outer(theta, delta, "-"))), N, L)
  colnames(X) <- paste0("I", seq_len(L))
  fit1 <- rasch(X)
  fit2 <- rasch(X)
  fit2$est$cov_tau <- NULL

  linked <- equate_tests(fit1, fit2, independent = TRUE)
  expect_false(linked$inferential)
  expect_true(all(is.na(linked$table$p)))
  expect_match(linked$note, "covariance is unavailable")

  # Marginal variances remain sufficient when the origin was fixed outside
  # the comparison.
  fixed <- equate_tests(fit1, fit2, shift = "none", independent = TRUE)
  expect_true(fixed$inferential)
  expect_true(any(is.finite(fixed$table$p)))
})

test_that("spread inference forwards repeated person IDs to the sandwich", {
  set.seed(4602)
  N <- 300; L <- 6
  theta <- rnorm(N)
  X <- matrix(rbinom(N * L, 1,
    plogis(outer(theta, seq(-1, 1, length.out = L), "-"))), N, L)
  colnames(X) <- paste0("I", seq_len(L))
  fit <- combine_items(rasch(X),
                       list(c("I1", "I2"), c("I3", "I4"), c("I5", "I6")))
  fit$person$id <- rep(sprintf("P%03d", seq_len(N / 2L)), each = 2L)
  pc <- rasch:::.pcml_pc_fit(fit$X)
  seen <- NULL
  out <- testthat::with_mocked_bindings(
    spread_test(fit),
    .pcml_pc_fit = function(X, n_components = 4, maxit = 60, tol = 1e-8,
                            cluster = NULL) {
      seen <<- cluster
      pc
    },
    .package = "rasch")

  expect_s3_class(out, "rasch_spread")
  expect_identical(seen, rasch:::.dif_ids(fit$person$id))
  expect_true(rasch:::.has_repeated_person_ids(seen))
})

test_that("frame bootstrap resamples people rather than response rows", {
  id <- rep(c("A1", "A2", "B1"), each = 2L)
  group <- rep(c("A", "A", "B"), each = 2L)
  set.seed(4603)
  nested <- rasch:::.frame_cluster_resample(id, group)
  expect_true(nested$stratified)
  expect_length(nested$rows, length(id))
  expect_true(all(vapply(split(nested$rows, nested$id), function(ii)
    length(unique(id[ii])) == 1L, logical(1))))
  expect_identical(as.integer(table(group[nested$rows])), c(4L, 2L))

  crossed_id <- rep(c("P1", "P2", "P3"), each = 2L)
  crossed_group <- rep(c("A", "B"), 3L)
  set.seed(4604)
  crossed <- rasch:::.frame_cluster_resample(crossed_id, crossed_group)
  expect_false(crossed$stratified)
  expect_true(rasch:::.frame_ids_cross_groups(crossed_id, crossed_group))
  expect_true(all(vapply(split(crossed$rows, crossed$id), function(ii)
    length(unique(crossed_id[ii])) == 1L, logical(1))))
  expect_true(all(lengths(split(crossed$rows, crossed$id)) == 2L))

  expect_identical(
    rasch:::.frame_person_support(rep(TRUE, length(id)), id), 3L)
})

test_that("conditional frame inference refuses cross-frame repeated persons", {
  fit <- list(
    est = list(converged = TRUE),
    frame_group = "group",
    factors = data.frame(group = c("A", "B")),
    person = data.frame(id = c("P1", "P1")))
  class(fit) <- c("rasch_efrm", "rasch")
  expect_error(
    frame_invariance(fit, se_method = "conditional"),
    "person appears in more than one frame.*cross-frame covariance")
})

test_that("DIF minimum support counts distinct persons rather than rows", {
  group <- factor(rep(c("A", "B"), each = 6L))
  id <- rep(sprintf("P%02d", seq_len(6L)), each = 2L)
  observed <- rep(TRUE, length(group))
  expect_identical(
    rasch:::.dif_cell_n(group, observed, id),
    c(A = 3L, B = 3L))
  expect_identical(
    rasch:::.dif_cell_n(group, group == "A", id),
    c(A = 3L, B = 0L))

  set.seed(4605)
  N <- 40L; L <- 6L
  theta <- rnorm(N)
  X <- matrix(rbinom(N * L, 1,
    plogis(outer(theta, seq(-1, 1, length.out = L), "-"))), N, L)
  colnames(X) <- paste0("I", seq_len(L))
  dat <- data.frame(id = rep(sprintf("P%03d", seq_len(N)), each = 2L),
                    grp = rep(rep(c("A", "B"), each = N / 2L), each = 2L),
                    X[rep(seq_len(N), each = 2L), , drop = FALSE],
                    check.names = FALSE)
  fit <- rasch(dat, id = "id", factors = "grp", items = colnames(X))
  compared <- compare_fits(first = fit, second = fit)
  expect_equal(compared$persons, rep(N, 2L))
  expect_error(
    dif_size(fit, "I1", by = "grp", min_n = 21L),
    "fewer than two usable levels")
})

test_that("tailored bootstrap accepts only a complete finite shift vector", {
  items <- c("I1", "I2")
  good <- list(table = data.frame(item = rev(items), shift = c(0.2, -0.1)))
  expect_identical(rasch:::.tailored_boot_shift(good, items), c(-0.1, 0.2))

  missing <- list(table = data.frame(item = "I1", shift = -0.1))
  nonfinite <- list(table = data.frame(item = items, shift = c(-0.1, NA)))
  duplicated <- list(table = data.frame(item = c("I1", "I1"),
                                        shift = c(-0.1, 0.2)))
  expect_null(rasch:::.tailored_boot_shift(missing, items))
  expect_null(rasch:::.tailored_boot_shift(nonfinite, items))
  expect_null(rasch:::.tailored_boot_shift(duplicated, items))
})

test_that("BTL fixed-origin equating does not re-estimate an anchored shift", {
  objects <- paste0("O", 1:3)
  fit <- structure(list(
    objects = data.frame(object = objects,
                         location = c(-0.8, 0, 0.8),
                         se = rep(0.2, 3)),
    converged = TRUE, m = 1L, categories = 0:1,
    thr_structure = "none", clustered = FALSE),
    class = "rasch_btl")
  bank <- data.frame(object = objects,
                     location = c(-1.0, -0.4, 0.2),
                     se = rep(0.3, 3))

  estimated <- btl_equate(fit, bank)
  fixed <- btl_equate(fit, bank, shift = "none")
  expect_false(estimated$inferential)
  expect_true(fixed$inferential)
  expect_identical(fixed$shift_method, "none")
  expect_identical(fixed$shift_setting, "none")
  expect_equal(fixed$shift, 0)
  expected_se <- rep(sqrt(0.2^2 + 0.3^2), 3)
  expect_equal(fixed$table$se_diff, expected_se)
  expect_equal(fixed$table$t,
               (fit$objects$location - bank$location) / expected_se)
  expect_equal(fixed$table$p_adj,
               p.adjust(fixed$table$p, method = "holm", n = 3L))

  one <- bank[1L, , drop = FALSE]
  one_link <- btl_equate(fit, one, shift = "none")
  expect_true(one_link$inferential)
  expect_identical(one_link$n_common, 1L)
  expect_true(is.finite(one_link$table$p_adj))
  expect_error(btl_equate(fit, one), "at least two common objects")
})
