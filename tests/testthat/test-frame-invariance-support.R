test_that("frame support agrees with conditional pairs for all response patterns", {
  m <- c(1L, 2L, 3L)
  X <- as.matrix(expand.grid(lapply(m, function(k) c(NA_integer_, 0:k))))
  pairs <- combn(seq_along(m), 2L)
  expected <- rep(FALSE, nrow(X))
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1L, k]; j <- pairs[2L, k]
    ok <- !is.na(X[, i]) & !is.na(X[, j]) &
      X[, i] + X[, j] > 0 & X[, i] + X[, j] < m[i] + m[j]
    ok[is.na(ok)] <- FALSE
    expected <- expected | ok
  }
  expect_identical(.frame_informative_rows(X, m), expected)
  expect_identical(.frame_informative_rows(X[, 3L, drop = FALSE], 3L),
                   rep(FALSE, nrow(X)))
  expect_identical(.frame_informative_rows(X[, FALSE, drop = FALSE], integer()),
                   rep(FALSE, nrow(X)))
  expect_identical(.frame_informative_rows(cbind(X, 0L), c(m, 0L)), expected)
})

test_that("frame support stays within the set and counts each person once", {
  # Each row is informative over both sets, but not necessarily within either.
  X <- rbind(c(0, 0, 1, 1), c(1, 1, 0, 0), c(0, 1, 1, 1),
             c(1, 0, 0, 0), c(NA, 1, 0, 1))
  id <- c("A", "B", "C", "C", "D")
  expect_true(all(.frame_informative_rows(X, rep(1L, 4L))))
  first <- .frame_informative_rows(X[, 1:2], c(1L, 1L))
  second <- .frame_informative_rows(X[, 3:4], c(1L, 1L))
  expect_identical(first, c(FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_identical(second, c(FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(.frame_person_support(first, id), 1L)
  expect_identical(.frame_person_support(second, id), 1L)
})

test_that("extreme padding cannot unlock frame-invariance inference", {
  for (categories in c(2L, 4L)) {
    d <- simulate_efrm(n_per_group = 30, items_per_set = 6, n_sets = 1,
                       n_groups = 2, n_categories = categories, seed = 909)
    sets <- attr(d, "truth")$item_sets
    items <- unlist(sets, use.names = FALSE)
    fit <- rasch_efrm(d, groups = "group", id = "id", item_sets = sets,
                      boot_reps = 0)
    extra <- d[rep(c(1L, 31L), each = 100L), ]
    extra$id <- paste0("extra", seq_len(nrow(extra)))
    extra[items] <- rep(rep(c(0L, categories - 1L), each = 50L), 2L)
    padded <- rasch_efrm(rbind(d, extra), groups = "group", id = "id",
                         item_sets = sets, boot_reps = 0)
    expect_equal(padded$thresholds$tau, fit$thresholds$tau, tolerance = 1e-10)
    for (method in c("conditional", "bootstrap")) {
      expect_error(frame_invariance(fit, se_method = method, boot_reps = 2),
                    "at least 50 persons contributing informative item pairs")
      expect_error(frame_invariance(padded, se_method = method, boot_reps = 2),
                    "at least 50 persons contributing informative item pairs")
    }
  }
})

test_that("supported frame-invariance comparisons remain available", {
  d <- simulate_efrm(n_per_group = 120, items_per_set = 6, n_sets = 1,
                     n_groups = 2, seed = 909)
  fit <- rasch_efrm(d, groups = "group", id = "id",
                    item_sets = attr(d, "truth")$item_sets, boot_reps = 0)
  z <- frame_invariance(fit)
  expect_s3_class(z, "rasch_frame_invariance")
  expect_true(all(is.finite(z$locations$p)))
})

test_that("frame support uses the separate calibration's category endpoints", {
  for (endpoint in c("upper", "lower")) {
    d <- simulate_efrm(n_per_group = 150, items_per_set = 4, n_sets = 1,
                       n_groups = 2, n_categories = 3, seed = 901)
    sets <- attr(d, "truth")$item_sets
    items <- unlist(sets, use.names = FALSE)
    d <- d[c(1:40, 151:300), ]
    g1 <- d$group == "g1"
    d[g1, items[1L]] <- if (endpoint == "upper")
      pmin(d[g1, items[1L]], 1L) else pmax(d[g1, items[1L]], 1L)
    extra <- d[rep(1L, 100L), ]
    extra$id <- paste0("extra", seq_len(nrow(extra)))
    extra[items] <- if (endpoint == "upper") 2L else 0L
    extra[items[1L]] <- 1L

    before <- rasch(d[g1, items], model = "PCM")
    after <- rasch(rbind(d[g1, items], extra[items]), model = "PCM")
    expect_equal(after$thresholds$tau, before$thresholds$tau)
    expect_equal(after$est$cov_tau, before$est$cov_tau)
    support <- function(f) .frame_person_support(
      .frame_informative_rows(f$X, f$m), f$person$id)
    expect_identical(support(after), support(before))
    expect_lt(support(after), 50L)

    fit <- rasch_efrm(d, groups = "group", id = "id", item_sets = sets,
                      boot_reps = 0)
    padded <- rasch_efrm(rbind(d, extra), groups = "group", id = "id",
                         item_sets = sets, boot_reps = 0)
    vm <- padded$virtual_map
    jj <- match(vm$vkey[vm$group == "g1"], colnames(padded$X))
    expect_gte(sum(.frame_informative_rows(padded$X[, jj], padded$m[jj])), 50L)
    for (method in c("conditional", "bootstrap")) {
      expect_error(frame_invariance(fit, se_method = method, boot_reps = 2),
                    "at least 50 persons contributing informative item pairs")
      expect_error(frame_invariance(padded, se_method = method, boot_reps = 2),
                    "at least 50 persons contributing informative item pairs")
    }
    # A bootstrap replicate does not repeat the observed minimum-person gate.
    # Resampling below that boundary must not itself select the null sample.
    if (endpoint == "upper")
      expect_type(.frame_invariance_conditional(padded, strict = FALSE), "list")
  }
})
