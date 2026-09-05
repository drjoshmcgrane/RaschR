test_that("equal external weights reproduce ordinary person estimates", {
  set.seed(912)
  d <- simulate_rasch(n_persons = 120, n_items = 6, seed = 912)
  f <- rasch(d)
  w <- stats::setNames(rep(3, ncol(f$X)), colnames(f$X))
  z <- weighted_person_estimates(f, w)
  w_pad <- w
  names(w_pad) <- paste0(" ", names(w_pad), " ")
  expect_equal(weighted_person_estimates(f, w_pad)$theta, z$theta)

  expect_equal(z$theta, f$person$theta, tolerance = 1e-7)
  expect_equal(z$se, f$person$se, tolerance = 1e-7)
  expect_equal(z$extreme, f$person$extreme)
  expect_equal(attr(z, "weighting")$normalised_weight,
               rep(1, ncol(f$X)))
})

test_that("external item weights use the weighted score sandwich", {
  set.seed(913)
  d <- simulate_rasch(n_persons = 100, n_items = 5, seed = 913)
  f <- rasch(d)
  w <- stats::setNames(c(0, 1, 2, 1, 0.5), colnames(f$X))
  z <- weighted_person_estimates(f, w)
  q <- unname(w[colnames(f$X)]); q <- q / mean(q)
  p <- which(rowSums(!is.na(f$X[, q > 0, drop = FALSE])) > 0)[1]
  cols <- which(q > 0 & !is.na(f$X[p, ]))
  mo <- lapply(cols, function(j)
    item_moments(z$theta[p], f$tau_list[[j]]))
  H <- sum(q[cols] * vapply(mo, `[[`, 0, "V"))
  J <- sum(q[cols]^2 * vapply(mo, `[[`, 0, "V"))

  expect_equal(z$se[p], sqrt(J) / H, tolerance = 1e-10)

  E <- vapply(mo, `[[`, 0, "E")
  m3 <- vapply(mo, `[[`, 0, "mu3")
  U <- sum(q[cols] * (f$X[p, cols] - E))
  corrected_score <- U + J * sum(q[cols] * m3) / (2 * H^2)
  old_frequency_score <- U + sum(q[cols] * m3) / (2 * H)
  expect_equal(corrected_score, 0, tolerance = 2e-7)
  expect_gt(abs(old_frequency_score), 1e-4)
  expect_equal(z$n_items, rowSums(!is.na(f$X[, q > 0, drop = FALSE])))
  expect_equal(weighted_person_estimates(f, 10 * w)$theta, z$theta,
               tolerance = 1e-9)

  # Relative scale must remain irrelevant at the floating-point boundary.
  tiny <- .Machine$double.xmin * .Machine$double.eps
  wt <- w; wt[wt > 0] <- tiny
  zt <- weighted_person_estimates(f, wt)
  wr <- stats::setNames(as.numeric(wt > 0), names(wt))
  zr <- weighted_person_estimates(f, wr)
  expect_equal(zt$theta, zr$theta, tolerance = 1e-9)
  expect_equal(zt$se, zr$se, tolerance = 1e-9)
})

test_that("set weights resolve lists and EFRM maps", {
  set.seed(914)
  d <- simulate_rasch(n_persons = 100, n_items = 6, seed = 914)
  f <- rasch(d)
  sets <- list(A = colnames(f$X)[1:3], B = colnames(f$X)[4:6])
  z <- weighted_person_estimates(f, c(A = 2, B = 1), by = "set",
                                 sets = sets)
  expect_equal(attr(z, "weighting")$set, rep(c("A", "B"), each = 3))
  z_pad <- weighted_person_estimates(
    f, c(" A " = 2, "B " = 1), by = "set",
    sets = list(" A " = paste0(" ", colnames(f$X)[1:3], " "),
                "B " = colnames(f$X)[4:6]))
  expect_equal(z_pad$theta, z$theta)
  set_vector <- stats::setNames(rep(c(" A ", "B "), each = 3),
                                paste0(" ", colnames(f$X), " "))
  expect_equal(weighted_person_estimates(
    f, c(A = 2, B = 1), by = "set", sets = set_vector)$theta, z$theta)
  expect_error(weighted_person_estimates(
    f, c(A = 2, " A " = 1, B = 1), by = "set", sets = sets),
    "non-empty named numeric")

  fe <- f
  class(fe) <- c("rasch_efrm", "rasch")
  fe$virtual_map <- data.frame(vkey = colnames(f$X), item = colnames(f$X))
  fe$set_of <- stats::setNames(rep(c("A", "B"), each = 3), colnames(f$X))
  fe$alpha_table <- data.frame(set = c("A", "B"))
  fe$linking <- list(alpha_edges = data.frame(converged = TRUE))
  expect_equal(weighted_person_estimates(fe, c(A = 2, B = 1), by = "set")$theta,
               z$theta)
})

test_that("item weighting preserves exact fitted names before trimming", {
  d <- simulate_rasch(n_persons = 100, n_items = 4, seed = 919)
  names(d)[2L] <- " I01 "
  f <- rasch(d)
  w <- stats::setNames(rep(1, ncol(f$X)), colnames(f$X))
  exact <- weighted_person_estimates(f, w)
  expect_equal(exact$theta, f$person$theta, tolerance = 1e-7)

  sets <- list(A = colnames(f$X)[1:2], B = colnames(f$X)[3:4])
  by_set <- weighted_person_estimates(f, c(A = 1, B = 1), by = "set",
                                      sets = sets)
  expect_equal(by_set$theta, f$person$theta, tolerance = 1e-7)

  # When no exact spaced name exists, padded selector syntax remains useful.
  names(d)[2L] <- "I01"
  f2 <- rasch(d)
  w2 <- stats::setNames(rep(1, ncol(f2$X)),
                        paste0(" ", colnames(f2$X), " "))
  expect_equal(weighted_person_estimates(f2, w2)$theta,
               f2$person$theta, tolerance = 1e-7)
})

test_that("equal weights reproduce MFRM and EFRM person estimates", {
  dm <- simulate_mfrm(n_persons = 50, n_items = 4, n_raters = 4, seed = 916)
  mf <- rasch_mfrm(dm, "person", "item", "score", facets = "rater")
  wm <- stats::setNames(rep(1, nrow(mf$item_effects)), mf$item_effects$item)
  zm <- weighted_person_estimates(mf, wm)
  expect_equal(zm$theta, mf$person$theta, tolerance = 1e-7)
  expect_equal(zm$se, mf$person$se, tolerance = 1e-7)

  de <- simulate_efrm(n_per_group = 100, items_per_set = 4, n_sets = 2,
                      n_groups = 2, set_unit_ratio = 1.25, seed = 917)
  tr <- attr(de, "truth")
  ef <- rasch_efrm(de, item_sets = tr$item_sets, groups = "group",
                   boot_reps = 0)
  we <- stats::setNames(rep(1, length(tr$item_sets)), names(tr$item_sets))
  ze <- weighted_person_estimates(ef, we, by = "set")
  expect_equal(ze$theta, ef$person$theta, tolerance = 1e-7)
  expect_equal(ze$se, ef$person$se, tolerance = 1e-7)
})

test_that("external weights are validated before estimation", {
  set.seed(915)
  f <- rasch(simulate_rasch(n_persons = 80, n_items = 5, seed = 915))
  nm <- colnames(f$X)
  expect_error(weighted_person_estimates(f, numeric()), "non-empty named")
  expect_error(weighted_person_estimates(f, c(a = -1)), "non-negative")
  expect_error(weighted_person_estimates(f, stats::setNames(rep(0, 5), nm)),
               "at least one")
  expect_error(weighted_person_estimates(f, c(I01 = 1)),
               "every fitted item")
  expect_error(weighted_person_estimates(f, stats::setNames(rep(1, 5), nm),
                                         by = "set"), "must map")
  bad_names <- stats::setNames(rep(1, 5), nm); names(bad_names)[1] <- "   "
  expect_error(weighted_person_estimates(f, bad_names), "named numeric")
  expect_error(weighted_person_estimates(
    f, c(A = 1, B = 1), by = "set",
    sets = list(A = c(nm[1:2], NA_character_), B = nm[3:5])),
    "item-name elements")
  failed <- f
  failed$est$converged <- FALSE
  expect_error(weighted_person_estimates(
    failed, stats::setNames(rep(1, 5), nm)), "did not converge")
})

test_that("external weights and set maps cannot be matrices", {
  fit <- rasch(simulate_rasch(120, 4, seed = 732))
  w <- matrix(1, 1L, 4L)
  names(w) <- colnames(fit$X)
  expect_error(weighted_person_estimates(fit, w), "named numeric vector")

  weights <- c(A = 1, B = 1)
  map <- matrix(c("A", "A", "B", "B"), 1L)
  names(map) <- colnames(fit$X)
  expect_error(weighted_person_estimates(fit, weights, by = "set", sets = map),
               "named item-to-set vector")

  sets <- list(A = matrix(colnames(fit$X)[1:2], 1L),
               B = colnames(fit$X)[3:4])
  expect_error(weighted_person_estimates(fit, weights, by = "set", sets = sets),
               "set list")
})

test_that("external-weight extremes follow responses rather than a score tolerance", {
  f <- rasch(simulate_rasch(n_persons = 100, n_items = 3, seed = 918))
  f$X[1, ] <- c(1L, 0L, 0L)
  w <- stats::setNames(c(1e-14, 1, 1), colnames(f$X))
  z <- weighted_person_estimates(f, w)
  expect_false(z$extreme[1])
  expect_gt(z$weighted_score[1], 0)
  expect_lt(z$weighted_score[1], 1e-12)
})
