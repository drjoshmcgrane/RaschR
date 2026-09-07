test_that("tiny weights on an observed subset retain person estimates", {
  set.seed(759)
  X <- matrix(rbinom(400 * 6, 1, .5), 400, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  X[1:100, 1:3] <- NA
  fit <- rasch(X)
  equal <- weighted_person_estimates(fit, setNames(rep(1, 6), colnames(X)))
  expect_true(all(is.finite(equal$theta[1:100])))
  expect_true(any(equal$extreme[1:100]))
  expect_true(any(!equal$extreme[1:100]))
  for (tiny in c(1e-100, 1e-160, 1e-200)) {
    weights <- setNames(c(rep(1, 3), rep(tiny, 3)), colnames(X))
    actual <- weighted_person_estimates(fit, weights)
    expect_equal(actual$theta[1:100], equal$theta[1:100], tolerance = 1e-8)
    expect_equal(actual$se[1:100], equal$se[1:100], tolerance = 1e-8)
    # Unequal weights within the same pattern must also retain their ratios.
    weights[4:6] <- tiny * c(1, 2, 3)
    actual <- weighted_person_estimates(fit, weights)
    reference <- weighted_person_estimates(fit,
      setNames(c(1, 1, 1, 1, 2, 3), colnames(X)))
    expect_equal(actual$theta[1:100], reference$theta[1:100], tolerance = 1e-8)
    expect_equal(actual$se[1:100], reference$se[1:100], tolerance = 1e-8)
  }
})

test_that("equal external weights reproduce ordinary person estimates", {
  set.seed(912)
  d <- simulate_rasch(n_persons = 120, n_items = 6, seed = 912)
  f <- rasch(d)
  w <- stats::setNames(rep(3, ncol(f$X)), colnames(f$X))
  z <- weighted_person_estimates(f, w)
  expect_identical(anyDuplicated(names(z)), 0L)
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
  expect_invisible(.validate_weighted_person_table(z, f))
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
  duplicated_vector <- set_vector
  names(duplicated_vector)[2L] <- paste0("  ", colnames(f$X)[1L], "  ")
  expect_error(weighted_person_estimates(
    f, c(A = 2, B = 1), by = "set", sets = duplicated_vector),
    "exactly once")
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
  extreme_range <- stats::setNames(
    c(.Machine$double.xmin * .Machine$double.eps,
      .Machine$double.xmax, rep(1, length(nm) - 2L)), nm)
  expect_error(weighted_person_estimates(f, extreme_range),
               "relative range")
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

test_that("person moments refuse finite inputs whose combined scale overflows", {
  expect_error(item_moments(.Machine$double.xmax,
                            c(-.Machine$double.xmax), disc = 2),
               "numerically representable")
  expect_error(item_moments(0, rep(.Machine$double.xmax, 2)),
               "numerically representable")
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

test_that("app weighted measures remain signed and enter reports and archives", {
  f <- rasch(simulate_rasch(n_persons = 90, n_items = 5, seed = 920))
  weights <- stats::setNames(c(2, 2, 1, 1, 0.5), colnames(f$X))
  state <- list(
    table = weighted_person_estimates(f, weights),
    weights = weights, by = "item", sets = NULL,
    filename = "weights.csv", fit_signature = .fit_boot_signature(f))
  state$result_signature <- .fit_boot_md5(state)
  attr(f, "report_person_weights") <- state
  plain_fit <- f
  attr(plain_fit, "report_person_weights") <- NULL

  html <- tempfile(fileext = ".html")
  explicit_html <- tempfile(fileext = ".html")
  out <- tempfile("rasch-weighted-export-")
  on.exit(unlink(c(html, explicit_html, out), recursive = TRUE, force = TRUE),
          add = TRUE)
  expect_identical(suppressWarnings(report_html(f, html)), html)
  page <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(page, "Externally weighted secondary person measures",
               fixed = TRUE)
  expect_match(page, "do not replace the ordinary person estimates",
               fixed = TRUE)
  expect_match(page, "Resolved external weights", fixed = TRUE)
  expect_identical(suppressWarnings(report_html(
    plain_fit, explicit_html, person_weights = state$table)), explicit_html)
  explicit_page <- paste(readLines(explicit_html, warn = FALSE),
                         collapse = "\n")
  expect_match(explicit_page, "Resolved external weights", fixed = TRUE)
  if (requireNamespace("rmarkdown", quietly = TRUE) &&
      rmarkdown::pandoc_available()) {
    rendered <- tempfile(fileext = ".html")
    on.exit(unlink(rendered, force = TRUE), add = TRUE)
    suppressWarnings(report_document(
      plain_fit, rendered, format = "html", person_weights = state$table))
    templated <- paste(readLines(rendered, warn = FALSE), collapse = "\n")
    expect_match(templated,
                 "Externally weighted secondary person measures",
                 fixed = TRUE)
    expect_match(templated, "Resolved external weights", fixed = TRUE)
  }

  suppressWarnings(save_outputs(
    plain_fit, out, formats = "png", item_plots = FALSE, dpi = 72,
    person_weights = state$table))
  csv <- file.path(out, "tables",
                   "person_estimates_externally_weighted.csv")
  weights_csv <- file.path(out, "tables",
                           "person_estimate_external_weights.csv")
  expect_true(file.exists(csv))
  expect_true(file.exists(weights_csv))
  saved <- utils::read.csv(csv, check.names = FALSE)
  saved_weights <- utils::read.csv(weights_csv, check.names = FALSE)
  expect_equal(saved$theta, state$table$theta, tolerance = 1e-12)
  expect_equal(saved$se, state$table$se, tolerance = 1e-12)
  expected_weights <- attr(state$table, "weighting", exact = TRUE)
  # read.csv() infers an all-missing character column as logical; compare its
  # values after restoring the declared type rather than weakening the check.
  saved_weights$set <- as.character(saved_weights$set)
  expect_equal(saved_weights, expected_weights, tolerance = 1e-12)

  changed_fit <- rasch(simulate_rasch(90, 5, seed = 921))
  attr(changed_fit, "report_person_weights") <- state
  expect_error(report_html(changed_fit, tempfile(fileext = ".html")),
               "different fitted model")

  changed_table <- f
  bad_state <- state
  bad_state$table$theta[1L] <- bad_state$table$theta[1L] + 0.01
  attr(changed_table, "report_person_weights") <- bad_state
  expect_error(report_html(changed_table, tempfile(fileext = ".html")),
               "changed since it was calculated")

  tampered <- state$table
  tampered$theta[1L] <- tampered$theta[1L] + 0.01
  expect_error(report_html(
    plain_fit, tempfile(fileext = ".html"), person_weights = tampered),
    "does not reproduce")
})
