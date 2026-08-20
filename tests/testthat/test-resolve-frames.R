
efrm_with_dif <- function(N = 400, K = 10, shift = 1.2, seed = 3) {
  set.seed(seed)
  u <- c(1.4^-0.5, 1.4^0.5)
  delta <- seq(-1.6, 1.6, length.out = K)
  items <- sprintf("I%02d", seq_len(K))
  sh <- rep(0, K); sh[4] <- shift
  mk <- function(g, s) {
    th <- stats::rnorm(N, 0, 1.3)
    vapply(seq_len(K), function(i)
      stats::rbinom(N, 1, stats::plogis(u[g] * (th - delta[i] - s[i]))),
      numeric(N))
  }
  X <- rbind(mk(1, rep(0, K)), mk(2, sh))
  colnames(X) <- items
  d <- data.frame(id = sprintf("P%04d", seq_len(2 * N)), X,
                  group = rep(c("g1", "g2"), each = N), check.names = FALSE)
  list(d = d, items = items,
       fit = rasch_efrm(d, items = items, item_sets = list(set1 = items),
                        groups = "group", id = "id", boot_reps = 0))
}

test_that("resolving keeps the item measuring inside each frame", {
  z <- efrm_with_dif()
  f <- z$fit
  r <- resolve_frames(f, "I04", boot_reps = 0)
  g <- drop_items(f, "I04", boot_reps = 0)

  expect_s3_class(r, "rasch_efrm")
  # one source item per frame, each answered by that frame's persons alone
  expect_true(all(c("I04 (g1)", "I04 (g2)") %in% names(r$set_of)))
  expect_false("I04" %in% names(r$set_of))
  expect_equal(length(names(r$set_of)), length(names(f$set_of)) + 1L)
  expect_setequal(unique(r$set_of), unique(f$set_of))

  # the point of resolving: the item still contributes to a person's score,
  # which dropping it cannot
  expect_equal(max(r$person$max_raw), max(f$person$max_raw))
  expect_lt(max(g$person$max_raw), max(f$person$max_raw))
  expect_equal(nrow(r$person), nrow(f$person))

  expect_match(paste(r$notes, collapse = " "), "I04 resolved by frame")
})

test_that("the two versions of a resolved item take separate locations", {
  z <- efrm_with_dif(shift = 1.2)
  r <- resolve_frames(z$fit, "I04", boot_reps = 0)
  loc <- r$item_arbitrary$location[
    match(c("I04 (g1)", "I04 (g2)"), r$item_arbitrary$item)]
  expect_true(all(is.finite(loc)))
  # the tie is broken, so the planted shift shows up as a location difference
  expect_gt(abs(diff(loc)), 0.5)
})

test_that("resolving refuses to leave the group units unidentified", {
  z <- efrm_with_dif()
  # the group units are identified by the items two frames share, and a
  # resolved item is shared by nobody
  expect_error(resolve_frames(z$fit, z$items, boot_reps = 0),
               "fewer than two items common to the frames")
  expect_error(resolve_frames(z$fit, z$items[1:9], boot_reps = 0),
               "unidentified")
  expect_error(resolve_frames(z$fit, "NOPE"), "not in the fit")
  expect_error(resolve_frames(z$fit, character()), "at least one item")
})

test_that("resolve_frames refuses the fits it cannot serve", {
  d <- simulate_rasch(200, 6, seed = 1)
  f <- rasch(d, id = "id")
  expect_error(resolve_frames(f, "I01"), "needs a frame model")
  expect_error(resolve_frames(f, "I01"), "split_items")
})

test_that("the unit-linking bootstrap accepts the replicate count it is given", {
  # a flat floor of 30 rejected every boot_reps below 30 without a single
  # replicate having failed, and reported it as a weak linking design
  skip_on_cran()
  d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 2,
                     n_groups = 2, set_unit_ratio = 1.3, seed = 11)
  tr <- attr(d, "truth")
  for (b in c(5, 20, 29)) {
    f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                    boot_reps = b)
    expect_s3_class(f, "rasch_efrm")
    expect_true(all(is.finite(f$alpha_table$alpha)))
  }
})
