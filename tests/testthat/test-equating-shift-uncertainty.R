test_that("a sole usable exact anchor fixes the Rasch and CJ shift", {
  set.seed(735)
  X <- matrix(rbinom(2400, 1, .5), 400, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  f <- rasch(X, anchors = data.frame(item = "I1", k = 1, tau = 0))
  bank <- f$items[, c("item", "location", "se", "max")]
  bank$location <- bank$location - c(1, rep(3, 5))
  bank$se <- c(0, rep(NA_real_, 5))
  z <- equate_tests(f, bank)
  expect_equal(z$shift, 1)
  expect_equal(z$table$adj_difference[1], 0)
  expect_false(z$inferential)

  d <- simulate_btl(6, 20, reps_per_pair = 5, seed = 736)
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge",
           anchors = c(O1 = 0))
  bank <- f$objects[, c("object", "location", "se")]
  anchor <- bank$object == "O1"
  bank$location <- bank$location - ifelse(anchor, 1, 3)
  bank$se <- ifelse(anchor, 0, NA_real_)
  z <- btl_equate(f, bank)
  expect_equal(z$shift, 1)
  expect_equal(z$shift_se, 0)
  expect_equal(z$equated$location[anchor], 0)
  expect_equal(z$equated$se[anchor], 0)
  expect_true(all(is.na(z$equated$se[!anchor])))
  expect_false(z$inferential)
})

test_that("a fixed CJ bank inherits estimated shift uncertainty", {
  d <- simulate_btl(6, 30, reps_per_pair = 5, seed = 737)
  tr <- attr(d, "truth")
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  bank <- data.frame(object = c("O1", "O2", "O3", "Onew"),
    location = c(unname(tr$location[c("O1", "O2", "O3")]), 2), se = 0)
  z <- btl_equate(f, bank)
  expect_true(z$inferential)
  expect_gt(z$shift_se, .1)
  expect_equal(z$equated$se, rep(z$shift_se, nrow(bank)))
  expect_equal(unname(attr(z$equated, "cov_location")),
               matrix(z$shift_se^2, nrow(bank), nrow(bank)))
  expect_equal(attr(z$equated, "df_location"), 29)
})

test_that("CJ equated covariance retains location-shift cross terms", {
  fit <- function(n, seed) {
    d <- simulate_btl(n, 30, reps_per_pair = 5, seed = seed)
    btl(d, "object_a", "object_b", "winner", judge = "judge")
  }
  a <- fit(6, 738)
  b <- fit(7, 739)
  z <- btl_equate(a, b, independent = TRUE)
  common <- z$table$object
  u <- 1 / (z$table$se_1^2 + z$table$se_2^2)
  u <- u / sum(u)
  A <- matrix(0, nrow(b$objects), nrow(a$objects))
  A[, match(common, a$objects$object)] <- matrix(u, nrow(A), length(u), byrow = TRUE)
  B <- diag(nrow(b$objects))
  B[, match(common, b$objects$object)] <-
    B[, match(common, b$objects$object)] - matrix(u, nrow(B), length(u), byrow = TRUE)
  expected <- A %*% a$cov_beta %*% t(A) + B %*% b$cov_beta %*% t(B)
  expect_equal(unname(attr(z$equated, "cov_location")), unname(expected))
  expect_equal(z$equated$se, unname(sqrt(diag(expected))))
  expect_gt(max(abs(z$equated$se^2 - (b$objects$se^2 + z$shift_se^2))), .001)

  unknown <- btl_equate(a, b)
  expect_true(all(is.na(unknown$equated$se)))
  expect_null(attr(unknown$equated, "cov_location"))
  expect_match(paste(unknown$notes, collapse = " "), "Equated location standard errors.*withheld")
  none <- btl_equate(a, b, shift = "none")
  expect_equal(none$equated$se, b$objects$se)

  # Marginal SEs alone do not identify the cross term. Sparse covariance
  # metadata must not produce a plausible but conditional shifted SE.
  bank <- b$objects[, c("object", "location", "se")]
  marginal <- btl_equate(a, bank)
  expect_true(all(is.na(marginal$equated$se)))
  attr(bank, "cov_location") <- b$cov_beta
  attr(bank, "df_location") <- 29
  joint <- btl_equate(a, bank)
  expect_equal(joint$equated$se, z$equated$se)
  expect_equal(attr(joint$equated, "cov_location"), attr(z$equated, "cov_location"))
  # Reordering a named bank must reorder, not change, its propagated covariance.
  perm <- c(7, 3, 1, 5, 2, 6, 4)
  bank_perm <- bank[perm, ]
  attr(bank_perm, "cov_location") <- b$cov_beta
  reordered <- btl_equate(a, bank_perm)
  expect_equal(reordered$equated$se, joint$equated$se[perm])
  expect_equal(attr(reordered$equated, "cov_location"),
               attr(joint$equated, "cov_location")[perm, perm])
})
