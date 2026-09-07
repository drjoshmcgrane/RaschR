test_that("dichotomous comparison information is symmetric in the tails", {
  gap <- c(0, 1, 20, 40, 100, 700)
  positive <- .btl_info_of_d(gap, 1L, NULL)
  negative <- .btl_info_of_d(-gap, 1L, NULL)
  expect_true(all(positive > 0))
  expect_identical(positive, negative)
  # The two-category model supplies the independent probability reference.
  reference <- vapply(gap, function(d) prod(item_moments(d, 0)$P), 0)
  expect_equal(positive, reference, tolerance = 1e-14)
  expect_equal(positive[gap == 40] / exp(-40), 1, tolerance = 1e-14)
  expect_identical(.btl_info_of_d(c(710, 740), 1L, NULL), exp(-c(710, 740)))
})

test_that("ordered comparison information need not peak at zero gap", {
  information <- .btl_info_of_d(c(-3, 0, 3), 2L, c(-3, 3))
  expect_equal(information[1L], information[3L], tolerance = 1e-14)
  expect_gt(information[3L], information[2L])
})

test_that("next-pair information retains a distant anchored comparison", {
  pairs <- t(utils::combn(LETTERS[1:3], 2))
  d <- data.frame(a = rep(pairs[, 1L], each = 20L),
                   b = rep(pairs[, 2L], each = 20L))
  d$winner <- ifelse(seq_len(nrow(d)) %% 2L, d$a, d$b)
  fit <- btl(d, "a", "b", winner = "winner", anchors = c(A = -20, C = 20))
  expect_true(fit$converged)
  result <- btl_next_pairs(fit, n = 3L, weight_se = FALSE)
  tail <- result[result$object_a == "C" & result$object_b == "A", ]
  expect_identical(nrow(tail), 1L)
  expect_equal(tail$expected_information / exp(-40), 1, tolerance = 1e-14)
})

test_that("targeting refuses an unrepresentable reference before changing layout", {
  d <- simulate_btl(n_objects = 5, n_judges = 15, reps_per_pair = 25,
                     seed = 9571)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(2, 1), mar = c(3, 2, 1, 1))
  before <- graphics::par(c("mfrow", "mar", "mfg"))
  for (grid in list(c(-1001, -1000), c(1000, 1001))) {
    expect_error(plot_btl_targeting(fit, grid = grid),
                  "grid contains no numerically measurable comparison information")
    expect_identical(graphics::par(c("mfrow", "mar", "mfg")), before)
  }
  # Information is tiny but representable here; normalising the curve before
  # scaling it avoids overflow in the right-axis multiplier.
  for (grid in list(c(-711, -710), c(710, 711)))
    expect_no_error(plot_btl_targeting(fit, grid = grid))
})
