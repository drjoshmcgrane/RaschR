# Paired-comparison targeting and information: the design-information
# analogue of the test-information function, its per-object accounting, the
# naive-vs-sandwich standard-error relationship, and the adaptive next-pair
# recommender (Pollitt 2012).

sim_btl_t <- function(beta, n_per_pair, seed = 1, judges = NULL) {
  set.seed(seed)
  pr <- t(utils::combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = n_per_pair),
                  b = rep(pr[, 2], each = n_per_pair),
                  stringsAsFactors = FALSE)
  p <- plogis(beta[d$a] - beta[d$b])
  d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
  if (!is.null(judges)) d$judge <- sample(judges, nrow(d), replace = TRUE)
  d
}

test_that("design information retains dichotomous and graded position effects", {
  for (m in c(1L, 3L)) {
    set.seed(991 + m)
    n <- 1200L
    beta <- setNames(seq(-1, 1, length.out = 5), LETTERS[1:5])
    a <- sample(names(beta), n, TRUE)
    b <- vapply(a, function(z) sample(setdiff(names(beta), z), 1L), "")
    tau <- if (m == 1L) 0 else c(-1.2, 0, 1.2)
    d <- data.frame(a, b)
    eta <- beta[a] - beta[b] + 1.3
    d$response <- vapply(eta, function(z)
      sample(0:m, 1L, prob = item_moments(z, tau)$P), 0L)
    fit <- btl(d, "a", "b", response = "response", position = TRUE)
    info <- btl_information(fit)
    P <- fit$fitted_prob
    E <- drop(P %*% (0:m))
    expected <- drop(P %*% (0:m)^2) - E^2
    expect_equal(info$comparisons$information, expected, tolerance = 1e-12)
    expect_equal(info$total, sum(fit$comparisons$weight * expected))
    expect_equal(sum(info$objects$information), 2 * info$total)
    equal_unit <- .btl_info_of_d(info$comparisons$gap, m, tau =
      if (m > 1L) fit$thresholds$tau else NULL)
    expect_gt(max(abs(expected - equal_unit)), 0.05)

    next_pairs <- btl_next_pairs(fit, weight_se = FALSE)
    position <- fit$dependence$estimate[fit$dependence$effect == "position"]
    expect_equal(next_pairs$expected_information,
      .btl_info_of_d(next_pairs$gap + position, m,
                    if (m > 1L) fit$thresholds$tau else NULL))
    broken <- fit; broken$fitted_prob <- P[-1L, , drop = FALSE]
    expect_error(btl_information(broken), "row-aligned fitted probabilities")
    if (m == 3L) {
      grid <- c(-1, 0, 1)
      curves <- list()
      grDevices::pdf(NULL)
      testthat::with_mocked_bindings(plot_btl_categories(fit, grid),
        lines = function(x, y, ...) curves[[length(curves) + 1L]] <<- y,
        .package = "rasch")
      grDevices::dev.off()
      want <- vapply(grid, function(z)
        item_moments(z + position, fit$thresholds$tau)$P, numeric(4L))
      expect_equal(do.call(rbind, curves), want)
    }
  }
})

test_that("history-dependent design information is conditional on fitted history", {
  d <- simulate_btl(5, 20, reps_per_pair = 50,
    dependence = list(exposure = 0.7, carry_over = 0.4), seed = 993)
  fit <- btl(d, "object_a", "object_b", "winner", judge = "judge", order = "order")
  expect_true(.btl_information_history(fit))
  P <- fit$fitted_prob[, 2L]
  info <- btl_information(fit)
  expect_equal(info$comparisons$information, P * (1 - P))
  expect_error(btl_next_pairs(fit, weight_se = FALSE), "judge and comparison history")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_btl_targeting(fit))

  object <- fit$objects$object[1L]
  cm <- fit$comparisons
  a <- cm$object_a == object; b <- cm$object_b == object
  opp <- c(cm$object_b[a], cm$object_a[b])
  weight <- c(cm$weight[a], cm$weight[b])
  expectation <- c(P[a], 1 - P[b])
  expected <- tapply(weight * expectation, opp, sum) / tapply(weight, opp, sum)
  captured <- list()
  testthat::with_mocked_bindings(
    plot_btl_icc(fit, object, min_n = 1),
    lines = function(x, y, ..., lwd = 1) {
      if (lwd == 3) captured[[length(captured) + 1L]] <<- list(x = x, y = y)
    }, .package = "rasch")
  loc <- fit$objects$location[match(names(expected), fit$objects$object)]
  expect_length(captured, 1L)
  expect_equal(as.numeric(captured[[1L]]$y), as.numeric(expected[order(loc)]))

  group <- setNames(rep(c("A", "B"), length.out = length(unique(cm$judge))),
                     unique(cm$judge))
  expect_no_error(plot_btl_icc(fit, object, group = group, min_n = 1))
})

test_that("dichotomous pair information falls with the location gap", {
  beta <- c(A = -2, B = -1, C = 0, D = 1, E = 2)
  d <- sim_btl_t(beta, 60, seed = 11)
  fit <- btl(d, "a", "b", "win")

  np <- btl_next_pairs(fit, n = choose(length(beta), 2), weight_se = FALSE)
  # gap is oriented to the stronger object, so it is non-negative; the
  # single-comparison information P(1-P) is a strictly decreasing function of
  # the gap, so ordering pairs by gap must order their information downward
  o <- np[order(np$gap), ]
  expect_true(all(diff(o$expected_information) <= 1e-12))
  # and the reference peaks at gap 0: max information at the smallest gap
  expect_equal(which.max(np$expected_information), which.min(np$gap))
})

test_that("per-object design information is the hand-computed sum", {
  # a tiny three-object design; each object central enough to win and lose
  beta <- c(A = -0.7, B = 0, C = 0.7)
  d <- sim_btl_t(beta, 20, seed = 5)
  fit <- btl(d, "a", "b", "win")
  info <- btl_information(fit)

  # rebuild the per-object information by hand from the fitted locations: a
  # comparison contributes w * P(1-P) to BOTH of its objects
  bl <- setNames(fit$objects$location, fit$objects$object)
  cmp <- fit$comparisons
  dd <- bl[cmp$object_a] - bl[cmp$object_b]
  p <- plogis(dd)
  Iw <- cmp$weight * p * (1 - p)
  hand <- tapply(c(Iw, Iw), c(cmp$object_a, cmp$object_b), sum)
  expect_equal(info$objects$information,
               as.numeric(hand[info$objects$object]))
  # total is the weighted sum over comparisons; se_naive inverts the design
  expect_equal(info$total, sum(Iw))
  expect_equal(info$objects$se_naive, 1 / sqrt(info$objects$information))
})

test_that("graded per-comparison information is the item_moments variance", {
  set.seed(3)
  beta <- c(A = -1.2, B = -0.5, C = 0, D = 0.6, E = 1.1)
  beta <- beta - mean(beta)
  tau <- c(-1.4, 0, 1.4)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 50), b = rep(pr[, 2], each = 50),
                  stringsAsFactors = FALSE)
  d$grade <- vapply(seq_len(nrow(d)), function(r) {
    p <- item_moments(beta[d$a[r]] - beta[d$b[r]], tau)$P
    sample(0:3, 1, prob = p)
  }, 0L)
  gfit <- btl(d, "a", "b", response = "grade")
  info <- btl_information(gfit)                     # runs on a graded fit

  # each comparison's information must equal the item_moments score variance
  # at its fitted location gap, to machine precision
  tt <- gfit$thresholds$tau
  cc <- info$comparisons
  expected <- vapply(cc$gap, function(z) item_moments(z, tt)$V, 0)
  expect_equal(cc$information, expected)
  # graded information is the response variance, bounded well above the
  # dichotomous ceiling of 0.25
  expect_true(max(cc$information) > 0.25)
})

test_that("naive and sandwich standard errors track each other", {
  # a clean, model-true design with judges: comparisons are effectively
  # independent, so the design-information SE and the clustered sandwich SE
  # should move together across objects
  beta <- setNames(seq(-2, 2, length.out = 8), LETTERS[1:8])
  d <- sim_btl_t(beta, 40, seed = 7, judges = sprintf("J%02d", 1:20))
  fit <- btl(d, "a", "b", "win", judge = "judge")
  info <- btl_information(fit)
  expect_true(fit$clustered)
  expect_gt(cor(info$objects$se, info$objects$se_naive), 0.7)
})

test_that("next-pair recommender adapts to poorly measured objects", {
  objs <- c("A", "B", "C", "D", "E", "F")
  # F sits in the middle of the scale, so even a handful of its comparisons
  # stay mixed (never extreme) once most are deleted
  beta <- c(A = -2, B = -1.2, C = -0.4, D = 0.6, E = 1.5, F = 0)
  set.seed(42)
  pr <- t(combn(objs, 2))
  dat <- data.frame(a = rep(pr[, 1], each = 25), b = rep(pr[, 2], each = 25),
                    stringsAsFactors = FALSE)
  dat$win <- ifelse(runif(nrow(dat)) < plogis(beta[dat$a] - beta[dat$b]),
                    dat$a, dat$b)
  fit_full <- btl(dat, "a", "b", "win")

  # weight_se = FALSE: the top pair is a closest-location pair
  all_np <- btl_next_pairs(fit_full, n = choose(length(objs), 2),
                           weight_se = FALSE)
  top <- btl_next_pairs(fit_full, n = 1, weight_se = FALSE)
  expect_equal(top$gap[1], min(all_np$gap))

  # n_existing is the count of comparisons already observed for a pair
  cnt_AB <- sum((dat$a == "A" & dat$b == "B") | (dat$a == "B" & dat$b == "A"))
  row_AB <- all_np[(all_np$object_a == "A" & all_np$object_b == "B") |
                   (all_np$object_a == "B" & all_np$object_b == "A"), ]
  expect_equal(row_AB$n_existing, cnt_AB)

  # delete most of F's comparisons, keeping three per pair, and refit: F is
  # now the least-measured object (largest se), so weight_se = TRUE must
  # promote its pairs into the recommendations
  isF <- dat$a == "F" | dat$b == "F"
  pk <- paste(pmin(dat$a, dat$b), pmax(dat$a, dat$b))
  keep <- rep(TRUE, nrow(dat))
  for (k in unique(pk[isF])) {
    rows <- which(pk == k)
    keep[rows[-(1:3)]] <- FALSE
  }
  fit2 <- btl(dat[keep, ], "a", "b", "win")
  ob <- fit2$objects
  expect_equal(ob$object[which.max(ob$se)], "F")     # F is worst measured

  wnp <- btl_next_pairs(fit2, n = 5, weight_se = TRUE)
  expect_true("F" %in% c(wnp$object_a, wnp$object_b))
  # and se-weighting reorders relative to raw information
  unw <- btl_next_pairs(fit2, n = 5, weight_se = FALSE)
  expect_false(identical(paste(wnp$object_a, wnp$object_b),
                         paste(unw$object_a, unw$object_b)))

  invalid <- fit2
  invalid$cov_beta[1, 1] <- -1e6
  expect_error(btl_next_pairs(invalid, weight_se = TRUE),
               "positive-semidefinite")
  expect_no_error(btl_next_pairs(invalid, weight_se = FALSE))

  misaligned <- fit2
  rownames(misaligned$cov_beta)[1] <- "not an object"
  expect_error(btl_next_pairs(misaligned, weight_se = TRUE),
               "cannot be aligned")

  no_cov <- fit2
  no_cov$cov_beta <- NULL
  no_cov$objects$se[1] <- NA_real_
  expect_error(btl_next_pairs(no_cov, weight_se = TRUE),
               "standard errors")
  expect_no_error(btl_next_pairs(no_cov, weight_se = FALSE))
  expect_error(btl_next_pairs(fit2, weight_se = matrix(TRUE)),
               "TRUE or FALSE")
})

test_that("the targeting plot draws without error", {
  beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
  d <- sim_btl_t(beta, 30, seed = 2)
  fit <- btl(d, "a", "b", "win")

  set.seed(3)
  bg <- c(A = -1.2, B = -0.4, C = 0.3, D = 0.5, E = 1)
  pr <- t(combn(names(bg), 2))
  dg <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40),
                   stringsAsFactors = FALSE)
  dg$grade <- vapply(seq_len(nrow(dg)), function(r) {
    p <- item_moments(bg[dg$a[r]] - bg[dg$b[r]], c(-1, 0, 1))$P
    sample(0:3, 1, prob = p)
  }, 0L)
  gfit <- btl(dg, "a", "b", response = "grade")

  pdf(NULL)
  on.exit(dev.off())
  expect_error(plot_btl_targeting(fit), NA)
  expect_error(plot_btl_targeting(gfit), NA)
})
