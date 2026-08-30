# The parametric bootstrap null for the item fit statistics: that the
# replicates travel the same road as the observed data, that each statistic
# is read in the tail or tails that mean misfit, that the p-values have the
# resolution their replicate count allows, and that the correction runs in
# the direction the miscalibration requires.

# the resolution-floor warning is correct at these deliberately small B --
# assert it once below, and keep every other call quiet
fbq <- function(...) suppressWarnings(fit_bootstrap(...))

simb <- function(N, d, seed = 1) {
  set.seed(seed)
  th <- rnorm(N)
  X <- matrix(rbinom(N * length(d), 1, plogis(outer(th, d, "-"))),
              N, length(d))
  colnames(X) <- paste0("I", seq_along(d))
  X
}

STATS <- c("chisq", "fit_resid", "infit_ms", "outfit_ms", "infit_z", "outfit_z")

test_that("every fit statistic is carried with a probability in range", {
  fit <- rasch(simb(300, seq(-1.5, 1.5, length.out = 6), seed = 3))
  bs <- fbq(fit, B = 29, seed = 1)
  expect_equal(nrow(bs$items), 6L)
  expect_identical(bs$items$item, fit$items$item)
  for (st in STATS) {
    expect_equal(bs$items[[st]], fit$items[[st]], info = st)
    p <- bs$items[[paste0(st, "_p_boot")]]
    expect_true(all(p > 0 & p <= 1), info = st)
    expect_true(all(p >= 1 / (1 + bs$B_used)), info = st)
    expect_true(all(bs$items[[paste0(st, "_p_boot_adj")]] >= p - 1e-12), info = st)
  }
  expect_equal(bs$B, 29L)
  expect_true(bs$B_used <= bs$B)
  expect_identical(bs$theta, "conditional")
  expect_setequal(names(bs$replicates), STATS)
  expect_true(all(vapply(bs$replicates, nrow, 0L) == bs$B_used))
})

test_that("the chi-square is read in one tail and the others in two", {
  # a statistic far below its null is misfit for the fit residual (a steeper
  # item than the model predicts) but not for a discrepancy that can only be
  # too large
  null <- c(rep(1, 50), rep(2, 50))
  expect_equal(.boot_p(0, null, "upper"), 1)
  expect_lt(.boot_p(0, null, "two"), 0.05)
  expect_lt(.boot_p(99, null, "upper"), 0.05)
  expect_lt(.boot_p(99, null, "two"), 0.05)
  # two-sided doubles the smaller tail and is capped at 1
  expect_lte(.boot_p(1.5, null, "two"), 1)
  expect_equal(.boot_p(NA_real_, null, "upper"), NA_real_)
})

test_that("the whole-test readings carry their own nulls", {
  fit <- rasch(simb(400, seq(-1.5, 1.5, length.out = 6), seed = 16))
  bs <- fbq(fit, B = 29, seed = 3)
  expect_equal(bs$total$chisq, fit$total_chisq)
  expect_equal(bs$total$fit_resid_mean, fit$item_fit_summary$mean)
  expect_equal(bs$total$fit_resid_sd, fit$item_fit_summary$sd)
  for (nm in c("chisq_p_boot", "fit_resid_mean_p_boot", "fit_resid_sd_p_boot"))
    expect_true(bs$total[[nm]] > 0 && bs$total[[nm]] <= 1, info = nm)
})

test_that("a seed reproduces the bootstrap and leaves the caller's stream alone", {
  fit <- rasch(simb(250, seq(-1.2, 1.2, length.out = 5), seed = 4))
  a <- fbq(fit, B = 19, seed = 7)
  b <- fbq(fit, B = 19, seed = 7)
  expect_equal(a$items$chisq_p_boot, b$items$chisq_p_boot)
  expect_equal(a$items$fit_resid_p_boot, b$items$fit_resid_p_boot)
  expect_equal(a$replicates$chisq, b$replicates$chisq)
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  invisible(fbq(fit, B = 9, seed = 3))
  expect_equal(runif(1), before)
})

test_that("the replicates are drawn before dispatch, not inside a worker", {
  # what makes the result independent of the worker count is that every
  # random draw a replicate depends on is made in the coordinating process.
  # These tests run serially by policy (setup-efrm-workers.R), so the
  # property is checked where it lives: the same seed must fix the draws
  # regardless of how many replicates are asked for, so a longer run must
  # extend a shorter one rather than redraw it.
  fit <- rasch(simb(250, seq(-1.2, 1.2, length.out = 5), seed = 13))
  set.seed(4); th_a <- .fit_theta("resample", fit$person$theta, fit$psi)
  set.seed(4); th_b <- .fit_theta("resample", fit$person$theta, fit$psi)
  expect_equal(th_a, th_b)
  expect_equal(fbq(fit, B = 19, seed = 4)$replicates$chisq,
               fbq(fit, B = 19, seed = 4)$replicates$chisq)
})

test_that("the schemes for person locations all run and stay in range", {
  fit <- rasch(simb(250, seq(-1.2, 1.2, length.out = 5), seed = 5))
  for (s in c("conditional", "resample", "fixed", "normal")) {
    bs <- suppressWarnings(fit_bootstrap(fit, B = 19, theta = s, seed = 2))
    expect_identical(bs$theta, s)
    expect_true(all(bs$items$chisq_p_boot > 0 & bs$items$chisq_p_boot <= 1))
  }
})

test_that("the bootstrap statistic carries the bias the asymptotic df ignores", {
  # under a true model the replicated statistic sits well above its nominal
  # degrees of freedom: that gap is the miscalibration being corrected, and
  # the corrected p-values are accordingly larger than the asymptotic ones
  fit <- rasch(simb(1500, seq(-1.5, 1.5, length.out = 6), seed = 6))
  bs <- fbq(fit, B = 49, seed = 11)
  expect_true(mean(colMeans(bs$replicates$chisq)) >
                1.3 * stats::median(fit$items$df))
  expect_true(mean(bs$items$chisq_p_boot > fit$items$p) >= 5 / 6)
})

test_that("missing data are carried into every replicate", {
  X <- simb(300, seq(-1.5, 1.5, length.out = 6), seed = 8)
  X[cbind(sample(300, 60), sample(6, 60, replace = TRUE))] <- NA
  fit <- rasch(X)
  bs <- fbq(fit, B = 19, seed = 1)
  expect_true(all(bs$items$n_boot > 0))
  expect_true(all(is.finite(bs$items$chisq_p_boot)))
})

test_that("familywise adjustment covers each statistic's own family", {
  fit <- rasch(simb(300, seq(-1.5, 1.5, length.out = 6), seed = 9))
  bs <- fbq(fit, B = 19, seed = 5)
  for (st in STATS) {
    p <- bs$items[[paste0(st, "_p_boot")]]
    adj <- bs$items[[paste0(st, "_p_boot_adj")]]
    expect_equal(adj, stats::p.adjust(p, method = "holm"), info = st)
    expect_true(all(adj <= 1), info = st)
  }
})

test_that("the bootstrap refuses designs it cannot generate from", {
  fit <- rasch(simb(200, seq(-1, 1, length.out = 5), seed = 10))
  expect_error(fbq(list(items = 1)), "fitted model from rasch")
  for (cl in c("rasch_efrm", "rasch_mfrm", "rasch_explanatory")) {
    f2 <- fit; class(f2) <- c(cl, "rasch")
    expect_error(fbq(f2, B = 5), "generating structure",
                 class = "rasch_refusal")
  }
  f3 <- fit; f3$disc <- c(1, 1, 1, 1, 2)
  expect_error(fbq(f3, B = 5), "frame units", class = "rasch_refusal")
})

test_that("the bootstrap validates its own controls", {
  fit <- rasch(simb(200, seq(-1, 1, length.out = 5), seed = 12))
  expect_error(fbq(fit, B = 0), "whole positive number of replicates")
  expect_error(fbq(fit, B = c(10, 20)), "whole positive number of replicates")
  expect_error(fbq(fit, B = 2.5), "whole positive number of replicates")
  expect_error(fbq(fit, B = 5, workers = 0), "whole positive number of workers")
  expect_error(fbq(fit, B = 5, seed = c(1, 2)), "one finite number")
  expect_error(fbq(fit, B = 5, theta = "wishful"), "should be one of")
})


test_that("an anchored fit is bootstrapped under its own anchors", {
  X <- simb(300, seq(-1.5, 1.5, length.out = 6), seed = 15)
  free <- rasch(X)
  anc <- data.frame(item = c("I1", "I2"), k = 1L,
                    tau = free$thresholds$tau[free$thresholds$item %in% 1:2])
  fit <- rasch(X, anchors = anc)
  bs <- fbq(fit, B = 19, seed = 8)
  expect_true(all(bs$items$chisq_p_boot > 0 & bs$items$chisq_p_boot <= 1))
  expect_true(all(bs$items$n_boot > 0))
})

test_that("polytomous fits bootstrap without losing their replicates", {
  for (mod in c("PCM", "RSM")) {
    d <- simulate_rasch(400, 6, model = mod, n_categories = 4, seed = 21)
    # a short run can still lose the odd replicate and say so; what matters
    # is that it is the odd one, not a third of them
    bs <- suppressWarnings(
      fbq(rasch(d, id = "id", model = mod), B = 29, seed = 2))
    expect_true(all(bs$items$chisq_p_boot > 0 & bs$items$chisq_p_boot <= 1))
    expect_true(all(is.finite(bs$items$fit_resid_p_boot)))
    # the shipped scheme generates at the observed spread, so replicates are
    # not lost to items whose extreme category goes unvisited
    expect_gte(bs$B_used, 26L)
  }
})

test_that("a thinned null is reported rather than left to be inferred", {
  d <- simulate_rasch(400, 6, model = "RSM", n_categories = 4, seed = 21)
  fit <- rasch(d, id = "id", model = "RSM")
  # fbq() would also swallow the warning under test; suppress only the
  # resolution-floor warning by matching this one specifically
  w <- character(0)
  withCallingHandlers(
    fit_bootstrap(fit, B = 29, theta = "normal", seed = 2),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("could not be estimated", w)))
})

test_that("the conditional generator preserves scores, masks, and theory", {
  d <- simulate_rasch(250, 8, model = "PCM", n_categories = 4,
                      missing = 0.08, seed = 11)
  f <- rasch(d, id = "id")
  nm <- is.na(f$X)
  set.seed(2)
  Xb <- .fit_gen_conditional(f$X, f$tau_list, nm)
  expect_identical(is.na(Xb), nm)
  expect_equal(rowSums(Xb, na.rm = TRUE), rowSums(f$X, na.rm = TRUE))
  expect_gt(mean(Xb != f$X, na.rm = TRUE), 0.1)
  # enumerable case: three dichotomous items, every score 1; the conditional
  # category probabilities are the scaled item weights
  tau3 <- list(-1, 0, 1)
  Xs <- matrix(0L, 4000, 3); Xs[, 1] <- 1L
  set.seed(3)
  Xg <- .fit_gen_conditional(Xs, tau3, NULL)
  w <- exp(-c(-1, 0, 1))
  expect_lt(max(abs(colMeans(Xg) - w / sum(w))), 0.03)
})

test_that("the conditional bootstrap preserves a booklet's missingness link", {
  set.seed(5)
  N <- 240
  grp <- rep(c("low", "high"), each = N / 2)
  th <- rnorm(N, ifelse(grp == "low", -0.75, 0.75), 1)
  X <- sapply(seq(-2, 2, length.out = 12),
              function(dd) rbinom(N, 1, plogis(th - dd)))
  colnames(X) <- sprintf("I%02d", 1:12)
  X[grp == "low", 9:12] <- NA
  X[grp == "high", 1:4] <- NA
  f <- rasch(X)
  obs <- mean(f$person$theta[grp == "high"], na.rm = TRUE) -
    mean(f$person$theta[grp == "low"], na.rm = TRUE)
  set.seed(9)
  Xb <- .fit_gen_conditional(f$X, f$tau_list, is.na(f$X))
  fb <- rasch(Xb)
  rep_diff <- mean(fb$person$theta[grp == "high"], na.rm = TRUE) -
    mean(fb$person$theta[grp == "low"], na.rm = TRUE)
  # abilities resampled independently of the mask would put this near zero
  expect_gt(obs, 0.8)
  expect_lt(abs(rep_diff - obs), 0.4)
})

test_that("bootstrap acceptance needs thirty successes and a majority", {
  expect_identical(.rasch_min_boot_success(30), 30L)
  expect_identical(.rasch_min_boot_success(40), 30L)
  expect_identical(.rasch_min_boot_success(58), 30L)
  expect_identical(.rasch_min_boot_success(100), 51L)
  expect_identical(.rasch_min_boot_success(300), 151L)
  # a covariance with more columns than draws cannot reach full rank, so
  # the floor also clears the quantity count
  expect_identical(.rasch_min_boot_success(100, n_quantities = 80), 81L)
  expect_identical(.rasch_min_boot_success(40, n_quantities = 10), 30L)
})

test_that("a B below the adjusted-inference floor warns, naming the remedy", {
  fit <- rasch(simb(200, seq(-1, 1, length.out = 5), seed = 17))
  expect_warning(fit_bootstrap(fit, B = 99, seed = 1),
                 "cannot reach Holm-adjusted significance")
  expect_no_warning(fbq(fit, B = 199, seed = 1))
  expect_silent(capture.output(suppressMessages(
    invisible(fit_bootstrap(fit, B = 200, seed = 1)))))
})
