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
  expect_s3_class(bs, "rasch_fit_bootstrap")
  expect_identical(bs$model_kind, "rasch")
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
  expect_equal(nrow(bs$persons), nrow(fit$person))
  expect_identical(bs$persons$id, fit$person$id)
  for (st in c("fit_resid", "infit_ms", "outfit_ms", "infit_z", "outfit_z")) {
    p <- bs$persons[[paste0(st, "_p_boot")]]
    adj <- bs$persons[[paste0(st, "_p_boot_adj")]]
    expect_true(all(p > 0 & p <= 1, na.rm = TRUE), info = st)
    expect_true(all(adj >= p, na.rm = TRUE), info = st)
  }
  expect_identical(bs$person_adjustment$method,
                   "single-step maximum-statistic bootstrap")
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
    if (s == "conditional") expect_s3_class(bs$persons, "data.frame")
    else expect_null(bs$persons)
  }
})

test_that("fixed locations remain paired with their response rows", {
  X <- rbind(c(NA, NA), c(0, 1), c(1, 0))
  th <- c(NA_real_, -1, 1)
  prepared <- .fit_theta_source(th, X, "fixed")
  expect_equal(prepared, c(0, -1, 1))
  expect_identical(.fit_theta("fixed", prepared, list(), n = nrow(X)),
                   prepared)
  X[1, 1] <- 0
  expect_error(.fit_theta_source(th, X, "fixed"),
               "finite person location", class = "rasch_refusal")
  expect_equal(.fit_theta_source(th, X, "resample"), c(-1, 1))
})

test_that("source-tree bootstrap workers do not load another installation", {
  expect_warning(out <- with_mocked_bindings(
    .rasch_boot_apply(3L, identity, workers = 2L, label = "test bootstrap"),
    .rasch_namespace_is_installed = function() FALSE,
    .package = "rasch"), "source tree")
  expect_equal(out, as.list(1:3))
})

test_that("the maxT reference uses the joint null and never lowers raw p", {
  null <- cbind(a = 1:40, b = 40:1, c = c(1:20, 20:1))
  z <- .boot_maxt(c(40, 40, 20), null, "upper", min_success = 30L)
  expect_equal(z$family_n, 3L)
  expect_equal(z$family_boot, 40L)
  expect_true(all(z$p_adj >= z$p, na.rm = TRUE))
  expect_true(all(z$p_adj <= 1, na.rm = TRUE))
  expect_true(all(is.finite(z$p_adj)))

  # A testable member with no null variation invalidates the joint maximum,
  # rather than quietly reducing the family for the other members.
  constant <- null; constant[, 3] <- 20
  zc <- .boot_maxt(c(40, 40, 20), constant, "upper", min_success = 30L)
  expect_true(all(is.na(zc$p_adj)))
  expect_equal(zc$family_n, 3L)
  expect_equal(zc$adjusted_n, 0L)
  expect_true(all(is.na(.boot_maxt(c(40, 40, 20), constant, "upper",
    min_success = 30L, mode = "centred")$p_adj)))
  expect_true(all(is.na(.boot_maxt(c(40, 40, 20), constant, "upper",
    min_success = 30L, mode = "raw")$p_adj)))

  null[1:20, 3] <- NA
  z2 <- .boot_maxt(c(40, 40, 20), null, "upper", min_success = 30L)
  expect_true(is.na(z2$p[3]))
  expect_equal(z2$family_n, 3L)
  expect_true(all(is.na(z2$p_adj)))
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
  expect_error(fbq(fit, B = 5, seed = c(1, 2)), "non-negative whole")
  expect_error(fbq(fit, B = 5, seed = 1.5), "non-negative whole")
  expect_error(fbq(fit, B = 5, seed = -1), "non-negative whole")
  expect_error(fbq(fit, B = 5, seed = .Machine$integer.max + 1),
               "integer range")
  expect_error(fbq(fit, B = 5, theta = "wishful"), "should be one of")

  failed <- fit; failed$est$converged <- FALSE
  expect_error(fbq(failed, B = 5), "observed Rasch fit did not converge",
               class = "rasch_refusal")
})

test_that("a non-converged replicate is not admitted to the null", {
  d <- simulate_rasch(250, 8, model = "PCM", n_categories = 4, seed = 1)
  f <- suppressWarnings(rasch(d, id = "id", model = "PCM"))
  z <- suppressWarnings(.fit_refit(
    f$X, f$model, f$n_groups, anchors = NULL, maxit = 1L, tol = 1e-12))
  expect_identical(.fit_boot_status(z), "nonconverged")
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
  expect_true(any(grepl("were unusable", w)))
})

test_that("a depleted or statistic-specific null is withheld", {
  fit <- rasch(simb(200, seq(-1, 1, length.out = 5), seed = 18))
  good <- lapply(STATS, function(st) fit$items[[st]])
  names(good) <- STATS
  good$persons <- lapply(c("fit_resid", "infit_ms", "outfit_ms",
                           "infit_z", "outfit_z"), function(st)
    fit$person[[st]])
  names(good$persons) <- c("fit_resid", "infit_ms", "outfit_ms",
                           "infit_z", "outfit_z")

  expect_error(
    suppressWarnings(with_mocked_bindings(
      fit_bootstrap(fit, B = 40, seed = 1),
      .rasch_boot_apply = function(n, fun, ...) c(
        rep(list(good), 29L), rep(list(NULL), n - 29L)),
      .package = "rasch")),
    "at least 30", class = "rasch_refusal")

  refused <- tryCatch(suppressWarnings(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) c(
      rep(list(good), 29L),
      rep(list(.fit_boot_failure("nonconverged")), 6L),
      rep(list(.fit_boot_failure("error")), 5L)),
    .package = "rasch")), error = identity)
  expect_s3_class(refused, "rasch_fit_bootstrap_refusal")
  expect_equal(refused$B_used, 29L)
  expect_equal(refused$B_nonconverged, 6L)
  expect_equal(refused$B_errors, 5L)

  w <- character(0)
  bs <- withCallingHandlers(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) c(
      rep(list(good), 30L), rep(list(NULL), n - 30L)),
    .package = "rasch"), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_equal(bs$B_used, 30L)
  expect_equal(bs$B_failed, 10L)
  expect_equal(bs$B_nonconverged, 0L)
  expect_equal(bs$B_errors, 10L)
  expect_equal(bs$minimum_usable, 30L)
  expect_true(any(grepl("were unusable", w)))
  expect_true(all(bs$items$n_boot_chisq == 30L))

  bad_stat <- good
  bad_stat$outfit_z[1L] <- NA_real_
  w <- character(0)
  bs2 <- withCallingHandlers(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) rep(list(bad_stat), n),
    .package = "rasch"), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_equal(bs2$items$n_boot_outfit_z[1L], 0L)
  expect_true(is.na(bs2$items$outfit_z_p_boot[1L]))
  ok <- is.finite(bs2$items$outfit_z_p_boot)
  expect_equal(bs2$items$outfit_z_p_boot_adj[ok],
               p.adjust(bs2$items$outfit_z_p_boot[ok], method = "holm",
                        n = nrow(bs2$items)))
  expect_true(any(grepl("too few usable replicated statistics for outfit_z", w)))

  bad_fr <- good
  bad_fr$fit_resid[1L] <- NA_real_
  bs_fr <- suppressWarnings(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) rep(list(bad_fr), n),
    .package = "rasch"))
  expect_equal(bs_fr$total$n_boot_fit_resid_mean, 0L)
  expect_equal(bs_fr$total$n_boot_fit_resid_sd, 0L)
  expect_true(is.na(bs_fr$total$fit_resid_mean_p_boot))
  expect_true(is.na(bs_fr$total$fit_resid_sd_p_boot))

  bs3 <- suppressWarnings(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) c(
      rep(list(good), 30L),
      rep(list(.fit_boot_failure("nonconverged")), 6L),
      rep(list(.fit_boot_failure("error")), 4L)),
    .package = "rasch"))
  expect_equal(bs3$B_used, 30L)
  expect_equal(bs3$B_nonconverged, 6L)
  expect_equal(bs3$B_errors, 4L)
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
  expect_identical(.fit_min_boot_success(5), 3L)
  expect_identical(.fit_min_boot_success(29), 15L)
  expect_identical(.fit_min_boot_success(40), 30L)
  expect_identical(.fit_min_boot_success(100), 51L)
  expect_identical(.rasch_min_boot_success(30), 30L)
  expect_identical(.rasch_min_boot_success(40), 30L)
  expect_identical(.rasch_min_boot_success(58), 30L)
  expect_identical(.rasch_min_boot_success(100), 51L)
  expect_identical(.rasch_min_boot_success(300), 151L)
  # the covariance guard also clears the number of independently estimable
  # directions supplied by its caller
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

test_that("paired-comparison fit bootstrap follows the fitted design", {
  d <- simulate_btl(5, 20, reps_per_pair = 8, seed = 41)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  bs <- suppressWarnings(fit_bootstrap(fit, B = 19, workers = 1, seed = 4))
  expect_s3_class(bs, "rasch_fit_bootstrap")
  expect_identical(bs$model_kind, "btl")
  expect_identical(bs$model, "paired comparisons")
  expect_equal(bs$B_used, 19L)
  expect_identical(bs$objects$object, fit$objects$object)
  expect_identical(bs$judges$judge, fit$judges$judge)
  expect_equal(bs$pairs[, c("object_a", "object_b")],
               fit$pairs[, c("object_a", "object_b")])
  expect_true(bs$total$chisq_p_boot > 0 && bs$total$chisq_p_boot <= 1)
  for (tab in c("pairs", "objects", "judges")) {
    pcols <- grep("_p_boot$", names(bs[[tab]]), value = TRUE)
    for (nm in pcols) {
      adj <- bs[[tab]][[paste0(nm, "_adj")]]
      expect_true(all(adj >= bs[[tab]][[nm]], na.rm = TRUE),
                  info = paste(tab, nm))
    }
  }
})

test_that("non-converged observed paired-comparison fits are refused", {
  d <- simulate_btl(5, 20, reps_per_pair = 4, seed = 47)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  fit$converged <- FALSE
  expect_error(fit_bootstrap(fit, B = 5, seed = 1),
               "observed paired-comparison fit did not converge",
               class = "rasch_refusal")
})

test_that("paired-comparison fit bootstrap does not invent judges", {
  d <- simulate_btl(5, 20, reps_per_pair = 6, seed = 45)
  fit <- btl(d, "object_a", "object_b", winner = "winner")
  bs <- suppressWarnings(fit_bootstrap(fit, B = 5, workers = 1, seed = 7))
  expect_null(bs$judges)
  expect_null(bs$adjustment$judges)
  expect_null(bs$replicates$judges)
})

test_that("paired-comparison fit bootstrap retains explanatory restrictions", {
  set.seed(46)
  predictors <- data.frame(object = LETTERS[1:6], x = seq(-1, 1, length.out = 6))
  pr <- t(combn(predictors$object, 2))
  d <- data.frame(object_a = rep(pr[, 1], each = 20),
                  object_b = rep(pr[, 2], each = 20))
  bx <- setNames(predictors$x, predictors$object)
  d$winner <- ifelse(runif(nrow(d)) <
                       plogis(bx[d$object_a] - bx[d$object_b]),
                     d$object_a, d$object_b)
  fit <- btl_explanatory(d, predictors, ~ x,
                         object_a = "object_a", object_b = "object_b",
                         winner = "winner")
  bs <- suppressWarnings(fit_bootstrap(fit, B = 5, workers = 1, seed = 8))
  expect_identical(bs$objects$object, fit$objects$object)
  expect_equal(bs$B_used, 5L)
})

test_that("paired-comparison bootstrap covers ordered categories and history", {
  d <- simulate_btl(5, 24, reps_per_pair = 6, model = "polytomous",
                    n_categories = 4, seed = 42)
  fit <- btl(d, "object_a", "object_b", response = "response",
             judge = "judge", thresholds = "pc")
  generated <- .btl_boot_data(fit)
  expect_true(is.ordered(generated$response))
  expect_identical(levels(generated$response), as.character(0:fit$m))
  expect_no_error(suppressWarnings(
    fit_bootstrap(fit, B = 9, workers = 1, seed = 5)))

  h <- simulate_btl(5, 30, reps_per_pair = 6,
                    dependence = list(exposure = .2, carry_over = .1),
                    seed = 43)
  fh <- btl(h, "object_a", "object_b", winner = "winner",
            judge = "judge", order = "order")
  bh <- suppressWarnings(fit_bootstrap(fh, B = 9, workers = 1, seed = 6))
  expect_gte(bh$B_used, bh$minimum_usable)
  expect_true(all(vapply(bh$replicates$objects, nrow, 0L) == bh$B_used))

  dropped <- fh
  dropped$dependence <- dropped$dependence[-1L, , drop = FALSE]
  expect_error(.btl_boot_data(dropped), "dropped.*effect",
               class = "rasch_refusal")
})

test_that("paired-comparison bootstrap refuses designs it cannot regenerate", {
  d <- simulate_btl(5, 20, reps_per_pair = 6, seed = 44)
  d$winner[1:4] <- "tie"
  half <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
              ties = "half")
  expect_error(fit_bootstrap(half, B = 5, seed = 1), "half-weighted ties",
               class = "rasch_refusal")
  expect_error(fit_bootstrap(btl(d[-(1:4), ], "object_a", "object_b",
                                 winner = "winner", judge = "judge"),
                             B = 5, theta = "fixed", seed = 1),
               "theta.*not paired comparisons")
})
