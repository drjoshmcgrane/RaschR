#!/usr/bin/env Rscript
# Validation of explanatory Rasch and comparative judgement models.
#
# Principal null cells assess the Kent comparison with a free calibration and
# coefficient bias, sandwich-SE calibration and coverage. Separate cells assess
# the familywise behaviour and power of Holm-adjusted fixed-departure
# diagnostics. Item/object structures are fixed within a scenario; only persons,
# judgements and responses are redrawn.
#
# Usage:
#   Rscript tools/simval/studies/explanatory-models.R        # full study
#   Rscript tools/simval/studies/explanatory-models.R quick  # short smoke run

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

quick <- identical(commandArgs(trailingOnly = TRUE)[1], "quick")
n_principal <- if (quick) 50L else 1000L
n_diag <- if (quick) 20L else 300L

draw_pcm <- function(theta, tau) {
  X <- sapply(tau, function(tt) vapply(theta, function(b)
    sample.int(length(tt) + 1L, 1L, prob = item_moments(b, tt)$P) - 1L,
    integer(1)))
  storage.mode(X) <- "integer"
  X
}

summarise_cell <- function(study, scenario, truth, estimates, ses, p_global,
                           n_attempted, n_refused, n_nonconv, notes) {
  ok <- is.finite(estimates) & is.finite(ses) & ses > 0
  rbind(
    sv_row(study, scenario, "coefficient", sum(ok),
      bias = mean(estimates[ok] - truth), emp_sd = stats::sd(estimates[ok]),
      mean_se = mean(ses[ok]), coverage95 = sv_coverage(estimates, ses, truth),
      effect = truth, n_attempted = n_attempted, n_refused = n_refused,
      n_nonconv = n_nonconv, notes = notes),
    sv_row(study, scenario, "global_kent_type1", sum(is.finite(p_global)),
      type1 = mean(p_global < .05, na.rm = TRUE),
      n_attempted = n_attempted, n_refused = n_refused,
      n_nonconv = n_nonconv,
      notes = paste(notes, "Kent comparison with free calibration; nominal alpha=.05"))
  )
}

run_rasch_cell <- function(kind, nrep, seed0) {
  item <- paste0("I", 1:10)
  format <- rep(c("A", "B"), 5)
  operation <- rep(0:1, each = 5)
  pred <- data.frame(item, operation, format, stringsAsFactors = FALSE)
  est <- se <- pg <- numeric(0); refused <- nonconv <- 0L
  for (r in seq_len(nrep)) {
    set.seed(seed0 + r); theta <- rnorm(600)
    if (kind == "LLTM") {
      delta <- .7 * operation + .35 * (format == "B")
      X <- sapply(delta, function(d) rbinom(length(theta), 1, plogis(theta - d)))
      form <- ~ operation + format; term <- "operation"; truth <- .7
    } else {
      tau <- lapply(seq_along(item), function(i) {
        fb <- as.integer(format[i] == "B")
        c(-.7 + .2 * fb, .7 + .6 * fb)
      })
      X <- draw_pcm(theta, tau)
      form <- ~ format + threshold + format:threshold
      term <- "formatB:threshold2"; truth <- .4
    }
    colnames(X) <- item
    f <- tryCatch(rasch_explanatory(X, pred, form), error = function(e) e)
    if (inherits(f, "error")) { refused <- refused + 1L; next }
    if (!isTRUE(f$est$converged) || !isTRUE(f$reference_fit$est$converged)) {
      nonconv <- nonconv + 1L; next
    }
    cf <- f$est$coefficients[term, ]
    z <- tryCatch(explanatory_test(f), error = function(e) e)
    if (inherits(z, "error") || !is.finite(z$p_kent)) {
      refused <- refused + 1L; next
    }
    est <- c(est, cf$estimate); se <- c(se, cf$se); pg <- c(pg, z$p_kent)
  }
  summarise_cell("explanatory-models", kind, truth, est, se, pg, nrep,
                 refused, nonconv,
                 "600 persons, 10 items; fixed balanced predictor design")
}

make_cj <- function(seed, graded = FALSE, judges = FALSE, departure = 0) {
  set.seed(seed)
  pred <- data.frame(object = LETTERS[1:8], domain = rep(0:1, each = 4),
                     format = rep(c("A", "B"), 4), stringsAsFactors = FALSE)
  beta <- .8 * pred$domain + .35 * (pred$format == "B")
  beta[1] <- beta[1] + departure
  names(beta) <- pred$object
  pr <- t(utils::combn(pred$object, 2L))
  reps <- if (judges) 20L else 30L
  d <- data.frame(a = rep(pr[, 1], each = reps),
                  b = rep(pr[, 2], each = reps), stringsAsFactors = FALSE)
  if (judges) d$judge <- rep(sprintf("J%02d", seq_len(reps)), nrow(pr))
  if (!graded) {
    p <- plogis(beta[d$a] - beta[d$b])
    d$winner <- ifelse(runif(nrow(d)) < p, d$a, d$b)
  } else {
    tau <- c(-.8, .8)
    d$response <- vapply(seq_len(nrow(d)), function(i)
      sample.int(3L, 1L,
        prob = item_moments(beta[d$a[i]] - beta[d$b[i]], tau)$P) - 1L,
      integer(1))
  }
  list(data = d, predictors = pred)
}

run_cj_cell <- function(graded, judges, nrep, seed0) {
  est <- se <- pg <- numeric(0); refused <- nonconv <- 0L
  for (r in seq_len(nrep)) {
    d <- make_cj(seed0 + r, graded, judges)
    args <- list(d$data, d$predictors, ~ domain + format,
                 object_a = "a", object_b = "b",
                 judge = if (judges) "judge" else NULL)
    args <- c(args, if (graded) list(response = "response")
              else list(winner = "winner"))
    f <- tryCatch(do.call(btl_explanatory, args), error = function(e) e)
    if (inherits(f, "error")) { refused <- refused + 1L; next }
    if (!isTRUE(f$converged) || !isTRUE(f$reference_fit$converged)) {
      nonconv <- nonconv + 1L; next
    }
    cf <- f$object_coefficients["domain", ]
    z <- tryCatch(explanatory_test(f), error = function(e) e)
    if (inherits(z, "error") || !is.finite(z$p_kent)) {
      refused <- refused + 1L; next
    }
    est <- c(est, cf$estimate); se <- c(se, cf$se); pg <- c(pg, z$p_kent)
  }
  scenario <- paste(if (graded) "ordered_CJ" else "dichotomous_CJ",
                    if (judges) "judge_clustered" else "independent",
                    sep = "_")
  summarise_cell("explanatory-models", scenario, .8, est, se, pg, nrep,
                 refused, nonconv,
                 "8 objects, all pairs; fixed balanced object-predictor design")
}

run_diagnostic_cell <- function(family = c("LLTM", "CJ"), effect, nrep,
                                seed0) {
  family <- match.arg(family)
  reject_any <- detect <- logical(0); refused <- nonconv <- 0L
  for (r in seq_len(nrep)) {
    if (family == "CJ") {
      d <- make_cj(seed0 + r, departure = effect)
      f <- tryCatch(btl_explanatory(d$data, d$predictors, ~ domain + format,
        "a", "b", winner = "winner"), error = function(e) e)
    } else {
      set.seed(seed0 + r)
      pred <- data.frame(item = paste0("I", 1:8),
        feature = rep(0:1, each = 4), stringsAsFactors = FALSE)
      delta <- .7 * pred$feature; delta[1] <- delta[1] + effect
      theta <- rnorm(500)
      X <- sapply(delta, function(d) rbinom(500, 1, plogis(theta - d)))
      colnames(X) <- pred$item
      f <- tryCatch(rasch_explanatory(X, pred, ~ feature),
                    error = function(e) e)
    }
    if (inherits(f, "error")) { refused <- refused + 1L; next }
    conv <- if (family == "CJ") f$converged else f$est$converged
    if (!isTRUE(conv)) { nonconv <- nonconv + 1L; next }
    dg <- tryCatch(explanatory_diagnostics(f), error = function(e) e)
    if (inherits(dg, "error")) { refused <- refused + 1L; next }
    nm <- if (family == "CJ") dg$object else dg$item
    target <- if (family == "CJ") "A" else "I1"
    reject_any <- c(reject_any, any(dg$p_adj < .05, na.rm = TRUE))
    detect <- c(detect, any(nm == target & dg$p_adj < .05, na.rm = TRUE))
  }
  n <- length(reject_any)
  scenario <- paste0(family, "_fixed_departure_", effect)
  if (effect == 0) {
    sv_row("explanatory-models", scenario, "diagnostic_familywise_type1", n,
      familywise = mean(reject_any), effect = effect, n_attempted = nrep,
      n_refused = refused, n_nonconv = nonconv,
      notes = "Holm family over every available fixed departure under the null")
  } else {
    rbind(
      sv_row("explanatory-models", scenario, "any_departure_power", n,
        power = mean(reject_any), effect = effect, n_attempted = nrep,
        n_refused = refused, n_nonconv = nonconv,
        notes = "probability that the Holm family detects any fixed departure"),
      sv_row("explanatory-models", scenario, "planted_departure_power", n,
        power = mean(detect), effect = effect, n_attempted = nrep,
        n_refused = refused, n_nonconv = nonconv,
        notes = "power for the planted departure on the first item or object")
    )
  }
}

rows <- rbind(
  run_rasch_cell("LLTM", n_principal, 7100000L),
  run_rasch_cell("LPCM", n_principal, 7200000L),
  run_cj_cell(FALSE, FALSE, n_principal, 7300000L),
  run_cj_cell(TRUE, FALSE, n_principal, 7400000L),
  run_cj_cell(FALSE, TRUE, n_principal, 7500000L),
  run_diagnostic_cell("LLTM", 0, n_diag, 7600000L),
  run_diagnostic_cell("LLTM", .8, n_diag, 7700000L),
  run_diagnostic_cell("CJ", 0, n_diag, 7800000L),
  run_diagnostic_cell("CJ", .8, n_diag, 7900000L))

sv_write(rows, if (quick) "explanatory-models-quick" else "explanatory-models")
