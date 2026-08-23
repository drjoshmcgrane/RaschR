#!/usr/bin/env Rscript
# Edge-case validation for explanatory Rasch models.
#
# Covers Kent calibration beyond the original balanced cells, sandwich-SE
# calibration when item maximum scores differ, and convergence of the free
# reference calibration at large N.

options(simval.script = "tools/simval/studies/explanatory-edge-cases.R")
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

quick <- identical(commandArgs(trailingOnly = TRUE)[1L], "quick")
n_rep <- if (quick) 100L else 1000L

draw_responses <- function(theta, tau) {
  X <- sapply(tau, function(tt) vapply(theta, function(th)
    sample.int(length(tt) + 1L, 1L,
      prob = item_moments(th, tt)$P) - 1L, integer(1)))
  storage.mode(X) <- "integer"
  X
}

run_cell <- function(scenario, n_persons, m, x, base, seed0) {
  items <- paste0("I", seq_along(m))
  predictors <- data.frame(item = items, x = x)
  estimates <- ses <- r_squared <- p_kent <- p_naive <- numeric(0)
  refused <- nonconv <- 0L
  errors <- character(0)
  for (r in seq_len(n_rep)) {
    set.seed(seed0 + r)
    theta <- rnorm(n_persons)
    tau <- lapply(seq_along(m), function(i)
      base[seq_len(m[i])] + 0.45 * x[i])
    X <- draw_responses(theta, tau)
    colnames(X) <- items
    fit <- tryCatch(rasch_explanatory(
      X, predictors, if (all(m == 1L)) ~ x else ~ x + threshold,
      maxit = 80L, tol = 1e-8), error = function(e) e)
    if (inherits(fit, "error")) {
      refused <- refused + 1L
      errors <- c(errors, conditionMessage(fit))
      next
    }
    if (!isTRUE(fit$est$converged) ||
        !isTRUE(fit$reference_fit$est$converged)) {
      nonconv <- nonconv + 1L
      next
    }
    test <- tryCatch(explanatory_test(fit), error = function(e) e)
    if (inherits(test, "error") || !is.finite(test$p)) {
      refused <- refused + 1L
      errors <- c(errors, if (inherits(test, "error"))
        conditionMessage(test) else "non-finite Kent probability")
      next
    }
    cf <- fit$est$coefficients["x", ]
    estimates <- c(estimates, cf$estimate)
    ses <- c(ses, cf$se)
    r_squared <- c(r_squared, test$r_squared)
    p_kent <- c(p_kent, test$p)
    p_naive <- c(p_naive, test$p_naive)
  }
  n <- length(estimates)
  note <- sprintf(
    "%d persons, %d items, maximum scores %s; first refusal: %s",
    n_persons, length(m), paste(sort(unique(m)), collapse = "/"),
    if (length(errors)) unique(errors)[1L] else "none")
  rbind(
    sv_row("explanatory-edge-cases", scenario, "coefficient_x", n,
      effect = 0.45, bias = mean(estimates - 0.45),
      emp_sd = stats::sd(estimates), mean_se = mean(ses),
      coverage95 = sv_coverage(estimates, ses, 0.45),
      n_attempted = n_rep, n_refused = refused, n_nonconv = nonconv,
      notes = note),
    sv_row("explanatory-edge-cases", scenario, "calibration_r2", n,
      effect = mean(r_squared), emp_sd = stats::sd(r_squared),
      n_attempted = n_rep, n_refused = refused, n_nonconv = nonconv,
      notes = paste(note, "The effect column records mean calibration R-squared")),
    sv_row("explanatory-edge-cases", scenario, "kent_type1", n,
      type1 = mean(p_kent < .05), n_attempted = n_rep,
      n_refused = refused, n_nonconv = nonconv, notes = note),
    sv_row("explanatory-edge-cases", scenario, "naive_type1", n,
      type1 = mean(p_naive < .05), n_attempted = n_rep,
      n_refused = refused, n_nonconv = nonconv,
      notes = paste(note, "Unscaled composite-likelihood reference; diagnostic only"))
  )
}

rows <- rbind(
  run_cell("LLTM_N300", 300L, rep(1L, 12L),
           rep(c(-1, -0.6, -0.2, 0.2, 0.6, 1), 2L), 0, 8100000L),
  run_cell("LPCM_four_category_N300", 300L, rep(3L, 12L),
           rep(c(-1, -0.3, 0.3, 1), 3L), c(-0.8, 0.2, 1), 8200000L),
  run_cell("LPCM_mixed_maximum_N800", 800L, rep(1:3, each = 4L),
           rep(c(-1, -0.3, 0.3, 1), 3L), c(-0.8, 0.2, 1), 8300000L),
  run_cell("LPCM_four_category_N2000", 2000L, rep(3L, 10L),
           rep(c(-1, -0.5, 0, 0.5, 1), 2L), c(-0.8, 0.2, 1), 8400000L)
)

sv_write(rows, if (quick) "explanatory-edge-cases-quick"
         else "explanatory-edge-cases")
