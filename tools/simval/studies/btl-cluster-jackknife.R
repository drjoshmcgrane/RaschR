# Comparison of the current CR1 judge-clustered sandwich with a delete-one-
# judge cluster jackknife. This is a validation study, not an implementation
# shortcut: row-level HC3 is invalid because comparisons remain dependent
# within judges.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "500"))

simulate_panel <- function(workloads, effect = 0, judge_sd = 0.5) {
  objects <- LETTERS[1:4]
  beta <- c(effect, rep(-effect / 3, 3)); names(beta) <- objects
  pairs <- t(utils::combn(objects, 2L))
  rows <- vector("list", length(workloads))
  for (j in seq_along(workloads)) {
    u <- rnorm(4L, 0, judge_sd); u <- u - mean(u); names(u) <- objects
    take <- rep(seq_len(nrow(pairs)), length.out = workloads[j])
    take <- sample(take)
    a <- pairs[take, 1L]; b <- pairs[take, 2L]
    y <- rbinom(length(take), 1, plogis(beta[a] + u[a] - beta[b] - u[b]))
    rows[[j]] <- data.frame(a = a, b = b, y = y,
                            judge = sprintf("J%02d", j))
  }
  do.call(rbind, rows)
}

fit_locations <- function(dat, clustered) {
  f <- btl(dat, "a", "b", response = "y",
           judge = if (clustered) "judge" else NULL)
  z <- stats::setNames(f$objects$location, f$objects$object)
  list(fit = f, estimate = z)
}

run_scenario <- function(label, workloads, effect = 0) {
  J <- length(workloads)
  est <- se_cr1 <- se_jk <- rep(NA_real_, NREP)
  refused <- nonconv <- 0L
  old <- options(rasch.btl_guard_override = TRUE)
  on.exit(options(old), add = TRUE)
  for (r in seq_len(NREP)) {
    dat <- simulate_panel(workloads, effect)
    full <- tryCatch(fit_locations(dat, TRUE), error = function(e) NULL)
    if (is.null(full)) { refused <- refused + 1L; next }
    if (!isTRUE(full$fit$converged)) { nonconv <- nonconv + 1L; next }
    est[r] <- full$estimate["A"]
    se_cr1[r] <- full$fit$objects$se[match("A", full$fit$objects$object)]
    loo <- matrix(NA_real_, J, 4L,
                  dimnames = list(NULL, LETTERS[1:4]))
    for (j in seq_len(J)) {
      fj <- tryCatch(fit_locations(dat[dat$judge != sprintf("J%02d", j), ],
                                   FALSE), error = function(e) NULL)
      if (!is.null(fj)) loo[j, ] <- fj$estimate[colnames(loo)]
    }
    if (all(is.finite(loo))) {
      lm <- colMeans(loo)
      V <- (J - 1) / J * crossprod(sweep(loo, 2L, lm, "-"))
      se_jk[r] <- sqrt(max(V["A", "A"], 0))
    }
  }
  okc <- is.finite(est) & is.finite(se_cr1) & se_cr1 > 0
  okj <- is.finite(est) & is.finite(se_jk) & se_jk > 0
  crit <- stats::qt(0.975, J - 1L)
  type_or_power <- function(se, ok)
    mean(abs(est[ok] / se[ok]) > crit)
  rows <- rbind(
    sv_row("btl-cluster-jackknife", label, "CR1 location uncertainty",
      sum(okc), bias = mean(est[okc]) - effect, emp_sd = stats::sd(est[okc]),
      mean_se = mean(se_cr1[okc]),
      coverage95 = mean(abs(est[okc] - effect) <= crit * se_cr1[okc]),
      type1 = if (effect == 0) type_or_power(se_cr1, okc) else NA_real_,
      power = if (effect != 0) type_or_power(se_cr1, okc) else NA_real_,
      effect = effect, n_attempted = NREP, n_refused = refused,
      n_nonconv = nonconv,
      notes = "current judge-clustered sandwich with CR1 scale and t(J-1)"),
    sv_row("btl-cluster-jackknife", label, "delete-one-judge uncertainty",
      sum(okj), bias = mean(est[okj]) - effect, emp_sd = stats::sd(est[okj]),
      mean_se = mean(se_jk[okj]),
      coverage95 = mean(abs(est[okj] - effect) <= crit * se_jk[okj]),
      type1 = if (effect == 0) type_or_power(se_jk, okj) else NA_real_,
      power = if (effect != 0) type_or_power(se_jk, okj) else NA_real_,
      effect = effect, n_attempted = NREP,
      n_refused = max(0L, NREP - sum(okj) - nonconv),
      n_nonconv = nonconv,
      notes = "delete-one-judge refits; t(J-1); comparison candidate only"))
  rows
}

set.seed(8.24e7)
rows <- rbind(
  run_scenario("10 judges; balanced", rep(48L, 10L)),
  run_scenario("20 judges; one carries 20 per cent", c(160L, rep(34L, 19L))),
  run_scenario("12 judges; concentrated below public guard",
               c(100L, rep(20L, 11L))),
  run_scenario("20 judges; one carries 20 per cent; effect 0.4",
               c(160L, rep(34L, 19L)), effect = 0.4))
sv_write(rows, "btl-cluster-jackknife")
