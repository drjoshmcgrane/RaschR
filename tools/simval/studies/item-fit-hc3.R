# Sensitivity of the conventional class-interval item-fit ANOVA to HC3
# covariance. This study does not presuppose that the established statistic
# should change: it tests whether a robust version improves null calibration
# without materially reducing detection of an item with a non-Rasch slope.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "1000"))

hc3_item_p <- function(fit) {
  out <- rep(NA_real_, ncol(fit$residuals))
  for (i in seq_along(out)) {
    ci <- if (is.null(fit$ci_item)) fit$person$class_interval else fit$ci_item[[i]]
    sel <- !is.na(fit$residuals[, i]) & !is.na(ci) & !fit$person$extreme
    d <- data.frame(z = fit$residuals[sel, i], interval = factor(ci[sel]))
    keep <- d$interval %in% names(which(table(d$interval) >= 2L))
    d <- droplevels(d[keep, , drop = FALSE])
    if (nlevels(d$interval) < 2L) next
    z <- .dif_type2(d, "interval", variance = "hc3")
    if (!is.null(z)) out[i] <- z$p[z$term == "interval"]
  }
  out
}

simulate_fit <- function(scenario, effect = 0) {
  standard <- scenario == "standard 15-item test"
  long <- scenario == "broad 30-item test"
  ten <- scenario == "10-item test"
  twelve <- scenario == "12-item test"
  is_missing <- scenario == "short test with MAR missingness"
  n <- if (standard) 500L else 600L
  L <- if (standard) 15L else if (long) 30L else if (ten) 10L else if (twelve)
    12L else 8L
  theta <- switch(scenario,
    `short targeted test` = rnorm(n),
    `short mistargeted test` = rnorm(n, -1.2, 1.25),
    `short test with MAR missingness` = rnorm(n),
    `10-item test` = rnorm(n),
    `12-item test` = rnorm(n),
    `standard 15-item test` = rnorm(n),
    `broad 30-item test` = rnorm(n))
  delta <- if (standard || ten || twelve) seq(-2, 2, length.out = L) else if (long)
    seq(-2.5, 2.5, length.out = L) else seq(-1.4, 1.4, length.out = L)
  disc <- rep(1, L); disc[3L] <- 1 + effect
  X <- matrix(rbinom(n * L, 1,
    plogis(sweep(outer(theta, delta, "-"), 2L, disc, "*"))), n, L)
  # Correct the simulation expression to a_i(theta-delta_i): sweep above
  # multiplies columns and therefore has exactly that form.
  if (is_missing) {
    pm <- plogis(-1.4 + 0.7 * theta)
    miss <- matrix(runif(n * L) < rep(pm, L), n, L)
    X[miss] <- NA_integer_
  }
  colnames(X) <- paste0("I", seq_len(L))
  rasch(X)
}

run_scenario <- function(scenario, effect = 0) {
  z <- matrix(NA_real_, NREP, 9L,
    dimnames = list(NULL, c("raw_classical", "raw_hc3", "fwer_classical",
                            "fwer_hc3", "power_classical", "power_hc3",
                            "raw_trait", "fwer_trait", "power_trait")))
  refused <- nonconv <- 0L
  for (r in seq_len(NREP)) {
    fit <- tryCatch(simulate_fit(scenario, effect), error = function(e) NULL)
    if (is.null(fit)) { refused <- refused + 1L; next }
    if (!isTRUE(fit$est$converged)) { nonconv <- nonconv + 1L; next }
    pc <- fit$item_anova$p
    ph <- hc3_item_p(fit)
    pt <- fit$item_trait$p
    z[r, "raw_classical"] <- mean(pc < 0.05, na.rm = TRUE)
    z[r, "raw_hc3"] <- mean(ph < 0.05, na.rm = TRUE)
    ac <- stats::p.adjust(pc, "holm"); ah <- stats::p.adjust(ph, "holm")
    z[r, "fwer_classical"] <- any(ac < 0.05, na.rm = TRUE)
    z[r, "fwer_hc3"] <- any(ah < 0.05, na.rm = TRUE)
    z[r, "raw_trait"] <- mean(pt < 0.05, na.rm = TRUE)
    at <- stats::p.adjust(pt, "holm")
    z[r, "fwer_trait"] <- any(at < 0.05, na.rm = TRUE)
    z[r, "power_trait"] <- at[3L] < 0.05
    z[r, "power_classical"] <- ac[3L] < 0.05
    z[r, "power_hc3"] <- ah[3L] < 0.05
  }
  mk <- function(quantity, col, field) {
    ok <- is.finite(z[, col]); n <- sum(ok)
    args <- list(study = "item-fit-hc3", scenario = scenario,
      quantity = quantity, n_reps = n, n_attempted = NREP,
      n_refused = NREP - n - nonconv, n_nonconv = nonconv, effect = effect,
      notes = paste("metric-specific denominator; unavailable diagnostic rows",
                    "count as refusals. HC3 is a sensitivity candidate;",
                    "conventional RMT ANOVA remains the comparator"))
    args[[field]] <- if (effect != 0 && field %in% c("type1", "familywise"))
      NA_real_ else mean(z[ok, col])
    if (field == "type1" && effect == 0) args$mc_override <- list(
      type1 = stats::sd(z[ok, col]) / sqrt(n))
    do.call(sv_row, args)
  }
  rbind(
    mk("classical item-wise Type I", "raw_classical", "type1"),
    mk("HC3 item-wise Type I", "raw_hc3", "type1"),
    mk("classical Holm FWER", "fwer_classical", "familywise"),
    mk("HC3 Holm FWER", "fwer_hc3", "familywise"),
    mk("item-trait item-wise Type I", "raw_trait", "type1"),
    mk("item-trait Holm FWER", "fwer_trait", "familywise"),
    mk("classical planted-item power", "power_classical", "power"),
    mk("HC3 planted-item power", "power_hc3", "power"),
    mk("item-trait planted-item power", "power_trait", "power"))
}

set.seed(8.23e7)
rows <- rbind(
  run_scenario("short targeted test"),
  run_scenario("short mistargeted test"),
  run_scenario("short test with MAR missingness"),
  run_scenario("10-item test"),
  run_scenario("12-item test"),
  run_scenario("standard 15-item test"),
  run_scenario("broad 30-item test"),
  run_scenario("short targeted test", effect = 0.5))
sv_write(rows, "item-fit-hc3")
