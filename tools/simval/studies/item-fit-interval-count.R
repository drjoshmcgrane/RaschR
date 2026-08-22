# Calibration and power of the conventional item-fit ANOVA across class-
# interval counts. Short tests have few attainable person scores, so allowing
# every score to define its own interval can make an item help determine the
# grouping used to test that same item.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "500"))

simulate_fit <- function(L, effect = 0) {
  n <- 600L
  theta <- rnorm(n)
  delta <- seq(if (L <= 8L) -1.4 else -2,
               if (L <= 8L) 1.4 else 2, length.out = L)
  disc <- rep(1, L); disc[3L] <- 1 + effect
  eta <- sweep(outer(theta, delta, "-"), 2L, disc, "*")
  X <- matrix(rbinom(n * L, 1, plogis(eta)), n, L)
  colnames(X) <- paste0("I", seq_len(L))
  rasch(X)
}

run_scenario <- function(label, L, effect = 0) {
  counts <- 2:10
  raw <- fwer <- power <- matrix(NA_real_, NREP, length(counts),
    dimnames = list(NULL, counts))
  refused <- nonconv <- 0L
  for (r in seq_len(NREP)) {
    fit <- tryCatch(simulate_fit(L, effect), error = function(e) NULL)
    if (is.null(fit)) { refused <- refused + 1L; next }
    if (!isTRUE(fit$est$converged)) { nonconv <- nonconv + 1L; next }
    for (ng in counts) {
      ci <- .class_intervals(fit$person$theta, fit$person$extreme, ng)
      ia <- .item_anova(fit$residuals, ci, fit$person$extreme)
      raw[r, as.character(ng)] <- mean(ia$p < 0.05, na.rm = TRUE)
      pa <- stats::p.adjust(ia$p, "holm")
      fwer[r, as.character(ng)] <- any(pa < 0.05, na.rm = TRUE)
      power[r, as.character(ng)] <- pa[3L] < 0.05
    }
  }
  rows <- list()
  for (ng in counts) {
    ok <- is.finite(raw[, as.character(ng)]) &
      is.finite(fwer[, as.character(ng)]) & is.finite(power[, as.character(ng)])
    n <- sum(ok)
    note <- sprintf("requested %d intervals; ties may reduce the realised count", ng)
    rows[[length(rows) + 1L]] <- sv_row(
      "item-fit-interval-count", label, "item-wise rejection", n,
      type1 = if (effect == 0) mean(raw[ok, as.character(ng)]) else NA_real_,
      power = if (effect != 0) mean(power[ok, as.character(ng)]) else NA_real_,
      mc_override = if (effect == 0) list(type1 =
        stats::sd(raw[ok, as.character(ng)]) / sqrt(n)) else list(),
      effect = effect, n_attempted = NREP, n_refused = refused,
      n_nonconv = nonconv, notes = note)
    rows[[length(rows) + 1L]] <- sv_row(
      "item-fit-interval-count", label, "Holm familywise rejection", n,
      familywise = if (effect == 0)
        mean(fwer[ok, as.character(ng)]) else NA_real_,
      power = if (effect != 0)
        mean(fwer[ok, as.character(ng)]) else NA_real_, effect = effect,
      n_attempted = NREP, n_refused = refused, n_nonconv = nonconv,
      notes = note)
  }
  do.call(rbind, rows)
}

set.seed(8.25e7)
rows <- rbind(
  run_scenario("8 items; model true", 8L),
  run_scenario("8 items; discrimination 1.5 on I3", 8L, effect = 0.5),
  run_scenario("15 items; model true", 15L),
  run_scenario("15 items; discrimination 1.5 on I3", 15L, effect = 0.5))
sv_write(rows, "item-fit-interval-count")
