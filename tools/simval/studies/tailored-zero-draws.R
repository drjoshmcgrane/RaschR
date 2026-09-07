# Check the bootstrap's no-tailoring boundary. This is an algorithm/accounting
# check, not an operating-characteristic study of Type I error or coverage.
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "tailored-zero-draws"
B <- 200L
rows <- lapply(c(35L, 150L), function(N) {
  d <- simulate_rasch(N, 4L, seed = 972000L + N)
  fit <- rasch(d, id = "id")
  observed <- fit$moments$E[!is.na(fit$X)]
  chance <- min(observed[is.finite(observed)]) + 1e-6
  stopifnot(chance > 0, chance < 1)
  set.seed(973000L + N)
  statuses <- character(B)
  zero <- logical(B)
  for (b in seq_len(B)) {
    sampled <- .tailored_boot_rows(fit$person$id)
    fb <- tryCatch(rasch(fit$X[sampled$rows, , drop = FALSE],
                        id = sampled$id), error = identity)
    if (inherits(fb, "error")) { statuses[b] <- "error"; next }
    if (!isTRUE(fb$est$converged)) {
      statuses[b] <- "nonconverged"; next
    }
    shift <- suppressWarnings(.tailored_boot_refit(
      fb, chance, NULL, fit$items$item, fit$m))
    statuses[b] <- .fit_boot_status(shift)
    if (statuses[b] == "ok") zero[b] <- all(shift == 0)
  }
  stopifnot(any(zero), sum(table(statuses)) == B)
  sv_row(STUDY, paste0("N=", N, "; four items; chance just above minimum fitted probability"),
    "fraction of attempted draws retained as exact zero changes", B,
    n_attempted = B, n_refused = 0, n_nonconv = sum(statuses == "nonconverged"),
    n_error = sum(statuses == "error"), n_withheld = 0,
    n_metric_unavailable = 0, n_boot_attempted = B,
    n_boot_used = sum(statuses == "ok"),
    n_boot_nonconv = sum(statuses == "nonconverged"),
    n_boot_errors = sum(statuses == "error"),
    notes = sprintf("zero draws=%d/%d; fraction=%.4f; MCSE=%.4f; remaining usable draws have nonzero shifts; not a null-calibration claim",
      sum(zero), B, mean(zero), mc_se_prop(mean(zero), B)))
})
sv_write(do.call(rbind, rows), STUDY)
