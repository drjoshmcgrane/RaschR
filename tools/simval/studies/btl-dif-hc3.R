# Null calibration of the between-judge BTL-DIF test when judge means have
# unequal precision. The residual model is the model tested by btl_dif() after
# its judge-level aggregation, isolated here so the covariance reference can
# be examined without repeatedly recalibrating the same object locations.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- 10000L
scenarios <- list(
  "8 high-variance, 16 low-variance judges" = c(8, 16, 2, 1),
  "16 high-variance, 8 low-variance judges" = c(16, 8, 2, 1),
  "8 judges per group, variance ratio 4" = c(8, 8, 2, 1)
)
rows <- list()
set.seed(4.8e7)
for (nm in names(scenarios)) {
  z <- scenarios[[nm]]; p_classical <- p_hc3 <- numeric(NREP)
  for (r in seq_len(NREP)) {
    d <- data.frame(
      residual = c(stats::rnorm(z[1L], sd = z[3L]),
                   stats::rnorm(z[2L], sd = z[4L])),
      group = factor(rep(c("A", "B"), z[1:2])))
    p_classical[r] <- .dif_type2(
      d, "group", "residual", variance = "classical")$p[1L]
    p_hc3[r] <- .dif_type2(
      d, "group", "residual", variance = "hc3")$p[1L]
  }
  rows[[length(rows) + 1L]] <- sv_row(
    "btl-dif-hc3", nm, "classical equal-variance Type I", NREP,
    type1 = mean(p_classical < 0.05), n_attempted = NREP,
    n_refused = 0L, n_nonconv = 0L,
    notes = "diagnostic comparison; superseded for BTL-DIF between-judge terms")
  rows[[length(rows) + 1L]] <- sv_row(
    "btl-dif-hc3", nm, "HC3 Type I", NREP,
    type1 = mean(p_hc3 < 0.05), n_attempted = NREP,
    n_refused = 0L, n_nonconv = 0L,
    notes = "public BTL-DIF between-judge covariance; residual F denominator")
}
sv_write(do.call(rbind, rows), "btl-dif-hc3")
