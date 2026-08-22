# Paired local-power comparison under exact balance and homoskedasticity.
# Effects shrink with sample size so power does not saturate, allowing the
# finite-sample efficiency of classical and HC3 uniform DIF tests to be
# compared throughout the sample-size range.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "5000"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "1")))
NITEM <- 8L

one_replicate <- function(n_cell, levels_group, effect) {
  group <- factor(rep(seq_len(levels_group), each = 4L * n_cell))
  ci <- factor(rep(rep(seq_len(4L), each = n_cell), levels_group))
  ci_mean <- c(-0.45, -0.15, 0.15, 0.45)[ci]
  Y <- matrix(rnorm(length(group) * NITEM), nrow = length(group),
              ncol = NITEM) + ci_mean
  Y[group == levels_group, 3L] <-
    Y[group == levels_group, 3L] + effect

  pc <- ph <- matrix(NA_real_, NITEM, 2L,
    dimnames = list(NULL, c("group", "group:ci")))
  for (i in seq_len(NITEM)) {
    d <- data.frame(z = Y[, i], group = group, ci = ci)
    fc <- .dif_type2(d, c("group", "ci", "group:ci"),
                     variance = "classical")
    fh <- .dif_type2(d, c("group", "ci", "group:ci"),
                     variance = "hc3", robust_terms = "group")
    pc[i, ] <- fc$p[match(colnames(pc), fc$term)]
    ph[i, ] <- fh$p[match(colnames(ph), fh$term)]
  }
  ac <- stats::p.adjust(as.vector(pc), "holm")[3L] < 0.05
  ah <- stats::p.adjust(as.vector(ph), "holm")[3L] < 0.05
  c(classical = ac, hybrid = ah, disagreement = ac != ah)
}

run_condition <- function(levels_group, n_cell, effect) {
  seeds <- sample.int(.Machine$integer.max, NREP)
  z <- do.call(rbind, parallel::mclapply(seq_len(NREP), function(r) {
    set.seed(seeds[r])
    one_replicate(n_cell, levels_group, effect)
  }, mc.cores = NCORE, mc.preschedule = TRUE))
  scenario <- sprintf("%d group levels; %d per group-by-interval cell; N=%d",
                      levels_group, n_cell,
                      levels_group * 4L * n_cell)
  mk <- function(quantity, value, note) sv_row(
    "dif-hc3-homoskedastic-local-power", scenario, quantity, NREP,
    power = value, effect = effect, n_attempted = NREP,
    n_refused = 0L, n_nonconv = 0L,
    notes = paste("balanced cells; iid N(0,1) errors; local alternative;", note))
  rbind(
    mk("classical planted-item power", mean(z[, "classical"]),
       "Holm-adjusted uniform term for I3"),
    mk("hybrid planted-item power", mean(z[, "hybrid"]),
       "Holm-adjusted HC3 uniform term for I3"),
    mk("adjusted-decision disagreement", mean(z[, "disagreement"]),
       "I3 decision differs between methods"))
}

set.seed(8.28e7)
designs <- list(
  c(2L, 10L, 0.50), c(2L, 30L, 0.30),
  c(2L, 75L, 0.20), c(2L, 150L, 0.14),
  c(3L, 10L, 0.50), c(3L, 50L, 0.224))
rows <- do.call(rbind, lapply(designs, function(x)
  run_condition(as.integer(x[1L]), as.integer(x[2L]), x[3L])))
sv_write(rows, "dif-hc3-homoskedastic-local-power")
