#!/usr/bin/env Rscript
# Direct finite-simulation calibration of the maximum-statistic adjustment.
# The observed vector and B null rows are iid under the global null. This
# isolates the standardisation algorithm from model fitting and makes the
# exchangeability requirement testable at high precision in seconds.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "20000"))
B <- as.integer(Sys.getenv("SV_B", "99"))
K <- as.integer(Sys.getenv("SV_FAMILY", "10"))
if (NREP < 1000L || B < 30L || K < 2L)
  stop("use at least 1,000 experiments, 30 null draws and two statistics")

one <- function(r, mode, side) {
  set.seed(930000L + r + 100000L * match(mode, c("centred", "studentised")) +
             10000L * (side == "two"))
  x <- matrix(stats::rnorm((B + 1L) * K), B + 1L, K)
  z <- .boot_maxt(x[1L, ], x[-1L, , drop = FALSE], side,
                  min_success = .fit_min_boot_success(B), mode = mode)
  any(z$p_adj <= .05, na.rm = TRUE)
}

rows <- list()
for (mode in c("centred", "studentised")) for (side in c("upper", "two")) {
  rejected <- vapply(seq_len(NREP), one, logical(1), mode = mode, side = side)
  rate <- mean(rejected)
  rows[[length(rows) + 1L]] <- sv_row(
    "maxT exchangeability", sprintf("%s, %s-sided, B=%d, family=%d",
      mode, if (side == "two") "two" else "upper", B, K),
    "global-null familywise error", n_reps = NREP, n_attempted = NREP,
    n_refused = 0L, n_nonconv = 0L, n_error = 0L, familywise = rate,
    notes = paste("observed vector and null rows iid normal; each null row",
                  "standardised against the other B-1 rows"))
}

out <- do.call(rbind, rows)
sv_write(out, "maxt-exchangeability")
print(out[, c("scenario", "n_reps", "familywise", "mc_se_familywise")],
      row.names = FALSE)
