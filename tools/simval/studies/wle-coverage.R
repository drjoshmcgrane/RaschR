# STUDY: wle-coverage
#
# Coverage study for person WLE (Warm, 1989) measures reported by rasch().
# Replaces the single-point WLE spot check (calibration.R block F) with a
# proper per-grid-point coverage study over a theta grid, test lengths, and
# targeting offsets.
#
# Design
# ------
# For each scenario (test length L in {5, 10, 25}; targeting offset in
# {0, +1}):
#   1. Draw ONE large calibration sample (N = 2000 persons ~ N(0, 1.3)) from
#      the FIXED item bank seq(-2, 2, length.out = L) (mean 0), fit rasch()
#      once, and hold the resulting tau_list ("as a user would have it":
#      estimated, not the true, item calibration) FIXED for every replicate
#      in the scenario.
#   2. Build the Warm WLE score-table (person_wle()) from that fixed,
#      estimated calibration -- also fixed across replicates.
#   3. For each nominal grid point g in seq(-3, 3, by = 0.75), simulate many
#      replicate response vectors from the SAME fixed item bank AT THE
#      KNOWN TRUE theta = g + offset, read off each replicate's raw score
#      and its WLE + reported SE from the fixed score-table, and compare
#      the estimate to the known truth theta.
#
# Why the offset is applied to the person grid, not to the item locations
# -------------------------------------------------------------------------
# The task brief frames "targeting offset" as shifting the item locations
# off the person distribution. That framing is analytically DEGENERATE
# under this package's identification convention: rasch()/pcml() explicitly
# recentre every calibration so the mean item location is exactly zero
# (R/estimation.R, "recentre so the mean item location is zero"; verified
# empirically to machine precision -- a calibration on items with true mean
# location = +1 comes back with mean estimated location exactly 0). A
# uniform additive shift of the TRUE item locations is therefore fully and
# exactly absorbed into that arbitrary recentring constant: scoring against
# a fixed external theta grid then shows an exact mechanical bias of
# -offset, not a genuine WLE miscalibration (confirmed: with dv_true =
# base + 1 the observed bias was -0.97..-1.04 at every grid point, i.e.
# the offset itself, not a function of theta or test length). That is an
# identification-indeterminacy artifact -- an unavoidable fact about
# unanchored Rasch calibration, not a property of the WLE estimator -- and
# reporting it as "WLE coverage" would conflate the two.
#
# The substantively meaningful version of "targeting offset" -- the person
# population's true ability sitting away from the item bank's difficulty
# range -- is captured, WITHOUT the origin confound, by keeping the item
# bank fixed (mean 0, identical across offset conditions, so the estimated
# calibration's origin exactly matches the true origin throughout) and
# shifting the TRUE theta at which replicate persons are probed. This is
# the same relative item-person mismatch the brief describes; it is immune
# to the recentring wash-out because only the item side of the model is
# subject to the identification constraint. Both the degenerate item-shift
# behaviour and this corrected design are exercised below: scenario "diag"
# reports the former as a documented calibration/identification check, and
# the main "L=.. offset=.." scenarios use the latter for the coverage
# analysis proper.
#
# Held fixed across replicates: the item bank AND the calibrations used to
# score persons -- ONE estimated calibration per test length, shared across
# both offset conditions so offset comparisons are not confounded with
# calibration noise. Every grid point is scored twice from the same
# replicate draws: once against the TRUE item locations (the conditional
# performance of the Warm WLE itself) and once against the estimated
# calibration (the user-facing composite of WLE and calibration error).
# Redrawn every replicate: the person's responses at the known true
# theta = grid + offset.
#
# Pass bands (informal): SE ratio in ~[0.85, 1.20], coverage95 close to the
# nominal 0.95 (allow ~+/-3 MC SE), no systematic bias trend with |theta|
# once inside the item range; degradation is expected once true theta runs
# past the item bank's [-2, 2] span, most visibly for offset = +1.
#
# Run from the package root:
#   Rscript tools/simval/studies/wle-coverage.R        # all scenarios
#   Rscript tools/simval/studies/wle-coverage.R 1       # chunk 1 (L in 5,10)
#   Rscript tools/simval/studies/wle-coverage.R 2       # chunk 2 (L = 25)

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

STUDY <- "wle-coverage"
grid  <- seq(-3, 3, by = 0.75)
Ncal  <- 2000                      # large calibration sample, as a user would have it

scenarios <- expand.grid(L = c(5, 10, 25), offset = c(0, 1))
scenarios <- scenarios[order(scenarios$L, scenarios$offset), ]

args  <- commandArgs(trailingOnly = TRUE)
chunk <- if (length(args) >= 1) suppressWarnings(as.integer(args[1])) else NA
if (!is.na(chunk)) {
  # chunk 1: L in {5, 10}; chunk 2: L = 25 -- roughly balances the (tiny)
  # per-scenario calibration-fit cost across two concurrent processes.
  keep <- if (chunk == 1) scenarios$L %in% c(5, 10) else scenarios$L == 25
  scenarios <- scenarios[keep, ]
}

rows <- list()
.cal_cache <- list()      # one calibration per test length, shared across offsets

t_start <- Sys.time()

# --- diagnostic: confirm the item-shift wash-out analytically, once -------
# (n_reps = 1: this is a single deterministic calibration fit, not a Monte
# Carlo claim; it documents WHY the offset is applied to the person grid
# below rather than to the item locations, per the design note above.)
{
  L <- 10; shift <- 1
  dv_shift <- seq(-2, 2, length.out = L) + shift
  set.seed(41000L)
  th_d <- rnorm(Ncal, 0, 1.3)
  Xd <- matrix(rbinom(Ncal * L, 1, plogis(outer(th_d, dv_shift, "-"))), Ncal, L)
  colnames(Xd) <- sprintf("I%02d", seq_len(L))
  fd <- rasch(Xd)
  mean_est_loc <- mean(vapply(fd$tau_list, mean, 0))
  rows[[length(rows) + 1]] <- sv_row(
    STUDY, "diag: item-shift identification", "mean estimated item location", n_reps = 1L,
    bias = mean_est_loc - mean(dv_shift),
    notes = sprintf(paste(
      "true item locations shifted by +%.1f (mean true loc = %.1f);",
      "rasch() recentres to mean estimated loc = %.6f (identification convention,",
      "R/estimation.R: 'recentre so the mean item location is zero'). A uniform",
      "item-location shift is therefore fully absorbed into the arbitrary origin",
      "-- confirms the item-shift framing of 'targeting offset' is degenerate for",
      "a coverage claim; the person-grid-shift design below avoids this."),
      shift, mean(dv_shift), mean_est_loc))
}

for (si in seq_len(nrow(scenarios))) {
  L <- scenarios$L[si]; offset <- scenarios$offset[si]
  principal <- (L == 10)                      # principal test length
  reps <- if (principal) 5000L else 1500L     # >=1000 principal, >=400 elsewhere

  # item bank is FIXED (mean 0) across offset conditions -- see design note
  # above for why the offset is applied to the person grid, not here.
  dv <- seq(-2, 2, length.out = L)

  # --- one large calibration sample PER TEST LENGTH, fit once, held fixed
  # across BOTH offset conditions (a per-offset calibration draw would
  # confound offset comparisons with calibration noise) --------------------
  cal_key <- as.character(L)
  if (is.null(.cal_cache[[cal_key]])) {
    seed_cal <- 42000L + L * 10L
    set.seed(seed_cal)
    th_cal <- rnorm(Ncal, 0, 1.3)
    Xcal <- matrix(rbinom(Ncal * L, 1, plogis(outer(th_cal, dv, "-"))), Ncal, L)
    colnames(Xcal) <- sprintf("I%02d", seq_len(L))
    fc <- rasch(Xcal)
    .cal_cache[[cal_key]] <- list(
      conv = isTRUE(fc$est$converged),
      st_hat = person_wle(fc$tau_list),
      st_true = person_wle(setNames(as.list(dv), sprintf("I%02d", seq_len(L)))))
  }
  cal <- .cal_cache[[cal_key]]
  conv <- cal$conv
  st <- cal$st_hat            # WLE score table from the ESTIMATED calibration
  st_true <- cal$st_true      # WLE score table from the TRUE item locations

  scen_lab <- sprintf("L=%d offset=%+d", L, offset)
  rows[[length(rows) + 1]] <- sv_row(
    STUDY, scen_lab, "calibration fit", n_reps = 1L,
    n_attempted = 1L, n_nonconv = if (conv) 0L else 1L,
    notes = sprintf("one calibration fit, N=%d persons ~ N(0,1.3), fixed item bank [-2,2]; converged=%s",
                     Ncal, conv))

  for (gi in seq_along(grid)) {
    g <- grid[gi]; theta_true <- g + offset      # true theta = nominal grid + targeting offset
    seed_rep <- 50000L + L * 1000L + offset * 100L + gi
    set.seed(seed_rep)
    p <- plogis(theta_true - dv)                        # true response probs, fixed item bank
    Xr <- matrix(rbinom(reps * L, 1, rep(p, each = reps)), reps, L)
    raw <- rowSums(Xr)
    key <- as.character(raw)
    ex_lo <- mean(raw == 0L); ex_hi <- mean(raw == L)
    extreme_rate <- ex_lo + ex_hi

    # two arms from the same replicate draws: TRUE calibration isolates the
    # conditional performance of the WLE itself; ESTIMATED calibration is
    # the user-facing composite (WLE + one fixed calibration's error).
    for (arm in c("true-calib", "est-calib")) {
      tab <- if (arm == "true-calib") st_true else st
      est <- unname(tab$theta[key]); se <- unname(tab$se[key])
      ok <- is.finite(est) & is.finite(se)               # Warm WLE: finite at every score

      rows[[length(rows) + 1]] <- sv_row(
        STUDY, scen_lab, sprintf("theta_true=%+.2f [%s]", theta_true, arm),
        n_reps = sum(ok), n_attempted = reps, n_refused = reps - sum(ok),
        bias = mean(est[ok] - theta_true), emp_sd = sd(est[ok]),
        mean_se = mean(se[ok]),
        coverage95 = sv_coverage(est[ok], se[ok], theta_true),
        notes = sprintf(
          "nominal_grid=%+.2f; extreme_rate=%.3f (raw=0: %.3f, raw=%d: %.3f); reps=%d; principal=%s; fixed item bank [-2,2]; one calibration per L shared across offsets",
          g, extreme_rate, ex_lo, L, ex_hi, reps, principal))
    }
  }
  cat(sprintf("scenario %s done (%.1fs elapsed)\n", scen_lab,
              as.numeric(Sys.time() - t_start, units = "secs")))
}

res <- do.call(rbind, rows)
suffix <- if (!is.na(chunk)) paste0("-chunk", chunk) else ""
sv_write(res, paste0("wle-coverage", suffix))

cat(sprintf("total wall time: %.1fs\n", as.numeric(Sys.time() - t_start, units = "secs")))
