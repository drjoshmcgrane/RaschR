# STUDY: btl-efrm-bias-sweep
#
# Does the small negative bias of the linked BTL-EFRM log set-unit estimate
# disappear as the within-set calibrations become more precise? Conditional
# fits are sufficient because this study concerns the reconciled point
# estimate, not the deliberately withheld conditional inference.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "btl-efrm-bias-sweep"

one_cell <- function(reps, R, seed0) {
  one <- function(r) {
    d <- simulate_btl_efrm(
      n_objects_per_set = 6, n_sets = 2, n_panels = 2,
      n_judges_per_panel = 6, reps_within = reps, reps_cross = reps,
      set_units = c(1, 1.4), set_origins = c(0, 0.5),
      seed = seed0 + r)
    tr <- attr(d, "truth")
    f <- tryCatch(btl_efrm(
      d, "object_a", "object_b", "winner", "judge", "panel",
      object_sets = tr$object_sets, se_method = "conditional"),
      error = function(e) NULL)
    if (is.null(f)) return(c(alpha = NA, kappa = NA, refused = 1,
                              nonconv = 0))
    if (!isTRUE(f$converged))
      return(c(alpha = NA, kappa = NA, refused = 0, nonconv = 1))
    c(alpha = log(f$alpha_table$alpha[f$alpha_table$set == "set2"]),
      kappa = f$kappa_table$kappa[f$kappa_table$set == "set2"],
      refused = 0, nonconv = 0)
  }
  cores <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  z <- if (cores > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(R), one, mc.cores = cores,
                       mc.set.seed = FALSE, mc.preschedule = FALSE)
  else lapply(seq_len(R), one)
  z <- do.call(rbind, z)
  ok <- is.finite(z[, "alpha"])
  refused <- sum(z[, "refused"]); nonconv <- sum(z[, "nonconv"])
  rows <- rbind(
    sv_row(STUDY, sprintf("%d repetitions per pair", reps),
           "log alpha[set2] bias", sum(ok),
           bias = mean(z[ok, "alpha"]) - log(1.4),
           emp_sd = sd(z[ok, "alpha"]), n_attempted = R,
           n_refused = refused, n_nonconv = nonconv),
    sv_row(STUDY, sprintf("%d repetitions per pair", reps),
           "kappa[set2] bias", sum(ok),
           bias = mean(z[ok, "kappa"]) - 0.5,
           emp_sd = sd(z[ok, "kappa"]), n_attempted = R,
           n_refused = refused, n_nonconv = nonconv))
  cat(sprintf("[%s] reps=%d: n=%d alpha bias %+.4f (MCSE %.4f)\n",
              format(Sys.time(), "%H:%M"), reps, sum(ok), rows$bias[1L],
              rows$emp_sd[1L] / sqrt(sum(ok))))
  rows
}

rows <- do.call(rbind, list(
  one_cell(10L, 500L, 510000L),
  one_cell(20L, 500L, 520000L),
  one_cell(50L, 500L, 530000L),
  one_cell(100L, 500L, 540000L)
))
sv_write(rows, STUDY)
