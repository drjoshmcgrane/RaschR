# STUDY: btl-efrm-supported-topup
#
# Null calibration of the BTL-EFRM set-unit inference outside the caution
# band used for panels with fewer than eight effective judges. The model is
# fitted with the public default of 200 judge-bootstrap replicates; outer
# replication is parallelised, so each fit uses one worker.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "btl-efrm-supported-topup"
N_REPS <- suppressWarnings(as.integer(Sys.getenv("SV_REPS", "500")))
CORES <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
if (!is.finite(N_REPS) || N_REPS < 1L)
  stop("SV_REPS must be a positive integer")
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L

one <- function(r) {
  blank <- c(estimate = NA, se = NA, p_raw = NA, p_omnibus_adj = NA,
             p_followup_adj = NA, refused = 0, nonconv = 0)
  d <- simulate_btl_efrm(
    n_objects_per_set = 6, n_sets = 2, n_panels = 2,
    n_judges_per_panel = 12, reps_within = 20, reps_cross = 20,
    panel_units = c(1, 1), set_units = c(1, 1), set_origins = c(0, 0),
    seed = 930000L + r)
  f <- tryCatch(suppressWarnings(btl_efrm(
    d, "object_a", "object_b", "winner", "judge", "panel",
    object_sets = attr(d, "truth")$object_sets,
    se_method = "judge_bootstrap", boot_reps = 200L, workers = 1L,
    seed = 940000L + r)), error = function(e) NULL)
  if (is.null(f)) { blank["refused"] <- 1; return(blank) }
  if (!isTRUE(f$converged)) { blank["nonconv"] <- 1; return(blank) }
  at <- f$alpha_table[f$alpha_table$set == "set2", , drop = FALSE]
  ot <- f$unit_omnibus[f$unit_omnibus$term == "set units (alpha)", ,
                       drop = FALSE]
  c(estimate = log(at$alpha), se = at$se_log_alpha, p_raw = at$p,
    p_omnibus_adj = ot$p_adj, p_followup_adj = at$p_adj,
    refused = 0, nonconv = 0)
}

z <- if (CORES > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(seq_len(N_REPS), one, mc.cores = CORES,
                     mc.set.seed = FALSE, mc.preschedule = FALSE)
} else lapply(seq_len(N_REPS), one)
z <- do.call(rbind, z)
refused <- sum(z[, "refused"]); nonconv <- sum(z[, "nonconv"])
ok <- is.finite(z[, "estimate"]) & is.finite(z[, "se"]) & z[, "se"] > 0
make_rate <- function(field, label) {
  kp <- ok & is.finite(z[, field])
  sv_row(STUDY, "2 panels, 12 judges/panel, 2 sets, 20 repeats/pair",
         label, sum(kp), type1 = mean(z[kp, field] < 0.05),
         n_attempted = N_REPS, n_refused = refused, n_nonconv = nonconv,
         notes = "null set-unit ratio; judge bootstrap B=200")
}
rows <- rbind(
  sv_row(STUDY, "2 panels, 12 judges/panel, 2 sets, 20 repeats/pair",
         "log alpha[set2] bias, SE calibration and coverage", sum(ok),
         bias = mean(z[ok, "estimate"]), emp_sd = stats::sd(z[ok, "estimate"]),
         mean_se = mean(z[ok, "se"]),
         coverage95 = mean(abs(z[ok, "estimate"]) <= 1.96 * z[ok, "se"]),
         n_attempted = N_REPS, n_refused = refused, n_nonconv = nonconv,
         notes = "null set-unit ratio; judge bootstrap B=200"),
  make_rate("p_raw", "raw set-unit Type I"),
  make_rate("p_omnibus_adj", "Holm-adjusted set-unit omnibus Type I"),
  make_rate("p_followup_adj", "Holm-adjusted set-unit follow-up Type I")
)
sv_write(rows, STUDY)
