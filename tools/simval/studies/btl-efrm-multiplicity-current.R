# STUDY: btl-efrm-multiplicity-current
#
# Null calibration of the complete BTL-EFRM decision families after the
# reconciled-panel refit. This is the BTL-EFRM arm of the older combined
# frame-unit-multiplicity study, with term-specific probabilities retained.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "btl-efrm-multiplicity-current"
N_REPS <- suppressWarnings(as.integer(Sys.getenv("SV_REPS", "1000")))
CORES <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
BOOT <- suppressWarnings(as.integer(Sys.getenv("SV_BOOT", "100")))
if (!is.finite(N_REPS) || N_REPS < 1L)
  stop("SV_REPS must be a positive integer")
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L
if (!is.finite(BOOT) || BOOT < 30L)
  stop("SV_BOOT must be an integer of at least 30")

one <- function(r) {
  nm <- c("phi_raw", "alpha_raw", "kappa_raw", "phi_omni_adj",
          "alpha_omni_adj", "kappa_omni_adj", "phi_follow_adj",
          "alpha_follow_adj", "kappa_follow_adj", "alpha_est", "alpha_se")
  blank <- setNames(rep(NA_real_, length(nm)), nm)
  d <- simulate_btl_efrm(
    n_objects_per_set = 6, n_sets = 2, n_panels = 2,
    n_judges_per_panel = 6, reps_within = 20, reps_cross = 20,
    panel_units = c(1, 1), set_units = c(1, 1), set_origins = c(0, 0),
    seed = 810000L + r)
  f <- tryCatch(suppressWarnings(btl_efrm(
    d, "object_a", "object_b", "winner", "judge", "panel",
    object_sets = attr(d, "truth")$object_sets,
    se_method = "judge_bootstrap", boot_reps = BOOT, workers = 1L,
    seed = 910000L + r)), error = function(e) NULL)
  if (is.null(f)) return(c(blank, refused = 1, nonconv = 0))
  if (!isTRUE(f$converged)) return(c(blank, refused = 0, nonconv = 1))
  get_omni <- function(term, column) {
    z <- f$unit_omnibus[f$unit_omnibus$term == term, column]
    if (length(z) == 1L) z else NA_real_
  }
  ph <- f$phi_table[f$phi_table$panel == "panel2", , drop = FALSE]
  al <- f$alpha_table[f$alpha_table$set == "set2", , drop = FALSE]
  ka <- f$kappa_table[f$kappa_table$set == "set2", , drop = FALSE]
  c(phi_raw = get_omni("panel units (phi)", "p"),
    alpha_raw = get_omni("set units (alpha)", "p"),
    kappa_raw = get_omni("set origins (kappa)", "p"),
    phi_omni_adj = get_omni("panel units (phi)", "p_adj"),
    alpha_omni_adj = get_omni("set units (alpha)", "p_adj"),
    kappa_omni_adj = get_omni("set origins (kappa)", "p_adj"),
    phi_follow_adj = ph$p_adj, alpha_follow_adj = al$p_adj,
    kappa_follow_adj = ka$p_adj, alpha_est = log(al$alpha),
    alpha_se = al$se_log_alpha, refused = 0, nonconv = 0)
}

z <- if (CORES > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(seq_len(N_REPS), one, mc.cores = CORES,
                     mc.set.seed = FALSE, mc.preschedule = FALSE)
} else lapply(seq_len(N_REPS), one)
z <- do.call(rbind, z)
refused <- sum(z[, "refused"]); nonconv <- sum(z[, "nonconv"])
base_ok <- rowSums(is.finite(z[, c("phi_raw", "alpha_raw", "kappa_raw",
                                    "phi_omni_adj", "alpha_omni_adj",
                                    "kappa_omni_adj", "phi_follow_adj",
                                    "alpha_follow_adj", "kappa_follow_adj")])) == 9L
scenario <- "2 panels, 6 judges/panel, 2 sets, 20 repeats/pair"
note <- sprintf("judge bootstrap B=%d; panels are in the documented caution band", BOOT)
rate <- function(field, label) {
  ok <- base_ok & is.finite(z[, field])
  sv_row(STUDY, scenario, label, sum(ok),
         type1 = mean(z[ok, field] < 0.05), n_attempted = N_REPS,
         n_refused = refused, n_nonconv = nonconv, notes = note)
}
rows <- rbind(
  sv_row(STUDY, scenario, "probability of any raw omnibus p below .05",
         sum(base_ok), familywise = mean(apply(z[base_ok,
           c("phi_raw", "alpha_raw", "kappa_raw"), drop = FALSE], 1L,
           function(x) any(x < 0.05))), n_attempted = N_REPS,
         n_refused = refused, n_nonconv = nonconv, notes = note),
  sv_row(STUDY, scenario, "Holm familywise error: omnibus decisions",
         sum(base_ok), familywise = mean(apply(z[base_ok,
           c("phi_omni_adj", "alpha_omni_adj", "kappa_omni_adj"),
           drop = FALSE], 1L, function(x) any(x < 0.05))),
         n_attempted = N_REPS, n_refused = refused, n_nonconv = nonconv,
         notes = note),
  sv_row(STUDY, scenario, "Holm familywise error: individual unit follow-ups",
         sum(base_ok), familywise = mean(apply(z[base_ok,
           c("phi_follow_adj", "alpha_follow_adj", "kappa_follow_adj"),
           drop = FALSE], 1L, function(x) any(x < 0.05))),
         n_attempted = N_REPS, n_refused = refused, n_nonconv = nonconv,
         notes = note),
  rate("phi_raw", "raw panel-unit Type I"),
  rate("alpha_raw", "raw set-unit Type I"),
  rate("kappa_raw", "raw set-origin Type I"),
  rate("phi_omni_adj", "Holm-adjusted panel-unit omnibus Type I"),
  rate("alpha_omni_adj", "Holm-adjusted set-unit omnibus Type I"),
  rate("kappa_omni_adj", "Holm-adjusted set-origin omnibus Type I")
)
ok_a <- is.finite(z[, "alpha_est"]) & is.finite(z[, "alpha_se"]) &
  z[, "alpha_se"] > 0
rows <- rbind(rows, sv_row(
  STUDY, scenario, "log alpha[set2] bias, SE calibration and coverage",
  sum(ok_a), bias = mean(z[ok_a, "alpha_est"]),
  emp_sd = stats::sd(z[ok_a, "alpha_est"]),
  mean_se = mean(z[ok_a, "alpha_se"]),
  coverage95 = mean(abs(z[ok_a, "alpha_est"]) <= 1.96 * z[ok_a, "alpha_se"]),
  n_attempted = N_REPS, n_refused = refused, n_nonconv = nonconv,
  notes = note))
sv_write(rows, STUDY)
