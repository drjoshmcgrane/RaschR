# STUDY: frame-unit-multiplicity
#
# Null calibration of the decision rule used for EFRM and BTL-EFRM unit
# families. The earlier studies report the marginal size of each raw omnibus
# test. This study records the probability that any test in the reported
# Holm-adjusted omnibus family rejects, and does the same for the separate
# family of individual unit follow-ups.
# Run from the package root. Set SV_CORES on Unix-like systems to parallelise.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "frame-unit-multiplicity"
ALPHA <- 0.05
R_EFRM <- suppressWarnings(as.integer(Sys.getenv("SV_EFRM_REPS", "1000")))
R_BTL <- suppressWarnings(as.integer(Sys.getenv("SV_BTL_REPS", "1000")))
if (!is.finite(R_EFRM) || R_EFRM < 1L)
  stop("SV_EFRM_REPS must be a positive integer")
if (!is.finite(R_BTL) || R_BTL < 1L)
  stop("SV_BTL_REPS must be a positive integer")
EFRM_BOOT <- suppressWarnings(as.integer(Sys.getenv("SV_EFRM_BOOT", "80")))
if (!is.finite(EFRM_BOOT) || EFRM_BOOT < 30L)
  stop("SV_EFRM_BOOT must be an integer of at least 30")
BTL_BOOT <- suppressWarnings(as.integer(Sys.getenv("SV_BTL_BOOT", "100")))
if (!is.finite(BTL_BOOT) || BTL_BOOT < 30L)
  stop("SV_BTL_BOOT must be an integer of at least 30")
CORES <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L

check_adjustment <- function(tab) {
  if (is.null(tab) || !nrow(tab)) return(FALSE)
  ok <- is.finite(tab$p)
  expected <- rep(NA_real_, nrow(tab))
  expected[ok] <- stats::p.adjust(tab$p[ok], method = "holm")
  isTRUE(all.equal(tab$p_adj, expected, tolerance = 1e-14,
                   check.attributes = FALSE))
}

family_result <- function(omnibus, followups) {
  if (!check_adjustment(omnibus) || !check_adjustment(followups))
    stop("reported adjusted probabilities do not match Holm adjustment")
  c(raw_omnibus_any = any(omnibus$p < ALPHA, na.rm = TRUE),
    adjusted_omnibus_any = any(omnibus$p_adj < ALPHA, na.rm = TRUE),
    raw_followup_any = any(followups$p < ALPHA, na.rm = TRUE),
    adjusted_followup_any = any(followups$p_adj < ALPHA, na.rm = TRUE))
}

one_efrm <- function(r) {
  blank <- c(raw_omnibus_any = NA, adjusted_omnibus_any = NA,
             raw_followup_any = NA, adjusted_followup_any = NA)
  d <- simulate_efrm(100, 4, n_sets = 2, n_groups = 2,
                     set_unit_ratio = 1, group_unit_ratio = 1,
                     seed = 610000L + r)
  f <- tryCatch(rasch_efrm(
    d, item_sets = attr(d, "truth")$item_sets, groups = "group",
    se_method = "hybrid", boot_reps = EFRM_BOOT, workers = 1L,
    seed = 710000L + r),
    error = function(e) NULL)
  if (is.null(f)) return(c(blank, refused = 1, nonconv = 0))
  if (!isTRUE(f$est$converged))
    return(c(blank, refused = 0, nonconv = 1))
  c(family_result(f$efrm_vs_rasch$unit_omnibus,
                  f$efrm_vs_rasch$unit_tests),
    refused = 0, nonconv = 0)
}

one_btl_efrm <- function(r) {
  blank <- c(raw_omnibus_any = NA, adjusted_omnibus_any = NA,
             raw_followup_any = NA, adjusted_followup_any = NA)
  d <- simulate_btl_efrm(
    n_objects_per_set = 6, n_sets = 2, n_panels = 2,
    n_judges_per_panel = 6, reps_within = 20, reps_cross = 20,
    panel_units = c(1, 1), set_units = c(1, 1),
    set_origins = c(0, 0), seed = 810000L + r)
  f <- tryCatch(btl_efrm(
    d, "object_a", "object_b", winner = "winner", judge = "judge",
    panels = "panel", object_sets = attr(d, "truth")$object_sets,
    se_method = "judge_bootstrap", boot_reps = BTL_BOOT,
    workers = 1L, seed = 910000L + r),
    error = function(e) NULL)
  if (is.null(f)) return(c(blank, refused = 1, nonconv = 0))
  if (!isTRUE(f$converged))
    return(c(blank, refused = 0, nonconv = 1))
  followups <- rbind(
    data.frame(p = f$phi_table$p, p_adj = f$phi_table$p_adj),
    data.frame(p = f$alpha_table$p, p_adj = f$alpha_table$p_adj),
    data.frame(p = f$kappa_table$p, p_adj = f$kappa_table$p_adj))
  c(family_result(f$unit_omnibus, followups),
    refused = 0, nonconv = 0)
}

run_reps <- function(fun, n_reps) {
  z <- if (CORES > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(n_reps), fun, mc.cores = CORES,
                       mc.set.seed = FALSE)
  else lapply(seq_len(n_reps), fun)
  bad <- vapply(z, inherits, logical(1), "try-error")
  if (any(bad)) stop("worker failure: ", as.character(z[[which(bad)[1L]]]))
  do.call(rbind, z)
}

summarise_family <- function(model, z, n_attempted, note) {
  refused <- sum(z[, "refused"])
  nonconv <- sum(z[, "nonconv"])
  ok <- is.finite(z[, "adjusted_omnibus_any"]) &
    is.finite(z[, "adjusted_followup_any"])
  make <- function(field, quantity) sv_row(
    STUDY, model, quantity, n_reps = sum(ok),
    familywise = mean(z[ok, field]),
    n_attempted = n_attempted, n_refused = refused, n_nonconv = nonconv,
    notes = note)
  rbind(
    make("raw_omnibus_any", "probability of any raw omnibus p below .05"),
    make("adjusted_omnibus_any", "Holm familywise error: omnibus decisions"),
    make("raw_followup_any", "probability of any raw follow-up p below .05"),
    make("adjusted_followup_any",
         "Holm familywise error: individual unit follow-ups"))
}

cat(sprintf("EFRM null family: %d replicates, B=%d, on %d core(s)\n",
            R_EFRM, EFRM_BOOT, CORES))
ze <- run_reps(one_efrm, R_EFRM)
rows_e <- summarise_family(
  "EFRM: 2 groups x 2 sets, 100 persons/group, 4 items/set", ze, R_EFRM,
  sprintf("hybrid covariance B=%d; Holm is applied across two omnibus tests and separately across all individual unit contrasts",
          EFRM_BOOT))
sv_write(rows_e, STUDY)

cat(sprintf("BTL-EFRM null family: %d replicates, B=%d, on %d core(s)\n",
            R_BTL, BTL_BOOT, CORES))
zb <- run_reps(one_btl_efrm, R_BTL)
rows_b <- summarise_family(
  "BTL-EFRM: 2 panels x 2 sets, 12 judges, 6 objects/set", zb, R_BTL,
  sprintf("judge bootstrap B=%d; Holm is applied across three omnibus tests and separately across all panel-unit, set-unit and set-origin contrasts",
          BTL_BOOT))
sv_write(rbind(rows_e, rows_b), STUDY)
