# STUDY: btl-efrm-current
#
# Current-code calibration of BTL-EFRM bootstrap inference. The judge
# bootstrap is assessed under the null and two planted alternatives. A second
# null cell checks the independent-outcome parametric bootstrap and its
# normal/chi-square references.
# Run from the package root; set SV_CORES on Unix-like systems to parallelise.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "btl-efrm-current"
BOOT <- 100L

one <- function(r, panel_ratio, set_ratio, origin, seed0, se_method) {
  blank <- c(phi = NA, phi_se = NA, alpha = NA, alpha_se = NA,
             kappa = NA, kappa_se = NA, p_phi = NA, p_alpha = NA,
             p_kappa = NA)
  pu <- c(panel_ratio^(-0.5), panel_ratio^0.5)
  d <- simulate_btl_efrm(
    n_objects_per_set = 6, n_sets = 2, n_panels = 2,
    n_judges_per_panel = 6, reps_within = 20, reps_cross = 20,
    panel_units = pu, set_units = c(1, set_ratio),
    set_origins = c(0, origin), seed = seed0 + r
  )
  tr <- attr(d, "truth")
  set.seed(seed0 + 100000L + r)
  f <- tryCatch(btl_efrm(
    d, "object_a", "object_b", "winner", "judge", "panel",
    object_sets = tr$object_sets, se_method = se_method,
    boot_reps = BOOT, workers = 1L), error = function(e) NULL)
  if (is.null(f)) return(c(blank, refused = 1, nonconv = 0))
  if (!isTRUE(f$converged))
    return(c(blank, refused = 0, nonconv = 1))
  getp <- function(term) {
    z <- f$unit_omnibus$p[f$unit_omnibus$term == term]
    if (length(z)) z else NA_real_
  }
  pt <- f$phi_table[f$phi_table$panel == "panel2", ]
  at <- f$alpha_table[f$alpha_table$set == "set2", ]
  kt <- f$kappa_table[f$kappa_table$set == "set2", ]
  c(phi = log(pt$phi), phi_se = pt$se_log_phi,
    alpha = log(at$alpha), alpha_se = at$se_log_alpha,
    kappa = kt$kappa, kappa_se = kt$se_kappa,
    p_phi = getp("panel units (phi)"),
    p_alpha = getp("set units (alpha)"),
    p_kappa = getp("set origins (kappa)"),
    refused = 0, nonconv = 0)
}

cell <- function(label, panel_ratio, set_ratio, origin, R, seed0,
                 se_method = "judge_bootstrap") {
  cores <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  z <- if (cores > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(R), one, panel_ratio = panel_ratio,
                       set_ratio = set_ratio, origin = origin, seed0 = seed0,
                       se_method = se_method,
                       mc.cores = cores, mc.set.seed = FALSE)
  else lapply(seq_len(R), one, panel_ratio = panel_ratio,
              set_ratio = set_ratio, origin = origin, seed0 = seed0,
              se_method = se_method)
  z <- do.call(rbind, z)
  refused <- sum(z[, "refused"]); nonconv <- sum(z[, "nonconv"])
  rows <- list()
  add_parameter <- function(name, se_name, truth, p_name, effect) {
    ok <- is.finite(z[, name]) & is.finite(z[, se_name]) & z[, se_name] > 0
    rows[[length(rows) + 1L]] <<- sv_row(
      STUDY, label,
      sprintf("%s bias, SE calibration, coverage and omnibus rejection", name),
      sum(ok), effect = effect, bias = mean(z[ok, name]) - truth,
      emp_sd = sd(z[ok, name]), mean_se = mean(z[ok, se_name]),
      coverage95 = mean(abs(z[ok, name] - truth) <= 1.96 * z[ok, se_name]),
      type1 = if (abs(truth) < 1e-12) mean(z[ok, p_name] < 0.05) else NA_real_,
      power = if (abs(truth) >= 1e-12) mean(z[ok, p_name] < 0.05) else NA_real_,
      n_attempted = R, n_refused = refused, n_nonconv = nonconv,
      notes = if (se_method == "judge_bootstrap")
        sprintf("judge bootstrap, B=%d; judges are resampled within panel", BOOT)
      else sprintf(paste("parametric bootstrap, B=%d; outcomes are drawn",
                         "independently conditional on the fitted design"), BOOT))
  }
  add_parameter("phi", "phi_se", log(panel_ratio) / 2, "p_phi", panel_ratio)
  add_parameter("alpha", "alpha_se", log(set_ratio), "p_alpha", set_ratio)
  add_parameter("kappa", "kappa_se", origin, "p_kappa", origin)
  cat(sprintf("[%s] %s: analysed %d/%d, refused %d, nonconverged %d\n",
              format(Sys.time(), "%H:%M"), label,
              R - refused - nonconv, R, refused, nonconv))
  do.call(rbind, rows)
}

rows <- rbind(
  cell("equal panel and set units; zero origin", 1, 1, 0,
       R = 300L, seed0 = 410000L),
  cell("panel ratio 1.3; set ratio 1.3; origin 0.5", 1.3, 1.3, 0.5,
       R = 180L, seed0 = 420000L),
  cell("panel ratio 1.6; set ratio 1.6; origin 0.8", 1.6, 1.6, 0.8,
       R = 120L, seed0 = 430000L),
  cell("equal units and origin; independent-outcome bootstrap", 1, 1, 0,
       R = 300L, seed0 = 440000L, se_method = "bootstrap")
)
sv_write(rows, STUDY)
