suppressWarnings(pkgload::load_all(".", quiet = TRUE))
t_start <- Sys.time()
alpha_level <- 0.05

extract_omnibus_p <- function(fit, term) {
  if (is.null(fit$unit_omnibus)) return(NA_real_)
  row <- fit$unit_omnibus[fit$unit_omnibus$term == term, ]
  if (!nrow(row)) return(NA_real_)
  row$p[1]
}

run_batch <- function(n_rep, panel_units, set_units, set_origins, se_method,
                       boot_reps = NA, seed0) {
  p_phi <- p_alpha <- p_kappa <- numeric(n_rep)
  n_fail <- 0
  fr_resid <- list()
  for (r in seq_len(n_rep)) {
    seed <- seed0 + r
    d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                            n_judges_per_panel = 6, reps_within = 20, reps_cross = 20,
                            panel_units = panel_units, set_units = set_units,
                            set_origins = set_origins, object_sd = 1, seed = seed)
    tr <- attr(d, "truth")
    args <- list(data = d, object_a = "object_a", object_b = "object_b",
                 winner = "winner", judge = "judge", panels = "panel",
                 object_sets = tr$object_sets, se_method = se_method)
    if (!is.na(boot_reps)) args$boot_reps <- boot_reps
    fit <- tryCatch(do.call(btl_efrm, args), error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$converged)) { n_fail <- n_fail + 1
      p_phi[r] <- p_alpha[r] <- p_kappa[r] <- NA; next }
    p_phi[r]   <- extract_omnibus_p(fit, "panel units (phi)")
    p_alpha[r] <- extract_omnibus_p(fit, "set units (alpha)")
    p_kappa[r] <- extract_omnibus_p(fit, "set origins (kappa)")
    if (!is.null(fit$frames)) fr_resid[[length(fr_resid) + 1]] <- fit$frames$fit_resid
  }
  list(p_phi = p_phi, p_alpha = p_alpha, p_kappa = p_kappa, n_fail = n_fail,
       n_rep = n_rep, fit_resid = unlist(fr_resid))
}

rate_row <- function(p, label) {
  ok <- is.finite(p)
  rate <- mean(p[ok] < alpha_level)
  n <- sum(ok)
  se <- sqrt(rate * (1 - rate) / n)
  data.frame(check = label, n = n, rate = rate, mc_se = se)
}

cat("---- NULL rate: conditional SE, 200 reps ----\n")
null_cond <- run_batch(200, c(1, 1), c(1, 1), c(0, 0), se_method = "conditional",
                        seed0 = 10000)
print(rbind(rate_row(null_cond$p_phi, "phi omnibus (conditional, null)"),
            rate_row(null_cond$p_alpha, "alpha omnibus (conditional, null)"),
            rate_row(null_cond$p_kappa, "kappa omnibus (conditional, null)")),
      row.names = FALSE)

cat("\n---- NULL rate: judge_bootstrap SE (boot_reps=50), 150 reps ----\n")
null_jb <- run_batch(150, c(1, 1), c(1, 1), c(0, 0), se_method = "judge_bootstrap",
                      boot_reps = 50, seed0 = 20000)
print(rbind(rate_row(null_jb$p_phi, "phi omnibus (judge_bootstrap, null)"),
            rate_row(null_jb$p_alpha, "alpha omnibus (judge_bootstrap, null)"),
            rate_row(null_jb$p_kappa, "kappa omnibus (judge_bootstrap, null)")),
      row.names = FALSE)

cat("\n---- Frame fit residual, null, pooled (conditional-SE null batch) ----\n")
fr <- null_cond$fit_resid
fr <- fr[is.finite(fr)]
cat("n cells:", length(fr), " mean:", mean(fr), " sd:", sd(fr),
    " reject rate |z|>1.96:", mean(abs(fr) > 1.96), "\n")

cat("\n---- POWER: planted panel-unit ratio (phi ratio ~1.56), conditional SE, 100 reps ----\n")
pow_phi <- run_batch(100, c(0.8, 1.25), c(1, 1), c(0, 0), se_method = "conditional",
                      seed0 = 30000)
print(rate_row(pow_phi$p_phi, "phi omnibus (conditional, phi planted)"), row.names = FALSE)

cat("\n---- POWER: planted set-unit ratio (alpha ratio 1.4), conditional SE, 100 reps ----\n")
pow_alpha <- run_batch(100, c(1, 1), c(1, 1.4), c(0, 0), se_method = "conditional",
                        seed0 = 40000)
print(rate_row(pow_alpha$p_alpha, "alpha omnibus (conditional, alpha planted)"), row.names = FALSE)

cat("\n---- POWER: planted set-origin difference (kappa=0.6), conditional SE, 100 reps ----\n")
pow_kappa <- run_batch(100, c(1, 1), c(1, 1), c(0, 0.6), se_method = "conditional",
                        seed0 = 50000)
print(rate_row(pow_kappa$p_kappa, "kappa omnibus (conditional, kappa planted)"), row.names = FALSE)

saveRDS(list(null_cond = null_cond, null_jb = null_jb, pow_phi = pow_phi,
             pow_alpha = pow_alpha, pow_kappa = pow_kappa),
        "tools/simval/round1/btl-efrm/03_null_power.rds")
cat("\ntotal time (s):", as.numeric(Sys.time() - t_start, units = "secs"), "\n")
