suppressWarnings(pkgload::load_all("."), quiet = TRUE))
t_start <- Sys.time()
alpha_level <- 0.05

extract_omnibus_p <- function(fit, term) {
  if (is.null(fit$unit_omnibus)) return(NA_real_)
  row <- fit$unit_omnibus[fit$unit_omnibus$term == term, ]
  if (!nrow(row)) return(NA_real_)
  row$p[1]
}

run_batch <- function(n_rep, panel_units, set_units, set_origins, seed0, boot_reps = 50) {
  p_phi <- p_alpha <- numeric(n_rep)
  n_fail <- 0
  for (r in seq_len(n_rep)) {
    seed <- seed0 + r
    d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                            n_judges_per_panel = 6, reps_within = 20, reps_cross = 20,
                            panel_units = panel_units, set_units = set_units,
                            set_origins = set_origins, object_sd = 1, seed = seed)
    tr <- attr(d, "truth")
    fit <- tryCatch(btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
                              object_sets = tr$object_sets, se_method = "judge_bootstrap",
                              boot_reps = boot_reps), error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$converged)) { n_fail <- n_fail + 1
      p_phi[r] <- p_alpha[r] <- NA; next }
    p_phi[r]   <- extract_omnibus_p(fit, "panel units (phi)")
    p_alpha[r] <- extract_omnibus_p(fit, "set units (alpha)")
  }
  list(p_phi = p_phi, p_alpha = p_alpha, n_fail = n_fail, n_rep = n_rep)
}

rate_row <- function(p, label) {
  ok <- is.finite(p)
  rate <- mean(p[ok] < alpha_level)
  n <- sum(ok)
  se <- sqrt(rate * (1 - rate) / n)
  data.frame(check = label, n = n, rate = rate, mc_se = se)
}

cat("---- POWER (judge_bootstrap default SE): planted phi ratio ~1.56, 40 reps ----\n")
pow_phi_jb <- run_batch(40, c(0.8, 1.25), c(1, 1), c(0, 0), seed0 = 60000)
print(rate_row(pow_phi_jb$p_phi, "phi omnibus (judge_bootstrap, phi planted)"), row.names = FALSE)

cat("\n---- POWER (judge_bootstrap default SE): planted alpha ratio 1.4, 40 reps ----\n")
pow_alpha_jb <- run_batch(40, c(1, 1), c(1, 1.4), c(0, 0), seed0 = 70000)
print(rate_row(pow_alpha_jb$p_alpha, "alpha omnibus (judge_bootstrap, alpha planted)"), row.names = FALSE)

saveRDS(list(pow_phi_jb = pow_phi_jb, pow_alpha_jb = pow_alpha_jb),
        "tools/simval/round1/btl-efrm/03b_power_jboot.rds")
cat("\ntotal time (s):", as.numeric(Sys.time() - t_start, units = "secs"), "\n")
