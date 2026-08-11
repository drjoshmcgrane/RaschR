suppressWarnings(pkgload::load_all(".", quiet = TRUE))
t_start <- Sys.time()

n_rep <- 40
panel_units <- c(0.8, 1.25); set_units <- c(1, 1.4); set_origins <- c(0, 0.5)

lphi2_est   <- numeric(n_rep)   # log phi[panel2] point estimate
lalpha2_est <- numeric(n_rep)   # log alpha[set2] point estimate

se_cond_lphi <- se_boot_lphi <- se_jboot_lphi <- numeric(n_rep)
se_cond_lalpha <- se_boot_lalpha <- se_jboot_lalpha <- numeric(n_rep)

n_fail <- 0
for (r in seq_len(n_rep)) {
  seed <- 5000 + r
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                          n_judges_per_panel = 6, reps_within = 20, reps_cross = 20,
                          panel_units = panel_units, set_units = set_units,
                          set_origins = set_origins, object_sd = 1, seed = seed)
  tr <- attr(d, "truth")

  fc <- tryCatch(btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
                           object_sets = tr$object_sets, se_method = "conditional"),
                 error = function(e) NULL)
  fj <- tryCatch(btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
                           object_sets = tr$object_sets, se_method = "judge_bootstrap",
                           boot_reps = 100), error = function(e) NULL)
  fb <- tryCatch(btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
                           object_sets = tr$object_sets, se_method = "bootstrap",
                           boot_reps = 100), error = function(e) NULL)
  if (is.null(fc) || is.null(fj) || is.null(fb)) { n_fail <- n_fail + 1
    lphi2_est[r] <- NA; lalpha2_est[r] <- NA
    se_cond_lphi[r] <- se_boot_lphi[r] <- se_jboot_lphi[r] <- NA
    se_cond_lalpha[r] <- se_boot_lalpha[r] <- se_jboot_lalpha[r] <- NA
    next
  }
  lphi2_est[r]   <- log(fc$phi_table$phi[fc$phi_table$panel == "panel2"])
  lalpha2_est[r] <- log(fc$alpha_table$alpha[fc$alpha_table$set == "set2"])

  se_cond_lphi[r]  <- fc$phi_table$se_log_phi[fc$phi_table$panel == "panel2"]
  se_jboot_lphi[r] <- fj$phi_table$se_log_phi[fj$phi_table$panel == "panel2"]
  se_boot_lphi[r]  <- fb$phi_table$se_log_phi[fb$phi_table$panel == "panel2"]

  se_cond_lalpha[r]  <- fc$alpha_table$se_log_alpha[fc$alpha_table$set == "set2"]
  se_jboot_lalpha[r] <- fj$alpha_table$se_log_alpha[fj$alpha_table$set == "set2"]
  se_boot_lalpha[r]  <- fb$alpha_table$se_log_alpha[fb$alpha_table$set == "set2"]
}

calib <- function(est, se_vec, method) {
  ok <- is.finite(est) & is.finite(se_vec)
  emp_sd <- sd(est[ok])
  mean_se <- mean(se_vec[ok])
  data.frame(method = method, n = sum(ok), emp_sd = emp_sd, mean_reported_se = mean_se,
             ratio = emp_sd / mean_se)
}

cat("=== log phi[panel2] SE calibration (n_rep=", n_rep, ", n_fail=", n_fail, ") ===\n", sep = "")
tab_phi <- rbind(calib(lphi2_est, se_cond_lphi, "conditional"),
                  calib(lphi2_est, se_boot_lphi, "parametric bootstrap"),
                  calib(lphi2_est, se_jboot_lphi, "judge bootstrap"))
print(tab_phi, row.names = FALSE)

cat("\n=== log alpha[set2] SE calibration (n_rep=", n_rep, ", n_fail=", n_fail, ") ===\n", sep = "")
tab_alpha <- rbind(calib(lalpha2_est, se_cond_lalpha, "conditional"),
                     calib(lalpha2_est, se_boot_lalpha, "parametric bootstrap"),
                     calib(lalpha2_est, se_jboot_lalpha, "judge bootstrap"))
print(tab_alpha, row.names = FALSE)

saveRDS(list(tab_phi = tab_phi, tab_alpha = tab_alpha, n_rep = n_rep, n_fail = n_fail),
        "tools/simval/round1/btl-efrm/02_se_calibration.rds")
cat("total time (s):", as.numeric(Sys.time() - t_start, units = "secs"), "\n")
