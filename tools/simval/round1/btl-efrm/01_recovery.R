suppressWarnings(pkgload::load_all("."), quiet = TRUE))
t_start <- Sys.time()

run_condition <- function(label, panel_units, set_units, set_origins, n_rep = 60, seed0 = 1000) {
  G <- length(panel_units); S <- length(set_units)
  free_s <- if (S > 1) 2:S else integer(0)
  free_g <- if (G > 1) 1:G else integer(0)
  # storage
  lphi_true <- if (G > 1) log(panel_units / exp(mean(log(panel_units)))) else numeric(0)
  lalpha_true <- if (S > 1) log(set_units[free_s] / set_units[1]) else numeric(0)
  kappa_true  <- if (S > 1) (set_origins[free_s] - set_origins[1]) else numeric(0)

  est_lphi <- matrix(NA_real_, n_rep, G)
  est_lalpha <- matrix(NA_real_, n_rep, max(S - 1, 0))
  est_kappa  <- matrix(NA_real_, n_rep, max(S - 1, 0))
  v_true_all <- v_est_all <- vector("list", n_rep)
  beta_true_all <- beta_est_all <- vector("list", n_rep)
  n_fail <- 0

  for (r in seq_len(n_rep)) {
    seed <- seed0 + r
    d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = S, n_panels = G,
                            n_judges_per_panel = 6, reps_within = 20, reps_cross = 20,
                            panel_units = panel_units, set_units = set_units,
                            set_origins = set_origins, object_sd = 1, seed = seed)
    tr <- attr(d, "truth")
    fit <- tryCatch(btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
                              object_sets = tr$object_sets, se_method = "conditional"),
                     error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$converged)) { n_fail <- n_fail + 1; next }
    if (G > 1) est_lphi[r, ] <- log(fit$phi_table$phi)[match(names(tr$phi), fit$phi_table$panel)]
    if (S > 1) {
      est_lalpha[r, ] <- log(fit$alpha_table$alpha)[match(names(tr$alpha)[free_s], fit$alpha_table$set)]
      est_kappa[r, ]  <- fit$kappa_table$kappa[match(names(tr$kappa)[free_s], fit$kappa_table$set)]
    }
    v_true_all[[r]] <- tr$v[fit$objects$object]
    v_est_all[[r]]  <- fit$objects$v
    beta_true_all[[r]] <- tr$beta[fit$objects$object]
    beta_est_all[[r]]  <- fit$objects$beta_set
  }

  ok <- !sapply(v_true_all, is.null)
  v_true <- unlist(v_true_all[ok]); v_est <- unlist(v_est_all[ok])
  b_true <- unlist(beta_true_all[ok]); b_est <- unlist(beta_est_all[ok])

  summarize_vec <- function(true, est) {
    ok2 <- is.finite(true) & is.finite(est)
    true <- true[ok2]; est <- est[ok2]
    data.frame(n = length(true),
               cor = if (length(true) > 2) cor(true, est) else NA,
               rmse = sqrt(mean((est - true)^2)),
               bias = mean(est - true))
  }

  list(label = label, n_rep = n_rep, n_fail = n_fail,
       phi = if (G > 1) summarize_vec(matrix(rep(lphi_true, each = n_rep), n_rep), est_lphi) else NULL,
       alpha = if (S > 1) summarize_vec(matrix(rep(lalpha_true, each = n_rep), n_rep), est_lalpha) else NULL,
       kappa = if (S > 1) summarize_vec(matrix(rep(kappa_true, each = n_rep), n_rep), est_kappa) else NULL,
       v = summarize_vec(v_true, v_est),
       beta = summarize_vec(b_true, b_est))
}

set.seed(1)
res_null <- run_condition("null (phi=alpha=1,kappa=0)",
                           panel_units = c(1, 1), set_units = c(1, 1), set_origins = c(0, 0))
res_mod  <- run_condition("moderate (phi ratio~1.56, alpha ratio 1.4, kappa 0.5)",
                           panel_units = c(0.8, 1.25), set_units = c(1, 1.4), set_origins = c(0, 0.5))
res_big  <- run_condition("large (phi ratio~2.25, alpha ratio 2.0, kappa 1.0)",
                           panel_units = c(0.667, 1.5), set_units = c(1, 2.0), set_origins = c(0, 1.0))

saveRDS(list(res_null, res_mod, res_big),
        "tools/simval/round1/btl-efrm/01_recovery.rds")

for (res in list(res_null, res_mod, res_big)) {
  cat("=== ", res$label, " (n_rep=", res$n_rep, ", n_fail=", res$n_fail, ") ===\n", sep = "")
  if (!is.null(res$phi))   { cat("  log phi:   "); print(res$phi, row.names = FALSE) }
  if (!is.null(res$alpha)) { cat("  log alpha: "); print(res$alpha, row.names = FALSE) }
  if (!is.null(res$kappa)) { cat("  kappa:     "); print(res$kappa, row.names = FALSE) }
  cat("  v (common scale): "); print(res$v, row.names = FALSE)
  cat("  beta (within-set): "); print(res$beta, row.names = FALSE)
}
cat("total time (s):", as.numeric(Sys.time() - t_start, units = "secs"), "\n")
