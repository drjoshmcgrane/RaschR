suppressWarnings(pkgload::load_all("."), quiet=TRUE))
t0 <- Sys.time()

n_per_group <- 200
items_per_set <- 8
n_sets <- 2
n_groups <- 2
set_unit_ratio <- 1.3
group_unit_ratio <- 1.2
n_reps <- 60

run_one <- function(se_method, boot_reps, seed) {
  d <- simulate_efrm(n_per_group, items_per_set, n_sets = n_sets, n_groups = n_groups,
                      set_unit_ratio = set_unit_ratio, group_unit_ratio = group_unit_ratio,
                      seed = seed)
  fit <- tryCatch(rasch_efrm(d, item_sets = attr(d, "truth")$item_sets, groups = "group",
                              se_method = se_method, boot_reps = boot_reps),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(log_alpha = log(fit$alpha_table$alpha), se_log_alpha = fit$alpha_table$se_log_alpha,
       log_phi = log(fit$phi_table$phi), se_log_phi = fit$phi_table$se_log_phi)
}

gspan <- function(ratio, n) { u <- exp(seq(0, log(ratio), length.out = n)); u / exp(mean(log(u))) }
true_log_alpha <- log(gspan(set_unit_ratio, n_sets))
true_log_phi <- log(gspan(group_unit_ratio, n_groups))

calib_summary <- function(method, boot_reps) {
  out <- vector("list", n_reps)
  for (r in seq_len(n_reps)) out[[r]] <- run_one(method, boot_reps, seed = 5000L + r)
  ok <- !vapply(out, is.null, TRUE)
  out <- out[ok]
  la <- do.call(rbind, lapply(out, `[[`, "log_alpha"))       # n_reps x S
  se_la <- do.call(rbind, lapply(out, `[[`, "se_log_alpha"))
  lp <- do.call(rbind, lapply(out, `[[`, "log_phi"))
  se_lp <- do.call(rbind, lapply(out, `[[`, "se_log_phi"))
  emp_sd_alpha <- apply(la, 2, sd)
  mean_se_alpha <- colMeans(se_la)
  emp_sd_phi <- apply(lp, 2, sd)
  mean_se_phi <- colMeans(se_lp)
  list(method = method, n_ok = sum(ok),
       alpha_ratio = emp_sd_alpha / mean_se_alpha,
       phi_ratio = emp_sd_phi / mean_se_phi,
       alpha_bias = colMeans(la) - true_log_alpha,
       phi_bias = colMeans(lp) - true_log_phi,
       emp_sd_alpha = emp_sd_alpha, mean_se_alpha = mean_se_alpha,
       emp_sd_phi = emp_sd_phi, mean_se_phi = mean_se_phi)
}

res_hybrid <- calib_summary("hybrid", 300)
cat("HYBRID  n_ok=", res_hybrid$n_ok, "\n")
cat("  alpha empSD/meanSE:", round(res_hybrid$alpha_ratio, 3), " bias(log):", round(res_hybrid$alpha_bias,3), "\n")
cat("  phi   empSD/meanSE:", round(res_hybrid$phi_ratio, 3), " bias(log):", round(res_hybrid$phi_bias,3), "\n")

t_mid <- Sys.time()
cat("elapsed after hybrid:", as.numeric(t_mid - t0, units = "secs"), "s\n")

res_boot <- calib_summary("bootstrap", 40)
cat("BOOTSTRAP(boot_reps=40)  n_ok=", res_boot$n_ok, "\n")
cat("  alpha empSD/meanSE:", round(res_boot$alpha_ratio, 3), " bias(log):", round(res_boot$alpha_bias,3), "\n")
cat("  phi   empSD/meanSE:", round(res_boot$phi_ratio, 3), " bias(log):", round(res_boot$phi_bias,3), "\n")

saveRDS(list(hybrid = res_hybrid, boot = res_boot),
        "tools/simval/round1/efrm/check2_results.rds")
cat("total elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
