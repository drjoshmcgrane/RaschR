suppressWarnings(pkgload::load_all("."), quiet=TRUE))
t0 <- Sys.time()

n_per_group <- 200
items_per_set <- 8
n_sets <- 2
n_groups <- 2
n_reps <- 200
alpha_level <- 0.05

apply_mcar <- function(d, p, seed) {
  set.seed(seed + 77777L)
  item_cols <- setdiff(names(d), c("id", "group"))
  X <- as.matrix(d[item_cols])
  mask <- matrix(runif(length(X)) < p, nrow(X), ncol(X))
  X[mask] <- NA
  d[item_cols] <- as.data.frame(X)
  d
}

run_condition <- function(set_unit_ratio, group_unit_ratio, n_reps, seed0, mcar = 0) {
  p_alpha <- rep(NA_real_, n_reps); p_phi <- rep(NA_real_, n_reps)
  n_ok <- 0L
  for (r in seq_len(n_reps)) {
    d <- simulate_efrm(n_per_group, items_per_set, n_sets = n_sets, n_groups = n_groups,
                        set_unit_ratio = set_unit_ratio, group_unit_ratio = group_unit_ratio,
                        seed = seed0 + r)
    if (mcar > 0) d <- apply_mcar(d, mcar, seed0 + r)
    fit <- tryCatch(rasch_efrm(d, item_sets = attr(d, "truth")$item_sets, groups = "group",
                                se_method = "hybrid"),
                     error = function(e) NULL)
    if (is.null(fit)) next
    n_ok <- n_ok + 1L
    uo <- fit$efrm_vs_rasch$unit_omnibus
    if (!is.null(uo)) {
      ra <- uo$p[uo$term == "set units (alpha)"]
      rp <- uo$p[uo$term == "group units (phi)"]
      if (length(ra)) p_alpha[r] <- ra
      if (length(rp)) p_phi[r] <- rp
    }
  }
  list(p_alpha = p_alpha, p_phi = p_phi, n_ok = n_ok)
}

mc_err <- function(p_hat, n) sqrt(p_hat * (1 - p_hat) / n)

report <- function(label, res) {
  fa <- mean(res$p_alpha < alpha_level, na.rm = TRUE)
  fp <- mean(res$p_phi < alpha_level, na.rm = TRUE)
  na <- sum(is.finite(res$p_alpha)); np <- sum(is.finite(res$p_phi))
  cat(sprintf("[%s] n_ok=%d  alpha-flag-rate=%.3f (n=%d, mcse=%.3f)  phi-flag-rate=%.3f (n=%d, mcse=%.3f)\n",
              label, res$n_ok, fa, na, mc_err(fa, na), fp, np, mc_err(fp, np)))
}

cat("== complete data ==\n")
r_null_alpha <- run_condition(1.0, 1.3, n_reps, seed0 = 10000)   # alpha null, phi planted
report("alpha-null/phi-planted(1.3)", r_null_alpha)

r_null_phi <- run_condition(1.3, 1.0, n_reps, seed0 = 20000)     # alpha planted, phi null
report("alpha-planted(1.3)/phi-null", r_null_phi)

t_mid <- Sys.time(); cat("elapsed after complete-data conditions:", as.numeric(t_mid - t0, units="secs"), "s\n")

cat("== 25% MCAR ==\n")
n_reps_mcar <- 120
r_mcar_alpha <- run_condition(1.0, 1.3, n_reps_mcar, seed0 = 30000, mcar = 0.25)
report("MCAR alpha-null/phi-planted(1.3)", r_mcar_alpha)
r_mcar_phi <- run_condition(1.3, 1.0, n_reps_mcar, seed0 = 40000, mcar = 0.25)
report("MCAR alpha-planted(1.3)/phi-null", r_mcar_phi)

saveRDS(list(r_null_alpha=r_null_alpha, r_null_phi=r_null_phi,
             r_mcar_alpha=r_mcar_alpha, r_mcar_phi=r_mcar_phi),
        "tools/simval/round1/efrm/check3_results.rds")
cat("total elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
