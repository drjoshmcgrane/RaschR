suppressWarnings(pkgload::load_all(".", quiet=TRUE))
source("tools/simval/round1/efrm/sim_crossed.R")
t0 <- Sys.time()

n_reps <- 100
alpha_level <- 0.05

run_condition <- function(region_effect, cohort_effect, n_reps, seed0) {
  p_region <- rep(NA_real_, n_reps); p_cohort <- rep(NA_real_, n_reps); p_int <- rep(NA_real_, n_reps)
  n_ok <- 0L
  for (r in seq_len(n_reps)) {
    d <- sim_crossed_phi(150, 10, region_effect = region_effect, cohort_effect = cohort_effect,
                          seed = seed0 + r)
    fit <- tryCatch(rasch_efrm(d, item_sets = list(only = attr(d, "item_names")),
                                groups = c("region", "cohort"),
                                items = attr(d, "item_names"), se_method = "hybrid"),
                     error = function(e) NULL)
    if (is.null(fit) || is.null(fit$phi_factorial_tests)) next
    n_ok <- n_ok + 1L
    tt <- fit$phi_factorial_tests
    p_region[r] <- tt$p[tt$term == "region"]
    p_cohort[r] <- tt$p[tt$term == "cohort"]
    p_int[r] <- tt$p[tt$term == "region:cohort"]
  }
  list(p_region = p_region, p_cohort = p_cohort, p_int = p_int, n_ok = n_ok)
}

mc_err <- function(p_hat, n) sqrt(p_hat * (1 - p_hat) / n)
report <- function(label, res) {
  fr <- mean(res$p_region < alpha_level, na.rm = TRUE); nr <- sum(is.finite(res$p_region))
  fc <- mean(res$p_cohort < alpha_level, na.rm = TRUE); nc <- sum(is.finite(res$p_cohort))
  fi <- mean(res$p_int < alpha_level, na.rm = TRUE); ni <- sum(is.finite(res$p_int))
  cat(sprintf("[%s] n_ok=%d region-flag=%.3f(mcse %.3f,n=%d) cohort-flag=%.3f(mcse %.3f,n=%d) interaction-flag=%.3f(n=%d)\n",
              label, res$n_ok, fr, mc_err(fr, nr), nr, fc, mc_err(fc, nc), nc, fi, ni))
}

r_null <- run_condition(0, 0, n_reps, seed0 = 50000)          # both null
report("both-null", r_null)

r_power <- run_condition(log(1.3), 0, n_reps, seed0 = 60000)  # region planted, cohort null
report("region-planted(1.3)/cohort-null", r_power)

saveRDS(list(null = r_null, power = r_power),
        "tools/simval/round1/efrm/check4_results.rds")
cat("total elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
