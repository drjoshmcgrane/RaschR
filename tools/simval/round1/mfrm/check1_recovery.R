source("helpers.R")

args <- commandArgs(trailingOnly = TRUE)

I <- 6; R <- 6; CAT <- 4
base_seed <- 1000

# MCAR fits are much slower per person (the person-WLE root-finder needs many
# more iterations with sparse rows), so that condition runs at a smaller N
# and fewer reps to stay in budget; this is called out in the report.
conditions <- c("complete", "mcar25", "connected_k3")
cond_n <- c(complete = 200, mcar25 = 50, connected_k3 = 200)
cond_reps <- c(complete = 50, mcar25 = 20, connected_k3 = 50)
if (length(args)) {
  cond_reps[["complete"]] <- as.integer(args[1])
  cond_reps[["mcar25"]] <- max(5, round(as.integer(args[1]) * 20 / 50))
  cond_reps[["connected_k3"]] <- as.integer(args[1])
}

run_one <- function(cond, seed) {
  N <- cond_n[[cond]]
  d <- simulate_mfrm(N, I, R, n_categories = CAT, theta_sd = 1.2, item_sd = 1,
                     rater_severity_sd = 0.7, seed = seed)
  truth <- attr(d, "truth")
  d_use <- switch(cond,
    complete = d,
    mcar25 = { set.seed(seed + 5e5); d$score[sample.int(nrow(d), round(0.25 * nrow(d)))] <- NA; d },
    connected_k3 = make_connected_incomplete(d, k = 3, seed = seed))
  mf <- tryCatch(rasch_mfrm(d_use, person = "person", item = "item", score = "score",
                            facets = "rater"), error = function(e) e)
  if (inherits(mf, "error")) return(list(ok = FALSE, msg = conditionMessage(mf)))

  fe <- mf$facet_effects$rater
  sev_true <- truth$severity[fe$level]
  ie <- mf$item_effects
  loc_true <- truth$difficulty[ie$item]
  it <- mf$item_thresholds
  tau_true <- true_base_tau(CAT)[it$k] + truth$difficulty[it$item]

  list(ok = TRUE,
       sev_est = fe$severity, sev_se = fe$se, sev_true = sev_true,
       loc_est = ie$location, loc_se = ie$se, loc_true = loc_true,
       tau_est = it$tau, tau_se = it$se, tau_true = tau_true)
}

results <- list()
for (cond in conditions) {
  REPS <- cond_reps[[cond]]
  cat("=== condition:", cond, sprintf(" (N=%d, reps=%d) ===\n", cond_n[[cond]], REPS))
  acc <- list(sev_est = c(), sev_se = c(), sev_true = c(),
             loc_est = c(), loc_se = c(), loc_true = c(),
             tau_est = c(), tau_se = c(), tau_true = c())
  n_fail <- 0; fail_msgs <- c()
  t0 <- Sys.time()
  for (i in seq_len(REPS)) {
    r <- run_one(cond, base_seed + i)
    if (!r$ok) { n_fail <- n_fail + 1; fail_msgs <- c(fail_msgs, r$msg); next }
    for (nm in names(acc)) acc[[nm]] <- c(acc[[nm]], r[[nm]])
  }
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("  %d reps in %.1fs, %d failures\n", REPS, dt, n_fail))
  if (n_fail) cat("  fail msgs (unique):", paste(unique(fail_msgs), collapse = " | "), "\n")

  sev_cor <- cor(acc$sev_est, acc$sev_true)
  sev_rmse <- sqrt(mean((acc$sev_est - acc$sev_true)^2))
  sev_bias <- mean(acc$sev_est - acc$sev_true)
  sev_z <- (acc$sev_est - acc$sev_true) / acc$sev_se
  sev_z_sd <- sd(sev_z, na.rm = TRUE)

  loc_cor <- cor(acc$loc_est, acc$loc_true)
  loc_rmse <- sqrt(mean((acc$loc_est - acc$loc_true)^2))
  loc_bias <- mean(acc$loc_est - acc$loc_true)

  tau_cor <- cor(acc$tau_est, acc$tau_true)
  tau_rmse <- sqrt(mean((acc$tau_est - acc$tau_true)^2))
  tau_bias <- mean(acc$tau_est - acc$tau_true)

  results[[cond]] <- list(n_fail = n_fail, n_ok = REPS - n_fail,
    sev_cor = sev_cor, sev_rmse = sev_rmse, sev_bias = sev_bias, sev_z_sd = sev_z_sd,
    loc_cor = loc_cor, loc_rmse = loc_rmse, loc_bias = loc_bias,
    tau_cor = tau_cor, tau_rmse = tau_rmse, tau_bias = tau_bias,
    n_sev = length(acc$sev_est), n_loc = length(acc$loc_est), n_tau = length(acc$tau_est))
  print(results[[cond]])
}

saveRDS(results, "check1_results.rds")
cat("\n=== DONE check1 ===\n")
