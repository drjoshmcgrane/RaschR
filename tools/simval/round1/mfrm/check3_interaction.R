source("tools/simval/round1/mfrm/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
REPS_NULL <- if (length(args) >= 1) as.integer(args[1]) else 100
REPS_POWER <- if (length(args) >= 2) as.integer(args[2]) else 60

N <- 150; I <- 6; R <- 6; CAT <- 4
base_seed <- 3000
ALPHA <- 0.05

run_one <- function(bias, seed) {
  interaction_arg <- if (!is.na(bias))
    list(item = "I3", rater = "R3", bias = bias) else NULL
  d <- simulate_mfrm(N, I, R, n_categories = CAT, theta_sd = 1.2, item_sd = 1,
                     rater_severity_sd = 0.7, interaction = interaction_arg, seed = seed)
  mf <- tryCatch(rasch_mfrm(d, person = "person", item = "item", score = "score",
                            facets = "rater", interaction = "rater"), error = function(e) e)
  if (inherits(mf, "error")) return(NULL)
  mf$interaction_test$p
}

cat("=== NULL: no interaction planted ===\n")
p_null <- c(); n_fail <- 0
t0 <- Sys.time()
for (i in seq_len(REPS_NULL)) {
  p <- run_one(NA, base_seed + i)
  if (is.null(p)) { n_fail <- n_fail + 1; next }
  p_null <- c(p_null, p)
}
cat(sprintf("  %d reps in %.1fs, %d failures\n", REPS_NULL,
            as.numeric(Sys.time() - t0, units = "secs"), n_fail))
null_rate <- mean(p_null < ALPHA, na.rm = TRUE)
cat(sprintf("  P(reject at alpha=.05) = %.4f (MC se %.4f), n=%d\n",
            null_rate, mc_se(null_rate, length(p_null)), length(p_null)))

cat("\n=== POWER: planted item x rater interaction (I3 x R3, bias=1.5) ===\n")
p_pow <- c(); n_fail2 <- 0
t0 <- Sys.time()
for (i in seq_len(REPS_POWER)) {
  p <- run_one(1.5, base_seed + 1e5 + i)
  if (is.null(p)) { n_fail2 <- n_fail2 + 1; next }
  p_pow <- c(p_pow, p)
}
cat(sprintf("  %d reps in %.1fs, %d failures\n", REPS_POWER,
            as.numeric(Sys.time() - t0, units = "secs"), n_fail2))
power_rate <- mean(p_pow < ALPHA, na.rm = TRUE)
cat(sprintf("  P(reject at alpha=.05) = %.4f (MC se %.4f), n=%d\n",
            power_rate, mc_se(power_rate, length(p_pow)), length(p_pow)))

saveRDS(list(null_rate = null_rate, n_null = length(p_null),
            power_rate = power_rate, n_power = length(p_pow)),
        "tools/simval/round1/mfrm/check3_results.rds")
cat("\n=== DONE check3 ===\n")
