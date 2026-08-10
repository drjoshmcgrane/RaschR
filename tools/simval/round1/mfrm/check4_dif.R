source("helpers.R")

args <- commandArgs(trailingOnly = TRUE)
REPS_NULL <- if (length(args) >= 1) as.integer(args[1]) else 80
REPS_POWER <- if (length(args) >= 2) as.integer(args[2]) else 60
REPS_SIZE <- if (length(args) >= 3) as.integer(args[3]) else 30

N <- 200; I <- 6; R <- 6; CAT <- 4
base_seed <- 4000
ALPHA <- 0.05
DIF_ITEM <- "I3"
DIF_SHIFT <- 0.7

run_dif <- function(shift, seed) {
  d <- sim_mfrm_dif(N, I, R, n_categories = CAT, theta_sd = 1.2, item_sd = 1,
                    rater_severity_sd = 0.7, dif_item = DIF_ITEM, dif_shift = shift, seed = seed)
  mf <- tryCatch(rasch_mfrm(d, person = "person", item = "item", score = "score",
                            facets = "rater", factors = "group"), error = function(e) e)
  if (inherits(mf, "error")) return(NULL)
  da <- tryCatch(dif_anova(mf), error = function(e) e)
  if (inherits(da, "error")) return(NULL)
  row <- da$summary[da$summary$item == DIF_ITEM, ]
  list(p_uniform = row$p_uniform_adj, flagged = row$uniform_DIF, mf = mf)
}

cat("=== NULL: no DIF planted (shift=0) ===\n")
flag_null <- c(); n_fail <- 0
t0 <- Sys.time()
for (i in seq_len(REPS_NULL)) {
  r <- run_dif(0, base_seed + i)
  if (is.null(r)) { n_fail <- n_fail + 1; next }
  flag_null <- c(flag_null, r$flagged)
}
cat(sprintf("  %d reps in %.1fs, %d failures\n", REPS_NULL,
            as.numeric(Sys.time() - t0, units = "secs"), n_fail))
null_rate <- mean(flag_null, na.rm = TRUE)
cat(sprintf("  P(dif_anova flags I3, adj alpha=.05) = %.4f (MC se %.4f), n=%d\n",
            null_rate, mc_se(null_rate, length(flag_null)), length(flag_null)))

cat(sprintf("\n=== POWER: planted person-group DIF on %s (shift=%.2f) ===\n", DIF_ITEM, DIF_SHIFT))
flag_pow <- c(); n_fail2 <- 0
t0 <- Sys.time()
for (i in seq_len(REPS_POWER)) {
  r <- run_dif(DIF_SHIFT, base_seed + 1e5 + i)
  if (is.null(r)) { n_fail2 <- n_fail2 + 1; next }
  flag_pow <- c(flag_pow, r$flagged)
}
cat(sprintf("  %d reps in %.1fs, %d failures\n", REPS_POWER,
            as.numeric(Sys.time() - t0, units = "secs"), n_fail2))
power_rate <- mean(flag_pow, na.rm = TRUE)
cat(sprintf("  P(dif_anova flags I3) = %.4f (MC se %.4f), n=%d\n",
            power_rate, mc_se(power_rate, length(flag_pow)), length(flag_pow)))

cat("\n=== dif_size magnitude recovery (facet severity should cancel) ===\n")
diffs <- c(); n_fail3 <- 0
t0 <- Sys.time()
for (i in seq_len(REPS_SIZE)) {
  d <- sim_mfrm_dif(N, I, R, n_categories = CAT, theta_sd = 1.2, item_sd = 1,
                    rater_severity_sd = 0.7, dif_item = DIF_ITEM, dif_shift = DIF_SHIFT,
                    seed = base_seed + 2e5 + i)
  mf <- tryCatch(rasch_mfrm(d, person = "person", item = "item", score = "score",
                            facets = "rater", factors = "group"), error = function(e) e)
  if (inherits(mf, "error")) { n_fail3 <- n_fail3 + 1; next }
  ds <- tryCatch(dif_size(mf, DIF_ITEM, by = "group"), error = function(e) e)
  if (inherits(ds, "error") || any(ds$levels$weak)) { n_fail3 <- n_fail3 + 1; next }
  # pairs$difference = level_a(A) - level_b(B); planted true B-A = DIF_SHIFT
  diffs <- c(diffs, -ds$pairs$difference)
}
cat(sprintf("  %d reps in %.1fs, %d dropped (weak/failed)\n", REPS_SIZE,
            as.numeric(Sys.time() - t0, units = "secs"), n_fail3))
mag_mean <- mean(diffs); mag_sd <- sd(diffs); mag_se <- mag_sd / sqrt(length(diffs))
cat(sprintf("  recovered B-A magnitude: mean=%.3f (true=%.2f), sd=%.3f, MC se of mean=%.3f, n=%d\n",
            mag_mean, DIF_SHIFT, mag_sd, mag_se, length(diffs)))

saveRDS(list(null_rate = null_rate, n_null = length(flag_null),
            power_rate = power_rate, n_power = length(flag_pow),
            mag_mean = mag_mean, mag_sd = mag_sd, mag_se = mag_se, n_mag = length(diffs)),
        "check4_results.rds")
cat("\n=== DONE check4 ===\n")
