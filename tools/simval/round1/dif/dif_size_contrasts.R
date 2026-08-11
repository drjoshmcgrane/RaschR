suppressWarnings(pkgload::load_all(".", quiet=TRUE))
t_start <- Sys.time()

N <- 500; n_items <- 20; reps <- 150; target_item <- "I06"
truth_uniform <- 0.7

est <- c(); se <- c(); lower <- c(); upper <- c(); weak_any <- c()
contrast_est <- c(); contrast_match_diff <- c(); fails <- 0

for (k in seq_len(reps)) {
  seed <- 90000 + k
  d <- simulate_rasch(N, n_items, dif = list(items = target_item, uniform = truth_uniform),
                       n_groups = 2, seed = seed)
  fit <- tryCatch(rasch(d, id = "id", factors = "group"), error = function(e) NULL)
  if (is.null(fit)) { fails <- fails + 1; next }
  ds <- tryCatch(dif_size(fit, target_item, by = "group"), error = function(e) NULL)
  if (is.null(ds)) { fails <- fails + 1; next }
  # pairs$difference is level_a - level_b = g1 - g2; truth is g2 = g1 + 0.7
  # (dif applies to the LAST group, i.e. g2), so g1 - g2 should be -0.7
  diff_g1_g2 <- ds$pairs$difference[1]
  est <- c(est, diff_g1_g2)
  se <- c(se, ds$pairs$se[1])
  lower <- c(lower, ds$pairs$lower[1])
  upper <- c(upper, ds$pairs$upper[1])
  weak_any <- c(weak_any, any(ds$levels$weak))

  dc <- tryCatch(dif_contrasts(fit, items = target_item), error = function(e) NULL)
  if (!is.null(dc)) {
    row <- dc$table[dc$table$item == target_item, ]
    # dif_contrasts label is "group: g2 - g1" -> should equal -(diff_g1_g2)
    contrast_est <- c(contrast_est, row$estimate[1])
    contrast_match_diff <- c(contrast_match_diff, row$estimate[1] - (-diff_g1_g2))
  }
}

true_diff <- -truth_uniform   # g1 - g2
bias <- mean(est, na.rm = TRUE) - true_diff
rmse <- sqrt(mean((est - true_diff)^2, na.rm = TRUE))
corr_with_truth <- NA  # single true value, so report bias/RMSE instead of correlation
coverage <- mean(lower <= true_diff & upper >= true_diff, na.rm = TRUE)
se_calib <- sd(est, na.rm = TRUE) / mean(se, na.rm = TRUE)

cat("=== dif_size recovery: planted uniform DIF", truth_uniform, "logits on", target_item, "===\n")
cat("n reps used:", length(est), " fails:", fails, "\n")
cat("mean estimate (g1-g2):", mean(est, na.rm=TRUE), " true:", true_diff, "\n")
cat("bias:", bias, "\n")
cat("rmse:", rmse, "\n")
cat("95% CI coverage:", coverage, "\n")
cat("SE calibration (empirical SD / mean reported SE):", se_calib, "\n")
cat("weak-category flagged (false-fire) rate on healthy data:", mean(weak_any), "\n")

cat("\n=== dif_contrasts vs dif_size agreement (2-level factor) ===\n")
cat("n compared:", length(contrast_match_diff), "\n")
cat("max |dif_contrasts - (-dif_size diff)|:", max(abs(contrast_match_diff), na.rm=TRUE), "\n")
cat("mean |diff|:", mean(abs(contrast_match_diff), na.rm=TRUE), "\n")

saveRDS(list(est=est, se=se, lower=lower, upper=upper, weak_any=weak_any,
             contrast_est=contrast_est, contrast_match_diff=contrast_match_diff,
             fails=fails, true_diff=true_diff),
        "tools/simval/round1/dif/dif_size_contrasts.rds")
cat("\ntotal time:", as.numeric(Sys.time() - t_start, units = "secs"), "s\n")
