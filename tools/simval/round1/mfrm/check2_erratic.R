source("tools/simval/round1/mfrm/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
REPS_NULL <- if (length(args) >= 1) as.integer(args[1]) else 100
REPS_POWER <- if (length(args) >= 2) as.integer(args[2]) else 100

N <- 150; I <- 6; R <- 6; CAT <- 4
base_seed <- 2000
ALPHA <- 0.05
ZCRIT <- qnorm(1 - ALPHA / 2)   # fit_resid ~ N(0,1) under the null

run_one <- function(erratic, seed) {
  d <- simulate_mfrm(N, I, R, n_categories = CAT, theta_sd = 1.2, item_sd = 1,
                     rater_severity_sd = 0.7, erratic_raters = erratic, seed = seed)
  truth <- attr(d, "truth")
  mf <- tryCatch(rasch_mfrm(d, person = "person", item = "item", score = "score",
                            facets = "rater"), error = function(e) e)
  if (inherits(mf, "error")) return(NULL)
  fe <- mf$facet_effects$rater
  list(level = fe$level, fit_resid = fe$fit_resid, fit_resid_pooled = fe$fit_resid_pooled,
       outfit_ms = fe$outfit_ms, infit_ms = fe$infit_ms, erratic = fe$level %in% truth$erratic)
}

## --- null: erratic_raters = 0, false-positive rate of clean raters ---------
cat("=== NULL: erratic_raters = 0 ===\n")
acc <- data.frame()
n_fail <- 0
t0 <- Sys.time()
for (i in seq_len(REPS_NULL)) {
  r <- run_one(0, base_seed + i)
  if (is.null(r)) { n_fail <- n_fail + 1; next }
  acc <- rbind(acc, data.frame(fit_resid = r$fit_resid, fit_resid_pooled = r$fit_resid_pooled,
                               outfit_ms = r$outfit_ms))
}
cat(sprintf("  %d reps in %.1fs, %d failures, %d rater-levels\n",
            REPS_NULL, as.numeric(Sys.time() - t0, units = "secs"), n_fail, nrow(acc)))
null_rate_fitresid <- mean(abs(acc$fit_resid) > ZCRIT, na.rm = TRUE)
null_rate_pooled <- mean(abs(acc$fit_resid_pooled) > ZCRIT, na.rm = TRUE)
cat(sprintf("  P(|fit_resid| > %.2f) = %.4f  (MC se %.4f)\n", ZCRIT, null_rate_fitresid,
            mc_se(null_rate_fitresid, nrow(acc))))
cat(sprintf("  P(|fit_resid_pooled| > %.2f) = %.4f  (MC se %.4f)\n", ZCRIT, null_rate_pooled,
            mc_se(null_rate_pooled, nrow(acc))))

## --- power: 1 erratic rater in 6 (erratic_raters ~ 1/6), detection --------
cat("\n=== POWER: erratic_raters = 1/6 ===\n")
accP <- data.frame()
n_fail2 <- 0
t0 <- Sys.time()
for (i in seq_len(REPS_POWER)) {
  r <- run_one(1/6, base_seed + 1e5 + i)
  if (is.null(r)) { n_fail2 <- n_fail2 + 1; next }
  accP <- rbind(accP, data.frame(fit_resid = r$fit_resid, fit_resid_pooled = r$fit_resid_pooled,
                                 outfit_ms = r$outfit_ms, erratic = r$erratic))
}
cat(sprintf("  %d reps in %.1fs, %d failures, %d rater-levels (%d erratic)\n",
            REPS_POWER, as.numeric(Sys.time() - t0, units = "secs"), n_fail2, nrow(accP),
            sum(accP$erratic)))
power_fitresid <- mean(abs(accP$fit_resid[accP$erratic]) > ZCRIT, na.rm = TRUE)
power_pooled <- mean(abs(accP$fit_resid_pooled[accP$erratic]) > ZCRIT, na.rm = TRUE)
clean_rate_fitresid <- mean(abs(accP$fit_resid[!accP$erratic]) > ZCRIT, na.rm = TRUE)
clean_rate_pooled <- mean(abs(accP$fit_resid_pooled[!accP$erratic]) > ZCRIT, na.rm = TRUE)
n_erratic <- sum(accP$erratic); n_clean <- sum(!accP$erratic)
cat(sprintf("  erratic raters: P(flag) fit_resid = %.4f (MC se %.4f), pooled = %.4f (MC se %.4f)  [n=%d]\n",
            power_fitresid, mc_se(power_fitresid, n_erratic),
            power_pooled, mc_se(power_pooled, n_erratic), n_erratic))
cat(sprintf("  clean raters (same batch): P(flag) fit_resid = %.4f (MC se %.4f), pooled = %.4f (MC se %.4f) [n=%d]\n",
            clean_rate_fitresid, mc_se(clean_rate_fitresid, n_clean),
            clean_rate_pooled, mc_se(clean_rate_pooled, n_clean), n_clean))

saveRDS(list(null_rate_fitresid = null_rate_fitresid, null_rate_pooled = null_rate_pooled,
            n_null = nrow(acc), power_fitresid = power_fitresid, power_pooled = power_pooled,
            n_erratic = n_erratic, clean_rate_fitresid = clean_rate_fitresid,
            clean_rate_pooled = clean_rate_pooled, n_clean = n_clean),
        "tools/simval/round1/mfrm/check2_results.rds")
cat("\n=== DONE check2 ===\n")
