source("tools/simval/round1/mfrm/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args)) as.integer(args[1]) else 40

N <- 150; I <- 6; R <- 8; CAT <- 4
base_seed <- 5000
ALPHA <- 0.05

run_one <- function(halo, seed) {
  d <- simulate_mfrm(N, I, R, n_categories = CAT, theta_sd = 1.2, item_sd = 1,
                     rater_severity_sd = 0.7, halo = halo, seed = seed)
  truth <- attr(d, "truth")
  mf <- tryCatch(rasch_mfrm(d, person = "person", item = "item", score = "score",
                            facets = "rater", interaction = "rater"), error = function(e) e)
  if (inherits(mf, "error")) return(NULL)
  fe <- mf$facet_effects$rater
  ie <- mf$interaction_effects
  list(fe_halo = fe$level %in% truth$halo, fe_fit_resid = fe$fit_resid,
       ie_halo = ie$level %in% truth$halo, ie_gamma = ie$gamma,
       ie_sig = ie$significant, omni_p = mf$interaction_test$p)
}

cat(sprintf("=== halo raters (2/8) vs clean, %d reps ===\n", REPS))
acc_fe <- data.frame(); acc_ie <- data.frame(); omni_p <- c(); n_fail <- 0
t0 <- Sys.time()
for (i in seq_len(REPS)) {
  r <- run_one(0.25, base_seed + i)
  if (is.null(r)) { n_fail <- n_fail + 1; next }
  acc_fe <- rbind(acc_fe, data.frame(halo = r$fe_halo, fit_resid = r$fe_fit_resid))
  acc_ie <- rbind(acc_ie, data.frame(halo = r$ie_halo, gamma = r$ie_gamma, sig = r$ie_sig))
  omni_p <- c(omni_p, r$omni_p)
}
cat(sprintf("  %d reps in %.1fs, %d failures\n", REPS,
            as.numeric(Sys.time() - t0, units = "secs"), n_fail))

cat("\n-- additive-model rater fit_resid (halo vs clean) --\n")
print(tapply(acc_fe$fit_resid, acc_fe$halo, function(z) c(mean = mean(z), sd = sd(z))))

cat("\n-- item x rater interaction gamma (halo vs clean) --\n")
print(tapply(abs(acc_ie$gamma), acc_ie$halo, mean))

sig_rate_halo <- mean(acc_ie$sig[acc_ie$halo])
sig_rate_clean <- mean(acc_ie$sig[!acc_ie$halo])
n_halo_cells <- sum(acc_ie$halo); n_clean_cells <- sum(!acc_ie$halo)
cat(sprintf("\nHolm-significant interaction cell rate: halo = %.4f (MC se %.4f, n=%d), clean = %.4f (MC se %.4f, n=%d)\n",
            sig_rate_halo, mc_se(sig_rate_halo, n_halo_cells), n_halo_cells,
            sig_rate_clean, mc_se(sig_rate_clean, n_clean_cells), n_clean_cells))

omni_reject <- mean(omni_p < ALPHA, na.rm = TRUE)
cat(sprintf("Omnibus interaction_test rejection rate (planted halo present): %.4f (MC se %.4f, n=%d)\n",
            omni_reject, mc_se(omni_reject, length(omni_p)), length(omni_p)))

## null comparison batch: no halo, same omnibus test
cat("\n=== NULL (no halo) reference batch ===\n")
omni_p0 <- c(); n_fail0 <- 0
t0 <- Sys.time()
for (i in seq_len(REPS)) {
  r <- run_one(0, base_seed + 1e5 + i)
  if (is.null(r)) { n_fail0 <- n_fail0 + 1; next }
  omni_p0 <- c(omni_p0, r$omni_p)
}
cat(sprintf("  %d reps in %.1fs, %d failures\n", REPS,
            as.numeric(Sys.time() - t0, units = "secs"), n_fail0))
omni_reject0 <- mean(omni_p0 < ALPHA, na.rm = TRUE)
cat(sprintf("Omnibus interaction_test rejection rate (no halo): %.4f (MC se %.4f, n=%d)\n",
            omni_reject0, mc_se(omni_reject0, length(omni_p0)), length(omni_p0)))

saveRDS(list(sig_rate_halo = sig_rate_halo, n_halo_cells = n_halo_cells,
            sig_rate_clean = sig_rate_clean, n_clean_cells = n_clean_cells,
            omni_reject = omni_reject, n_omni = length(omni_p),
            omni_reject0 = omni_reject0, n_omni0 = length(omni_p0)),
        "tools/simval/round1/mfrm/check5_results.rds")
cat("\n=== DONE check5 ===\n")
