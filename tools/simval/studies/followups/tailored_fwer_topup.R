# FWER top-up for the tailored bootstrap principal claim: fresh seeds,
# same design as tools/simval/studies/tailored-bootstrap.R family (a).
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
N <- 300L; I <- 8L; diff <- seq(-1.5, 1.8, length.out = I)
BUDGET_S <- 10800; t0 <- Sys.time()
any_sig <- logical(0); floor_w <- 0L; n_ref <- 0L
r <- 0L
while (as.numeric(Sys.time() - t0, units = "secs") < BUDGET_S) {
  r <- r + 1L
  dat <- simulate_rasch(N, I, difficulty = diff, seed = 31e6 + r)
  fit0 <- tryCatch(rasch(dat, id = "id"), error = function(e) NULL)
  if (is.null(fit0) || !isTRUE(fit0$est$converged)) { n_ref <- n_ref + 1L; next }
  ta <- withCallingHandlers(
    tryCatch(tailored_analysis(fit0, chance = 0.25, se_method = "bootstrap",
                               boot_reps = 399), error = function(e) NULL),
    warning = function(w) {
      if (grepl("smallest achievable", conditionMessage(w))) floor_w <<- floor_w + 1L
      invokeRestart("muffleWarning")
    })
  if (is.null(ta)) { n_ref <- n_ref + 1L; next }
  any_sig <- c(any_sig, any(ta$table$significant %in% TRUE))
}
n <- length(any_sig); fw <- mean(any_sig)
cat(sprintf("TOPUP FWER: %d additional reps, familywise=%.4f (mc %.4f), refusals %d, floor-warned %d\n",
    n, fw, sqrt(max(fw*(1-fw), 0)/n), n_ref, floor_w))
