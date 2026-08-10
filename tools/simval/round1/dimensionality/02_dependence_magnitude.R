suppressWarnings(pkgload::load_all("."), quiet=TRUE))

N <- 800; I <- 10

# ---- NULL: no planted dependence, test false-positive rate at alpha=0.05 ----
NREPS <- 150
null_d <- numeric(NREPS); null_p <- numeric(NREPS)
for (r in seq_len(NREPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, seed = 3000 + r)
  f <- rasch(s)
  dm <- dependence_magnitude(f, dependent = "I06", independent = "I05")
  null_d[r] <- dm$d; null_p[r] <- dm$p
}
cat("=== NULL (no planted dependence) ===\n")
cat(sprintf("mean d-hat = %.4f (expect ~0), reject rate at alpha=0.05 = %.3f (n=%d, MC err ~ %.3f)\n",
            mean(null_d), mean(null_p < 0.05), NREPS,
            sqrt(0.05*0.95/NREPS)))

# ---- RECOVERY across planted strengths: bias/RMSE + SE calibration ----
strengths <- c(1, 3, 6)
RREPS <- 60
recov <- lapply(strengths, function(str) {
  # first get the true simulated d-in-logits at this strength via a large
  # reference sample (truth is defined by the generative shift, not a closed
  # form -- estimate it with a very large N as ground truth)
  truth_s <- simulate_rasch(n_persons = 20000, n_items = I,
                             dependence = list(pairs = list(c(5, 6)), strength = str),
                             seed = 9000 + str)
  truth_f <- rasch(truth_s)
  truth_dm <- dependence_magnitude(truth_f, dependent = "I06", independent = "I05")
  truth_d <- truth_dm$d

  ds <- numeric(RREPS); ses <- numeric(RREPS)
  for (r in seq_len(RREPS)) {
    s <- simulate_rasch(n_persons = N, n_items = I,
                         dependence = list(pairs = list(c(5, 6)), strength = str),
                         seed = 4000 + str * 1000 + r)
    f <- rasch(s)
    dm <- dependence_magnitude(f, dependent = "I06", independent = "I05")
    ds[r] <- dm$d; ses[r] <- dm$se
  }
  list(strength = str, truth_d = truth_d, mean_d = mean(ds), bias = mean(ds) - truth_d,
       rmse = sqrt(mean((ds - truth_d)^2)),
       empirical_sd = sd(ds), mean_se = mean(ses),
       se_calibration = sd(ds) / mean(ses))
})

cat("\n=== RECOVERY across planted strengths (N=800, I=10, dep pair I05->I06) ===\n")
for (rr in recov) {
  cat(sprintf("strength=%.0f: truth_d(N=20000)=%.3f, mean_d(N=800)=%.3f, bias=%.3f, RMSE=%.3f, empSD=%.3f, meanSE=%.3f, SE-calib=%.3f\n",
              rr$strength, rr$truth_d, rr$mean_d, rr$bias, rr$rmse,
              rr$empirical_sd, rr$mean_se, rr$se_calibration))
}
