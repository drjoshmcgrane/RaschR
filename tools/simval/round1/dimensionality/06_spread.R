suppressWarnings(pkgload::load_all("."), quiet=TRUE))

N <- 600; I <- 8
REPS <- 100

# ---- NULL: combine two genuinely independent items -> spread above bound, not flagged ----
null_flag <- logical(REPS); null_spread <- numeric(REPS)
for (r in seq_len(REPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, seed = 11000 + r)
  f <- rasch(s)
  fc <- combine_items(f, list(c("I01", "I02")))
  st <- spread_test(fc)
  null_flag[r] <- st$dependent
  null_spread[r] <- st$spread
}
cat("=== NULL: combined pair is genuinely independent ===\n")
cat(sprintf("false-positive (flagged dependent) rate = %.3f (n=%d, MC err=%.3f); mean spread = %.3f (bound = 0.690)\n",
            mean(null_flag), REPS, sqrt(mean(null_flag)*(1-mean(null_flag))/REPS), mean(null_spread)))

# ---- POWER: combine a genuinely dependent pair (planted via simulate_rasch dependence=) ----
pow_flag <- logical(REPS); pow_spread <- numeric(REPS)
for (r in seq_len(REPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I,
                       dependence = list(pairs = list(c(4, 5)), strength = 3), seed = 12000 + r)
  f <- rasch(s)
  fc <- combine_items(f, list(c("I04", "I05")))
  st <- spread_test(fc)
  pow_flag[r] <- st$dependent
  pow_spread[r] <- st$spread
}
cat("\n=== POWER: combined pair carries planted dependence (strength=3) ===\n")
cat(sprintf("true-positive (flagged dependent) rate = %.3f (n=%d, MC err=%.3f); mean spread = %.3f (bound = 0.690)\n",
            mean(pow_flag), REPS, sqrt(mean(pow_flag)*(1-mean(pow_flag))/REPS), mean(pow_spread)))
