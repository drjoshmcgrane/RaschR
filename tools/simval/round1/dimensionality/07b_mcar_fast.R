suppressWarnings(pkgload::load_all(".", quiet=TRUE))
N <- 700; I <- 10; MISS <- 0.25

# dependence_magnitude null + power (reduced reps)
NREPS <- 60
null_p <- numeric(NREPS); null_d <- numeric(NREPS)
for (r in seq_len(NREPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS, seed = 43000 + r)
  f <- rasch(s)
  dm <- dependence_magnitude(f, dependent = "I06", independent = "I05")
  null_d[r] <- dm$d; null_p[r] <- dm$p
}
cat(sprintf("[dependence_magnitude null, MCAR] mean d=%.4f, reject@0.05=%.3f (n=%d, MCerr=%.3f)\n",
            mean(null_d), mean(null_p < 0.05), NREPS, sqrt(0.05*0.95/NREPS)))

POWREPS <- 30
pow_d <- numeric(POWREPS)
for (r in seq_len(POWREPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS,
                       dependence = list(pairs = list(c(5, 6)), strength = 3), seed = 44000 + r)
  f <- rasch(s)
  dm <- dependence_magnitude(f, dependent = "I06", independent = "I05")
  pow_d[r] <- dm$d
}
cat(sprintf("[dependence_magnitude power strength=3, MCAR] mean d=%.3f (n=%d)\n", mean(pow_d), POWREPS))

# dimensionality_test null + power (reduced reps)
NREPS2 <- 100
null_flag <- logical(NREPS2)
for (r in seq_len(NREPS2)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS, seed = 45000 + r)
  f <- rasch(s)
  dt <- dimensionality_test(f)
  null_flag[r] <- isTRUE(dt$multidimensional)
}
cat(sprintf("[dimensionality_test null, MCAR] flag rate = %.3f (n=%d, MCerr=%.3f)\n",
            mean(null_flag), NREPS2, sqrt(mean(null_flag)*(1-mean(null_flag))/NREPS2)))

POWREPS2 <- 60
pow_flag2 <- logical(POWREPS2)
for (r in seq_len(POWREPS2)) {
  s <- simulate_rasch(n_persons = N, n_items = 16, missing = MISS,
                       second_dim = list(items = 9:16, rho = 0.2), seed = 46000 + r)
  f <- rasch(s)
  dt <- dimensionality_test(f)
  pow_flag2[r] <- isTRUE(dt$multidimensional)
}
cat(sprintf("[dimensionality_test power (2nd dim rho=0.2), MCAR] flag rate = %.3f (n=%d, MCerr=%.3f)\n",
            mean(pow_flag2), POWREPS2, sqrt(mean(pow_flag2)*(1-mean(pow_flag2))/POWREPS2)))

# spread_test + combine_items null + power
REPS6 <- 40
null_flag6 <- logical(REPS6)
for (r in seq_len(REPS6)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS, seed = 47000 + r)
  f <- rasch(s)
  fc <- combine_items(f, list(c("I01", "I02")))
  st <- spread_test(fc)
  null_flag6[r] <- st$dependent
}
cat(sprintf("[spread_test null, MCAR] flag rate = %.3f (n=%d)\n", mean(null_flag6), REPS6))

pow_flag6 <- logical(REPS6)
for (r in seq_len(REPS6)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS,
                       dependence = list(pairs = list(c(4, 5)), strength = 3), seed = 48000 + r)
  f <- rasch(s)
  fc <- combine_items(f, list(c("I04", "I05")))
  st <- spread_test(fc)
  pow_flag6[r] <- st$dependent
}
cat(sprintf("[spread_test power (planted I04-I05 strength=3), MCAR] flag rate = %.3f (n=%d)\n", mean(pow_flag6), REPS6))
