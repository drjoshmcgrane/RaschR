suppressWarnings(pkgload::load_all(".", quiet=TRUE))

N <- 700; I <- 10
MISS <- 0.25

cat("############ MCAR (25%) condition ############\n\n")

# ---- (1) Q3 null + power ----
REPS <- 150
null_heur <- logical(REPS)
for (r in seq_len(REPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS, seed = 21000 + r)
  f <- rasch(s)
  rc <- residual_correlations(f, flag = 0.2)
  null_heur[r] <- nrow(rc$flagged) > 0
}
cat(sprintf("[Q3 null] any-pair false-alarm rate = %.3f (n=%d, MCerr=%.3f)\n",
            mean(null_heur), REPS, sqrt(mean(null_heur)*(1-mean(null_heur))/REPS)))

POWREPS <- 100
pow_flag <- logical(POWREPS); target_q3star <- numeric(POWREPS)
for (r in seq_len(POWREPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS,
                       dependence = list(pairs = list(c(5, 6)), strength = 3), seed = 22000 + r)
  f <- rasch(s)
  rc <- residual_correlations(f, flag = 0.2)
  pr <- rc$pairs
  target <- pr[(pr$item_a == "I05" & pr$item_b == "I06") | (pr$item_a == "I06" & pr$item_b == "I05"), ]
  pow_flag[r] <- isTRUE(target$flagged); target_q3star[r] <- target$q3_star
}
cat(sprintf("[Q3 power] target-pair flag rate = %.3f (n=%d, MCerr=%.3f), mean q3* = %.3f\n",
            mean(pow_flag), POWREPS, sqrt(mean(pow_flag)*(1-mean(pow_flag))/POWREPS), mean(target_q3star)))

# ---- (2) dependence_magnitude null + recovery ----
NREPS <- 100
null_d <- numeric(NREPS); null_p <- numeric(NREPS)
for (r in seq_len(NREPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS, seed = 23000 + r)
  f <- rasch(s)
  dm <- dependence_magnitude(f, dependent = "I06", independent = "I05")
  null_d[r] <- dm$d; null_p[r] <- dm$p
}
cat(sprintf("\n[dependence_magnitude null] mean d=%.4f, reject@0.05=%.3f (n=%d, MCerr=%.3f)\n",
            mean(null_d), mean(null_p < 0.05), NREPS, sqrt(0.05*0.95/NREPS)))

RREPS <- 50
ds <- numeric(RREPS); ses <- numeric(RREPS)
for (r in seq_len(RREPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS,
                       dependence = list(pairs = list(c(5, 6)), strength = 3), seed = 24000 + r)
  f <- rasch(s)
  dm <- dependence_magnitude(f, dependent = "I06", independent = "I05")
  ds[r] <- dm$d; ses[r] <- dm$se
}
cat(sprintf("[dependence_magnitude recovery, strength=3] mean d=%.3f, empSD=%.3f, meanSE=%.3f, SE-calib=%.3f\n",
            mean(ds), sd(ds), mean(ses), sd(ds)/mean(ses)))

# ---- (4) dimensionality_test null + power ----
NREPS2 <- 150
null_flag <- logical(NREPS2)
for (r in seq_len(NREPS2)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS, seed = 25000 + r)
  f <- rasch(s)
  dt <- dimensionality_test(f)
  null_flag[r] <- isTRUE(dt$multidimensional)
}
cat(sprintf("\n[dimensionality_test null] flag rate = %.3f (n=%d, MCerr=%.3f)\n",
            mean(null_flag), NREPS2, sqrt(mean(null_flag)*(1-mean(null_flag))/NREPS2)))

POWREPS2 <- 100
pow_flag2 <- logical(POWREPS2)
for (r in seq_len(POWREPS2)) {
  s <- simulate_rasch(n_persons = N, n_items = 16, missing = MISS,
                       second_dim = list(items = 9:16, rho = 0.2), seed = 26000 + r)
  f <- rasch(s)
  dt <- dimensionality_test(f)
  pow_flag2[r] <- isTRUE(dt$multidimensional)
}
cat(sprintf("[dimensionality_test power, 2nd dim rho=0.2] flag rate = %.3f (n=%d, MCerr=%.3f)\n",
            mean(pow_flag2), POWREPS2, sqrt(mean(pow_flag2)*(1-mean(pow_flag2))/POWREPS2)))

# ---- (6) spread_test + combine_items under MCAR ----
REPS6 <- 60
null_flag6 <- logical(REPS6); pow_flag6 <- logical(REPS6)
refusal6 <- 0L
for (r in seq_len(REPS6)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS, seed = 27000 + r)
  f <- rasch(s)
  fc <- tryCatch(combine_items(f, list(c("I01", "I02"))), error = function(e) e)
  if (inherits(fc, "error")) { refusal6 <- refusal6 + 1L; next }
  st <- spread_test(fc)
  null_flag6[r] <- st$dependent
}
cat(sprintf("\n[spread_test null under MCAR] flag rate = %.3f (n=%d, refusals=%d)\n",
            mean(null_flag6[seq_len(REPS6 - refusal6)]), REPS6 - refusal6, refusal6))

for (r in seq_len(REPS6)) {
  s <- simulate_rasch(n_persons = N, n_items = I, missing = MISS,
                       dependence = list(pairs = list(c(4, 5)), strength = 3), seed = 28000 + r)
  f <- rasch(s)
  fc <- tryCatch(combine_items(f, list(c("I04", "I05"))), error = function(e) e)
  if (inherits(fc, "error")) next
  st <- spread_test(fc)
  pow_flag6[r] <- st$dependent
}
cat(sprintf("[spread_test power under MCAR, planted dep I04-I05] flag rate = %.3f (n=%d)\n",
            mean(pow_flag6), REPS6))
