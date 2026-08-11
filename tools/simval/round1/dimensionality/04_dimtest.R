suppressWarnings(pkgload::load_all(".", quiet=TRUE))

# ---- NULL: unidimensional, short test (8 dichotomous items -> caution expected) ----
NREPS <- 300
N <- 500; I <- 8
null_flag <- logical(NREPS); null_caution <- logical(NREPS); null_na <- logical(NREPS)
for (r in seq_len(NREPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, seed = 8000 + r)
  f <- rasch(s)
  dt <- dimensionality_test(f)
  null_na[r] <- is.na(dt$multidimensional)
  null_flag[r] <- isTRUE(dt$multidimensional)
  null_caution[r] <- !is.null(dt$caution)
}
cat("=== NULL (unidimensional, I=8 -> short subtests, alpha=0.05) ===\n")
cat(sprintf("flag rate (multidimensional=TRUE) = %.3f (n=%d, MC err=%.3f); NA/withheld rate = %.3f; caution fired = %.3f\n",
            mean(null_flag), NREPS, sqrt(mean(null_flag)*(1-mean(null_flag))/NREPS),
            mean(null_na), mean(null_caution)))

# ---- POWER: planted second dimension, longer test so subtests exceed 15 points (no caution) ----
N2 <- 500; I2 <- 16
POWREPS <- 200
pow_flag <- logical(POWREPS); pow_caution <- logical(POWREPS); pow_na <- logical(POWREPS)
for (r in seq_len(POWREPS)) {
  s <- simulate_rasch(n_persons = N2, n_items = I2,
                       second_dim = list(items = 9:16, rho = 0.2), seed = 8500 + r)
  f <- rasch(s)
  dt <- dimensionality_test(f)
  pow_na[r] <- is.na(dt$multidimensional %||% NA)
  pow_flag[r] <- isTRUE(dt$multidimensional)
  pow_caution[r] <- !is.null(dt$caution)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
cat("\n=== POWER (second_dim rho=0.2 on half the items, I=16 -> longer subtests) ===\n")
cat(sprintf("flag rate (multidimensional=TRUE) = %.3f (n=%d, MC err=%.3f); NA/withheld rate = %.3f; caution fired = %.3f\n",
            mean(pow_flag), POWREPS, sqrt(mean(pow_flag)*(1-mean(pow_flag))/POWREPS),
            mean(pow_na), mean(pow_caution)))

# ---- caution fires for short subtests, verdict still returned (single illustrative case) ----
s <- simulate_rasch(n_persons = 500, n_items = 8, seed = 1)
f <- rasch(s)
dt <- dimensionality_test(f)
cat("\n=== Caution-firing illustrative check (I=8 dichotomous, ~4 score points/subtest) ===\n")
cat(sprintf("caution present: %s ; verdict (multidimensional) still returned: %s (value=%s)\n",
            !is.null(dt$caution), !is.na(dt$multidimensional), dt$multidimensional))
cat("caution text:", dt$caution, "\n")
