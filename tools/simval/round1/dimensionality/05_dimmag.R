suppressWarnings(pkgload::load_all(".", quiet=TRUE))

N <- 600
d <- rep(seq(-1.5, 1.5, length.out = 6), 2)
subs <- list(paste0("I", 1:6), paste0("I", 7:12))

gen <- function(load, seed) {
  set.seed(seed)
  common <- rnorm(N); u1 <- rnorm(N); u2 <- rnorm(N)
  X <- sapply(1:12, function(i) rbinom(N, 1,
    plogis(common + load * (if (i <= 6) u1 else u2) - d[i])))
  colnames(X) <- paste0("I", 1:12); X
}

REPS <- 40
rho_b <- numeric(REPS); A_b <- numeric(REPS)
rho_u <- numeric(REPS); A_u <- numeric(REPS)
for (r in seq_len(REPS)) {
  fb <- rasch(gen(0.9, 10000 + r))
  fu <- rasch(gen(0.0, 20000 + r))
  dmb <- dimensionality_magnitude(fb, subs)
  dmu <- dimensionality_magnitude(fu, subs)
  rho_b[r] <- dmb$table$rho[1]; A_b[r] <- dmb$table$A[1]  # PSI row
  rho_u[r] <- dmu$table$rho[1]; A_u[r] <- dmu$table$A[1]
}
cat(sprintf("Bifactor (load=0.9): mean rho=%.3f (sd %.3f), mean A=%.3f (sd %.3f)\n",
            mean(rho_b), sd(rho_b), mean(A_b), sd(A_b)))
cat(sprintf("Unidimensional control (load=0): mean rho=%.3f (sd %.3f), mean A=%.3f (sd %.3f)\n",
            mean(rho_u), sd(rho_u), mean(A_u), sd(A_u)))
cat(sprintf("Direction check: rho lower under bifactor in %d/%d reps; A lower under bifactor in %d/%d reps\n",
            sum(rho_b < rho_u), REPS, sum(A_b < A_u), REPS))
