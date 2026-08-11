suppressWarnings(pkgload::load_all(".", quiet = TRUE))
I <- 10L; m <- 3L; N <- 500L; NR <- 300L
delta <- seq(-1.2, 1.2, length.out = I); names(delta) <- sprintf("I%02d", 1:I)
taus <- lapply(seq_len(I), function(i) .sim_thresholds(delta[i], m, 1, FALSE,
                                                       pattern = c(-1.8, 0.2, 4.2)))
est <- se <- array(NA_real_, c(NR, I, m)); withheld <- matrix(0, I, m)
set.seed(22e6)
for (r in seq_len(NR)) {
  theta <- rnorm(N)
  X <- sapply(seq_len(I), function(j) .sim_item(theta, taus[[j]]))
  colnames(X) <- names(delta)
  X[matrix(runif(N * I) < 0.30, N, I)] <- NA
  f <- tryCatch(suppressWarnings(rasch(X, model = "PCM")), error = function(e) NULL)
  if (is.null(f) || !isTRUE(f$est$converged)) next
  th <- f$thresholds
  for (q in seq_len(nrow(th))) {
    est[r, th$item[q], th$k[q]] <- th$tau[q]
    se[r, th$item[q], th$k[q]] <- th$se[q]
    if (is.na(th$se[q])) withheld[th$item[q], th$k[q]] <- withheld[th$item[q], th$k[q]] + 1
  }
}
cat("post-fix: item, k2 withheld-rate, k2 se_ratio among REPORTED cells\n")
for (j in c(1, 2, 3, 10)) {
  ok <- is.finite(est[, j, 2]) & is.finite(se[, j, 2])
  rr <- if (sum(ok) >= 50) sd(est[ok, j, 2]) / mean(se[ok, j, 2]) else NA
  cat(sprintf("%-5s withheld=%.2f  reported se_ratio=%.2f (n=%d)\n",
      names(delta)[j], withheld[j, 2] / NR, rr, sum(ok)))
}
