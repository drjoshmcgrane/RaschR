suppressWarnings(pkgload::load_all(".", quiet = TRUE))
I <- 10L; m <- 3L; N <- 500L; NR <- 300L
delta <- seq(-1.2, 1.2, length.out = I); names(delta) <- sprintf("I%02d", 1:I)
step_rare <- c(-1.8, 0.2, 4.2)
taus <- lapply(seq_len(I), function(i) .sim_thresholds(delta[i], m, 1, FALSE, pattern = step_rare))
est <- se <- array(NA_real_, c(NR, I, m))
set.seed(21e6)
for (r in seq_len(NR)) {
  theta <- rnorm(N)
  X <- sapply(seq_len(I), function(j) .sim_item(theta, taus[[j]]))
  colnames(X) <- names(delta)
  X[matrix(runif(N * I) < 0.30, N, I)] <- NA
  f <- tryCatch(suppressWarnings(rasch(X, model = "PCM")), error = function(e) NULL)
  if (is.null(f) || !isTRUE(f$est$converged)) next
  th <- f$thresholds
  for (q in seq_len(nrow(th))) { est[r, th$item[q], th$k[q]] <- th$tau[q]
                                 se[r, th$item[q], th$k[q]] <- th$se[q] }
}
cat("per-item se_ratio (empirical SD / mean reported SE), fresh seeds, their exact design:\n")
cat("item   k1     k2     k3\n")
pool_sd <- pool_se <- matrix(NA_real_, I, m)
for (j in seq_len(I)) {
  rr <- vapply(seq_len(m), function(k) {
    ok <- is.finite(est[, j, k]) & is.finite(se[, j, k])
    if (sum(ok) < 50) return(NA_real_)
    pool_sd[j, k] <<- sd(est[ok, j, k]); pool_se[j, k] <<- mean(se[ok, j, k])
    sd(est[ok, j, k]) / mean(se[ok, j, k])
  }, 0)
  cat(sprintf("%-5s %5.2f  %5.2f  %5.2f\n", names(delta)[j], rr[1], rr[2], rr[3]))
}
for (k in seq_len(m)) {
  ok <- is.finite(pool_sd[, k])
  cat(sprintf("k=%d POOLED-style ratio sqrt(mean sd^2)/mean se = %.2f  vs per-item median ratio = %.2f\n",
      k, sqrt(mean(pool_sd[ok, k]^2)) / mean(pool_se[ok, k]),
      median(pool_sd[ok, k] / pool_se[ok, k])))
}
