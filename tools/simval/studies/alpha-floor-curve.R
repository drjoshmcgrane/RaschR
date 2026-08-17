# The fixed-length bias floor of the corrected (and raw) EFRM unit
# estimator as a function of items per set. N = 20,000 per replicate
# suppresses sampling noise so the mean over 25 replicates estimates the
# probability-limit offset to ~0.0016. Ratio 1.4 throughout.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
al <- c(1.4^-0.5, 1.4^0.5); lt <- log(1.4)
run_len <- function(ips, R, gen_theta, seedf) {
  res <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("corr","raw")))
  delta <- seq(-1.5, 1.5, length.out = ips)
  for (r in seq_len(R)) {
    set.seed(seedf + r); n <- 20000
    th <- gen_theta(n)
    v_c <- v_r <- numeric(2); okm <- rep(TRUE, n); fits <- vector("list", 2)
    for (k in 1:2) {
      X <- sapply(delta, function(d) rbinom(n, 1, plogis(al[k] * (th - d))))
      colnames(X) <- sprintf("I%02d", seq_len(ips))
      f <- rasch(X); fits[[k]] <- f
      okm <- okm & !f$person$extreme
    }
    for (k in 1:2) {
      f <- fits[[k]]
      uh <- f$person$theta[okm]; seh <- f$person$se[okm]
      lm <- .person_link_moments(as.matrix(f$X), f$tau_list)
      v_c[k] <- (var(uh) - mean(lm$w[okm])) / mean(lm$g[okm])^2
      v_r[k] <- var(uh) - mean(seh^2)
    }
    res[r, ] <- c(0.5 * (log(v_c[2]) - log(v_c[1])) - lt,
                  0.5 * (log(v_r[2]) - log(v_r[1])) - lt)
  }
  res
}
cat("== normal persons ==\n")
for (ips in c(4, 6, 8, 12, 16, 24, 32)) {
  z <- run_len(ips, 25, function(n) rnorm(n, 0, 1.3), 70e3 + ips * 100)
  cat(sprintf("items %2d: corrected %+.4f (SE %.4f) | raw %+.4f\n", ips,
      mean(z[,1]), sd(z[,1])/sqrt(nrow(z)), mean(z[,2])))
}
cat("== skewed persons (chi-sq 3) ==\n")
for (ips in c(4, 8, 16, 32)) {
  z <- run_len(ips, 25, function(n) (rchisq(n, 3) - 3)/sqrt(6) * 1.3,
               75e3 + ips * 100)
  cat(sprintf("items %2d: corrected %+.4f (SE %.4f) | raw %+.4f\n", ips,
      mean(z[,1]), sd(z[,1])/sqrt(nrow(z)), mean(z[,2])))
}
