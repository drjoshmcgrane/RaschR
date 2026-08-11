# At what adjacent-category count does the joint threshold block's sandwich
# SE degrade for the item's OTHER thresholds? PCM, 6 items, m=3; item 3's
# top category made progressively rare by raising its last threshold.
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
N <- 500
simP <- function(t, tau) { x <- 0:length(tau); p <- exp(x*t - c(0, cumsum(tau))); p / sum(p) }
base <- list(c(-1.0, 0, 1.0), c(-1.2, -0.1, 0.9), c(-0.8, 0.1, 1.1),
             c(-1.1, 0, 1.2), c(-0.9, -0.2, 0.8), c(-1.0, 0.2, 1.0))
for (miss in c(0, 0.3)) for (shift3 in c(2.2, 2.8, 3.4, 4.0)) {
  tt <- base; tt[[3]][3] <- tt[[3]][3] + shift3   # rarify item 3's top category
  NR <- 400
  est <- ses <- array(NA_real_, c(NR, 3))         # item 3's three thresholds
  cnt_top <- rep(NA_real_, NR); flagged <- dropped <- rep(FALSE, NR)
  set.seed(9.9e6 + shift3 * 1e3)
  for (r in seq_len(NR)) {
    th <- rnorm(N, 0, 1.2)
    X <- sapply(seq_along(tt), function(i)
      vapply(th, function(t) sample(0:3, 1, prob = simP(t, tt[[i]])), 0L))
    colnames(X) <- sprintf("P%d", 1:6)
    if (miss > 0) X[matrix(runif(length(X)) < miss, nrow(X))] <- NA
    f <- tryCatch(rasch(X, model = "PCM"), error = function(e) NULL)
    if (is.null(f)) next
    thr <- f$thresholds
    r3 <- thr[thr$item == 3, ]
    for (k in 1:3) {                       # match by threshold index: a
      row <- r3[r3$k == k, ]               # zero-count top category drops
      if (nrow(row) == 1) {                # threshold 3 entirely (the item
        est[r, k] <- row$tau               # re-scores at its observed max)
        ses[r, k] <- row$se
      }
    }
    dropped[r] <- nrow(r3) < 3
    cnt_top[r] <- sum(X[, 3] == 3, na.rm = TRUE)
    flagged[r] <- any(is.na(r3$se))
  }
  ratios <- vapply(1:3, function(k) {
    fin <- is.finite(est[, k]) & is.finite(ses[, k])
    if (sum(fin) < 30) return(NA_real_)
    sd(est[fin, k]) / mean(ses[fin, k])
  }, 0)
  cat(sprintf("miss=%.1f shift=%.1f: top-cat count=%.1f, guard-flag=%.2f, dropped=%.2f, se_ratio k1/k2/k3 = %.2f / %.2f / %.2f (n=%d)\n",
      miss, shift3, mean(cnt_top, na.rm = TRUE), mean(flagged), mean(dropped),
      ratios[1], ratios[2], ratios[3], sum(is.finite(est[,1]))))
}
