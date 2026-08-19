# STUDY: humphry-item-side
#
# Is the item-side variance-ratio argument (Humphry 2005, eq. 2.27) biased,
# and by how much at the sample sizes it is applied at?
#
# The argument: the same items are calibrated in two frames of reference;
# the ratio of the standard deviations of the item location estimates
# estimates the ratio of the frames' units. The estimates carry error, so
# the observed SD is inflated in both frames, and the ratio is attenuated
# toward one. The inflation does not cancel between frames, because the
# frame with the larger unit produces more deterministic responses and so
# more precise estimates.
#
# Three estimators of the unit ratio per replicate:
#   raw   sd(delta_hat_2) / sd(delta_hat_1)                  (as applied)
#   naive corrected by subtracting mean(se^2) per frame       (the person-side
#         construction this package used to use, transplanted)
#   cov   corrected by subtracting tr(V)/(K-1), the right term under the
#         sum-zero constraint that makes the item errors correlated
#
# Design sweep: persons per frame 200 to 5,000, items 8 to 40, planted unit
# ratio 1.30, independent person samples per frame (as in the Year 5 vs
# Year 7 application). Serial.
#   Rscript tools/simval/studies/humphry-item-side.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "humphry-item-side"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

R <- 30L
ratio <- 1.30
rho <- c(ratio^-0.5, ratio^0.5)
lt <- log(ratio)

for (K in c(8L, 20L, 40L)) {
  delta <- seq(-2, 2, length.out = K)
  for (N in c(200L, 500L, 1000L, 2000L, 5000L)) {
    est <- matrix(NA_real_, R, 3,
                  dimnames = list(NULL, c("raw", "naive", "cov")))
    n_ref <- 0L
    for (r in seq_len(R)) {
      set.seed(120000 + K * 10000 + N + r)
      v <- lapply(1:2, function(g) {
        th <- rnorm(N, 0, 1.3)
        X <- sapply(delta, function(d)
          rbinom(N, 1, plogis(rho[g] * (th - d))))
        colnames(X) <- sprintf("I%02d", seq_len(K))
        f <- tryCatch(rasch(X), error = function(e) NULL)
        if (is.null(f)) return(NULL)
        V <- f$est$cov_tau
        list(v_obs = var(f$items$location),
             e_naive = mean(f$items$se^2),
             e_cov = sum(diag(V)) / (K - 1))
      })
      if (any(vapply(v, is.null, TRUE))) { n_ref <- n_ref + 1L; next }
      pos <- function(z) if (z > 0) z else NA_real_
      est[r, "raw"] <- sqrt(v[[2]]$v_obs / v[[1]]$v_obs)
      est[r, "naive"] <- sqrt(pos(v[[2]]$v_obs - v[[2]]$e_naive) /
                              pos(v[[1]]$v_obs - v[[1]]$e_naive))
      est[r, "cov"] <- sqrt(pos(v[[2]]$v_obs - v[[2]]$e_cov) /
                            pos(v[[1]]$v_obs - v[[1]]$e_cov))
    }
    for (m in colnames(est)) {
      ok <- is.finite(est[, m])
      add(sprintf("%d items/frame, N = %d per frame, unit ratio 1.30", K, N),
          sprintf("log unit ratio, %s SD ratio", m), sum(ok),
          bias = mean(log(est[ok, m])) - lt,
          emp_sd = sd(log(est[ok, m])),
          effect = N, n_attempted = R, n_refused = n_ref,
          notes = sprintf("items %d", K))
    }
    cat(sprintf("[%s] K = %2d, N = %4d: raw %+.4f  naive %+.4f  cov %+.4f\n",
                format(Sys.time(), "%H:%M"), K, N,
                mean(log(est[, "raw"]), na.rm = TRUE) - lt,
                mean(log(est[, "naive"]), na.rm = TRUE) - lt,
                mean(log(est[, "cov"]), na.rm = TRUE) - lt))
  }
}

sv_write(do.call(rbind, rows), "humphry-item-side")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
