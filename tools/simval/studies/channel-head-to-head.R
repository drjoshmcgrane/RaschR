# STUDY: channel-head-to-head
#
# Humphry's estimator and ours on identical data.
#
#   item-side   ratio of the SDs of the item locations calibrated in each
#               frame (Humphry 2005, eq. 2.27)
#   person-side ratio of the true-score variances of the person estimates
#               in each frame, with the truncated-score-moment correction
#               this package uses for item-set units
#
# The two are never both available in a real design: his needs the same
# items in both frames, ours needs the same persons, and item sets
# partition the items. Simulating a design with both -- same persons, same
# items, two frames differing only in unit, responses conditionally
# independent given theta -- puts the channels on the same data so the
# comparison is of information, not of design.
#
# Planted unit ratio 1.30, dichotomous.
# Serial. Rscript tools/simval/studies/channel-head-to-head.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "channel-head-to-head"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

ratio <- 1.30
rho <- c(ratio^-0.5, ratio^0.5)
lt <- log(ratio)
R <- 100L

for (K in c(8L, 12L, 20L, 40L)) {
  delta <- seq(-2, 2, length.out = K)
  for (N in c(250L, 980L)) {
    m <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("item", "person")))
    n_ref <- 0L
    for (r in seq_len(R)) {
      set.seed(170000 + K * 5000 + N + r)
      th <- rnorm(N, 0, 1.3)
      fits <- lapply(1:2, function(g) {
        X <- vapply(delta, function(d)
          rbinom(N, 1, plogis(rho[g] * (th - d))), numeric(N))
        colnames(X) <- sprintf("I%02d", seq_len(K))
        tryCatch(rasch(X), error = function(e) NULL)
      })
      if (any(vapply(fits, is.null, TRUE))) { n_ref <- n_ref + 1L; next }
      # item-side: dispersion of item locations, one calibration per frame
      m[r, "item"] <- log(sd(fits[[2]]$items$location) /
                          sd(fits[[1]]$items$location))
      # person-side: corrected true-score variance ratio over common persons
      ok <- !fits[[1]]$person$extreme & !fits[[2]]$person$extreme
      if (sum(ok) < 30) { n_ref <- n_ref + 1L; next }
      v <- vapply(1:2, function(g) {
        f <- fits[[g]]
        uh <- f$person$theta[ok]
        lm <- .person_link_moments(as.matrix(f$X), f$tau_list)
        (var(uh) - mean(lm$w[ok])) / mean(lm$g[ok])^2
      }, 0)
      if (any(!is.finite(v)) || any(v <= 0)) { n_ref <- n_ref + 1L; next }
      m[r, "person"] <- 0.5 * (log(v[2]) - log(v[1]))
    }
    for (ch in colnames(m)) {
      ok <- is.finite(m[, ch])
      add(sprintf("%d items, N = %d, unit ratio 1.30", K, N),
          sprintf("log unit ratio, %s channel", ch), sum(ok),
          bias = mean(m[ok, ch]) - lt, emp_sd = sd(m[ok, ch]),
          effect = N, n_attempted = R, n_refused = n_ref,
          notes = sprintf("items %d", K))
    }
    oi <- is.finite(m[, "item"]); op <- is.finite(m[, "person"])
    cat(sprintf("[%s] K = %2d, N = %3d: item %+.4f (sd %.4f) | person %+.4f (sd %.4f) | SD ratio %.2fx\n",
                format(Sys.time(), "%H:%M"), K, N,
                mean(m[oi, "item"]) - lt, sd(m[oi, "item"]),
                mean(m[op, "person"]) - lt, sd(m[op, "person"]),
                sd(m[op, "person"]) / sd(m[oi, "item"])))
  }
}

sv_write(do.call(rbind, rows), "channel-head-to-head")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
