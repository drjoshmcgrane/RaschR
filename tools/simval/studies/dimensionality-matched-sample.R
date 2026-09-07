# Numerical conformance of the matched-sample reliability decomposition.
# This checks sample handling, not unbiased recovery or Type I error.
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "dimensionality-matched-sample"
B <- 20L
scenarios <- c("complete_null", "complete_bifactor", "missing_at_random",
               "ability_related_missing", "pcm_missing")
rows <- lapply(seq_along(scenarios), function(s) {
  scenario <- scenarios[s]
  ans <- lapply(seq_len(B), function(b) {
    set.seed(917000L + 100L * s + b)
    N <- 1200L
    theta <- if (scenario == "ability_related_missing")
      c(rnorm(N / 2, sd = .6), rnorm(N / 2, sd = 3)) else rnorm(N)
    u <- if (scenario == "complete_bifactor")
      matrix(rnorm(N * 2), N, 2) * .8 else matrix(0, N, 2)
    delta <- rep(seq(-1, 1, length.out = 4), 2)
    X <- sapply(seq_len(8L), function(i) {
      th <- theta + u[, ceiling(i / 4)]
      if (scenario != "pcm_missing")
        return(rbinom(N, 1, plogis(th - delta[i])))
      tau <- delta[i] + c(-.5, .5)
      eta <- outer(th, 0:2) - matrix(c(0, cumsum(tau)), N, 3, byrow = TRUE)
      prob <- exp(eta - apply(eta, 1, max))
      vapply(seq_len(N), function(n) sample(0:2, 1, prob = prob[n, ]), 0L)
    })
    colnames(X) <- paste0("I", 1:8)
    if (scenario == "ability_related_missing") X[601:1200, c(1, 5)] <- NA
    if (scenario %in% c("missing_at_random", "pcm_missing"))
      X[matrix(runif(N * 8) < .05, N, 8)] <- NA
    fit <- rasch(X)
    z <- dimensionality_magnitude(fit, list(paste0("I", 1:4), paste0("I", 5:8)))
    p <- fit$person; q <- z$refit$person
    keep <- complete.cases(X) & complete.cases(z$refit$X) &
      is.finite(p$theta) & is.finite(p$se) & p$se >= 0 &
      is.finite(q$theta) & is.finite(q$se) & q$se >= 0
    reliability <- function(th, se) max(1 - mean(se^2) / var(th), 0)
    r <- c(reliability(p$theta[keep], p$se[keep]),
           reliability(q$theta[keep], q$se[keep]))
    expected <- if (all(r > 0)) max(2 * (r[1] / r[2] - 1) * 7 / 6, 0) else NA_real_
    stopifnot(isTRUE(all.equal(z$table$c2[1], expected)),
              z$table$n[1] == sum(keep),
              z$table$n_excluded[1] == sum(!keep),
              z$table$n[2] == sum(complete.cases(X)))
    error <- max(abs(c(z$table$run1[1], z$table$subtest[1]) - r))
    c(error = error, n = sum(keep), unavailable = is.na(expected),
      old_psi = fit$psi$PSI, matched_psi = r[1])
  })
  a <- do.call(rbind, ans)
  stopifnot(max(a[, "error"]) < 1e-12)
  row <- sv_row(STUDY, scenario, "matched-sample arithmetic", B,
    n_attempted = B, n_refused = 0L, n_nonconv = 0L, n_error = 0L,
    notes = paste("1200 response rows; eight items in two equal subscales;",
      "20 fixed seeds; checks use an independent reliability calculation;",
      "no claim about estimator bias, coverage or rejection rates"))
  row$max_absolute_reliability_error <- max(a[, "error"])
  row$minimum_matched_n <- min(a[, "n"])
  row$maximum_matched_n <- max(a[, "n"])
  row$n_magnitude_unavailable <- sum(a[, "unavailable"])
  row$mean_original_psi <- mean(a[, "old_psi"])
  row$mean_matched_psi <- mean(a[, "matched_psi"])
  row
})
sv_write(do.call(rbind, rows), STUDY)
