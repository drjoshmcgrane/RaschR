# Verify design information against numerical likelihood curvature.
# This is an algorithm-conformance study, not an error-rate calibration.
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "btl-information-context"
beta <- setNames(seq(-1, 1, length.out = 5L), LETTERS[1:5])
run_cell <- function(m, position) {
  ans <- lapply(seq_len(20L), function(r) {
    set.seed(994000L + r + 100L * m)
    a <- sample(names(beta), 600L, TRUE)
    b <- vapply(a, function(z) sample(setdiff(names(beta), z), 1L), "")
    tau <- if (m == 1L) 0 else c(-1.2, 0, 1.2)
    d <- data.frame(a, b)
    d$response <- vapply(beta[a] - beta[b] + position, function(z)
      sample(0:m, 1L, prob = item_moments(z, tau)$P), 0L)
    fit <- btl(d, "a", "b", response = "response", position = TRUE)
    if (!isTRUE(fit$converged)) return(c(error = NA, old_relative = NA))
    info <- btl_information(fit)
    P <- fit$fitted_prob
    h <- 1e-4
    # Exponential tilting changes the location gap by +/-h. The observed
    # response term is linear and drops out of this second derivative.
    curvature <- (log(drop(P %*% exp(h * (0:m)))) +
                   log(drop(P %*% exp(-h * (0:m))))) / h^2
    error <- max(abs(curvature - info$comparisons$information))
    old <- .btl_info_of_d(info$comparisons$gap, m,
                          if (m > 1L) fit$thresholds$tau else NULL)
    c(error = error, old_relative = sum(old) / info$total - 1)
  })
  z <- do.call(rbind, ans); ok <- complete.cases(z)
  stopifnot(all(ok), max(z[, "error"]) < 1e-6)
  row <- sv_row(STUDY, sprintf("%d categories; position=%.1f", m + 1L, position),
    "agreement with numerical per-comparison likelihood curvature", sum(ok),
    n_attempted = length(ans), n_refused = 0L, n_nonconv = sum(!ok), n_error = 0L,
    notes = "fixed object locations and thresholds; 20 response draws per condition; no inference-calibration claim")
  row$max_absolute_curvature_error <- max(z[, "error"])
  row$mean_relative_error_without_context <- mean(z[, "old_relative"])
  row
}
rows <- lapply(c(1L, 3L), function(m)
  do.call(rbind, lapply(c(-1.2, 0, 1.2), function(position) run_cell(m, position))))
sv_write(do.call(rbind, rows), STUDY)
