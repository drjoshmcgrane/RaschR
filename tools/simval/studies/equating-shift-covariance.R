# Fixed-weight covariance conformance, not coverage or Type I error.
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "equating-shift-covariance"
B <- 20L
fit <- function(n, seed) {
  d <- simulate_btl(n, 30, reps_per_pair = 5, seed = seed)
  btl(d, "object_a", "object_b", "winner", judge = "judge")
}
checks <- lapply(seq_len(B), function(r) {
  a <- fit(6, 740000L + r)
  b <- fit(7, 741000L + r)
  stopifnot(a$converged, b$converged)
  z <- btl_equate(a, b, independent = TRUE)
  stopifnot(z$inferential)
  common <- z$table$object
  u <- 1 / (z$table$se_1^2 + z$table$se_2^2)
  u <- u / sum(u)
  A <- matrix(0, nrow(b$objects), nrow(a$objects))
  A[, match(common, a$objects$object)] <- matrix(u, nrow(A), length(u), byrow = TRUE)
  H <- diag(nrow(b$objects))
  H[, match(common, b$objects$object)] <-
    H[, match(common, b$objects$object)] - matrix(u, nrow(H), length(u), byrow = TRUE)
  expected <- A %*% a$cov_beta %*% t(A) + H %*% b$cov_beta %*% t(H)
  err_joint <- max(abs(attr(z$equated, "cov_location") - expected),
                   abs(z$equated$se^2 - diag(expected)))
  bank <- data.frame(object = c("O1", "O2", "O3", "Onew"),
                     location = c(-1, 0, 1, 2), se = 0)
  fixed <- btl_equate(a, bank)
  stopifnot(fixed$inferential, fixed$shift_se > 0)
  err_fixed <- max(abs(fixed$equated$se - fixed$shift_se),
                   abs(attr(fixed$equated, "cov_location") - fixed$shift_se^2))
  c(joint = err_joint, fixed = err_fixed)
})
checks <- do.call(rbind, checks)
stopifnot(all(is.finite(checks)), max(checks) < 1e-10)
rows <- lapply(colnames(checks), function(s) {
  row <- sv_row(STUDY, s, "equated covariance linear-transform agreement", B,
    n_attempted = B, n_refused = 0L, n_nonconv = 0L, n_error = 0L,
    notes = paste("30 judges per calibration; six and seven objects;",
      "joint case includes a non-common object and negative location-shift",
      "covariance; fixed bank inherits shared shift uncertainty;",
      "weights held at their fitted values; conformance only"))
  row$maximum_absolute_error <- max(checks[, s])
  row
})
sv_write(do.call(rbind, rows), STUDY)
