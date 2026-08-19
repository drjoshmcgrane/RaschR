# STUDY: humphry-pgd-misfit
#
# The unit ratio in Humphry (2005, ch. 4) is read off the dispersion of the
# common items' locations, so anything that moves those locations in one
# frame and not the other enters the estimate as if it were a unit
# difference. His own data show this is not hypothetical: certain items
# suggest "a combination of PGD and ISD", and after transformation the
# RMSD of 0.24 against an RMSE of 0.12 is item-level departure well beyond
# error.
#
# Two departures are planted on a subset of the 12 common items:
#   dif   a location shift in Year 7 only (uniform DIF), signs alternating
#   isd   a discrimination multiplier in Year 7 only (item-specific
#         discrimination on top of the person-group unit)
#
# Three dispersion ratios are compared, because a robust measure should
# resist a minority of aberrant items where the standard deviation cannot:
#   sd      ratio of standard deviations (the estimator as applied)
#   mad     ratio of median absolute deviations
#   trim    ratio of 20% trimmed standard deviations
#
# Planted unit ratio 1.306, N = 980 per year, ability gap 0.5 logits.
# Serial. Rscript tools/simval/studies/humphry-pgd-misfit.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "humphry-pgd-misfit"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

phi <- c(0.875, 1.143)
truth <- phi[2] / phi[1]
lt <- log(truth)
yr5_loc <- c(-0.60, -0.28, -0.68, -0.01, -0.12, -0.10,
              0.78,  0.33, -0.09, -1.30,  1.43,  0.65)
delta_common <- yr5_loc / phi[1]
K_common <- length(delta_common)
K_unique <- 18L
N <- 980L
gap <- 0.5
R <- 200L

trim_sd <- function(x, p = 0.2) {
  q <- quantile(x, c(p / 2, 1 - p / 2), names = FALSE)
  sd(x[x >= q[1] & x <= q[2]])
}

run_cell <- function(label, n_aff, kind, size) {
  m <- matrix(NA_real_, R, 3, dimnames = list(NULL, c("sd", "mad", "trim")))
  n_ref <- 0L
  aff <- if (n_aff > 0) seq_len(n_aff) else integer(0)
  sgn <- rep(c(1, -1), length.out = max(n_aff, 1))
  for (r in seq_len(R)) {
    set.seed(160000 + n_aff * 1000 + round(size * 100) + r)
    loc <- vector("list", 2)
    for (g in 1:2) {
      th <- rnorm(N, if (g == 1) 0 else gap, 1.2)
      d_uniq <- seq(-1.5, 1.5, length.out = K_unique) +
        if (g == 1) 0 else gap
      d_com <- delta_common
      disc <- rep(1, K_common)
      if (g == 2 && n_aff > 0) {
        if (kind == "dif") d_com[aff] <- d_com[aff] + sgn * size
        if (kind == "isd") disc[aff] <- size
      }
      X_com <- vapply(seq_len(K_common), function(i)
        rbinom(N, 1, plogis(phi[g] * disc[i] * (th - d_com[i]))), numeric(N))
      X_uni <- vapply(d_uniq, function(d)
        rbinom(N, 1, plogis(phi[g] * (th - d))), numeric(N))
      X <- cbind(X_com, X_uni)
      colnames(X) <- c(sprintf("C%02d", seq_len(K_common)),
                       sprintf("U%02d", seq_len(K_unique)))
      f <- tryCatch(rasch(X), error = function(e) NULL)
      if (is.null(f)) { loc[[g]] <- NULL; break }
      loc[[g]] <- f$items$location[seq_len(K_common)]
    }
    if (any(vapply(loc, is.null, TRUE))) { n_ref <- n_ref + 1L; next }
    a <- loc[[1]]; b <- loc[[2]]
    m[r, "sd"] <- log(sd(b) / sd(a))
    m[r, "mad"] <- log(mad(b) / mad(a))
    m[r, "trim"] <- log(trim_sd(b) / trim_sd(a))
  }
  for (est in colnames(m)) {
    ok <- is.finite(m[, est])
    add(label, sprintf("log unit ratio, %s dispersion ratio", est), sum(ok),
        bias = mean(m[ok, est]) - lt, emp_sd = sd(m[ok, est]),
        effect = size, n_attempted = R, n_refused = n_ref,
        notes = sprintf("planted ratio %.3f; recovered %.3f", truth,
                        exp(mean(m[ok, est]))))
  }
  cat(sprintf("[%s] %-44s sd %+.4f (%.3f) | mad %+.4f (%.3f) | trim %+.4f (%.3f)\n",
              format(Sys.time(), "%H:%M"), label,
              mean(m[, "sd"], na.rm = TRUE) - lt, exp(mean(m[, "sd"], na.rm = TRUE)),
              mean(m[, "mad"], na.rm = TRUE) - lt, exp(mean(m[, "mad"], na.rm = TRUE)),
              mean(m[, "trim"], na.rm = TRUE) - lt, exp(mean(m[, "trim"], na.rm = TRUE))))
}

run_cell("clean (no item-level departure)", 0L, "none", 0)
run_cell("2 of 12 items, uniform DIF 0.5", 2L, "dif", 0.5)
run_cell("4 of 12 items, uniform DIF 0.5", 4L, "dif", 0.5)
run_cell("2 of 12 items, uniform DIF 1.0", 2L, "dif", 1.0)
run_cell("2 of 12 items, discrimination x1.5", 2L, "isd", 1.5)
run_cell("4 of 12 items, discrimination x1.5", 4L, "isd", 1.5)

sv_write(do.call(rbind, rows), "humphry-pgd-misfit")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
