# STUDY: humphry-isd-replication
#
# Humphry's ITEM-SET discrimination study is estimated on the person side,
# with the construction this package replaced. From the thesis: "ISD for
# each item set was estimated using a matrix of log ratios of standard
# deviations for common persons across the sets", with "the effects of
# measurement error on variance removed using Equation (2.29)" -- that is,
# var(WLE) minus the mean squared standard error (Andrich 1982).
#
# His reported recovery is good (mean absolute difference 0.010 against
# planted ISDs), which does not obviously square with the +5% bias this
# package measured for the same construction at 8 dichotomous items per
# set. His design is replicated exactly to settle it.
#
# Design (Tables 3.9a, 3.9b, 3.10): 40 items from -4.0 to 4.0 in steps of
# 0.2 omitting 0, interleaved into 4 sets of 10 so each set spans the full
# range; one person group, N = 1000, ability SD 1.51; ISDs 0.604, 0.906,
# 1.209, 1.511 (ratios 2:3:4:5, product 1).
#
# Reported per-set behaviour worth reproducing: his correction
# under-recovers each set's ability SD by about 7 per cent (set 1 expected
# 0.91, estimated 0.84; set 4 expected 2.29, estimated 2.14), yet the
# ratios survive because the per-set distortions largely cancel.
#
# Three estimators, all on the log-ratio matrix he describes:
#   humphry    var(WLE) - mean(se^2)                  (equation 2.29)
#   raw        var(WLE), uncorrected
#   corrected  truncated-score-moment correction      (what this package ships)
# Serial. Rscript tools/simval/studies/humphry-isd-replication.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "humphry-isd-replication"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

alpha <- c(0.604, 0.906, 1.209, 1.511)
S <- length(alpha)
la_true <- log(alpha)
loc_all <- setdiff(seq(-4, 4, by = 0.2), 0)          # 40 locations
set_of <- rep(seq_len(S), length.out = length(loc_all))   # interleaved
N <- 1000L
theta_sd <- 1.51
R <- 100L

# solve the full matrix of pairwise log-ratios by least squares under
# sum(log alpha) = 0, as the thesis describes
solve_matrix <- function(v) {
  A <- matrix(0, 0, S); y <- numeric(0)
  for (a in 1:(S - 1)) for (b in (a + 1):S) {
    row <- numeric(S); row[b] <- 1; row[a] <- -1
    A <- rbind(A, row); y <- c(y, 0.5 * (log(v[b]) - log(v[a])))
  }
  A <- rbind(A, rep(1, S)); y <- c(y, 0)
  drop(qr.solve(A, y))
}

est <- array(NA_real_, c(R, S, 3),
             dimnames = list(NULL, NULL, c("humphry", "raw", "corrected")))
per_set <- matrix(NA_real_, R, S)      # his Table 3.10 check
n_ref <- 0L
for (r in seq_len(R)) {
  set.seed(190000 + r)
  th <- rnorm(N, 0, theta_sd)
  v <- matrix(NA_real_, S, 3)
  okall <- rep(TRUE, N)
  fits <- vector("list", S)
  for (s in seq_len(S)) {
    d <- loc_all[set_of == s]
    X <- vapply(d, function(dd)
      rbinom(N, 1, plogis(alpha[s] * (th - dd))), numeric(N))
    colnames(X) <- sprintf("S%dI%02d", s, seq_along(d))
    f <- tryCatch(rasch(X), error = function(e) NULL)
    if (is.null(f)) break
    fits[[s]] <- f
    okall <- okall & !f$person$extreme
  }
  if (any(vapply(fits, is.null, TRUE))) { n_ref <- n_ref + 1L; next }
  for (s in seq_len(S)) {
    f <- fits[[s]]
    uh <- f$person$theta[okall]; seh <- f$person$se[okall]
    lm <- .person_link_moments(as.matrix(f$X), f$tau_list)
    v[s, 1] <- var(uh) - mean(seh^2)
    v[s, 2] <- var(uh)
    v[s, 3] <- (var(uh) - mean(lm$w[okall])) / mean(lm$g[okall])^2
    per_set[r, s] <- sqrt(max(v[s, 1], 0))
  }
  if (any(v[, 1] <= 0)) { n_ref <- n_ref + 1L; next }
  for (k in 1:3) est[r, , k] <- exp(solve_matrix(v[, k]))
}

for (k in dimnames(est)[[3]]) {
  m <- est[, , k]
  ok <- stats::complete.cases(m)
  bias_log <- colMeans(log(m[ok, , drop = FALSE])) - la_true
  add(sprintf("Humphry ISD design: 4 sets x 10 items, N = 1000"),
      sprintf("log ISD bias, %s correction", k), sum(ok),
      bias = mean(bias_log), emp_sd = mean(apply(log(m[ok, , drop = FALSE]), 2, sd)),
      n_attempted = R, n_refused = n_ref,
      notes = sprintf("recovered %s vs planted %s; end-to-end ratio %.3f (planted %.3f)",
        paste(sprintf("%.3f", colMeans(m[ok, , drop = FALSE])), collapse = "/"),
        paste(sprintf("%.3f", alpha), collapse = "/"),
        mean(m[ok, S] / m[ok, 1]), alpha[S] / alpha[1]))
  cat(sprintf("[%s] %-10s recovered %s | mean log bias %+.4f | end ratio %.3f vs %.3f\n",
      format(Sys.time(), "%H:%M"), k,
      paste(sprintf("%.3f", colMeans(m[ok, , drop = FALSE])), collapse = " "),
      mean(bias_log), mean(m[ok, S] / m[ok, 1]), alpha[S] / alpha[1]))
}
cat(sprintf("[check] per-set SD under equation 2.29: %s (his Table 3.10: 0.84 1.30 1.74 2.14)\n",
    paste(sprintf("%.2f", colMeans(per_set, na.rm = TRUE)), collapse = " ")))

sv_write(do.call(rbind, rows), "humphry-isd-replication")
cat(sprintf("TOTAL elapsed: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
