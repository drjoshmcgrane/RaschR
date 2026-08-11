suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
# the guard this sweep calibrates would withhold SEs exactly in the 4-7
# effective-judge region whose inflation is the evidence; the
# simulation-only override keeps the covariance observable (see the guard
# comment in R/btl.R)
old_opt <- options(rasch.btl_guard_override = TRUE)
on.exit(options(old_opt), add = TRUE)
rows <- list()
K <- 8; RPP <- 40
objs <- sprintf("O%d", seq_len(K))
beta <- setNames(seq(-1, 1, length.out = K), objs)
pr <- t(utils::combn(objs, 2))
gen <- function(J, share) {
  jids <- sprintf("J%d", seq_len(J))
  d <- data.frame(object_a = rep(pr[, 1], each = RPP),
                  object_b = rep(pr[, 2], each = RPP), stringsAsFactors = FALSE)
  n <- nrow(d)
  prob <- if (is.na(share)) rep(1 / J, J) else c(share, rep((1 - share) / (J - 1), J - 1))
  d$judge <- sample(jids, n, replace = TRUE, prob = prob)
  lp <- beta[d$object_a] - beta[d$object_b]
  d$winner <- ifelse(rbinom(n, 1, plogis(lp)) == 1, d$object_a, d$object_b)
  d
}
cells <- expand.grid(J = c(10, 20), share = c(NA, 0.15, 0.25, 0.35, 0.5))
NR <- 500
for (ci in seq_len(nrow(cells))) {
  J <- cells$J[ci]; share <- cells$share[ci]
  rej <- eff <- rep(NA_real_, NR)
  set.seed(8.8e6 + ci * 1e4)
  for (r in seq_len(NR)) {
    d <- gen(J, share)
    f <- tryCatch(btl(d, "object_a", "object_b", winner = "winner", judge = "judge"),
                  error = function(e) NULL)
    if (is.null(f) || anyNA(f$objects$se)) next   # override keeps SEs finite
    est <- setNames(f$objects$location, f$objects$object)[objs]
    V <- f$cov_beta; rownames(V) <- colnames(V) <- f$objects$object
    z <- vapply(seq_len(nrow(pr)), function(e) {
      a <- pr[e, 1]; b <- pr[e, 2]
      dd <- (est[a] - est[b]) - (beta[a] - beta[b])
      dd / sqrt(V[a, a] + V[b, b] - 2 * V[a, b])
    }, 0)
    rej[r] <- mean(abs(z) > qt(0.975, J - 1))   # the package's t reference, df = clusters - 1
    sh <- tapply(rep(1, nrow(d)), d$judge, sum); sh <- sh / sum(sh)
    eff[r] <- 1 / sum(sh^2)
  }
  ok <- is.finite(rej)
  lab <- sprintf("J=%d share=%s", J, ifelse(is.na(share), "balanced", sprintf("%.2f", share)))
  cat(sprintf("%s: type1=%.4f (per-rep MCSE %.4f) mean nc_eff=%.1f n=%d\n",
      lab, mean(rej[ok]), sd(rej[ok]) / sqrt(sum(ok)), mean(eff[ok]), sum(ok)))
  rows[[length(rows) + 1]] <- sv_row("btl-share-sweep", lab,
    "type1 of pairwise contrast t-tests (guard overridden to observe the withheld region)",
    sum(ok), type1 = mean(rej[ok]),
    mc_override = list(type1 = sd(rej[ok]) / sqrt(sum(ok))),
    n_attempted = NR, n_refused = NR - sum(ok),
    notes = sprintf("mean effective clusters %.1f; per-replicate proportion over 28 contrasts; reference qt(.975, J-1)", mean(eff[ok])))
}
sv_write(do.call(rbind, rows), "btl-share-sweep")
