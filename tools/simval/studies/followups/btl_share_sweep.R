suppressWarnings(pkgload::load_all(".", quiet = TRUE))
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
    if (is.null(f) || !all(is.finite(f$objects$se))) next
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
  cat(sprintf("J=%2d share=%s: type1=%.4f (per-rep MCSE %.4f) mean nc_eff=%.1f n=%d\n",
      J, ifelse(is.na(share), "bal ", sprintf("%.2f", share)),
      mean(rej[ok]), sd(rej[ok]) / sqrt(sum(ok)), mean(eff[ok]), sum(ok)))
}
