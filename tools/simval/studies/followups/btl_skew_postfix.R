suppressWarnings(pkgload::load_all(".", quiet = TRUE))
K <- 8; RPP <- 40
objs <- sprintf("O%d", seq_len(K)); beta <- setNames(seq(-1, 1, length.out = K), objs)
pr <- t(utils::combn(objs, 2))
gen <- function(J, share) {
  jids <- sprintf("J%d", seq_len(J))
  d <- data.frame(object_a = rep(pr[,1], each = RPP), object_b = rep(pr[,2], each = RPP))
  prob <- if (is.na(share)) rep(1/J, J) else c(share, rep((1-share)/(J-1), J-1))
  d$judge <- sample(jids, nrow(d), TRUE, prob)
  lp <- beta[d$object_a] - beta[d$object_b]
  d$winner <- ifelse(rbinom(nrow(d), 1, plogis(lp)) == 1, d$object_a, d$object_b)
  d
}
for (cfg in list(c(10, 0.5), c(50, 0.5), c(10, 0.35), c(10, NA))) {
  J <- cfg[1]; sh <- cfg[2]; NR <- 200
  wh <- note <- 0L
  set.seed(23e6 + J * 100)
  for (r in seq_len(NR)) {
    f <- btl(gen(J, sh), "object_a", "object_b", winner = "winner", judge = "judge")
    if (!f$cl$inference_available) wh <- wh + 1L
    else if (any(grepl("uneven", f$notes))) note <- note + 1L
  }
  cat(sprintf("J=%2d share=%s: withheld=%.2f cautioned=%.2f (n=%d)\n",
      J, ifelse(is.na(sh), "bal", sprintf("%.2f", sh)), wh/NR, note/NR, NR))
}
