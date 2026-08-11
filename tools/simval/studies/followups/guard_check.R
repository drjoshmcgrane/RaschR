suppressWarnings(pkgload::load_all(".", quiet = TRUE))
K <- 8; RPP <- 40
objs <- sprintf("O%d", seq_len(K))
beta <- setNames(seq(-1, 1, length.out = K), objs)
pr <- t(utils::combn(objs, 2))
gen <- function(J, share, seed) {
  set.seed(seed)
  jids <- sprintf("J%d", seq_len(J))
  d <- data.frame(object_a = rep(pr[, 1], each = RPP),
                  object_b = rep(pr[, 2], each = RPP), stringsAsFactors = FALSE)
  prob <- if (is.na(share)) rep(1/J, J) else c(share, rep((1 - share)/(J - 1), J - 1))
  d$judge <- sample(jids, nrow(d), replace = TRUE, prob = prob)
  lp <- beta[d$object_a] - beta[d$object_b]
  d$winner <- ifelse(rbinom(nrow(d), 1, plogis(lp)) == 1, d$object_a, d$object_b)
  d
}
for (cfg in list(list("balanced J=10", 10, NA), list("share=0.25 J=10", 10, 0.25),
                 list("share=0.50 J=20", 20, 0.5), list("J=9 balanced", 9, NA))) {
  f <- btl(gen(cfg[[2]], cfg[[3]], 42), "object_a", "object_b",
           winner = "winner", judge = "judge")
  cat(sprintf("%-16s: inference=%s nc_eff=%.1f se[1]=%s\n", cfg[[1]],
      f$cl$inference_available, f$cl$n_units_effective,
      format(round(f$objects$se[1], 3))))
  cn <- grep("effective|uneven|concentrat", f$notes, value = TRUE)
  if (length(cn)) cat("   note: ", substr(cn[1], 1, 100), "\n")
}
