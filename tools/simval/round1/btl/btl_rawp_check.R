suppressWarnings(pkgload::load_all("."), quiet = TRUE))
sim_btl_beta <- function(beta, n_judges = 12, reps_per_pair = 25, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  objs <- names(beta); K <- length(objs)
  jids <- sprintf("J%d", seq_len(n_judges))
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(object_a = rep(pr[, 1], each = reps_per_pair),
                   object_b = rep(pr[, 2], each = reps_per_pair),
                   stringsAsFactors = FALSE)
  d$judge <- sample(jids, nrow(d), TRUE)
  lp <- beta[d$object_a] - beta[d$object_b]
  resp <- as.integer(stats::runif(nrow(d)) < stats::plogis(lp))
  d$winner <- ifelse(resp == 1L, d$object_a, d$object_b)
  rownames(d) <- NULL
  d
}
K <- 8
beta0 <- setNames(as.numeric(scale(seq_len(K))) * 1, sprintf("O%d", seq_len(K)))
set.seed(52)
n <- 80
raw_flags <- list(); adj_flags <- list()
for (r in seq_len(n)) {
  dA <- sim_btl_beta(beta0, seed = 700000 + r); dB <- sim_btl_beta(beta0, seed = 800000 + r)
  f1 <- btl(dA, "object_a", "object_b", winner = "winner", judge = "judge")
  f2 <- btl(dB, "object_a", "object_b", winner = "winner", judge = "judge")
  eq <- btl_equate(f1, f2, independent = TRUE)
  raw_flags[[r]] <- eq$table$p < 0.05
  adj_flags[[r]] <- eq$table$p_adj < 0.05
}
rf <- unlist(raw_flags); rf <- rf[!is.na(rf)]
af <- unlist(adj_flags); af <- af[!is.na(af)]
cat(sprintf("raw p<0.05 rate: %.4f (n=%d)\n", mean(rf), length(rf)))
cat(sprintf("Holm-adjusted p_adj<0.05 rate: %.4f (n=%d)\n", mean(af), length(af)))
