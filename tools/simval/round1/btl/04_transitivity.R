suppressWarnings(pkgload::load_all("."), quiet=TRUE))

n_rep <- 150
K <- 7; J <- 10; reps_pp <- 20

## --- model-true (consistent) full round-robin: circular rate reference ----
rate_null <- rep(NA_real_, n_rep)
zeta_null <- rep(NA_real_, n_rep)
for (i in seq_len(n_rep)) {
  d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = reps_pp, seed = 7000 + i)
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  tr <- btl_transitivity(f)
  rate_null[i] <- tr$summary$circular_rate
  zeta_null[i] <- tr$summary$zeta
}

## --- planted intransitivity: force a rock-paper-scissors triple among 3 --
## objects (O1 beats O2, O2 beats O3, O3 beats O1 at high, model-defying
## probability) layered on an otherwise model-true round robin among the
## remaining objects (custom generator: simulate_btl has no cyclic-triple
## option).
make_rps <- function(seed, K = 7, J = 10, reps = 20, p_cyc = 0.85) {
  d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = reps, seed = seed)
  set.seed(seed + 555001)
  is12 <- (d$object_a == "O1" & d$object_b == "O2") | (d$object_a == "O2" & d$object_b == "O1")
  is23 <- (d$object_a == "O2" & d$object_b == "O3") | (d$object_a == "O3" & d$object_b == "O2")
  is31 <- (d$object_a == "O3" & d$object_b == "O1") | (d$object_a == "O1" & d$object_b == "O3")
  win_cyc <- function(rows, winner_name, loser_name) {
    if (!any(rows)) return(NULL)
    hit_a <- d$object_a[rows] == winner_name
    ifelse(runif(sum(rows)) < p_cyc,
           ifelse(hit_a, d$object_a[rows], d$object_b[rows]),
           ifelse(hit_a, d$object_b[rows], d$object_a[rows]))
  }
  d$winner[is12] <- win_cyc(is12, "O1", "O2")
  d$winner[is23] <- win_cyc(is23, "O2", "O3")
  d$winner[is31] <- win_cyc(is31, "O3", "O1")
  d
}
rate_plant <- rep(NA_real_, n_rep)
zeta_plant <- rep(NA_real_, n_rep)
for (i in seq_len(n_rep)) {
  d <- make_rps(8000 + i)
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  tr <- btl_transitivity(f)
  rate_plant[i] <- tr$summary$circular_rate
  zeta_plant[i] <- tr$summary$zeta
}

cat("== Transitivity: model-true vs planted rock-paper-scissors (K=7 full round-robin) ==\n")
cat("model-true: mean circular rate", round(mean(rate_null),4),
    " sd", round(sd(rate_null),4), " MC err(mean)~", round(sd(rate_null)/sqrt(n_rep),4), "\n")
cat("model-true: mean zeta (consistency)", round(mean(zeta_null),4), "\n")
cat("planted RPS: mean circular rate", round(mean(rate_plant),4),
    " sd", round(sd(rate_plant),4), " MC err(mean)~", round(sd(rate_plant)/sqrt(n_rep),4), "\n")
cat("planted RPS: mean zeta (consistency)", round(mean(zeta_plant),4), "\n")
cat("t-test model-true vs planted circular rate: t=",
    round(t.test(rate_plant, rate_null)$statistic,2),
    " p=", format.pval(t.test(rate_plant, rate_null)$p.value, digits=3), "\n")
