suppressWarnings(pkgload::load_all(".", quiet=TRUE))

## --- (9a) ~25% MCAR: recovery + SE calibration + pairwise chi-sq null ------
n_rep <- 150
K <- 10; J <- 14; reps_pp <- 25
cor_v <- rep(NA_real_, n_rep); rmse_v <- rep(NA_real_, n_rep)
rej <- rep(NA, n_rep)
LOC <- list(); SE <- list(); OBJ <- list()
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = reps_pp, seed = 50000 + i)
    set.seed(50000 + i + 777)
    d <- d[runif(nrow(d)) >= 0.25, ]   # ~25% MCAR
    f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
    rec <- sim_recovery(f, d)
    cor_v[i] <- rec$summary$correlation
    rmse_v[i] <- rec$summary$rmse
    rej[i] <- f$total_p < 0.05
    LOC[[i]] <- f$objects$location; SE[[i]] <- f$objects$se; OBJ[[i]] <- f$objects$object
    TRUE
  }, error = function(e) FALSE)
}
cat("== 25% MCAR: recovery + null pairwise chi-square (K=10,J=14) ==\n")
cat("mean correlation:", round(mean(cor_v, na.rm=TRUE), 4),
    " mean rmse:", round(mean(rmse_v, na.rm=TRUE), 4),
    " n=", sum(!is.na(cor_v)), "\n")
cat("pairwise chi-sq rejection rate at alpha=0.05:", round(mean(rej, na.rm=TRUE), 4),
    " MC err~", round(sqrt(0.05*0.95/sum(!is.na(rej))), 4), "\n")
okrows <- !vapply(LOC, is.null, TRUE)
objs <- OBJ[[which(okrows)[1]]]
LOCm <- t(sapply(which(okrows), function(i) LOC[[i]][match(objs, OBJ[[i]])]))
SEm  <- t(sapply(which(okrows), function(i) SE[[i]][match(objs, OBJ[[i]])]))
ratio <- apply(LOCm, 2, sd) / colMeans(SEm)
cat("SE calibration under MCAR -- mean(emp_sd/mean_se) across objects:",
    round(mean(ratio), 3), " range:", round(range(ratio), 3), "\n")

## --- (9b) sparse-but-connected structural design (a "chain" plus a few
## cross-links, far from complete round-robin) ------------------------------
set.seed(60000)
K2 <- 12
objs2 <- sprintf("O%d", seq_len(K2))
beta2 <- setNames(as.numeric(scale(seq_len(K2))), objs2)
# chain: adjacent objects only, plus wraparound + a couple of chords to keep
# it well connected without being a full round robin
chain_pairs <- rbind(
  cbind(objs2[-K2], objs2[-1]),
  cbind(objs2[c(1,4,7)], objs2[c(6,9,12)])  # a few chords for redundancy
)
reps_each <- 60
d3 <- data.frame(object_a = rep(chain_pairs[,1], each = reps_each),
                  object_b = rep(chain_pairs[,2], each = reps_each),
                  stringsAsFactors = FALSE)
d3$judge <- sample(sprintf("J%d", 1:14), nrow(d3), TRUE)
p3 <- plogis(beta2[d3$object_a] - beta2[d3$object_b])
d3$winner <- ifelse(runif(nrow(d3)) < p3, d3$object_a, d3$object_b)
f3 <- btl(d3, "object_a", "object_b", winner = "winner", judge = "judge")
cat("\n== Sparse-but-connected structural design (chain + chords, K=12,",
    nrow(chain_pairs), "of", choose(K2,2), "possible pairs used) ==\n")
cat("converged:", f3$converged, "\n")
loc_est <- setNames(f3$objects$location, f3$objects$object)
cat("recovery correlation:", round(cor(beta2[names(loc_est)], loc_est), 4),
    " rmse:", round(sqrt(mean((loc_est - (beta2[names(loc_est)] -
      mean(beta2[names(loc_est)])))^2)), 4), "\n")
cat("mean se (sparse design):", round(mean(f3$objects$se), 4),
    " (compare to a full round-robin at similar N, se in scripts 01/03 ~0.15-0.19)\n")

## --- (9c) separated design: Ford refusal is a PASS -------------------------
set.seed(70000)
K3 <- 8
objs3 <- sprintf("O%d", seq_len(K3))
beta3 <- setNames(as.numeric(scale(seq_len(K3))), objs3)
cluster1 <- objs3[1:4]; cluster2 <- objs3[5:8]
pr1 <- t(utils::combn(cluster1, 2)); pr2 <- t(utils::combn(cluster2, 2))
pr <- rbind(pr1, pr2)   # NOTE: no cross-cluster pairs at all -- two disjoint tournaments
d4 <- data.frame(object_a = rep(pr[,1], each = 30), object_b = rep(pr[,2], each = 30),
                  stringsAsFactors = FALSE)
p4 <- plogis(beta3[d4$object_a] - beta3[d4$object_b])
d4$winner <- ifelse(runif(nrow(d4)) < p4, d4$object_a, d4$object_b)
res <- tryCatch(btl(d4, "object_a", "object_b", winner = "winner"),
                 error = function(e) e)
cat("\n== Separated design (two disjoint clusters, no cross-cluster comparisons) ==\n")
if (inherits(res, "error")) {
  cat("btl() REFUSED as expected. Error message:\n  ", conditionMessage(res), "\n")
} else {
  cat("btl() did NOT refuse -- returned a fit (this would be a FAIL for the Ford guard)\n")
}
