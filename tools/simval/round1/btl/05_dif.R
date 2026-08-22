suppressWarnings(pkgload::load_all(".", quiet=TRUE))

n_rep <- 200
K <- 6; J <- 12  # >= 10 judges, clears btl_dif's cluster-inference guard
grp <- setNames(rep(c("g1", "g2"), each = J/2), sprintf("J%d", 1:J))

## --- null: common object locations, group arbitrary/unrelated -------------
## simulate_btl() does not add judge-specific object deviations in this call.
## The separate source regression in test-statistical-validity.R constructs
## those deviations explicitly; this historical battery cell checks the
## equal-location null with random allocation of comparisons to judges.
any_uniform_flag <- logical(0)
per_object_flag <- logical(0)
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = 40, seed = 9000 + i)
    f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
    dif <- btl_dif(f, grp)
    any_uniform_flag <- c(any_uniform_flag, any(dif$summary$uniform_DIF))
    per_object_flag <- c(per_object_flag, dif$summary$uniform_DIF)
    TRUE
  }, error = function(e) FALSE)
}
cat("== btl_dif null calibration (common locations, arbitrary group, K=6,J=12) ==\n")
cat("replicates completed:", sum(!is.na(any_uniform_flag)), "/", n_rep, "\n")
cat("per-object false uniform-DIF flag rate (BH-adjusted, alpha=0.05):",
    round(mean(per_object_flag, na.rm=TRUE), 4),
    " n_object_tests=", length(per_object_flag),
    " MC err~", round(sqrt(0.05*0.95/length(per_object_flag)), 4), "\n")
cat("replicate-level 'any object flagged' false rate:",
    round(mean(any_uniform_flag, na.rm=TRUE), 4), "\n")

## --- power: planted judge-group location shift on one object -------------
pow_flag <- logical(0)
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    beta <- setNames(as.numeric(scale(seq_len(K))), sprintf("O%d", seq_len(K)))
    pr <- t(utils::combn(names(beta), 2))
    d <- data.frame(object_a = rep(pr[,1], each = 40), object_b = rep(pr[,2], each = 40),
                     stringsAsFactors = FALSE)
    set.seed(9500 + i)
    d$judge <- sample(sprintf("J%d", 1:J), nrow(d), TRUE)
    shift <- ifelse(d$judge %in% sprintf("J%d", 1:(J/2)) & d$object_a == "O3", 1.0,
              ifelse(d$judge %in% sprintf("J%d", 1:(J/2)) & d$object_b == "O3", -1.0, 0))
    p <- plogis(beta[d$object_a] - beta[d$object_b] + shift)
    d$winner <- ifelse(runif(nrow(d)) < p, d$object_a, d$object_b)
    f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
    dif <- btl_dif(f, grp, objects = "O3")
    pow_flag <- c(pow_flag, isTRUE(dif$summary$uniform_DIF[dif$summary$object == "O3"]))
    TRUE
  }, error = function(e) FALSE)
}
cat("\n== btl_dif power (planted 1.0-logit judge-group shift on O3, K=6,J=12) ==\n")
cat("replicates completed:", sum(!is.na(pow_flag)), "/", n_rep, "\n")
cat("uniform-DIF detection rate on O3:", round(mean(pow_flag, na.rm=TRUE), 4), "\n")
