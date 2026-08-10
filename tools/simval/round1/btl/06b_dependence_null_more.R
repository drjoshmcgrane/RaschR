suppressWarnings(pkgload::load_all("."), quiet=TRUE))
n_rep <- 400
K <- 8; J <- 14; reps_pp <- 25
alpha <- 0.05
p_exp <- rep(NA_real_, n_rep); p_cry <- rep(NA_real_, n_rep)
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = reps_pp,
                       dependence = list(exposure = 0, carry_over = 0),
                       seed = 21000 + i)
    f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge", order = "order")
    p_exp[i] <- f$dependence$p[f$dependence$effect == "exposure"]
    p_cry[i] <- f$dependence$p[f$dependence$effect == "carry_over"]
    TRUE
  }, error = function(e) FALSE)
}
cat("== Dependence null calibration, n_rep=", n_rep, "==\n")
cat("exposure rejection rate:", round(mean(p_exp<alpha,na.rm=TRUE),4),
    " n=", sum(!is.na(p_exp)), " MC err~", round(sqrt(0.05*0.95/sum(!is.na(p_exp))),4), "\n")
cat("carry_over rejection rate:", round(mean(p_cry<alpha,na.rm=TRUE),4),
    " n=", sum(!is.na(p_cry)), "\n")
