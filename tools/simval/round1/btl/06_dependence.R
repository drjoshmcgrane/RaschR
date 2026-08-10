suppressWarnings(pkgload::load_all("."), quiet=TRUE))

n_rep <- 150
K <- 8; J <- 14; reps_pp <- 25

## --- planted exposure + carry-over: recovery -------------------------------
exq_true <- 0.6; cry_true <- 0.5
est_exp <- rep(NA_real_, n_rep); se_exp <- rep(NA_real_, n_rep)
est_cry <- rep(NA_real_, n_rep); se_cry <- rep(NA_real_, n_rep)
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = reps_pp,
                       dependence = list(exposure = exq_true, carry_over = cry_true),
                       seed = 10000 + i)
    f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge", order = "order")
    est_exp[i] <- f$dependence$estimate[f$dependence$effect == "exposure"]
    se_exp[i]  <- f$dependence$se[f$dependence$effect == "exposure"]
    est_cry[i] <- f$dependence$estimate[f$dependence$effect == "carry_over"]
    se_cry[i]  <- f$dependence$se[f$dependence$effect == "carry_over"]
    TRUE
  }, error = function(e) FALSE)
}
cat("== Dependence (order effect) recovery, planted exposure=", exq_true,
    "carry_over=", cry_true, "==\n")
cat("exposure:   mean est", round(mean(est_exp,na.rm=TRUE),4),
    " bias", round(mean(est_exp,na.rm=TRUE)-exq_true,4),
    " rmse", round(sqrt(mean((est_exp-exq_true)^2,na.rm=TRUE)),4),
    " emp_sd/mean_se", round(sd(est_exp,na.rm=TRUE)/mean(se_exp,na.rm=TRUE),3), "\n")
cat("carry_over: mean est", round(mean(est_cry,na.rm=TRUE),4),
    " bias", round(mean(est_cry,na.rm=TRUE)-cry_true,4),
    " rmse", round(sqrt(mean((est_cry-cry_true)^2,na.rm=TRUE)),4),
    " emp_sd/mean_se", round(sd(est_cry,na.rm=TRUE)/mean(se_cry,na.rm=TRUE),3), "\n")

## --- null: order present but no true dependence, rejection rate ----------
alpha <- 0.05
p_exp <- rep(NA_real_, n_rep); p_cry <- rep(NA_real_, n_rep)
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = reps_pp,
                       dependence = list(exposure = 0, carry_over = 0),
                       seed = 11000 + i)
    f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge", order = "order")
    p_exp[i] <- f$dependence$p[f$dependence$effect == "exposure"]
    p_cry[i] <- f$dependence$p[f$dependence$effect == "carry_over"]
    TRUE
  }, error = function(e) FALSE)
}
cat("\n== Dependence null calibration (order present, no true effect) ==\n")
cat("exposure rejection rate at alpha=0.05:", round(mean(p_exp<alpha,na.rm=TRUE),4),
    " n=", sum(!is.na(p_exp)),
    " MC err~", round(sqrt(0.05*0.95/sum(!is.na(p_exp))),4), "\n")
cat("carry_over rejection rate at alpha=0.05:", round(mean(p_cry<alpha,na.rm=TRUE),4),
    " n=", sum(!is.na(p_cry)), "\n")
