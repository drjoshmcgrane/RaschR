suppressWarnings(pkgload::load_all(".", quiet=TRUE))

flag_z <- 2.5   # same threshold the package's own print.rasch_btl uses

## --- (3a) judge fit: null rate vs erratic-judge power ----------------------
n_rep <- 200
J <- 12
null_judge_flags <- logical(0)
pow_judge_flags <- logical(0)
pow_reg_flags <- logical(0)
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    dn <- simulate_btl(n_objects = 8, n_judges = J, reps_per_pair = 25, seed = 2000 + i)
    fn <- btl(dn, "object_a", "object_b", winner = "winner", judge = "judge")
    null_judge_flags <- c(null_judge_flags, abs(fn$judges$fit_resid) > flag_z)

    de <- simulate_btl(n_objects = 8, n_judges = J, reps_per_pair = 25,
                        erratic_judges = 1/J, seed = 3000 + i)  # exactly 1 erratic judge
    fe <- btl(de, "object_a", "object_b", winner = "winner", judge = "judge")
    erratic_id <- attr(de, "truth")$erratic
    is_err <- fe$judges$judge %in% erratic_id
    pow_judge_flags <- c(pow_judge_flags, abs(fe$judges$fit_resid[is_err]) > flag_z)
    pow_reg_flags   <- c(pow_reg_flags,   abs(fe$judges$fit_resid[!is_err]) > flag_z)
    TRUE
  }, error = function(e) FALSE)
}
cat("== Judge fit: null vs power (flag |fit_resid| >", flag_z, ") ==\n")
cat("null false-flag rate (regular judges, no erratic present):",
    round(mean(null_judge_flags, na.rm=TRUE), 4),
    " n=", sum(!is.na(null_judge_flags)),
    " MC err~", round(sqrt(0.05*0.95/length(null_judge_flags)),4), "\n")
cat("power: erratic-judge flag rate:", round(mean(pow_judge_flags, na.rm=TRUE), 4),
    " n=", sum(!is.na(pow_judge_flags)), "\n")
cat("power run, regular-judge false-flag rate (same replicates):",
    round(mean(pow_reg_flags, na.rm=TRUE), 4),
    " n=", sum(!is.na(pow_reg_flags)), "\n")

## --- (3b) object fit: null rate vs an inconsistent object ------------------
# simulate_btl has no "erratic object" knob; build one by corrupting all rows
# involving a target object to a coin-flip outcome, independent of location,
# on top of an otherwise-model-true draw (a custom generator layered on the
# package simulator, per the brief's "build custom generators only where the
# package's own simulators do not support the design").
make_inconsistent_object <- function(seed, K = 8, J = 12, reps = 25, target = "O1") {
  d <- simulate_btl(n_objects = K, n_judges = J, reps_per_pair = reps, seed = seed)
  hit <- d$object_a == target | d$object_b == target
  set.seed(seed + 999983)
  coin <- sample(c(TRUE, FALSE), sum(hit), TRUE)
  d$winner[hit] <- ifelse(coin, d$object_a[hit], d$object_b[hit])
  d
}
null_obj_flags <- logical(0)
pow_obj_flags <- logical(0)
pow_obj_reg_flags <- logical(0)
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    dn <- simulate_btl(n_objects = 8, n_judges = 12, reps_per_pair = 25, seed = 4000 + i)
    fn <- btl(dn, "object_a", "object_b", winner = "winner", judge = "judge")
    null_obj_flags <- c(null_obj_flags, abs(fn$objects$fit_resid) > flag_z)

    di <- make_inconsistent_object(5000 + i)
    fi <- btl(di, "object_a", "object_b", winner = "winner", judge = "judge")
    is_bad <- fi$objects$object == "O1"
    pow_obj_flags <- c(pow_obj_flags, abs(fi$objects$fit_resid[is_bad]) > flag_z)
    pow_obj_reg_flags <- c(pow_obj_reg_flags, abs(fi$objects$fit_resid[!is_bad]) > flag_z)
    TRUE
  }, error = function(e) FALSE)
}
cat("\n== Object fit: null vs power (flag |fit_resid| >", flag_z, ") ==\n")
cat("null false-flag rate (regular objects):", round(mean(null_obj_flags, na.rm=TRUE), 4),
    " n=", sum(!is.na(null_obj_flags)), "\n")
cat("power: inconsistent-object flag rate:", round(mean(pow_obj_flags, na.rm=TRUE), 4),
    " n=", sum(!is.na(pow_obj_flags)), "\n")
cat("power run, regular-object false-flag rate (same replicates):",
    round(mean(pow_obj_reg_flags, na.rm=TRUE), 4),
    " n=", sum(!is.na(pow_obj_reg_flags)), "\n")

## --- (3c) pairwise chi-square: null rejection rate + df check -------------
alpha <- 0.05
rej <- logical(0); df_ok <- logical(0)
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    dn <- simulate_btl(n_objects = 8, n_judges = 12, reps_per_pair = 25, seed = 6000 + i)
    fn <- btl(dn, "object_a", "object_b", winner = "winner", judge = "judge")
    rej <- c(rej, fn$total_p < alpha)
    expected_df <- sum(fn$pairs$n >= 2) - nrow(fn$objects)  # np for dichotomous, no judge/thr params beyond locations
    df_ok <- c(df_ok, fn$total_df == expected_df)
    TRUE
  }, error = function(e) FALSE)
}
cat("\n== Pairwise chi-square (dichotomous, full round-robin, K=8,J=12) ==\n")
cat("rejection rate at alpha=0.05:", round(mean(rej, na.rm=TRUE), 4),
    " n=", sum(!is.na(rej)),
    " MC err~", round(sqrt(0.05*0.95/length(rej)), 4), "\n")
cat("df formula (pairs_used - n_objects) matches fit$total_df in all reps:",
    all(df_ok, na.rm = TRUE), " (", sum(!df_ok, na.rm=TRUE), "mismatches )\n")
