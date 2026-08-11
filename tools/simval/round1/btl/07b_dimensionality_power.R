suppressWarnings(pkgload::load_all(".", quiet=TRUE))
n_rep <- 60
flag_camp <- rep(NA, n_rep)
for (i in seq_len(n_rep)) {
  ok <- tryCatch({
    d <- simulate_btl(n_objects = 15, n_judges = 20, reps_per_pair = 60,
                       object_sd = 1.5, second_attribute = list(rho = 0.0),
                       seed = 33000 + i)
    f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
    bd <- btl_dimensionality(f, reps = 100)
    flag_camp[i] <- isTRUE(bd$leading_structured)
    TRUE
  }, error = function(e) FALSE)
}
cat("== btl_dimensionality power, stronger planted judge-camp second attribute ==\n")
cat("(K=15,J=20,reps_pp=60,object_sd=1.5,rho=0)\n")
cat("flag rate:", round(mean(flag_camp, na.rm = TRUE), 4), " n=", sum(!is.na(flag_camp)), "\n")
