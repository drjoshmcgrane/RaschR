suppressWarnings(pkgload::load_all("."), quiet=TRUE))
batch <- sim_replicate(simulate_btl, 5, n_objects = 10, n_judges = 14,
                        reps_per_pair = 20, seed = 100)
fits <- sim_apply(batch, function(dd) {
  fit <- btl(dd, "object_a", "object_b", winner = "winner", judge = "judge")
  list(loc = fit$objects$location, se = fit$objects$se,
       obj = fit$objects$object, converged = fit$converged)
})
str(fits[[1]])
cat("n_failed:", attr(fits, "n_failed"), "\n")
