suppressWarnings(pkgload::load_all(".", quiet=TRUE))
batch <- sim_replicate(simulate_btl, 5, n_objects = 10, n_judges = 14,
                        reps_per_pair = 20, seed = 100)
dd <- batch[[1]]
fit <- tryCatch(btl(dd, "object_a", "object_b", winner = "winner", judge = "judge"),
                 error=function(e) e)
print(fit)
