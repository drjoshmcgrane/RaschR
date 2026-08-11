suppressWarnings(pkgload::load_all(".", quiet=TRUE))
d <- simulate_btl(n_objects=8, n_judges=12, reps_per_pair=25, model="graded", n_categories=4, seed=500)
ff <- btl(d, "object_a","object_b", response="response", judge="judge", thresholds="free")
print(ff$components)
v <- ff$components$estimate[ff$components$component=="spread"]
print(v)
