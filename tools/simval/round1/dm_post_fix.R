suppressWarnings(pkgload::load_all(".", quiet=TRUE))
run <- function(nreps, gen, dep, ind, tag) {
  p <- numeric(nreps); d <- numeric(nreps); se <- numeric(nreps)
  for (r in seq_len(nreps)) {
    dm <- dependence_magnitude(gen(r), dependent = dep, independent = ind)
    p[r] <- dm$p; d[r] <- dm$d; se[r] <- dm$se
  }
  cat(sprintf("%s: reject@.05 = %.4f (n=%d, MC %.4f); SD(d)=%.4f meanSE=%.4f ratio=%.3f\n",
      tag, mean(p < .05), nreps, sqrt(.05*.95/nreps), sd(d), mean(se), sd(d)/mean(se)))
}
run(400, function(r) rasch(simulate_rasch(800, 10, seed = 910000 + r)),
    "I06", "I05", "dichotomous null (post-fix)")
run(200, function(r) rasch(simulate_rasch(600, 8, model = "PCM",
    n_categories = 3, seed = 920000 + r), model = "PCM"),
    "I05", "I04", "PCM null (post-fix)")
