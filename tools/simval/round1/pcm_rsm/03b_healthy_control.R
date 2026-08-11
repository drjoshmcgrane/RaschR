suppressWarnings(pkgload::load_all(".", quiet=TRUE))
healthy_rate <- c()
for (r in seq_len(60)) {
  d <- simulate_rasch(n_persons = 1000, n_items = 5, model = "PCM",
                       n_categories = 5, difficulty = c(-1, 1),
                       threshold_spread = 1.0, seed = 82000 + r)
  f <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
  if (!inherits(f, "error")) healthy_rate <- c(healthy_rate, any(f$thresholds$weak))
}
cat(sprintf("healthy design (min category counts comfortably >3): weak-flag-fires rate over %d reps = %.4f\n",
            length(healthy_rate), mean(healthy_rate)))
