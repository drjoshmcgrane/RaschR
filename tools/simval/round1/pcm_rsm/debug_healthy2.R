suppressWarnings(pkgload::load_all("."), quiet=TRUE))
mins <- c(); weak <- c()
for (r in 1:40) {
  d <- simulate_rasch(n_persons = 1000, n_items = 5, model = "PCM",
                       n_categories = 5, difficulty = c(-1, 1),
                       threshold_spread = 1.0, seed = 81000 + r)
  f <- rasch(d, model="PCM")
  cc <- sapply(f$thresholds_diag, function(x) x$category_counts)
  mins <- c(mins, min(cc))
  weak <- c(weak, any(f$thresholds$weak))
}
cat("min category count across reps: min=", min(mins), " median=", median(mins), "\n")
cat("weak rate:", mean(weak), "\n")
