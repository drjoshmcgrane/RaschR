suppressWarnings(pkgload::load_all(".", quiet=TRUE))
d <- simulate_rasch(n_persons = 600, n_items = 5, model = "PCM",
                     n_categories = 5, difficulty = c(-2, 2),
                     threshold_spread = 1.4, seed = 80201)
f <- rasch(d, model="PCM")
print(sapply(f$thresholds_diag, function(x) x$category_counts))
cat("weak:", any(f$thresholds$weak), "\n")
print(f$notes)
