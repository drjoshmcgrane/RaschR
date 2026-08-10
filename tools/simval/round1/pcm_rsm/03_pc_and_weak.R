suppressWarnings(pkgload::load_all("."), quiet=TRUE))
options(width=140)

## ---------------------------------------------------------------
## (4) pc_components on sparse categories: pc route stays ordered;
## on rich, well-separated data pc_components (saturated, n_components
## >= m) matches free estimation closely.
## ---------------------------------------------------------------
cat("=== pc_components: sparse categories ===\n")
sparse_ordered <- c(); sparse_free_disordered <- c(); sparse_free_weak <- c()
nrep <- 60
for (r in seq_len(nrep)) {
  d <- simulate_rasch(n_persons = 150, n_items = 6, model = "PCM",
                       n_categories = 5, difficulty = c(-1.5, 1.5),
                       threshold_spread = 0.55, seed = 60000 + r)
  fit_pc <- tryCatch(rasch(d, model = "PCM", pc_components = 2),
                      error = function(e) e)
  fit_free <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
  if (!inherits(fit_pc, "error")) {
    ord <- vapply(fit_pc$thresholds_diag, function(x) x$ordered, TRUE)
    sparse_ordered <- c(sparse_ordered, mean(ord))
  }
  if (!inherits(fit_free, "error")) {
    ordf <- vapply(fit_free$thresholds_diag, function(x) x$ordered, TRUE)
    sparse_free_disordered <- c(sparse_free_disordered, mean(!ordf))
    sparse_free_weak <- c(sparse_free_weak, any(fit_free$thresholds$weak, na.rm=TRUE))
  }
}
cat(sprintf("pc_components=2 route: item-level ordered rate = %.4f (n reps=%d)\n",
            mean(sparse_ordered), length(sparse_ordered)))
cat(sprintf("free-threshold route on same sparse design: item-level disordered rate = %.4f; any-weak-flag rate = %.4f\n\n",
            mean(sparse_free_disordered), mean(sparse_free_weak)))

cat("=== pc_components vs free on RICH, well-separated data ===\n")
cor_tau <- c(); rmse_tau <- c()
nrep2 <- 40
for (r in seq_len(nrep2)) {
  d <- simulate_rasch(n_persons = 2000, n_items = 6, model = "PCM",
                       n_categories = 4, difficulty = c(-2, 2),
                       threshold_spread = 1.5, seed = 70000 + r)
  fit_pc <- tryCatch(rasch(d, model = "PCM", pc_components = 4), error = function(e) e)
  fit_free <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
  if (!inherits(fit_pc, "error") && !inherits(fit_free, "error")) {
    cor_tau <- c(cor_tau, cor(fit_pc$thresholds$tau, fit_free$thresholds$tau))
    rmse_tau <- c(rmse_tau, sqrt(mean((fit_pc$thresholds$tau - fit_free$thresholds$tau)^2)))
  }
}
cat(sprintf("cor(pc4 tau, free tau) = %.5f (mean over %d reps, MCse=%.5f); rmse = %.5f\n\n",
            mean(cor_tau), length(cor_tau), sd(cor_tau)/sqrt(length(cor_tau)), mean(rmse_tau)))

## ---------------------------------------------------------------
## (5) Weak-category guard: fires (se NA + note) when a category has
## <3 responses; does NOT fire on healthy data.
## ---------------------------------------------------------------
cat("=== weak-category guard ===\n")
# trigger design: small N, many categories, one item pushed to have a
# near-empty middle category via extreme per-item threshold spacing
set.seed(80001)
d_bad <- simulate_rasch(n_persons = 60, n_items = 5, model = "PCM",
                         n_categories = 5, difficulty = c(-1, 1),
                         threshold_spread = 3.5, seed = 80001)
fit_bad <- rasch(d_bad, model = "PCM")
cat("category counts per item (trigger design):\n")
print(sapply(fit_bad$thresholds_diag, function(x) x$category_counts))
cat("any weak flagged:", any(fit_bad$thresholds$weak), "\n")
cat("any se NA among weak rows all NA:",
    all(is.na(fit_bad$thresholds$se[fit_bad$thresholds$weak])), "\n")
cat("weak-category note present:",
    any(grepl("weakly determined", fit_bad$notes)), "\n\n")

# repeat over reps to get a rate for the trigger design
trig_rate <- c()
for (r in seq_len(40)) {
  d <- simulate_rasch(n_persons = 60, n_items = 5, model = "PCM",
                       n_categories = 5, difficulty = c(-1, 1),
                       threshold_spread = 3.5, seed = 80100 + r)
  f <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
  if (!inherits(f, "error")) trig_rate <- c(trig_rate, any(f$thresholds$weak))
}
cat(sprintf("trigger design: weak-flag-fires rate over %d reps = %.3f\n\n", length(trig_rate), mean(trig_rate)))

# healthy control: ample N, well-separated categories -> guard should not fire
healthy_rate <- c()
for (r in seq_len(40)) {
  d <- simulate_rasch(n_persons = 600, n_items = 5, model = "PCM",
                       n_categories = 5, difficulty = c(-2, 2),
                       threshold_spread = 1.4, seed = 80200 + r)
  f <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
  if (!inherits(f, "error")) healthy_rate <- c(healthy_rate, any(f$thresholds$weak))
}
cat(sprintf("healthy design: weak-flag-fires rate over %d reps = %.3f\n\n", length(healthy_rate), mean(healthy_rate)))
