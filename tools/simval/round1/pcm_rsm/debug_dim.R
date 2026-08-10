suppressWarnings(pkgload::load_all("."), quiet=TRUE))
make_booklet <- function(d) {
  N <- nrow(d)
  grp <- sample(c("A","B"), N, replace = TRUE)
  d2 <- d
  d2[grp == "A", c("I07","I08","I09","I10")] <- NA
  d2[grp == "B", c("I01","I02","I03","I04")] <- NA
  d2
}
for (r in 1:60) {
  set.seed(3000+r)
  d <- simulate_rasch(n_persons = 600, n_items = 10, model = "PCM",
                       n_categories = 4, difficulty = c(-2, 2),
                       threshold_spread = 1.3, seed = 3000 + r)
  d <- make_booklet(d)
  fit <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
  if (inherits(fit, "error")) { cat("rep",r,"ERROR:", conditionMessage(fit), "\n"); next }
  if (length(fit$thresholds$tau) != 30) cat("rep", r, "tau length", length(fit$thresholds$tau), "notes:", paste(fit$notes,collapse="; "), "\n")
}
