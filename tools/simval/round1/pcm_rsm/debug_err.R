suppressWarnings(pkgload::load_all("."), quiet=TRUE))
make_booklet <- function(d) {
  N <- nrow(d)
  grp <- sample(c("A","B"), N, replace = TRUE)
  d2 <- d
  d2[grp == "A", c("I07","I08","I09","I10")] <- NA
  d2[grp == "B", c("I01","I02","I03","I04")] <- NA
  d2
}
for (r in 1:5) {
  d <- simulate_rasch(n_persons = 500, n_items = 8, model = "PCM",
                       n_categories = 4, difficulty = c(-2, 2),
                       threshold_spread = 1.3, seed = 30000 + r)
  d <- make_booklet(d)
  fit <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
  if (inherits(fit, "error")) cat("rep", r, "ERROR:", conditionMessage(fit), "\n")
  else cat("rep", r, "ok\n")
}
