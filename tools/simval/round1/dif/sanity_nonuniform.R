suppressWarnings(pkgload::load_all(".", quiet=TRUE))
target_item <- "I06"
for (shift in c(0.6, 1.0, 1.5, 2.0)) {
  flags <- c()
  for (k in 1:40) {
    d <- simulate_rasch(500, 20, dif = list(items = target_item, nonuniform = shift),
                         n_groups = 2, seed = 500000 + k)
    fit <- rasch(d, id="id", factors="group")
    da <- dif_anova(fit)
    s <- da$summary
    flags <- c(flags, isTRUE(s$nonuniform_DIF[s$item==target_item]))
  }
  cat("nonuniform shift=", shift, " detection rate=", mean(flags), "\n")
}
