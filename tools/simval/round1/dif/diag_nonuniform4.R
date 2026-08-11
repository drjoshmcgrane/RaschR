suppressWarnings(pkgload::load_all(".", quiet=TRUE))
target_item <- "I10"
for (N in c(1500)) {
 for (shift in c(1.5, 2.5)) {
  flags <- c()
  for (k in 1:25) {
    d <- simulate_rasch(N, 20, dif = list(items = target_item, nonuniform = shift),
                         n_groups = 2, seed = 700000 + k)
    fit <- rasch(d, id="id", factors="group")
    da <- dif_anova(fit)
    s <- da$summary
    flags <- c(flags, isTRUE(s$nonuniform_DIF[s$item==target_item]))
  }
  cat("N=",N," shift=", shift, " detection rate=", mean(flags), "\n")
 }
}
