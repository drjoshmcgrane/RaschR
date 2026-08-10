suppressWarnings(pkgload::load_all("."), quiet=TRUE))
n_items_b <- 20; skip_idx <- 16:20; N_b <- 600
d <- simulate_rasch(N_b, n_items_b, n_groups = 2, seed = 1)
d[d$group == "g1", sprintf("I%02d", skip_idx)] <- NA
fit <- rasch(d, id = "id", factors = "group")
da <- dif_anova(fit)
print(table(da$summary$item))
print(da$notes)
