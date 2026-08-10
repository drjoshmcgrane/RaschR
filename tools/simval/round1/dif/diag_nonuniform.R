suppressWarnings(pkgload::load_all("."), quiet=TRUE))
target_item <- "I06"
d <- simulate_rasch(3000, 20, dif = list(items = target_item, nonuniform = 2.5),
                     n_groups = 2, seed = 1)
fit <- rasch(d, id="id", factors="group")
da <- dif_anova(fit)
print(da$summary[da$summary$item==target_item,])
cat("n_groups (class intervals) used:", da$n_groups, "\n")

# check raw residual pattern by class interval x group for I06
Z <- fit$residuals[, target_item]
grp <- fit$factors$group
ci <- .rasch_dummy <- NULL
