suppressWarnings(pkgload::load_all(".", quiet=TRUE))
target_item <- "I06"
d <- simulate_rasch(5000, 20, dif = list(items = target_item, nonuniform = 2.5),
                     n_groups = 2, seed = 1)
fit <- rasch(d, id="id", factors="group")
Z <- fit$residuals[, target_item]
theta <- fit$person$theta
grp <- fit$factors$group
ci <- cut(theta, quantile(theta, seq(0,1,0.2), na.rm=TRUE), include.lowest=TRUE)
tab <- aggregate(Z, list(ci=ci, grp=grp), mean, na.rm=TRUE)
print(tab[order(tab$ci, tab$grp),])
cat("\nraw p(x=1) by ci x group:\n")
X <- fit$X[,target_item]
tab2 <- aggregate(X, list(ci=ci, grp=grp), mean, na.rm=TRUE)
print(tab2[order(tab2$ci, tab2$grp),])
cat("\nitem location est (pooled):", fit$items$location[fit$items$item==target_item], "\n")
truth <- attr(d,"truth")
cat("planted: delta=", truth$difficulty[target_item], " disc base=1, nonuniform shift=2.5 -> disc g2=3.5\n")
