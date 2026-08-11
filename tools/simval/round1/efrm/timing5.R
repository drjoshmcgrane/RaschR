suppressWarnings(pkgload::load_all(".", quiet=TRUE))
source("tools/simval/round1/efrm/sim_crossed.R")
d <- sim_crossed_phi(150, 10, region_effect = log(1.3), cohort_effect = 0, seed = 1)
t1 <- system.time(fit <- rasch_efrm(d, item_sets = list(only = attr(d,"item_names")),
                                     groups = c("region","cohort"),
                                     items = attr(d,"item_names"), se_method = "hybrid"))
print(t1)
print(fit$phi_table)
print(fit$phi_factorial)
print(fit$phi_factorial_tests)
