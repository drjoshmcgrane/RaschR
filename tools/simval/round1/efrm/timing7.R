suppressWarnings(pkgload::load_all(".", quiet=TRUE))
set.seed(7)
n_per_group <- 250
d <- simulate_efrm(n_per_group, 8, n_sets=2, n_groups=2, set_unit_ratio=1.3, group_unit_ratio=1, seed=7)
tr <- attr(d,"truth")
N <- nrow(d)
cohort <- factor(sample(c("C1","C2"), N, replace=TRUE))
d$cohort <- cohort

# inject uniform DIF on item S1I01 for cohort==C2: regenerate that column
shift <- 0.9
alpha_true <- setNames(tr$alpha, names(tr$item_sets)); phi_true <- tr$phi
set_of_item <- setNames(rep(names(tr$item_sets), each=8), unlist(tr$item_sets))
theta <- tr$theta
grp_idx <- as.integer(d$group)
rho <- unname(alpha_true[set_of_item["S1I01"]]) * phi_true[grp_idx]
delta_i <- tr$difficulty["S1I01"]
shift_vec <- ifelse(cohort=="C2", shift, 0)
set.seed(701)
d$S1I01 <- rbinom(N, 1, plogis(rho*(theta - delta_i - shift_vec)))

fit <- rasch_efrm(d, item_sets=tr$item_sets, groups="group", factors="cohort", se_method="hybrid")
cat("frame_group:", fit$frame_group, "\n")
da <- dif_anova(fit, factors=c("group","cohort"))
cat(da$notes, sep="\n")
print(names(da))
print(head(da$summary))
