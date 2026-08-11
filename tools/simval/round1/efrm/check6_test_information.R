suppressWarnings(pkgload::load_all(".", quiet=TRUE))
t0 <- Sys.time()

# ---- build a linking-design EFRM dataset: most persons take one set only,
# a linking subsample per group also takes the other set. simulate_efrm()
# itself administers every item to every person (fully crossed), so the
# booklet/linking structure is created here by hand, post-hoc, via NA.
set.seed(42)
n_per_group <- 300
d <- simulate_efrm(n_per_group, 8, n_sets = 2, n_groups = 2,
                    set_unit_ratio = 1.3, group_unit_ratio = 1.2, seed = 42)
tr <- attr(d, "truth")
set1_items <- tr$item_sets$set1; set2_items <- tr$item_sets$set2

g1_idx <- which(d$group == "g1"); g2_idx <- which(d$group == "g2")
link_frac <- 0.15
n_link1 <- round(length(g1_idx) * link_frac)
n_link2 <- round(length(g2_idx) * link_frac)
set.seed(43)
link1 <- sample(g1_idx, n_link1)   # g1 persons who ALSO take set2
link2 <- sample(g2_idx, n_link2)   # g2 persons who ALSO take set1

# default: g1 takes only set1, g2 takes only set2
d[setdiff(g1_idx, link1), set2_items] <- NA
d[setdiff(g2_idx, link2), set1_items] <- NA
# link1/link2 already have both sets from simulate_efrm -- leave as is

fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", se_method = "hybrid")

cat("design blocks (administration patterns) found by the package:\n")
blocks <- rasch:::.design_blocks(fit)
print(names(blocks))

ti <- test_information(fit)
cat("\ntest_information() design levels:\n")
print(unique(ti$design))

# hand-compute information for each block at a few theta grid points and
# compare against test_information()'s reported curve for that design
grid_check <- c(-2, -1, 0, 1, 2)
ok_all <- TRUE
for (nm in names(blocks)) {
  ii <- blocks[[nm]]
  hand <- vapply(grid_check, function(th)
    sum(vapply(ii, function(i)
      fit$disc[i]^2 * item_moments(th, fit$tau_list[[i]], disc = fit$disc[i])$V, 0)), 0)
  pkg <- ti$info[ti$design == nm][match(grid_check, ti$theta[ti$design == nm])]
  diff <- max(abs(hand - pkg))
  cat(sprintf("  [%s] n_items=%d  max|hand-pkg| = %.10f\n", nm, length(ii), diff))
  if (!isTRUE(diff < 1e-8)) ok_all <- FALSE
}
cat("ALL BLOCKS MATCH HAND COMPUTATION:", ok_all, "\n")

saveRDS(list(blocks = names(blocks), ok_all = ok_all,
             design_levels = unique(ti$design)),
        "tools/simval/round1/efrm/check6_results.rds")

# ---- unlinked design: no person takes more than one set -> must be refused
cat("\n== unlinked design refusal check ==\n")
set.seed(44)
d2 <- simulate_efrm(n_per_group, 8, n_sets = 2, n_groups = 2,
                     set_unit_ratio = 1.3, group_unit_ratio = 1.0, seed = 44)
tr2 <- attr(d2, "truth")
g1_idx2 <- which(d2$group == "g1"); g2_idx2 <- which(d2$group == "g2")
d2[g1_idx2, tr2$item_sets$set2] <- NA   # g1: only set1
d2[g2_idx2, tr2$item_sets$set1] <- NA   # g2: only set2
refused <- tryCatch({
  rasch_efrm(d2, item_sets = tr2$item_sets, groups = "group", se_method = "hybrid")
  FALSE
}, error = function(e) { cat("refused with message:\n  ", conditionMessage(e), "\n"); TRUE })
cat("REFUSED (expected TRUE):", refused, "\n")

saveRDS(list(refused = refused), "tools/simval/round1/efrm/check6b_unlinked.rds")
cat("total elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
