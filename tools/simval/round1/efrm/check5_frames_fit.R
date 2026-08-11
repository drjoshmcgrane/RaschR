suppressWarnings(pkgload::load_all(".", quiet=TRUE))
t0 <- Sys.time()

n_reps <- 20
n_per_group <- 300
items_per_set <- 8

infit_all <- c(); outfit_all <- c()
for (r in seq_len(n_reps)) {
  seed <- 200000L + r
  d <- simulate_efrm(n_per_group, items_per_set, n_sets = 3, n_groups = 3,
                      set_unit_ratio = 1.0, group_unit_ratio = 1.0, seed = seed)
  fit <- tryCatch(rasch_efrm(d, item_sets = attr(d, "truth")$item_sets, groups = "group",
                              se_method = "hybrid"), error = function(e) NULL)
  if (is.null(fit)) next
  infit_all <- c(infit_all, fit$frames$infit_ms)
  outfit_all <- c(outfit_all, fit$frames$outfit_ms)
}
cat(sprintf("n_reps_ok=%d  n_frames=%d\n", n_reps, length(infit_all)))
cat(sprintf("infit_ms:  mean=%.4f sd=%.4f range=[%.3f, %.3f]\n",
            mean(infit_all), sd(infit_all), min(infit_all), max(infit_all)))
cat(sprintf("outfit_ms: mean=%.4f sd=%.4f range=[%.3f, %.3f]\n",
            mean(outfit_all), sd(outfit_all), min(outfit_all), max(outfit_all)))
saveRDS(list(infit = infit_all, outfit = outfit_all),
        "tools/simval/round1/efrm/check5_results.rds")
cat("elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
