suppressWarnings(pkgload::load_all("."), quiet=TRUE))
t_start <- Sys.time()

# The first nonuniform-DIF power run (item I06, delta ~ -1.2, shift 0.6,
# N=500) turned out badly underpowered: diagnostics showed a discrimination
# shift on an item whose difficulty sits near the edge of the observed
# trait range produces almost no crossover within the class intervals used
# (it reads mostly as uniform DIF, not an interaction) -- confirmed with a
# targeted check (N=3000, shift=2.5 on I06: eta2_nonuniform stayed ~0.009).
# Re-run with a item near the CENTRE of the trait range (delta ~ 0, so the
# crossover point falls inside the observed class intervals) and a shift
# demonstrated (in a smaller diagnostic run) to give strong power at
# N=1500 (0.88-1.0 detection at shift 1.5-2.5).
N <- 1500; n_items <- 20; reps <- 150; target_item <- "I10"; shift <- 1.5

detect_right <- c(); detect_any <- c(); other_false <- c(); fails <- 0
for (k in seq_len(reps)) {
  seed <- 720000 + k
  d <- simulate_rasch(N, n_items, dif = list(items = target_item, nonuniform = shift),
                       n_groups = 2, seed = seed)
  fit <- tryCatch(rasch(d, id = "id", factors = "group"), error = function(e) NULL)
  if (is.null(fit)) { fails <- fails + 1; next }
  da <- tryCatch(dif_anova(fit), error = function(e) NULL)
  if (is.null(da)) { fails <- fails + 1; next }
  s <- da$summary
  right <- s$nonuniform_DIF[s$item == target_item]
  detect_right <- c(detect_right, isTRUE(right))
  detect_any <- c(detect_any, any(s$nonuniform_DIF))
  other <- s$nonuniform_DIF[s$item != target_item]
  other_false <- c(other_false, other)
}
cat("=== power (corrected): non-uniform DIF (discrimination +", shift, ") on", target_item,
    "(centred item, N=", N, ") ===\n")
cat("detection rate (right item, nonuniform term):", mean(detect_right), " n=", length(detect_right), "\n")
cat("any-nonuniform-flag rate:", mean(detect_any), "\n")
cat("false nonuniform flag rate on OTHER items:", mean(other_false), " n=", length(other_false), "\n")
cat("fails:", fails, "\n")
saveRDS(list(detect_right=detect_right, other_false=other_false, fails=fails,
             N=N, shift=shift, target_item=target_item),
        "tools/simval/round1/dif/power_nonuniform_corrected.rds")
cat("total time:", as.numeric(Sys.time() - t_start, units = "secs"), "s\n")
