suppressWarnings(pkgload::load_all("."), quiet=TRUE))
source("tools/simval/round1/dif/gen_repeated.R")
t_start <- Sys.time()

N <- 300; n_items <- 8; reps <- 250; target_item_idx <- 4
target_item <- sprintf("I%02d", target_item_idx)

detect_right <- c(); detect_any <- c(); other_false <- c(); fails <- 0
for (k in seq_len(reps)) {
  seed <- 80000 + k
  # long-format data: one row per person per occasion, then the documented
  # stack_data() workflow -> rasch(id, factors, items) -> dif_anova(within=)
  long <- gen_repeated(N = N, n_items = n_items, g_levels = 1, occ_levels = 2,
                        occ_dif_item = target_item_idx, occ_dif_shift = 0.7,
                        seed = seed)
  item_cols <- sprintf("I%02d", seq_len(n_items))
  stacked <- stack_data(long, person = "person", time = "occasion",
                         items = item_cols)
  fit <- tryCatch(rasch(stacked, id = "id", factors = "time"),
                   error = function(e) NULL)
  if (is.null(fit)) { fails <- fails + 1; next }
  da <- tryCatch(dif_anova(fit, within = "time"), error = function(e) NULL)
  if (is.null(da)) { fails <- fails + 1; next }
  s <- da$summary
  right <- s$uniform_DIF[s$item == target_item]
  detect_right <- c(detect_right, isTRUE(right))
  detect_any <- c(detect_any, any(s$uniform_DIF))
  other <- s$uniform_DIF[s$item != target_item]
  other_false <- c(other_false, other)
}
cat("=== power: within-person occasion DIF 0.7 logits on", target_item, "via stack_data workflow ===\n")
cat("detection rate (right item, uniform 'time' term):", mean(detect_right), " n=", length(detect_right), "\n")
cat("any-uniform-flag rate:", mean(detect_any), "\n")
cat("false uniform flag rate on OTHER items:", mean(other_false), " n=", length(other_false), "\n")
cat("fails:", fails, "\n")
saveRDS(list(detect_right=detect_right, other_false=other_false, fails=fails),
        "tools/simval/round1/dif/power_within.rds")
cat("total time:", as.numeric(Sys.time() - t_start, units = "secs"), "s\n")
