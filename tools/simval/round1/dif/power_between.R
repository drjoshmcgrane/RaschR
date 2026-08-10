suppressWarnings(pkgload::load_all("."), quiet=TRUE))
t_start <- Sys.time()

N <- 500; n_items <- 20; reps <- 250
target_item <- "I06"

# --- planted uniform DIF (0.7 logits) ---------------------------------
detect_right_u <- c(); detect_any_u <- c(); other_false_u <- c(); fails_u <- 0
for (k in seq_len(reps)) {
  seed <- 60000 + k
  d <- simulate_rasch(N, n_items, dif = list(items = target_item, uniform = 0.7),
                       n_groups = 2, seed = seed)
  fit <- tryCatch(rasch(d, id = "id", factors = "group"), error = function(e) NULL)
  if (is.null(fit)) { fails_u <- fails_u + 1; next }
  da <- tryCatch(dif_anova(fit), error = function(e) NULL)
  if (is.null(da)) { fails_u <- fails_u + 1; next }
  s <- da$summary
  right <- s$uniform_DIF[s$item == target_item]
  detect_right_u <- c(detect_right_u, isTRUE(right))
  detect_any_u <- c(detect_any_u, any(s$uniform_DIF))
  other <- s$uniform_DIF[s$item != target_item]
  other_false_u <- c(other_false_u, other)
}
cat("=== power: uniform DIF 0.7 logits on", target_item, "(between, 2-level) ===\n")
cat("detection rate (right item, uniform term):", mean(detect_right_u), " n=", length(detect_right_u), "\n")
cat("any-uniform-flag rate:", mean(detect_any_u), "\n")
cat("false uniform flag rate on OTHER items:", mean(other_false_u), " n=", length(other_false_u), "\n")
cat("fails:", fails_u, "\n\n")

# --- planted non-uniform DIF (discrimination shift) --------------------
detect_right_n <- c(); detect_any_n <- c(); other_false_n <- c(); fails_n <- 0
for (k in seq_len(reps)) {
  seed <- 70000 + k
  d <- simulate_rasch(N, n_items, dif = list(items = target_item, nonuniform = 0.6),
                       n_groups = 2, seed = seed)
  fit <- tryCatch(rasch(d, id = "id", factors = "group"), error = function(e) NULL)
  if (is.null(fit)) { fails_n <- fails_n + 1; next }
  da <- tryCatch(dif_anova(fit), error = function(e) NULL)
  if (is.null(da)) { fails_n <- fails_n + 1; next }
  s <- da$summary
  right <- s$nonuniform_DIF[s$item == target_item]
  detect_right_n <- c(detect_right_n, isTRUE(right))
  detect_any_n <- c(detect_any_n, any(s$nonuniform_DIF))
  other <- s$nonuniform_DIF[s$item != target_item]
  other_false_n <- c(other_false_n, other)
}
cat("=== power: non-uniform DIF (discrimination +0.6) on", target_item, "(between, 2-level) ===\n")
cat("detection rate (right item, nonuniform term):", mean(detect_right_n), " n=", length(detect_right_n), "\n")
cat("any-nonuniform-flag rate:", mean(detect_any_n), "\n")
cat("false nonuniform flag rate on OTHER items:", mean(other_false_n), " n=", length(other_false_n), "\n")
cat("fails:", fails_n, "\n")

saveRDS(list(detect_right_u=detect_right_u, other_false_u=other_false_u,
             detect_right_n=detect_right_n, other_false_n=other_false_n,
             fails_u=fails_u, fails_n=fails_n),
        "tools/simval/round1/dif/power_between.rds")
cat("\ntotal time:", as.numeric(Sys.time() - t_start, units = "secs"), "s\n")
