suppressWarnings(pkgload::load_all(".", quiet=TRUE))
t_start <- Sys.time()

N <- 500; n_items <- 20; target_item <- "I06"

# --- power: strong single planted item, resolves first, minimal cascade ---
reps_power <- 80
first_right <- c(); n_splits <- c(); fails <- 0
for (k in seq_len(reps_power)) {
  seed <- 100000 + k
  d <- simulate_rasch(N, n_items, dif = list(items = target_item, uniform = 1.0),
                       n_groups = 2, seed = seed)
  fit <- tryCatch(rasch(d, id = "id", factors = "group"), error = function(e) NULL)
  if (is.null(fit)) { fails <- fails + 1; next }
  rd <- tryCatch(resolve_dif(fit), error = function(e) NULL)
  if (is.null(rd)) { fails <- fails + 1; next }
  sp <- rd$splits
  n_splits <- c(n_splits, if (is.null(sp)) 0 else nrow(sp))
  first_right <- c(first_right, !is.null(sp) && nrow(sp) >= 1 && sp$item[1] == target_item)
}
cat("=== resolve_dif power: planted 1.0 logit uniform DIF on", target_item, "===\n")
cat("n reps used:", length(first_right), " fails:", fails, "\n")
cat("rate planted item resolved FIRST:", mean(first_right), "\n")
cat("mean number of splits per run:", mean(n_splits), " (1 = no cascade)\n")
cat("distribution of n_splits:\n"); print(table(n_splits))

# --- null: no planted DIF at all, should stop with zero splits (no cascade) ---
reps_null <- 150
n_splits0 <- c(); fails0 <- 0
for (k in seq_len(reps_null)) {
  seed <- 110000 + k
  d <- simulate_rasch(N, n_items, n_groups = 2, seed = seed)
  fit <- tryCatch(rasch(d, id = "id", factors = "group"), error = function(e) NULL)
  if (is.null(fit)) { fails0 <- fails0 + 1; next }
  rd <- tryCatch(resolve_dif(fit), error = function(e) NULL)
  if (is.null(rd)) { fails0 <- fails0 + 1; next }
  sp <- rd$splits
  n_splits0 <- c(n_splits0, if (is.null(sp)) 0 else nrow(sp))
}
cat("\n=== resolve_dif null calibration: no planted DIF ===\n")
cat("n reps used:", length(n_splits0), " fails:", fails0, "\n")
cat("rate of zero splits (correctly stops immediately):", mean(n_splits0 == 0), "\n")
cat("mean number of splits:", mean(n_splits0), "\n")
cat("distribution of n_splits:\n"); print(table(n_splits0))

saveRDS(list(first_right=first_right, n_splits=n_splits, fails=fails,
             n_splits0=n_splits0, fails0=fails0),
        "tools/simval/round1/dif/resolve_dif.rds")
cat("\ntotal time:", as.numeric(Sys.time() - t_start, units = "secs"), "s\n")
