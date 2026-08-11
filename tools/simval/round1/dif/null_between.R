suppressWarnings(pkgload::load_all(".", quiet=TRUE))
t_start <- Sys.time()

run_between_null <- function(n_groups, n_persons, n_items, reps, seed0) {
  uni_flags <- c(); non_flags <- c(); any_uni <- c(); any_non <- c()
  fails <- 0
  for (k in seq_len(reps)) {
    seed <- seed0 + k
    d <- tryCatch(simulate_rasch(n_persons, n_items, n_groups = n_groups,
                                  seed = seed), error = function(e) NULL)
    if (is.null(d)) { fails <- fails + 1; next }
    fit <- tryCatch(rasch(d, id = "id", factors = "group"),
                     error = function(e) NULL)
    if (is.null(fit)) { fails <- fails + 1; next }
    da <- tryCatch(dif_anova(fit), error = function(e) NULL)
    if (is.null(da)) { fails <- fails + 1; next }
    s <- da$summary
    uni_flags <- c(uni_flags, s$uniform_DIF)
    non_flags <- c(non_flags, s$nonuniform_DIF)
    any_uni <- c(any_uni, any(s$uniform_DIF))
    any_non <- c(any_non, any(s$nonuniform_DIF))
  }
  list(uni_flags = uni_flags, non_flags = non_flags,
       any_uni = any_uni, any_non = any_non, fails = fails,
       n_reps = reps - fails)
}

cat("=== between 2-level, N=500, I=20, reps=250 ===\n")
r2 <- run_between_null(n_groups = 2, n_persons = 500, n_items = 20,
                        reps = 250, seed0 = 10000)
cat("item-rep uniform flag rate:", mean(r2$uni_flags), " n=", length(r2$uni_flags), "\n")
cat("item-rep nonuniform flag rate:", mean(r2$non_flags), " n=", length(r2$non_flags), "\n")
cat("any-flag(uniform) per rep rate:", mean(r2$any_uni), "\n")
cat("any-flag(nonuniform) per rep rate:", mean(r2$any_non), "\n")
cat("fails:", r2$fails, "/", 250, "\n")
saveRDS(r2, "tools/simval/round1/dif/null_between2.rds")

cat("\n=== between 3-level, N=750, I=20, reps=250 ===\n")
r3 <- run_between_null(n_groups = 3, n_persons = 750, n_items = 20,
                        reps = 250, seed0 = 20000)
cat("item-rep uniform flag rate:", mean(r3$uni_flags), " n=", length(r3$uni_flags), "\n")
cat("item-rep nonuniform flag rate:", mean(r3$non_flags), " n=", length(r3$non_flags), "\n")
cat("any-flag(uniform) per rep rate:", mean(r3$any_uni), "\n")
cat("any-flag(nonuniform) per rep rate:", mean(r3$any_non), "\n")
cat("fails:", r3$fails, "/", 250, "\n")
saveRDS(r3, "tools/simval/round1/dif/null_between3.rds")

cat("\ntotal time:", as.numeric(Sys.time() - t_start, units = "secs"), "s\n")
