suppressWarnings(pkgload::load_all("."), quiet=TRUE))
source("tools/simval/round1/dif/gen_repeated.R")
t_start <- Sys.time()

run_within_null <- function(occ_levels, N, n_items, reps, seed0) {
  uni_flags <- c(); non_flags <- c(); any_uni <- c(); any_non <- c(); fails <- 0
  for (k in seq_len(reps)) {
    seed <- seed0 + k
    d <- gen_repeated(N = N, n_items = n_items, g_levels = 1,
                       occ_levels = occ_levels, seed = seed)
    fit <- tryCatch(rasch(d, id = "person", factors = "occasion"),
                     error = function(e) NULL)
    if (is.null(fit)) { fails <- fails + 1; next }
    da <- tryCatch(dif_anova(fit, within = "occasion"), error = function(e) NULL)
    if (is.null(da)) { fails <- fails + 1; next }
    s <- da$summary
    uni_flags <- c(uni_flags, s$uniform_DIF)
    non_flags <- c(non_flags, s$nonuniform_DIF)
    any_uni <- c(any_uni, any(s$uniform_DIF))
    any_non <- c(any_non, any(s$nonuniform_DIF))
  }
  list(uni_flags = uni_flags, non_flags = non_flags, any_uni = any_uni,
       any_non = any_non, fails = fails)
}

cat("=== within occasion, 2 levels, N=300 persons, I=8, reps=250 ===\n")
w2 <- run_within_null(2, 300, 8, 250, 30000)
cat("item-rep uniform flag rate:", mean(w2$uni_flags), "n=", length(w2$uni_flags), "\n")
cat("item-rep nonuniform flag rate:", mean(w2$non_flags), "n=", length(w2$non_flags), "\n")
cat("any-flag(uniform) per rep:", mean(w2$any_uni), "\n")
cat("any-flag(nonuniform) per rep:", mean(w2$any_non), "\n")
cat("fails:", w2$fails, "\n")
saveRDS(w2, "tools/simval/round1/dif/null_within2.rds")

cat("\n=== within occasion, 3 levels (GG), N=300 persons, I=8, reps=200 ===\n")
w3 <- run_within_null(3, 300, 8, 200, 40000)
cat("item-rep uniform flag rate:", mean(w3$uni_flags), "n=", length(w3$uni_flags), "\n")
cat("item-rep nonuniform flag rate:", mean(w3$non_flags), "n=", length(w3$non_flags), "\n")
cat("any-flag(uniform) per rep:", mean(w3$any_uni), "\n")
cat("any-flag(nonuniform) per rep:", mean(w3$any_non), "\n")
cat("fails:", w3$fails, "\n")
saveRDS(w3, "tools/simval/round1/dif/null_within3.rds")

cat("\ntotal time:", as.numeric(Sys.time() - t_start, units = "secs"), "s\n")
