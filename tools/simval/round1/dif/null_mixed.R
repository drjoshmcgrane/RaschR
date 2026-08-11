suppressWarnings(pkgload::load_all(".", quiet=TRUE))
source("tools/simval/round1/dif/gen_repeated.R")
t_start <- Sys.time()

N <- 300; n_items <- 8; reps <- 250
terms_of_interest <- c("group", "occasion", "group:occasion")
flag_u <- setNames(vector("list", length(terms_of_interest)), terms_of_interest)
flag_n <- setNames(vector("list", length(terms_of_interest)), terms_of_interest)
any_u <- setNames(vector("list", length(terms_of_interest)), terms_of_interest)
any_n <- setNames(vector("list", length(terms_of_interest)), terms_of_interest)
fails <- 0

for (k in seq_len(reps)) {
  seed <- 50000 + k
  d <- gen_repeated(N = N, n_items = n_items, g_levels = 2, occ_levels = 2,
                     seed = seed)
  fit <- tryCatch(rasch(d, id = "person", factors = c("group", "occasion")),
                   error = function(e) NULL)
  if (is.null(fit)) { fails <- fails + 1; next }
  da <- tryCatch(dif_anova(fit, within = "occasion", effects = "factorial"),
                  error = function(e) NULL)
  if (is.null(da)) { fails <- fails + 1; next }
  s <- da$summary
  for (tt in terms_of_interest) {
    sel <- s$term == tt
    if (!any(sel)) next
    flag_u[[tt]] <- c(flag_u[[tt]], s$uniform_DIF[sel])
    flag_n[[tt]] <- c(flag_n[[tt]], s$nonuniform_DIF[sel])
    any_u[[tt]] <- c(any_u[[tt]], any(s$uniform_DIF[sel]))
    any_n[[tt]] <- c(any_n[[tt]], any(s$nonuniform_DIF[sel]))
  }
}

for (tt in terms_of_interest) {
  cat(sprintf("term %s: uniform item-rep rate=%.4f (n=%d), nonuniform item-rep rate=%.4f (n=%d), any-uniform/rep=%.4f, any-nonuniform/rep=%.4f\n",
              tt, mean(flag_u[[tt]]), length(flag_u[[tt]]),
              mean(flag_n[[tt]]), length(flag_n[[tt]]),
              mean(any_u[[tt]]), mean(any_n[[tt]])))
}
cat("fails:", fails, "/", reps, "\n")
saveRDS(list(flag_u=flag_u, flag_n=flag_n, any_u=any_u, any_n=any_n, fails=fails),
        "tools/simval/round1/dif/null_mixed.rds")
cat("total time:", as.numeric(Sys.time() - t_start, units = "secs"), "s\n")
