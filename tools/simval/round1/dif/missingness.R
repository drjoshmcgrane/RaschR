suppressWarnings(pkgload::load_all(".", quiet=TRUE))
source("tools/simval/round1/dif/gen_repeated.R")
t_start <- Sys.time()

# ---------------------------------------------------------------------
# (a) 25% MCAR, between 2-level factor, null (no planted DIF)
# NOTE: simulate_rasch(..., missing = 0.25) makes rasch() ~7s/fit (vs
# ~0.25s complete) -- a real fixed cost of the missing-data estimation
# path, not a bug (converges in 4 iterations either way; see
# bench_missing2.R). Reps cut to 40 to fit the time budget; MC error
# reported accordingly.
# ---------------------------------------------------------------------
N <- 500; n_items <- 20; reps <- 40
uni_flags <- c(); non_flags <- c(); any_uni <- c(); any_non <- c(); fails <- 0
for (k in seq_len(reps)) {
  seed <- 120000 + k
  d <- simulate_rasch(N, n_items, n_groups = 2, missing = 0.25, seed = seed)
  fit <- tryCatch(rasch(d, id = "id", factors = "group"), error = function(e) NULL)
  if (is.null(fit)) { fails <- fails + 1; next }
  da <- tryCatch(dif_anova(fit), error = function(e) NULL)
  if (is.null(da)) { fails <- fails + 1; next }
  s <- da$summary
  uni_flags <- c(uni_flags, s$uniform_DIF)
  non_flags <- c(non_flags, s$nonuniform_DIF)
  any_uni <- c(any_uni, any(s$uniform_DIF))
  any_non <- c(any_non, any(s$nonuniform_DIF))
}
cat("=== (a) 25% MCAR, between 2-level, null, N=500 I=20, reps=", reps, "===\n")
cat("item-rep uniform flag rate:", mean(uni_flags), "n=", length(uni_flags), "\n")
cat("item-rep nonuniform flag rate:", mean(non_flags), "n=", length(non_flags), "\n")
cat("any-flag(uniform)/rep:", mean(any_uni), " any-flag(nonuniform)/rep:", mean(any_non), "\n")
cat("fails:", fails, "\n\n")
res_a <- list(uni_flags=uni_flags, non_flags=non_flags, any_uni=any_uni, any_non=any_non, fails=fails)

# ---------------------------------------------------------------------
# (b) Differential missingness by group (structural): group g1 skips the
# last 5 items entirely (never sees them); group g2 answers all items.
# No planted DIF. Check the UNAFFECTED items (1:15, answered by both
# groups) show no spurious DIF; confirm the skipped items are excluded
# (not falsely flagged) rather than crashing.
# ---------------------------------------------------------------------
n_items_b <- 20; skip_idx <- 16:20; N_b <- 600; reps_b <- 200
uni_flags_b <- c(); non_flags_b <- c(); skipped_present <- c(); fails_b <- 0
for (k in seq_len(reps_b)) {
  seed <- 130000 + k
  d <- simulate_rasch(N_b, n_items_b, n_groups = 2, seed = seed)
  d[d$group == "g1", sprintf("I%02d", skip_idx)] <- NA
  fit <- tryCatch(rasch(d, id = "id", factors = "group"), error = function(e) NULL)
  if (is.null(fit)) { fails_b <- fails_b + 1; next }
  da <- tryCatch(dif_anova(fit), error = function(e) NULL)
  if (is.null(da)) { fails_b <- fails_b + 1; next }
  s <- da$summary
  unaff <- s[!s$item %in% sprintf("I%02d", skip_idx), ]
  uni_flags_b <- c(uni_flags_b, unaff$uniform_DIF)
  non_flags_b <- c(non_flags_b, unaff$nonuniform_DIF)
  skipped_present <- c(skipped_present,
                        any(s$item %in% sprintf("I%02d", skip_idx)))
}
cat("=== (b) differential item-missingness by group (structural), null, reps=", reps_b, "===\n")
cat("n reps used:", reps_b - fails_b, " fails:", fails_b, "\n")
cat("unaffected-item uniform flag rate:", mean(uni_flags_b), "n=", length(uni_flags_b), "\n")
cat("unaffected-item nonuniform flag rate:", mean(non_flags_b), "n=", length(non_flags_b), "\n")
cat("rate group-skipped items appear (falsely) in DIF table (should be 0/near-0, single-level items are dropped):", mean(skipped_present), "\n\n")
res_b <- list(uni_flags_b=uni_flags_b, non_flags_b=non_flags_b,
              skipped_present=skipped_present, fails_b=fails_b)

# ---------------------------------------------------------------------
# (c) Structural / differential missingness in the mixed (within-person)
# design -- the "round-12 guarantee": a COMMON (non-DIF) occasion shift on
# one item at T2, plus group g2 disproportionately missing the T2 wave
# (dropout). No planted uniform/non-uniform GROUP DIF. Confirm this does
# not manufacture spurious group (between) or group:occasion DIF.
# ---------------------------------------------------------------------
N_c <- 320; n_items_c <- 8; reps_c <- 200
common_shift_item <- 5
uni_g <- c(); non_g <- c(); uni_occ <- c(); non_occ_g_occ <- c(); fails_c <- 0
for (k in seq_len(reps_c)) {
  seed <- 140000 + k
  set.seed(seed)
  long <- gen_repeated(N = N_c, n_items = n_items_c, g_levels = 2, occ_levels = 2,
                        occ_dif_item = common_shift_item, occ_dif_shift = 1.0,
                        seed = seed)
  # structural dropout: group g2 misses T2 with prob 0.55; group g1 with
  # prob 0.10 (differential, unrelated to the group's actual responses)
  is_g2_t2 <- long$group == "g2" & long$occasion == "T2"
  is_g1_t2 <- long$group == "g1" & long$occasion == "T2"
  drop_g2 <- is_g2_t2 & (runif(nrow(long)) < 0.55)
  drop_g1 <- is_g1_t2 & (runif(nrow(long)) < 0.10)
  item_cols <- sprintf("I%02d", seq_len(n_items_c))
  long[drop_g2 | drop_g1, item_cols] <- NA
  # drop rows that are now entirely NA (person not observed at that
  # occasion at all) -- rasch() needs at least one observed item per row
  keep <- rowSums(!is.na(long[, item_cols])) > 0
  long <- long[keep, ]
  fit <- tryCatch(rasch(long, id = "person", factors = c("group", "occasion")),
                   error = function(e) NULL)
  if (is.null(fit)) { fails_c <- fails_c + 1; next }
  da <- tryCatch(dif_anova(fit, within = "occasion", effects = "factorial"),
                  error = function(e) NULL)
  if (is.null(da)) { fails_c <- fails_c + 1; next }
  s <- da$summary
  g_row <- s[s$term == "group", ]
  occ_row <- s[s$term == "occasion", ]
  int_row <- s[s$term == "group:occasion", ]
  uni_g <- c(uni_g, g_row$uniform_DIF)
  non_g <- c(non_g, g_row$nonuniform_DIF)
  uni_occ <- c(uni_occ, occ_row$uniform_DIF[occ_row$item == sprintf("I%02d", common_shift_item)])
  non_occ_g_occ <- c(non_occ_g_occ, int_row$uniform_DIF)  # group:occasion "uniform" IS the interaction test
}
cat("=== (c) mixed design, differential occasion dropout by group (round-12 guarantee), reps=", reps_c, "===\n")
cat("n reps used:", reps_c - fails_c, " fails:", fails_c, "\n")
cat("spurious GROUP uniform-DIF flag rate (should be ~alpha, not inflated):", mean(uni_g), "n=", length(uni_g), "\n")
cat("spurious GROUP nonuniform-DIF flag rate:", mean(non_g), "n=", length(non_g), "\n")
cat("common-shift item correctly flagged on OCCASION term (sanity, expect high):", mean(uni_occ, na.rm=TRUE), "n=", length(uni_occ), "\n")
cat("spurious GROUP:OCCASION interaction flag rate:", mean(non_occ_g_occ), "n=", length(non_occ_g_occ), "\n")

saveRDS(list(a=res_a, b=res_b,
             c=list(uni_g=uni_g, non_g=non_g, uni_occ=uni_occ,
                    non_occ_g_occ=non_occ_g_occ, fails_c=fails_c)),
        "tools/simval/round1/dif/missingness.rds")
cat("\ntotal time:", as.numeric(Sys.time() - t_start, units = "secs"), "s\n")
