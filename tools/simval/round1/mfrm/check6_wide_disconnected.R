source("tools/simval/round1/mfrm/helpers.R")

cat("=== Wide vs long input equivalence ===\n")
d_long <- simulate_mfrm(120, 5, 6, n_categories = 4, rater_severity_sd = 0.7, seed = 6001)
mf_long <- rasch_mfrm(d_long, person = "person", item = "item", score = "score", facets = "rater")

# reshape to wide: one row per person x rater, one column per item
d_wide <- reshape(d_long, idvar = c("person", "rater"), timevar = "item",
                  direction = "wide")
item_cols <- paste0("score.", sort(unique(d_long$item)))
names(d_wide)[match(item_cols, names(d_wide))] <- sort(unique(d_long$item))
mf_wide <- rasch_mfrm(d_wide, person = "person", facets = "rater",
                      items = sort(unique(d_long$item)))

sev_long <- mf_long$facet_effects$rater[order(mf_long$facet_effects$rater$level), ]
sev_wide <- mf_wide$facet_effects$rater[order(mf_wide$facet_effects$rater$level), ]
max_abs_diff_sev <- max(abs(sev_long$severity - sev_wide$severity))
max_abs_diff_se <- max(abs(sev_long$se - sev_wide$se))

loc_long <- mf_long$item_effects[order(mf_long$item_effects$item), ]
loc_wide <- mf_wide$item_effects[order(mf_wide$item_effects$item), ]
max_abs_diff_loc <- max(abs(loc_long$location - loc_wide$location))

cat(sprintf("  max|severity(long) - severity(wide)| = %.2e\n", max_abs_diff_sev))
cat(sprintf("  max|se(long) - se(wide)| = %.2e\n", max_abs_diff_se))
cat(sprintf("  max|item location(long) - item location(wide)| = %.2e\n", max_abs_diff_loc))
cat(sprintf("  loglik long = %.6f, wide = %.6f\n", mf_long$est$loglik, mf_wide$est$loglik))

cat("\n=== Disconnected judging plan: identification guard must refuse ===\n")
d_disc <- make_disconnected(d_long)
cat(sprintf("  disconnected design: %d rows (of %d), raters/persons split into 2 non-overlapping blocks\n",
            nrow(d_disc), nrow(d_long)))
res <- tryCatch(rasch_mfrm(d_disc, person = "person", item = "item", score = "score",
                           facets = "rater"), error = function(e) e)
refused <- inherits(res, "error")
cat(sprintf("  refused with error: %s\n", refused))
if (refused) cat("  message:", conditionMessage(res), "\n")

## sanity: the SAME two-block partition, but with one linking person added
## back who saw both rater blocks, should NOT be refused (control check)
cat("\n=== Control: same split + one linking person -> should connect and fit ===\n")
raters <- sort(unique(d_long$rater)); persons <- sort(unique(d_long$person))
half_r <- ceiling(length(raters) / 2)
blockA_r <- raters[seq_len(half_r)]; blockB_r <- raters[(half_r + 1):length(raters)]
linking_person <- persons[1]
d_link <- d_disc
extra <- d_long[d_long$person == linking_person & d_long$rater %in% blockB_r, ]
d_link <- rbind(d_link, extra)
res2 <- tryCatch(rasch_mfrm(d_link, person = "person", item = "item", score = "score",
                            facets = "rater"), error = function(e) e)
connected_ok <- !inherits(res2, "error")
cat(sprintf("  fits without error once linked: %s\n", connected_ok))

saveRDS(list(max_abs_diff_sev = max_abs_diff_sev, max_abs_diff_se = max_abs_diff_se,
            max_abs_diff_loc = max_abs_diff_loc, refused = refused,
            connected_ok = connected_ok), "tools/simval/round1/mfrm/check6_results.rds")
cat("\n=== DONE check6 ===\n")
