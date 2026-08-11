suppressWarnings(pkgload::load_all(".", quiet=TRUE))

set.seed(1)
d <- simulate_btl(n_objects = 10, n_judges = 14, reps_per_pair = 25, seed = 40001)
f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
info <- btl_information(f)

cat("== Information sanity: peaks at gap zero (dichotomous P(1-P)) ==\n")
cmp <- info$comparisons
# information should be a strictly decreasing function of |gap|
ord <- order(abs(cmp$gap))
mono_viol <- mean(diff(cmp$information[ord]) > 1e-9)  # fraction of "increases" as |gap| grows
cat("fraction of adjacent (sorted by |gap|) pairs where information INCREASES:",
    round(mono_viol, 4), "(should be ~0; P(1-P) is exactly monotone decreasing in |gap|)\n")
cat("correlation(information, -|gap|):", round(cor(cmp$information, -abs(cmp$gap)), 4), "\n")
cat("max information (should be at gap=0, value 0.25 for dichotomous P(1-P)):",
    round(max(cmp$information), 4), "\n")

cat("\n== Design information additivity ==\n")
manual_total <- sum(cmp$weight * cmp$information)
cat("total from btl_information:", round(info$total, 6),
    " manual sum(weight*information):", round(manual_total, 6),
    " match:", isTRUE(all.equal(info$total, manual_total)), "\n")
# each object's design information should equal the weighted sum of its own comparisons' info
obj_check <- sapply(info$objects$object, function(o) {
  sel <- cmp$object_a == o | cmp$object_b == o
  sum(cmp$weight[sel] * cmp$information[sel])
})
cat("per-object information matches manual recomputation:",
    isTRUE(all.equal(unname(obj_check), info$objects$information)), "\n")

cat("\n== se_naive vs fitted se: descriptive relationship ==\n")
cat("correlation(se, se_naive):", round(cor(info$objects$se, info$objects$se_naive), 4), "\n")
cat("mean se:", round(mean(info$objects$se),4), " mean se_naive:", round(mean(info$objects$se_naive),4), "\n")

cat("\n== Targeting sanity: btl_next_pairs favours near-neighbour / poorly-measured pairs ==\n")
nxt <- btl_next_pairs(f, n = 8)
all_pairs_gap <- info$pairs$gap
cat("mean |gap| of top-8 recommended pairs:", round(mean(abs(nxt$gap)), 4),
    " vs mean |gap| over ALL", nrow(info$pairs), "observed pairs:",
    round(mean(abs(all_pairs_gap)), 4), "\n")
cat("recommended pairs' mean gap smaller (more informative/near-neighbour) than average:",
    mean(abs(nxt$gap)) < mean(abs(all_pairs_gap)), "\n")

## sanity on an UNBALANCED design: one object heavily under-compared should
## rank near the top of next_pairs (its se/info is comparatively poor)
set.seed(2)
d2 <- simulate_btl(n_objects = 8, n_judges = 12, reps_per_pair = 25, seed = 40002)
# drop 90% of comparisons involving O1 to starve it of information
drop <- (d2$object_a == "O1" | d2$object_b == "O1") & runif(nrow(d2)) < 0.9
d2 <- d2[!drop, ]
f2 <- btl(d2, "object_a", "object_b", winner = "winner", judge = "judge")
info2 <- btl_information(f2)
nxt2 <- btl_next_pairs(f2, n = 10)
cat("\nunder-compared object O1: n_comparisons =", info2$objects$n_comparisons[info2$objects$object=="O1"],
    " (others average", round(mean(info2$objects$n_comparisons[info2$objects$object!="O1"]),1), ")\n")
cat("O1 appears in top-10 recommended next pairs:",
    sum(nxt2$object_a == "O1" | nxt2$object_b == "O1"), "times (of 10)\n")
cat("O1 has the lowest design information among objects:",
    info2$objects$object[which.min(info2$objects$information)] == "O1", "\n")
