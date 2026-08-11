suppressWarnings(pkgload::load_all(".", quiet=TRUE))

## --- (2a) graded threshold recovery: free thresholds -----------------------
set.seed(21)
dg <- simulate_btl(n_objects = 10, n_judges = 14, reps_per_pair = 30,
                    model = "graded", n_categories = 4, seed = 21)
fg_free <- btl(dg, "object_a", "object_b", response = "response",
               judge = "judge", thresholds = "free")
rec_free <- sim_recovery(fg_free, dg)
cat("== Graded object-location recovery, thresholds='free' ==\n")
print(rec_free$summary)
cat("thresholds (tau) estimated (m=3, symmetric around 0 by construction):\n")
print(fg_free$thresholds)

fg_pc <- btl(dg, "object_a", "object_b", response = "response",
             judge = "judge", thresholds = "pc")
rec_pc <- sim_recovery(fg_pc, dg)
cat("\n== Graded object-location recovery, thresholds='pc' ==\n")
print(rec_pc$summary)
cat("free vs pc location correlation:",
    round(cor(fg_free$objects$location, fg_pc$objects$location), 6), "\n")

## Threshold recovery across replicates (spread component ~ planted spread)
n_rep <- 120
batch <- sim_replicate(simulate_btl, n_rep, n_objects = 8, n_judges = 12,
                        reps_per_pair = 25, model = "graded", n_categories = 4,
                        seed = 500)
spread_free <- rep(NA_real_, n_rep); spread_pc <- rep(NA_real_, n_rep)
for (i in seq_len(n_rep)) {
  dd <- batch[[i]]
  ok <- tryCatch({
    ff <- btl(dd, "object_a", "object_b", response = "response", judge = "judge", thresholds = "free")
    fp <- btl(dd, "object_a", "object_b", response = "response", judge = "judge", thresholds = "pc")
    spread_free[i] <- ff$components$estimate[ff$components$component == "spread"]
    spread_pc[i]   <- fp$components$estimate[fp$components$component == "spread"]
    TRUE
  }, error = function(e) FALSE)
}
# planted thresholds: .sim_thresholds(0, m, 1.2) -- recover the "true" spread
# by rebuilding the same generator with delta=0, m=3, spread=1.2 and taking
# its own spread component definition (v1-projection) for an apples-to-apples target
true_tau <- rasch:::.sim_thresholds(0, 3, 1.2)
v1 <- seq_len(3) - (3+1)/2; v1 <- v1/sqrt(sum(v1^2))
true_spread <- sum(v1 * true_tau)
cat("\n== Threshold spread-component recovery across", n_rep, "replicates ==\n")
cat("planted spread:", round(true_spread,4), "\n")
cat("free: mean est", round(mean(spread_free,na.rm=TRUE),4),
    " bias", round(mean(spread_free,na.rm=TRUE)-true_spread,4),
    " rmse", round(sqrt(mean((spread_free-true_spread)^2,na.rm=TRUE)),4), "\n")
cat("pc:   mean est", round(mean(spread_pc,na.rm=TRUE),4),
    " bias", round(mean(spread_pc,na.rm=TRUE)-true_spread,4),
    " rmse", round(sqrt(mean((spread_pc-true_spread)^2,na.rm=TRUE)),4), "\n")

## --- (2b) winner+margin entry == response entry ---------------------------
set.seed(22)
beta <- c(A = -1.2, B = -0.4, C = 0.3, D = 0.9, E = 1.6)
pr <- t(utils::combn(names(beta), 2))
n_each <- 60
d <- data.frame(a = rep(pr[,1], each = n_each), b = rep(pr[,2], each = n_each),
                stringsAsFactors = FALSE)
tau <- c(-1.3, 0, 1.3)  # 4-category symmetric
P <- t(vapply(seq_len(nrow(d)), function(r)
  item_moments(beta[d$a[r]] - beta[d$b[r]], tau)$P, numeric(4)))
resp <- vapply(seq_len(nrow(d)), function(r) sample(0:3, 1, prob = P[r,]), 0L)
d$response <- resp
# response 0,1 = b wins (margins "big","small"); 2,3 = a wins (small,big)
d$winner <- ifelse(resp >= 2, d$a, d$b)
d$margin <- factor(ifelse(resp %in% c(0,3), "big", "small"),
                    levels = c("small","big"), ordered = TRUE)

f_resp <- btl(d, "a", "b", response = "response")
f_marg <- btl(d, "a", "b", winner = "winner", margin = "margin")
cat("\n== winner+margin entry vs response entry ==\n")
cat("max abs diff in object locations:",
    round(max(abs(f_resp$objects$location - f_marg$objects$location[
      match(f_resp$objects$object, f_marg$objects$object)])), 10), "\n")
cat("max abs diff in thresholds (spread only comparable term):",
    round(max(abs(f_resp$components$estimate - f_marg$components$estimate)), 10), "\n")
cat("loglik resp:", round(f_resp$loglik,6), " loglik margin:", round(f_marg$loglik,6), "\n")
