suppressWarnings(pkgload::load_all(".", quiet=TRUE))
set.seed(42)
d <- simulate_rasch(n_persons=600, n_items=10, model="PCM", n_categories=4,
                     difficulty=c(-2,2), threshold_spread=1.3, seed=42)
tr <- attr(d,"truth")
fit <- rasch(d, model="PCM")
shift <- mean(unlist(tr$thresholds))
tau_truth <- unlist(tr$thresholds) - shift
tau_est <- fit$thresholds$tau
cat("cor tau:", cor(tau_truth, tau_est), " rmse:", sqrt(mean((tau_truth-tau_est)^2)), "\n")
loc_truth <- tr$difficulty - shift
loc_est <- fit$items$location
cat("cor loc:", cor(loc_truth, loc_est), " rmse:", sqrt(mean((loc_truth-loc_est)^2)), "\n")
cat("any weak flagged:", any(fit$thresholds$weak), "\n")
print(fit$items[,c("item","location","se")])
