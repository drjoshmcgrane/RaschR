suppressWarnings(pkgload::load_all(".", quiet=TRUE))
set.seed(1)

## --- (1a) Point recovery: one largeish dataset, full round robin ----------
d <- simulate_btl(n_objects = 12, n_judges = 16, reps_per_pair = 30, seed = 11)
f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
rec <- sim_recovery(f, d)
cat("== Recovery (single draw, K=12,J=16,reps=30) ==\n")
print(rec$summary)

## --- (1b) SE calibration across replicates, judge-clustered SE ------------
n_rep <- 200
batch <- sim_replicate(simulate_btl, n_rep, n_objects = 10, n_judges = 14,
                        reps_per_pair = 20, seed = 100)
# sim_apply is built for scalar-per-replicate statistics; here we need each
# replicate's whole location/se vector, so we do the resilient lapply by
# hand on top of the sim_replicate batch (same tryCatch discipline).
fits <- lapply(batch, function(dd) tryCatch({
  fit <- btl(dd, "object_a", "object_b", winner = "winner", judge = "judge")
  list(loc = fit$objects$location, se = fit$objects$se,
       obj = fit$objects$object, converged = fit$converged)
}, error = function(e) NULL))
ok <- vapply(fits, function(x) is.list(x) && isTRUE(x$converged), TRUE)
cat("\n== SE calibration: converged", sum(ok), "/", n_rep, "==\n")
objs <- fits[[which(ok)[1]]]$obj
LOC <- t(vapply(fits[ok], function(x) x$loc[match(objs, x$obj)], numeric(length(objs))))
SE  <- t(vapply(fits[ok], function(x) x$se[match(objs, x$obj)], numeric(length(objs))))
emp_sd <- apply(LOC, 2, sd)
mean_se <- colMeans(SE)
ratio <- emp_sd / mean_se
cat("per-object empirical SD / mean reported clustered SE:\n")
print(data.frame(object = objs, emp_sd = round(emp_sd,4), mean_se = round(mean_se,4),
                  ratio = round(ratio,3)))
cat("overall mean ratio:", round(mean(ratio), 3),
    " MC error of mean ratio (~1/sqrt(2*(n-1))):", round(1/sqrt(2*(sum(ok)-1)), 3), "\n")

## --- (1c) btl <-> pcml equivalence spot check ------------------------------
# A dichotomous L-item Rasch dataset: pcml()'s pairwise-conditional
# likelihood sums, over every item pair, the conditional likelihood of who
# (within that pair) got the item right, given exactly one of the two was
# right -- which IS a BTL comparison between the two items. Build the
# equivalent long-format BTL "tournament" (one winner row per person per
# informative item pair) and confirm btl()'s conditional MLE reproduces
# pcml()'s item locations to numerical precision.
set.seed(42)
N <- 3000
delta <- c(I1 = -0.9, I2 = -0.1, I3 = 0.3, I4 = 0.7)
theta <- rnorm(N, 0, 1.4)
Y <- sapply(delta, function(dd) rbinom(N, 1, plogis(theta - dd)))
colnames(Y) <- names(delta)
pc <- pcml(Y)

items <- names(delta); L <- length(items)
pairs <- t(utils::combn(items, 2))
rows <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(r) {
  ia <- pairs[r, 1]; ib <- pairs[r, 2]
  inf <- Y[, ia] != Y[, ib]
  if (!any(inf)) return(NULL)
  data.frame(a = ia, b = ib, person = which(inf),
             win = ifelse(Y[inf, ia] == 1, ia, ib),
             stringsAsFactors = FALSE)
}))
# unclustered: each row treated as an independent trial (loses the fact that
# the same person supplies rows for several item pairs)
fb <- btl(rows, "a", "b", winner = "win")
# person-clustered: acknowledges that a person's several pairwise verdicts
# are not independent trials -- this is the estimating-equation structure
# pcml's own sandwich is built on (it sums each person's contribution once)
fbj <- btl(rows, "a", "b", winner = "win", judge = "person")

# pcml's tau is a DIFFICULTY (p_correct = plogis(theta - tau)): given exactly
# one of a pair correct, P(item i correct) = plogis(tau_j - tau_i), i.e. the
# "wins" convention of BTL (P(a beats b) = plogis(beta_a - beta_b)) has
# beta = -tau. thr rows are already in item order (id 1..L = I1..I4).
loc_pcml_beta <- setNames(-pc$thr$tau, items)
loc_btl  <- setNames(fb$objects$location[match(items, fb$objects$object)], items)
cat("\n== btl <-> pcml equivalence spot check (L=4 items) ==\n")
print(data.frame(item = items, pcml_as_beta = round(loc_pcml_beta, 6),
                  btl = round(loc_btl, 6),
                  abs_diff = round(abs(loc_pcml_beta - loc_btl), 8)))
cat("max abs difference in locations:", round(max(abs(loc_pcml_beta - loc_btl)), 8), "\n")
se_pcml <- pc$thr$se
se_btl_unclustered <- fb$objects$se[match(items, fb$objects$object)]
se_btl_clustered <- fbj$objects$se[match(items, fbj$objects$object)]
cat("correlation(pcml se, btl se unclustered):", round(cor(se_pcml, se_btl_unclustered), 4), "\n")
cat("correlation(pcml se, btl se person-clustered):", round(cor(se_pcml, se_btl_clustered), 4), "\n")
print(data.frame(item = items, se_pcml = round(se_pcml, 4),
                  se_btl_unclustered = round(se_btl_unclustered, 4),
                  se_btl_clustered = round(se_btl_clustered, 4)))
