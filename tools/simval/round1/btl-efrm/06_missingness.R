suppressWarnings(pkgload::load_all(".", quiet = TRUE))
t_start <- Sys.time()

## (a) ~25% MCAR: randomly delete 25% of comparison rows, refit, check recovery
n_rep <- 40
panel_units <- c(0.8, 1.25); set_units <- c(1, 1.4); set_origins <- c(0, 0.5)
lphi2 <- lalpha2 <- kappa2 <- numeric(n_rep)
n_fail <- 0
for (r in seq_len(n_rep)) {
  seed <- 8000 + r
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                          n_judges_per_panel = 6, reps_within = 20, reps_cross = 20,
                          panel_units = panel_units, set_units = set_units,
                          set_origins = set_origins, object_sd = 1, seed = seed)
  tr <- attr(d, "truth")
  set.seed(seed + 99999)
  drop <- runif(nrow(d)) < 0.25            # ~25% MCAR
  d_mcar <- d[!drop, ]
  fit <- tryCatch(btl_efrm(d_mcar, "object_a", "object_b", "winner", "judge", "panel",
                            object_sets = tr$object_sets, se_method = "conditional"),
                   error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged)) { n_fail <- n_fail + 1
    lphi2[r] <- lalpha2[r] <- kappa2[r] <- NA; next }
  lphi2[r]   <- log(fit$phi_table$phi[fit$phi_table$panel == "panel2"])
  lalpha2[r] <- log(fit$alpha_table$alpha[fit$alpha_table$set == "set2"])
  kappa2[r]  <- fit$kappa_table$kappa[fit$kappa_table$set == "set2"]
}
true_lphi2 <- log(panel_units[2] / exp(mean(log(panel_units))))
true_lalpha2 <- log(set_units[2] / set_units[1])
true_kappa2 <- set_origins[2] - set_origins[1]

sumr <- function(est, true, label) {
  ok <- is.finite(est)
  data.frame(parameter = label, n = sum(ok), true = true,
             mean_est = mean(est[ok]), bias = mean(est[ok]) - true,
             rmse = sqrt(mean((est[ok] - true)^2)),
             mc_se_bias = sd(est[ok]) / sqrt(sum(ok)))
}
cat("---- 25% MCAR recovery (n_rep=", n_rep, ", n_fail=", n_fail, ") ----\n", sep = "")
print(rbind(sumr(lphi2, true_lphi2, "log phi[panel2]"),
            sumr(lalpha2, true_lalpha2, "log alpha[set2]"),
            sumr(kappa2, true_kappa2, "kappa[set2]")), row.names = FALSE)

## (b) structural missing / booklet design: a "linking-only" panel that NEVER
## judges within-set pairs (only cross-set) -- it has no representation in the
## stage-1 phi-reconciliation network and should be correctly refused.
set.seed(42)
mk_pairs <- function(objs, n) {
  pr <- t(utils::combn(objs, 2))
  data.frame(object_a = rep(pr[, 1], n), object_b = rep(pr[, 2], n), stringsAsFactors = FALSE)
}
withinA <- mk_pairs(c("A1", "A2", "A3", "A4"), 20); withinA$judge_pool <- "panelA"
withinB <- mk_pairs(c("B1", "B2", "B3", "B4"), 20); withinB$judge_pool <- "panelB"
crossC <- expand.grid(object_a = c("A1", "A2", "A3", "A4"),
                       object_b = c("B1", "B2", "B3", "B4"), stringsAsFactors = FALSE)
crossC <- crossC[rep(seq_len(nrow(crossC)), 10), ]; crossC$judge_pool <- "panelC"
allbk <- rbind(withinA, withinB, crossC)
allbk$winner <- ifelse(runif(nrow(allbk)) < 0.5, allbk$object_a, allbk$object_b)
# judges nested in their pool (booklet design: each judge sees only their pool's items)
allbk$judge <- paste0(allbk$judge_pool, "_J", sample(1:8, nrow(allbk), replace = TRUE))
panel_map <- setNames(sub("_J.*", "", unique(allbk$judge)), unique(allbk$judge))
out <- tryCatch(
  btl_efrm(allbk, "object_a", "object_b", "winner", "judge", panels = panel_map,
           object_sets = list(set1 = c("A1", "A2", "A3", "A4"),
                               set2 = c("B1", "B2", "B3", "B4")),
           se_method = "conditional"),
  error = function(e) e)
cat("\n---- structural/booklet: linking-only panel (never judges within-set) ----\n")
if (inherits(out, "error")) {
  cat("Correctly refused:", grepl("linked|connect", conditionMessage(out), ignore.case = TRUE), "\n")
  cat("Message:", conditionMessage(out), "\n")
} else {
  cat("UNEXPECTEDLY FIT (no error raised) -- converged:", out$converged, "\n")
}

saveRDS(list(mcar = list(lphi2 = lphi2, lalpha2 = lalpha2, kappa2 = kappa2, n_fail = n_fail),
             booklet_error = if (inherits(out, "error")) conditionMessage(out) else NULL),
        "tools/simval/round1/btl-efrm/06_missingness.rds")
cat("\ntotal time (s):", as.numeric(Sys.time() - t_start, units = "secs"), "\n")
