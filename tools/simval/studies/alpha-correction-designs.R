# STUDY: alpha-correction-designs
#
# Design-robustness of the corrected EFRM item-set unit estimator
# (companion to alpha-correction.csv, which carries the primary
# calibration): longer sets, pairwise-only person overlap over a 3-set
# linking graph, booklet-style within-set missingness (per-pattern
# correction moments), and a strongly skewed person distribution where
# the estimator's population-free construction is put against TAM's
# normal-population MML anchor.
# Serial. Rscript tools/simval/studies/alpha-correction-designs.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "alpha-correction-designs"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M")),
                          sprintf(...), "\n")

fit_alpha <- function(d, item_sets) {
  f <- tryCatch(rasch_efrm(d, item_sets = item_sets, groups = "group",
                           id = "id", boot_reps = 300),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  om <- f$efrm_vs_rasch$unit_omnibus
  list(alpha = setNames(f$alpha_table$alpha, f$alpha_table$set),
       se = setNames(f$alpha_table$se_log_alpha, f$alpha_table$set),
       p_om = om$p[grepl("alpha", om$term)][1])
}

## ---- A: 15 items per set --------------------------------------------------
R <- 100L; la <- se <- rep(NA_real_, R); n_ref <- 0L
for (r in seq_len(R)) {
  s <- simulate_efrm(n_per_group = 500, items_per_set = 15, n_sets = 2,
                     n_groups = 1, set_unit_ratio = 1.4, seed = 30e3 + r)
  z <- fit_alpha(s, attr(s, "truth")$item_sets)
  if (is.null(z)) { n_ref <- n_ref + 1L; next }
  la[r] <- log(z$alpha["set2"]); se[r] <- z$se["set2"]
}
ok <- is.finite(la); tr <- log(1.4)/2
add("15 items/set, ratio 1.4, N=500", "log alpha[set2] bias/coverage", sum(ok),
    bias = mean(la[ok]) - tr, emp_sd = sd(la[ok]), mean_se = mean(se[ok]),
    coverage95 = mean(abs(la[ok] - tr) <= 1.96 * se[ok]),
    n_attempted = R, n_refused = n_ref,
    notes = "raw construction at this length runs -0.036 in the log ratio (scratch harness); corrected expected ~0")
tick("A: bias %+.4f cover %.3f", mean(la[ok]) - tr,
     mean(abs(la[ok] - tr) <= 1.96 * se[ok]))

## ---- B: 3 sets, pairwise-only person overlap ------------------------------
pairwise_blank <- function(s, item_sets, seed) {
  set.seed(seed + 5e8)
  arms <- sample(rep(1:3, length.out = nrow(s)))
  drop_set <- c(3L, 1L, 2L)[arms]           # arm k lacks set drop_set[k]
  for (k in 1:3) s[drop_set == k, item_sets[[k]]] <- NA
  s
}
runB <- function(ratio, R, seed0, what) {
  lr <- pom <- rep(NA_real_, R); n_ref <- 0L
  for (r in seq_len(R)) {
    s <- simulate_efrm(n_per_group = 600, items_per_set = 8, n_sets = 3,
                       n_groups = 1, set_unit_ratio = ratio, seed = seed0 + r)
    its <- attr(s, "truth")$item_sets
    s <- pairwise_blank(s, its, seed0 + r)
    z <- fit_alpha(s, its)
    if (is.null(z)) { n_ref <- n_ref + 1L; next }
    lr[r] <- log(z$alpha["set3"] / z$alpha["set1"]); pom[r] <- z$p_om
  }
  ok <- is.finite(lr)
  add(sprintf("3 sets, pairwise person overlap, ratio %.1f", ratio),
      sprintf("log alpha[set3/set1] %s", what), sum(ok),
      bias = mean(lr[ok]) - log(ratio), emp_sd = sd(lr[ok]),
      type1 = if (ratio == 1) mean(pom[ok] < 0.05) else NA_real_,
      power = if (ratio != 1) mean(pom[ok] < 0.05) else NA_real_,
      n_attempted = R, n_refused = n_ref,
      notes = "each person takes 2 of 3 sets; every linking edge rests on a distinct third of the sample")
  tick("B ratio %.1f: bias %+.4f rej %.3f", ratio,
       mean(lr[ok]) - log(ratio), mean(pom[ok] < 0.05))
}
runB(1,   150L, 31e3, "null")
runB(1.4, 100L, 32e3, "bias/power")

## ---- C: booklet missingness within sets -----------------------------------
booklet_blank <- function(s, item_sets, seed) {
  set.seed(seed + 6e8)
  forms <- list(1:5, 4:8, c(1:3, 6:8))
  for (k in seq_along(item_sets)) {
    fm <- sample(rep(1:3, length.out = nrow(s)))
    for (b in 1:3) {
      drop_cols <- item_sets[[k]][-forms[[b]]]
      s[fm == b, drop_cols] <- NA
    }
  }
  s
}
runC <- function(ratio, R, seed0, what) {
  la <- se <- pom <- rep(NA_real_, R); n_ref <- 0L
  for (r in seq_len(R)) {
    s <- simulate_efrm(n_per_group = 600, items_per_set = 8, n_sets = 2,
                       n_groups = 1, set_unit_ratio = ratio, seed = seed0 + r)
    its <- attr(s, "truth")$item_sets
    s <- booklet_blank(s, its, seed0 + r)
    z <- fit_alpha(s, its)
    if (is.null(z)) { n_ref <- n_ref + 1L; next }
    la[r] <- log(z$alpha["set2"]); se[r] <- z$se["set2"]; pom[r] <- z$p_om
  }
  ok <- is.finite(la); tr <- log(ratio)/2
  add(sprintf("booklet forms (5-6 of 8 items), ratio %.1f", ratio),
      sprintf("log alpha[set2] %s", what), sum(ok),
      bias = mean(la[ok]) - tr, emp_sd = sd(la[ok]), mean_se = mean(se[ok]),
      coverage95 = mean(abs(la[ok] - tr) <= 1.96 * se[ok]),
      type1 = if (ratio == 1) mean(pom[ok] < 0.05) else NA_real_,
      power = if (ratio != 1) mean(pom[ok] < 0.05) else NA_real_,
      n_attempted = R, n_refused = n_ref,
      notes = "three rotated forms per set; per-pattern correction moments exercised")
  tick("C ratio %.1f: bias %+.4f cover %.3f rej %.3f", ratio,
       mean(la[ok]) - tr, mean(abs(la[ok] - tr) <= 1.96 * se[ok]),
       mean(pom[ok] < 0.05))
}
runC(1,   150L, 33e3, "null")
runC(1.4, 100L, 34e3, "bias/power")

## ---- R: ratio sweep -------------------------------------------------------
# bias, coverage, and rejection across the range of unit ratios; larger
# ratios make the two sets' information increasingly asymmetric
for (ratio in c(1, 1.15, 1.3, 1.5, 1.7, 2)) {
  R <- 80L; la <- se <- pom <- rep(NA_real_, R); n_ref <- 0L
  for (r in seq_len(R)) {
    s <- simulate_efrm(n_per_group = 500, items_per_set = 8, n_sets = 2,
                       n_groups = 1, set_unit_ratio = ratio,
                       seed = round(40e3 + 1e3 * ratio * 10) + r)
    z <- fit_alpha(s, attr(s, "truth")$item_sets)
    if (is.null(z)) { n_ref <- n_ref + 1L; next }
    la[r] <- log(z$alpha["set2"]); se[r] <- z$se["set2"]; pom[r] <- z$p_om
  }
  ok <- is.finite(la); tr <- log(ratio)/2
  add(sprintf("ratio sweep %.2f, 8 items/set, N=500", ratio),
      "log alpha[set2] bias/coverage", sum(ok),
      bias = mean(la[ok]) - tr, emp_sd = sd(la[ok]), mean_se = mean(se[ok]),
      coverage95 = mean(abs(la[ok] - tr) <= 1.96 * se[ok]),
      type1 = if (ratio == 1) mean(pom[ok] < 0.05) else NA_real_,
      power = if (ratio != 1) mean(pom[ok] < 0.05) else NA_real_,
      effect = ratio, n_attempted = R, n_refused = n_ref)
  tick("sweep %.2f: bias %+.4f cover %.3f rej %.3f", ratio,
       mean(la[ok]) - tr, mean(abs(la[ok] - tr) <= 1.96 * se[ok]),
       mean(pom[ok] < 0.05))
}

## ---- P: polytomous items (4 categories) -----------------------------------
runP <- function(ratio, R, seed0, what) {
  la <- se <- pom <- rep(NA_real_, R); n_ref <- 0L
  for (r in seq_len(R)) {
    d <- simulate_efrm(n_per_group = 500, items_per_set = 8, n_sets = 2,
                       n_groups = 1, set_unit_ratio = ratio,
                       n_categories = 4, seed = seed0 + r)
    z <- fit_alpha(d, attr(d, "truth")$item_sets)
    if (is.null(z)) { n_ref <- n_ref + 1L; next }
    la[r] <- log(z$alpha["set2"]); se[r] <- z$se["set2"]; pom[r] <- z$p_om
  }
  ok <- is.finite(la); tr <- log(ratio)/2
  add(sprintf("polytomous (4 categories), ratio %.1f, 8 items/set, N=500", ratio),
      sprintf("log alpha[set2] %s", what), sum(ok),
      bias = mean(la[ok]) - tr, emp_sd = sd(la[ok]), mean_se = mean(se[ok]),
      coverage95 = mean(abs(la[ok] - tr) <= 1.96 * se[ok]),
      type1 = if (ratio == 1) mean(pom[ok] < 0.05) else NA_real_,
      power = if (ratio != 1) mean(pom[ok] < 0.05) else NA_real_,
      n_attempted = R, n_refused = n_ref)
  tick("P ratio %.1f: bias %+.4f cover %.3f rej %.3f", ratio,
       mean(la[ok]) - tr, mean(abs(la[ok] - tr) <= 1.96 * se[ok]),
       mean(pom[ok] < 0.05))
}
runP(1,   150L, 36e3, "null")
runP(1.4, 100L, 37e3, "bias/power")

## ---- D: skew sweep, corrected vs TAM --------------------------------------
# theta from centred scaled chi-square(df): skewness sqrt(8/df)
delta <- seq(-1.5, 1.5, length.out = 8)
al_t <- c(1.4^-0.5, 1.4^0.5)
items <- paste0("S", rep(1:2, each = 8), "I", sprintf("%02d", rep(1:8, 2)))
isets <- list(set1 = items[1:8], set2 = items[9:16])
for (df in c(30, 8, 3, 1)) {
  R <- 40L; lr <- matrix(NA_real_, R, 2); n_ref <- 0L
  for (r in seq_len(R)) {
    set.seed(round(35e3 + df * 100) + r); n <- 1000
    th <- (rchisq(n, df) - df) / sqrt(2 * df) * 1.3
    X <- do.call(cbind, lapply(1:2, function(k)
      sapply(delta, function(dd) rbinom(n, 1, plogis(al_t[k] * (th - dd))))))
    colnames(X) <- items
    d <- data.frame(id = sprintf("P%04d", 1:n), X, group = "g1",
                    check.names = FALSE)
    z <- fit_alpha(d, isets)
    mt <- tryCatch(TAM::tam.mml.2pl(resp = X, irtmodel = "2PL.groups",
                                    est.slopegroups = rep(1:2, each = 8),
                                    control = list(progress = FALSE)),
                   error = function(e) NULL)
    if (is.null(z) || is.null(mt)) { n_ref <- n_ref + 1L; next }
    a_t <- mt$B[items, 2, 1]
    lr[r, ] <- c(log(z$alpha["set2"] / z$alpha["set1"]),
                 log(mean(a_t[9:16]) / mean(a_t[1:8])))
  }
  ok <- stats::complete.cases(lr)
  add(sprintf("skewed persons chi-sq(%d) (skew %.2f), ratio 1.4, N=1000",
              df, sqrt(8 / df)),
      "log unit ratio: corrected vs TAM anchor", sum(ok),
      bias = mean(lr[ok, 1]) - log(1.4), emp_sd = sd(lr[ok, 1]),
      effect = sqrt(8 / df), n_attempted = R, n_refused = n_ref,
      notes = sprintf("TAM normal-population anchor bias %+.4f (sd %.4f) on the same data; raw pre-correction construction ran +0.084 at skew 1.63",
        mean(lr[ok, 2]) - log(1.4), sd(lr[ok, 2])))
  tick("D skew %.2f: corrected %+.4f, TAM %+.4f", sqrt(8 / df),
       mean(lr[ok, 1]) - log(1.4), mean(lr[ok, 2]) - log(1.4))
}

sv_write(do.call(rbind, rows), "alpha-correction-designs")
cat(sprintf("TOTAL elapsed: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
