# STUDY: comparison-validation
#
# Validates the app's automatic model-comparison cards under their null
# model and at least two departure magnitudes (the release requirement).
# Card surfaces already validated elsewhere are recorded as citation rows
# pointing at their provenance rather than re-run:
#   lr_test size/power                -> results/lr-test.csv
#   lr_test small-sample edge        -> results/lr-smalln-topup.csv
#   EFRM unit omnibus                 -> results/efrm-fix-sweep.csv
#   BTL-EFRM unit omnibus (F ref)    -> results/round2-followups.csv
#   MFRM interaction omnibus          -> results/custom-wald-tests.csv
# The four NEW inferential surfaces validated here:
#   (C) CL-AIC/CL-BIC selection, partial credit vs rating scale
#   (B) CL-AIC/CL-BIC selection, free vs principal-component thresholds
#   (D) CL-AIC/CL-BIC selection, comparative judgement free vs PC thresholds
#   (E) paired-comparison effect tests: position, exposure, carry-over
#
# Single-CPU serial by design. Truths are FIXED per scenario; responses
# are redrawn each replicate. Run from the package root:
#   Rscript tools/simval/studies/comparison-validation.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "comparison-validation"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

## ---- citation rows -------------------------------------------------------
cite <- list(
  c("PCM vs RSM formal test (lr_test)", "size 4.7% @2,000; power 25-74% by effect", "results/lr-test.csv"),
  c("lr_test small-sample edge", "6.1% among 1,927/2,000 admissible (N=300, 12 four-category items)", "results/lr-smalln-topup.csv"),
  c("EFRM vs equal-frame unit omnibus", "size 4.9% @1,200 across eight designs; power 53-100%", "results/efrm-fix-sweep.csv"),
  c("BTL-EFRM vs equal-unit omnibus", "kappa 5.5% @400 post F-reference; power cells in custom-wald", "results/round2-followups.csv"),
  c("MFRM additive vs interactive omnibus", "size 4.9% @1,000 (q=10), 4.3-5.2% @600 (q=25); power 44-99%", "results/custom-wald-tests.csv"))
for (cc in cite)
  add(cc[1], "validated elsewhere (citation row)", n_reps = 0L,
      notes = sprintf("%s; see %s", cc[2], cc[3]))

## ---- shared samplers -----------------------------------------------------
simP <- function(t, tau) { x <- 0:length(tau); p <- exp(x * t - c(0, cumsum(tau))); p / sum(p) }
draw_poly <- function(theta, taus) {
  X <- sapply(taus, function(tt) vapply(theta, function(b)
    sample(0:length(tt), 1, prob = simP(b, tt)), 0L))
  colnames(X) <- sprintf("I%02d", seq_along(taus)); X
}
sel_rates <- function(scen, truth_lab, fit_a, fit_b, lab_a, lab_b, n_reps, gen) {
  pick_aic <- pick_bic <- rep(NA, n_reps); n_ref <- 0L
  for (r in seq_len(n_reps)) {
    X <- gen(r)
    fa <- tryCatch(fit_a(X), error = function(e) NULL)
    fb <- tryCatch(fit_b(X), error = function(e) NULL)
    if (is.null(fa) || is.null(fb)) { n_ref <- n_ref + 1L; next }
    cmp <- tryCatch(compare_fits(a = fa, b = fb), error = function(e) NULL)
    if (is.null(cmp) || any(!is.finite(cmp$cl_aic))) { n_ref <- n_ref + 1L; next }
    pick_aic[r] <- cmp$label[which.min(cmp$cl_aic)]
    pick_bic[r] <- cmp$label[which.min(cmp$cl_bic)]
  }
  ok <- !is.na(pick_aic)
  for (crit in c("aic", "bic")) {
    pk <- if (crit == "aic") pick_aic else pick_bic
    add(scen, sprintf("cl_%s selects '%s'", crit, lab_a), sum(ok),
        power = mean(pk[ok] == "a"),
        n_attempted = n_reps, n_refused = n_ref,
        notes = sprintf("truth: %s; alternatives '%s' vs '%s'", truth_lab, lab_a, lab_b))
  }
  cat(sprintf("[%s] %s: aic->%s %.3f (n=%d)\n", format(Sys.time(), "%H:%M"),
              scen, lab_a, mean(pick_aic[ok] == "a"), sum(ok)))
}

## ---- (C) PCM vs RSM selection -------------------------------------------
N <- 500L; I <- 8L
delta <- seq(-1.5, 1.5, length.out = I)
common <- c(-0.9, 0.9)
for (eff in c(0, 0.2, 0.45)) {
  set.seed(31e6 + eff * 100)
  jit <- if (eff == 0) matrix(0, I, 2) else
    cbind(runif(I, -eff, eff), runif(I, -eff, eff))
  taus <- lapply(seq_len(I), function(i) delta[i] + common + jit[i, ] - mean(jit[i, ]))
  truth_lab <- if (eff == 0) "RSM-true (common pattern)" else
    sprintf("PCM-true, pattern jitter +/-%.2f", eff)
  sel_rates(sprintf("PCM vs RSM selection, eff=%.2f", eff), truth_lab,
            function(X) rasch(X, model = "PCM"),
            function(X) rasch(X, model = "RSM"),
            "PCM", "RSM", 400L,
            function(r) { set.seed(32e6 + eff * 1e5 + r)
                          draw_poly(rnorm(N, 0, 1.2), taus) })
}

## ---- (B) free vs principal-component thresholds (item side) --------------
# 4-category items so the alternatives differ in dimension: free = 3
# thresholds per item, pc_components = 2 = location + spread (equally
# spaced). PC-true is a pure linear trend; the departure loads the
# skewness (quadratic) contrast, orthogonal to location and spread.
for (skew in c(0, 0.3, 0.6)) {
  set.seed(33e6 + skew * 100)
  gap <- runif(I, 0.8, 1.2)
  taus <- lapply(seq_len(I), function(i)
    delta[i] + gap[i] * c(-1, 0, 1) + skew * c(1, -2, 1) / 2)
  truth_lab <- if (skew == 0) "PC-true (location + spread only)" else
    sprintf("free-true, skewness contrast %.2f", skew)
  sel_rates(sprintf("free vs PC thresholds, skew=%.2f", skew), truth_lab,
            function(X) rasch(X, model = "PCM"),
            function(X) rasch(X, model = "PCM", pc_components = 2),
            "free", "pc2", 400L,
            function(r) { set.seed(34e6 + skew * 1e5 + r)
                          draw_poly(rnorm(N, 0, 1.2), taus) })
}

## ---- (D) comparative judgement free vs PC thresholds ---------------------
K <- 8L; objs <- sprintf("O%d", seq_len(K))
beta <- setNames(seq(-1, 1, length.out = K), objs)
pr <- t(utils::combn(objs, 2))
cj_gen <- function(taus_cj, seed) {
  set.seed(seed)
  d <- data.frame(object_a = rep(pr[, 1], each = 12L),
                  object_b = rep(pr[, 2], each = 12L), stringsAsFactors = FALSE)
  lp <- beta[d$object_a] - beta[d$object_b]
  m <- length(taus_cj)
  d$response <- vapply(lp, function(l) sample(0:m, 1, prob = simP(l, taus_cj)), 0L)
  d$judge <- sample(sprintf("J%d", 1:14), nrow(d), replace = TRUE)
  d
}
# btl's "pc" structure is the normalised linear contrast: exactly equally
# spaced symmetric thresholds (1 parameter; "free" has floor(m/2) = 2).
# The null must therefore be equally spaced; the departure perturbs the
# symmetric pairs away from equal spacing.
for (dev in c(0, 0.35, 0.7)) {
  taus_cj <- c(-1.5 - dev, -0.5 + dev, 0.5 - dev, 1.5 + dev)
  truth_lab <- if (dev == 0) "PC-true (equally spaced symmetric)" else
    sprintf("free-true symmetric, pair deviation %.2f", dev)
  sel_rates(sprintf("CJ free vs PC thresholds, dev=%.2f", dev), truth_lab,
            function(d) btl(d, "object_a", "object_b", response = "response",
                            judge = "judge", thresholds = "free"),
            function(d) btl(d, "object_a", "object_b", response = "response",
                            judge = "judge", thresholds = "pc"),
            "free", "pc", 400L,
            function(r) cj_gen(taus_cj, 35e6 + dev * 1e5 + r))
}

## ---- (E) paired-comparison effect tests ----------------------------------
# simulate_btl only emits a per-judge order column when dependence is
# planted; the null and position designs need one too (null-consistent:
# no order effect is planted by construction)
add_order <- function(d) {
  d <- d[order(d$judge), ]
  d$order <- stats::ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  d
}
eff_null <- function(n_reps) {
  p_pos <- p_exp <- p_cry <- rep(NA_real_, n_reps); n_ref <- 0L
  for (r in seq_len(n_reps)) {
    d <- add_order(simulate_btl(n_objects = 8, n_judges = 14,
                                reps_per_pair = 10, seed = 36e6 + r))
    f <- tryCatch(btl(d, "object_a", "object_b", winner = "winner",
                      judge = "judge", order = "order", position = TRUE),
                  error = function(e) NULL)
    if (is.null(f) || is.null(f$dependence)) { n_ref <- n_ref + 1L; next }
    dp <- f$dependence
    p_pos[r] <- dp$p[dp$effect == "position"][1]
    p_exp[r] <- dp$p[dp$effect == "exposure"][1]
    p_cry[r] <- dp$p[dp$effect == "carry_over"][1]
  }
  for (nm in c("pos", "exp", "cry")) {
    pv <- get(paste0("p_", nm)); ok <- is.finite(pv)
    eff_lab <- c(pos = "position", exp = "exposure", cry = "carry_over")[nm]
    add(sprintf("effect test null: %s", eff_lab), "type1 at alpha=.05",
        sum(ok), type1 = mean(pv[ok] < 0.05),
        n_attempted = n_reps, n_refused = n_ref,
        notes = "no effects planted; one fit supplies all three tests")
    cat(sprintf("[%s] null %s: %.4f (n=%d)\n", format(Sys.time(), "%H:%M"),
                eff_lab, mean(pv[ok] < 0.05), sum(ok)))
  }
}
eff_null(800L)

eff_power <- function(effect, mag, n_reps) {
  pv <- rep(NA_real_, n_reps); n_ref <- 0L
  for (r in seq_len(n_reps)) {
    seed <- 37e6 + match(effect, c("position", "exposure", "carry_over")) * 1e6 +
      mag * 1e5 + r
    if (effect == "position") {
      d <- add_order(simulate_btl(n_objects = 8, n_judges = 14,
                                  reps_per_pair = 10, seed = seed))
      tr <- attr(d, "truth")$location
      set.seed(seed + 5e8)
      lp <- tr[d$object_a] - tr[d$object_b] + mag   # first-listed advantage
      win_a <- rbinom(nrow(d), 1, plogis(lp)) == 1
      d$winner <- ifelse(win_a, d$object_a, d$object_b)
    } else {
      dep <- setNames(list(mag), effect)
      d <- simulate_btl(n_objects = 8, n_judges = 14, reps_per_pair = 10,
                        dependence = dep, seed = seed)
    }
    f <- tryCatch(btl(d, "object_a", "object_b", winner = "winner",
                      judge = "judge", order = "order", position = TRUE),
                  error = function(e) NULL)
    if (is.null(f) || is.null(f$dependence)) { n_ref <- n_ref + 1L; next }
    pv[r] <- f$dependence$p[f$dependence$effect == effect][1]
  }
  ok <- is.finite(pv)
  add(sprintf("effect test power: %s, magnitude %.1f", effect, mag),
      "power at alpha=.05", sum(ok), power = mean(pv[ok] < 0.05), effect = mag,
      n_attempted = n_reps, n_refused = n_ref)
  cat(sprintf("[%s] power %s %.1f: %.3f (n=%d)\n", format(Sys.time(), "%H:%M"),
              effect, mag, mean(pv[ok] < 0.05), sum(ok)))
}
for (effect in c("position", "exposure", "carry_over"))
  for (mag in c(0.3, 0.6)) eff_power(effect, mag, 300L)

# carry-over null top-up at more clusters: round 1 saw 7.5% (400 reps)
# and this study ~8% at 14 judges; check the elevation is the documented
# few-cluster effect by rerunning the null with 30 judges and a per-judge
# history length comparable to round 1's design
cry_null30 <- function(n_reps) {
  pv <- rep(NA_real_, n_reps); n_ref <- 0L
  for (r in seq_len(n_reps)) {
    d <- add_order(simulate_btl(n_objects = 8, n_judges = 30,
                                reps_per_pair = 25, seed = 38e6 + r))
    f <- tryCatch(btl(d, "object_a", "object_b", winner = "winner",
                      judge = "judge", order = "order"),
                  error = function(e) NULL)
    if (is.null(f) || is.null(f$dependence)) { n_ref <- n_ref + 1L; next }
    pv[r] <- f$dependence$p[f$dependence$effect == "carry_over"][1]
  }
  ok <- is.finite(pv)
  add("effect test null: carry_over, 30 judges", "type1 at alpha=.05",
      sum(ok), type1 = mean(pv[ok] < 0.05),
      n_attempted = n_reps, n_refused = n_ref,
      notes = "few-cluster check: 14-judge nulls run ~8%; 30 judges, 700 comparisons")
  cat(sprintf("[%s] null carry_over J=30: %.4f (n=%d)\n",
              format(Sys.time(), "%H:%M"), mean(pv[ok] < 0.05), sum(ok)))
}
cry_null30(400L)

sv_write(do.call(rbind, rows), "comparison-validation")
cat(sprintf("TOTAL elapsed: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
