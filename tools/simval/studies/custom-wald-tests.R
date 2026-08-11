# Simulation validation: null size and power of the package's newer
# custom test procedures.
#   (a) dif_contrasts() on repeated/multifactor designs (R/dif.R)
#   (b) rasch_mfrm(..., interaction=) omnibus interaction_test F (R/mfrm.R)
#   (c) rasch_efrm() unit_omnibus Wald chi-square tests (R/efrm.R)
#   (d) btl_efrm() unit_omnibus Wald chi-square tests, judge bootstrap (R/btl-efrm.R)
#
# Run from the package root:
#   Rscript tools/simval/studies/custom-wald-tests.R          # everything, serially
#   Rscript tools/simval/studies/custom-wald-tests.R 1        # chunk 1: family (a) only
#   Rscript tools/simval/studies/custom-wald-tests.R 2        # chunk 2: families (b)+(c)+(d)
#   Rscript tools/simval/studies/custom-wald-tests.R combine  # stitch the two chunk CSVs
#
# Chunks 1 and 2 were run as two concurrent Rscript processes (the package's
# "at most 2 concurrent" budget) and then stitched with `combine`; running
# with no argument reproduces the identical combined CSV in one process
# (slower, ~85 CPU-minutes, but deterministic and self-contained).
#
# METHODOLOGY: true generating parameters are fixed within a scenario; only
# responses (and persons, where the scenario says so) are redrawn across
# replicates. A replicate that errors during simulation/fit is a REFUSAL; one
# that fits but is flagged not-converged is a NONCONV; both are excluded from
# n_reps (the analysed replicates each rate is computed on) and reported
# separately (n_attempted, n_refused, n_nonconv -> refusal_rate, nonconv_rate
# in the harness).
#
# BENCHMARKED replicate costs (this machine, single core), used to size reps
# to the ~70-minute total wall-clock budget at <= 2 concurrent processes:
#   (a) dif design 1 (within only)   N=200: 0.55s  N=500: 0.62s
#       dif design 2 (2x2 factorial) N=200: 0.61-0.72s  N=500: 0.84-0.95s
#   (b) MFRM interaction_test        R=3: 0.31-0.37s   R=6: 1.08-1.36s
#   (c) EFRM unit_omnibus (hybrid SE)                  0.68s
#   (d) BTL-EFRM unit_omnibus (judge_bootstrap, boot_reps=150)  2.2s
# Chunk 1 (family a) ~36 CPU-minutes; chunk 2 (b+c+d) ~47 CPU-minutes; run
# concurrently the observed wall clock was well inside the 70-minute budget.
#
# The PRINCIPAL null-calibration claim of the whole study is family (a),
# design 1 (within-only occasion DIF), N = 200, at 1,000 replicates (cost
# 0.55s/rep, comfortably under the 1.5s/rep threshold for the >=1000 rule).

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

SCALE <- as.numeric(Sys.getenv("SV_SCALE", "1"))   # replicate-count scale factor, for smoke tests
nrep <- function(n) max(1L, round(n * SCALE))

args <- commandArgs(trailingOnly = TRUE)
chunk <- if (length(args) >= 1) args[1] else "all"

STUDY <- "custom-wald-tests"
ALPHA <- 0.05
t_study <- Sys.time()

# ---------------------------------------------------------------------------
# Family (a): dif_contrasts() -- repeated-measures and 2x2 between x within
# ---------------------------------------------------------------------------

# Generator follows the construction of the dif-repeated-measures vignette
# and the package's own dif_anova roxygen example: dichotomous Rasch, 8
# items, an optional between-person "group" factor (imbalanced split) and a
# within-person "occasion" factor (2 waves, stacked). A uniform DIF shift can
# be planted on one item for the non-reference group level (between) and/or
# for the second occasion (within).
gen_dif <- function(N, design = c("within", "factorial"), n_items = 8,
                    occ_dif_item = 6, occ_dif_shift = 0,
                    group_dif_item = 3, group_dif_shift = 0,
                    group_split = 0.5, theta_sd = 1, seed = NULL) {
  design <- match.arg(design)
  if (!is.null(seed)) set.seed(seed)
  d <- seq(-2, 2, length.out = n_items)
  theta <- rnorm(N, 0, theta_sd)
  if (design == "factorial") {
    n1 <- round(N * group_split)
    group <- factor(c(rep("g1", n1), rep("g2", N - n1)))
  } else group <- factor(rep("g1", N))
  occ_lv <- c("T1", "T2")
  make_wave <- function(occ) {
    shift <- matrix(0, N, n_items)
    if (design == "factorial")
      shift[group == "g2", group_dif_item] <- group_dif_shift
    if (occ == "T2") shift[, occ_dif_item] <- occ_dif_shift
    matrix(rbinom(N * n_items, 1,
                  plogis(outer(theta, d, "-") - shift)), N, n_items)
  }
  X <- rbind(make_wave("T1"), make_wave("T2"))
  colnames(X) <- sprintf("I%02d", seq_len(n_items))
  data.frame(X, person = rep(seq_len(N), 2), group = rep(as.character(group), 2),
             occasion = rep(occ_lv, each = N), stringsAsFactors = FALSE)
}

# One replicate: fit + dif_contrasts, returns the result table or NULL.
dif_rep <- function(design, N, group_split = 0.5,
                    occ_dif_shift = 0, group_dif_shift = 0, seed) {
  d <- gen_dif(N, design = design, occ_dif_item = 6, occ_dif_shift = occ_dif_shift,
              group_dif_item = 3, group_dif_shift = group_dif_shift,
              group_split = group_split, seed = seed)
  factors <- if (design == "within") "occasion" else c("group", "occasion")
  fit <- tryCatch(rasch(d, id = "person", factors = factors), error = function(e) NULL)
  if (is.null(fit)) return(list(table = NULL, refused = TRUE, nonconv = FALSE))
  if (!isTRUE(fit$est$converged))
    return(list(table = NULL, refused = FALSE, nonconv = TRUE))
  dc <- if (design == "within")
    tryCatch(dif_contrasts(fit, within = "occasion", id = fit$person$id),
             error = function(e) NULL)
  else
    tryCatch(dif_contrasts(fit, id = fit$person$id), error = function(e) NULL)
  if (is.null(dc)) return(list(table = NULL, refused = TRUE, nonconv = FALSE))
  list(table = dc$table, refused = FALSE, nonconv = FALSE)
}

dif_batch <- function(design, N, n_reps, seed0, group_split = 0.5,
                      occ_dif_shift = 0, group_dif_shift = 0) {
  tabs <- vector("list", n_reps)
  n_refused <- 0L; n_nonconv <- 0L
  for (r in seq_len(n_reps)) {
    res <- dif_rep(design, N, group_split, occ_dif_shift, group_dif_shift, seed0 + r)
    if (res$refused) { n_refused <- n_refused + 1L; next }
    if (res$nonconv) { n_nonconv <- n_nonconv + 1L; next }
    tabs[[r]] <- res$table
  }
  list(tabs = Filter(Negate(is.null), tabs), n_attempted = n_reps,
       n_refused = n_refused, n_nonconv = n_nonconv)
}

# Null-size rows: pool every item x contrast cell (per-cell rate = "type1";
# any cell significant in a replicate = "familywise", the quantity the
# family's Holm adjustment nominally controls).
dif_null_rows <- function(scenario, batch, notes = "") {
  tabs <- batch$tabs
  n_reps <- length(tabs)
  # per-replicate proportion of significant cells: the pooled rate is its
  # mean, and the MC SE is sd(per-rep)/sqrt(n_reps) -- cells within one
  # replicate share a fit, so the plug-in binomial SE over all cells
  # would understate the Monte Carlo error
  per_rep <- vapply(tabs, function(t) mean(t$significant), 0)
  any_sig <- vapply(tabs, function(t) any(t$significant), logical(1))
  list(
    sv_row(STUDY, scenario, "type1_cell (per item x contrast, Holm-adjusted)",
           n_reps = n_reps, type1 = mean(per_rep),
           mc_override = list(type1 = stats::sd(per_rep) / sqrt(n_reps)),
           effect = 0,
           n_attempted = batch$n_attempted, n_refused = batch$n_refused,
           n_nonconv = batch$n_nonconv, notes = notes),
    sv_row(STUDY, scenario, "familywise (any cell significant per replicate)",
           n_reps = n_reps, familywise = mean(any_sig), effect = 0,
           n_attempted = batch$n_attempted, n_refused = batch$n_refused,
           n_nonconv = batch$n_nonconv, notes = notes)
  )
}

# Power row: the flag rate on the specific planted item x contrast cell.
dif_power_row <- function(scenario, batch, item, contrast, effect, notes = "") {
  tabs <- batch$tabs
  vals <- vapply(tabs, function(t) {
    row <- t[t$item == item & t$contrast == contrast, ]
    if (!nrow(row)) NA else row$significant[1]
  }, logical(1))
  ok <- !is.na(vals)
  sv_row(STUDY, scenario, sprintf("power (item %s, contrast '%s')", item, contrast),
         n_reps = sum(ok), power = mean(vals[ok]), effect = effect,
         n_attempted = batch$n_attempted, n_refused = batch$n_refused,
         n_nonconv = batch$n_nonconv, notes = notes)
}

run_family_a <- function() {
  rows <- list()
  t0 <- Sys.time()

  ## Design 1: within-only (single occasion factor, 2 waves, stacked)
  cat("[a] design 1 (within-only), N=200, null, ", nrep(1000), "reps (PRINCIPAL null claim)\n")
  b <- dif_batch("within", 200, nrep(1000), seed0 = 1e5,
                occ_dif_shift = 0)
  rows <- c(rows, dif_null_rows("dif design1(within) N=200 null", b,
    "PRINCIPAL null-calibration claim of the whole study"))

  cat("[a] design 1, N=500, null, ", nrep(300), "reps\n")
  b <- dif_batch("within", 500, nrep(300), seed0 = 2e5, occ_dif_shift = 0)
  rows <- c(rows, dif_null_rows("dif design1(within) N=500 null", b))

  for (N in c(200, 500)) for (eff in c(0.4, 0.8)) {
    seed0 <- 3e5 + N * 10 + eff * 1000
    cat(sprintf("[a] design 1, N=%d, power eff=%.1f, %d reps\n", N, eff, nrep(150)))
    b <- dif_batch("within", N, nrep(150), seed0 = seed0, occ_dif_shift = eff)
    rows[[length(rows) + 1]] <- dif_power_row(
      sprintf("dif design1(within) N=%d power occasion-DIF %.1f", N, eff),
      b, "I06", "occasion: T2 - T1", eff)
  }
  cat("design 1 elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")

  ## Design 2: 2x2 factorial (between "group" x within "occasion")
  t1 <- Sys.time()
  for (N in c(200, 500)) for (split in c(0.5, 0.15)) {
    seed0 <- 6e5 + N * 100 + split * 1000
    cat(sprintf("[a] design 2, N=%d, split=%.2f, null, %d reps\n", N, split, nrep(150)))
    b <- dif_batch("factorial", N, nrep(150), seed0 = seed0, group_split = split,
                  occ_dif_shift = 0, group_dif_shift = 0)
    rows <- c(rows, dif_null_rows(
      sprintf("dif design2(2x2) N=%d split=%.2f null", N, split), b))
  }
  for (N in c(200, 500)) for (split in c(0.5, 0.15)) for (eff in c(0.4, 0.8)) {
    seed0 <- 9e5 + N * 1000 + split * 10000 + eff * 100000
    cat(sprintf("[a] design 2, N=%d, split=%.2f, power eff=%.1f, %d reps\n",
                N, split, eff, nrep(100)))
    b <- dif_batch("factorial", N, nrep(100), seed0 = seed0, group_split = split,
                  occ_dif_shift = 0, group_dif_shift = eff)
    rows[[length(rows) + 1]] <- dif_power_row(
      sprintf("dif design2(2x2) N=%d split=%.2f power group-DIF %.1f", N, split, eff),
      b, "I03", "group: g2 - g1", eff)
  }
  cat("design 2 elapsed:", as.numeric(Sys.time() - t1, units = "secs"), "s\n")
  cat("family (a) total elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Family (b): MFRM interaction_test omnibus F
# ---------------------------------------------------------------------------

# FIXED TRUTH: simulate_mfrm() redraws the rater severities from every
# seed, which would vary the truth across replicates and fold that
# variation into the replicate distribution. Here the truth -- severities
# (the only random component of simulate_mfrm's truth; its item locations
# and thresholds are deterministic grids), items, thresholds, and the
# planted interaction -- is drawn ONCE per rater count R (7.7e6 + R,
# shared across every N, null, and effect cell), and only
# the persons and their responses are redrawn each replicate, using the
# package's own categorical sampler (.sim_item, exposed by load_all).
mfrm_truth <- function(R, seed_truth) {
  set.seed(seed_truth)
  I <- 6L
  list(lambda = setNames(as.numeric(scale(stats::rnorm(R))) * 0.7,
                         sprintf("R%d", seq_len(R))),
       delta = setNames(seq(-1, 1, length.out = I), sprintf("I%d", seq_len(I))),
       base_tau = c(-1.2, 0, 1.2))          # .sim_thresholds(0, 3, 1.2)
}

mfrm_rep <- function(N, truth, bias = NA, item = "I3", rater = "R3", seed) {
  lambda <- truth$lambda; delta <- truth$delta; base_tau <- truth$base_tau
  I <- length(delta); R <- length(lambda)
  set.seed(seed)
  theta <- stats::rnorm(N, 0, 1.2)
  grid <- expand.grid(p = seq_len(N), i = seq_len(I), r = seq_len(R))
  score <- integer(nrow(grid))
  for (i in seq_len(I)) for (r in seq_len(R)) {
    sel <- grid$i == i & grid$r == r
    tau_ir <- base_tau + delta[i] + lambda[r] +
      if (!is.na(bias) && names(delta)[i] == item && names(lambda)[r] == rater)
        bias else 0
    score[sel] <- .sim_item(theta[grid$p[sel]], tau_ir)
  }
  d <- data.frame(person = sprintf("P%03d", grid$p), item = names(delta)[grid$i],
                  rater = names(lambda)[grid$r], score = score,
                  stringsAsFactors = FALSE)
  fit <- tryCatch(rasch_mfrm(d, person = "person", item = "item", score = "score",
                             facets = "rater", interaction = "rater"),
                   error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA, df = NA, refused = TRUE, nonconv = FALSE))
  if (!isTRUE(fit$est$converged))
    return(list(p = NA, df = NA, refused = FALSE, nonconv = TRUE))
  list(p = fit$interaction_test$p, df = fit$interaction_test$df,
       refused = FALSE, nonconv = FALSE)
}

mfrm_batch <- function(N, R, bias, n_reps, seed0) {
  # ONE truth per rater count R, shared across every N, null, and effect
  # cell (a per-scenario truth would let the planted rater's severity
  # differ between effect sizes, confounding the power comparison);
  # response seeds still differ per cell via seed0
  truth <- mfrm_truth(R, seed_truth = 7.7e6 + R)
  p <- rep(NA_real_, n_reps); n_refused <- 0L; n_nonconv <- 0L; dfq <- NA
  for (r in seq_len(n_reps)) {
    res <- mfrm_rep(N, truth, bias, seed = seed0 + r)
    if (res$refused) { n_refused <- n_refused + 1L; next }
    if (res$nonconv) { n_nonconv <- n_nonconv + 1L; next }
    p[r] <- res$p; dfq <- res$df
  }
  list(p = p[!is.na(p)], n_attempted = n_reps, n_refused = n_refused,
       n_nonconv = n_nonconv, df = dfq)
}

run_family_b <- function() {
  rows <- list()
  t0 <- Sys.time()

  cat("[b] null n=50 R=3,", nrep(1000), "reps (>=800 confirmation cell)\n")
  b <- mfrm_batch(50, 3, NA, nrep(1000), seed0 = 1.1e6)
  rows[[length(rows) + 1]] <- sv_row(STUDY, sprintf("mfrm null n=50 R=3 (q=%s)", b$df),
    "type1 (interaction_test omnibus F)", n_reps = length(b$p), type1 = mean(b$p < ALPHA),
    effect = 0, n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv,
    notes = ">=800-rep confirmation of NEWS-documented ~5.7% at n=50 (there: n_items/raters unspecified in NEWS; here I=6, R=3, q=10)")

  cat("[b] null n=200 R=3,", nrep(300), "reps\n")
  b <- mfrm_batch(200, 3, NA, nrep(300), seed0 = 1.2e6)
  rows[[length(rows) + 1]] <- sv_row(STUDY, sprintf("mfrm null n=200 R=3 (q=%s)", b$df),
    "type1 (interaction_test omnibus F)", n_reps = length(b$p), type1 = mean(b$p < ALPHA),
    effect = 0, n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv)

  cat("[b] null n=50 R=6,", nrep(150), "reps\n")
  b <- mfrm_batch(50, 6, NA, nrep(150), seed0 = 1.3e6)
  rows[[length(rows) + 1]] <- sv_row(STUDY, sprintf("mfrm null n=50 R=6 (q=%s)", b$df),
    "type1 (interaction_test omnibus F)", n_reps = length(b$p), type1 = mean(b$p < ALPHA),
    effect = 0, n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv)

  cat("[b] null n=200 R=6,", nrep(150), "reps\n")
  b <- mfrm_batch(200, 6, NA, nrep(150), seed0 = 1.4e6)
  rows[[length(rows) + 1]] <- sv_row(STUDY, sprintf("mfrm null n=200 R=6 (q=%s)", b$df),
    "type1 (interaction_test omnibus F)", n_reps = length(b$p), type1 = mean(b$p < ALPHA),
    effect = 0, n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv)

  for (R in c(3, 6)) for (bias in c(0.4, 0.8)) {
    n_reps <- if (R == 3) nrep(150) else nrep(120)
    seed0 <- 1.5e6 + R * 1e4 + bias * 1e5
    cat(sprintf("[b] power N=200 R=%d bias=%.1f, %d reps\n", R, bias, n_reps))
    b <- mfrm_batch(200, R, bias, n_reps, seed0 = seed0)
    rows[[length(rows) + 1]] <- sv_row(STUDY,
      sprintf("mfrm power N=200 R=%d bias=%.1f (q=%s)", R, bias, b$df),
      "power (interaction_test omnibus F)", n_reps = length(b$p), power = mean(b$p < ALPHA),
      effect = bias, n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv)
  }
  cat("family (b) elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Family (c): EFRM unit_omnibus Wald tests
# ---------------------------------------------------------------------------

efrm_rep <- function(set_ratio, group_ratio, seed) {
  d <- simulate_efrm(200, 8, n_sets = 2, n_groups = 2, set_unit_ratio = set_ratio,
                     group_unit_ratio = group_ratio, seed = seed)
  fit <- tryCatch(rasch_efrm(d, item_sets = attr(d, "truth")$item_sets, groups = "group",
                             se_method = "hybrid"), error = function(e) NULL)
  if (is.null(fit)) return(list(p_alpha = NA, p_phi = NA, refused = TRUE, nonconv = FALSE))
  if (!isTRUE(fit$est$converged))
    return(list(p_alpha = NA, p_phi = NA, refused = FALSE, nonconv = TRUE))
  uo <- fit$efrm_vs_rasch$unit_omnibus
  pa <- uo$p[uo$term == "set units (alpha)"]; pg <- uo$p[uo$term == "group units (phi)"]
  list(p_alpha = if (length(pa)) pa else NA, p_phi = if (length(pg)) pg else NA,
       refused = FALSE, nonconv = FALSE)
}

efrm_batch <- function(set_ratio, group_ratio, n_reps, seed0) {
  pa <- pg <- rep(NA_real_, n_reps); n_refused <- 0L; n_nonconv <- 0L
  for (r in seq_len(n_reps)) {
    res <- efrm_rep(set_ratio, group_ratio, seed0 + r)
    if (res$refused) { n_refused <- n_refused + 1L; next }
    if (res$nonconv) { n_nonconv <- n_nonconv + 1L; next }
    pa[r] <- res$p_alpha; pg[r] <- res$p_phi
  }
  list(p_alpha = pa, p_phi = pg, n_attempted = n_reps, n_refused = n_refused,
       n_nonconv = n_nonconv)
}

run_family_c <- function() {
  rows <- list()
  t0 <- Sys.time()

  cat("[c] null set=1 group=1,", nrep(800), "reps\n")
  b <- efrm_batch(1, 1, nrep(800), seed0 = 2.1e6)
  ok_a <- !is.na(b$p_alpha); ok_g <- !is.na(b$p_phi)
  rows[[length(rows) + 1]] <- sv_row(STUDY, "efrm null (set=1,group=1) N=200/grp",
    "type1 (set units alpha omnibus)", n_reps = sum(ok_a), type1 = mean(b$p_alpha[ok_a] < ALPHA),
    effect = 0, n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv,
    notes = "same replicates give the group-units null below (one fit per rep)")
  rows[[length(rows) + 1]] <- sv_row(STUDY, "efrm null (set=1,group=1) N=200/grp",
    "type1 (group units phi omnibus)", n_reps = sum(ok_g), type1 = mean(b$p_phi[ok_g] < ALPHA),
    effect = 0, n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv)

  for (set_ratio in c(1.2, 1.5)) {
    cat(sprintf("[c] power set_ratio=%.1f (group=1), %d reps\n", set_ratio, nrep(150)))
    b <- efrm_batch(set_ratio, 1, nrep(150), seed0 = 2.2e6 + set_ratio * 1e5)
    ok_a <- !is.na(b$p_alpha); ok_g <- !is.na(b$p_phi)
    rows[[length(rows) + 1]] <- sv_row(STUDY,
      sprintf("efrm power set_ratio=%.1f group=1 N=200/grp", set_ratio),
      "power (set units alpha omnibus, target)", n_reps = sum(ok_a),
      power = mean(b$p_alpha[ok_a] < ALPHA), effect = set_ratio,
      n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv)
    rows[[length(rows) + 1]] <- sv_row(STUDY,
      sprintf("efrm power set_ratio=%.1f group=1 N=200/grp", set_ratio),
      "off-target flag rate (group units phi omnibus)", n_reps = sum(ok_g),
      power = mean(b$p_phi[ok_g] < ALPHA), effect = set_ratio,
      n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv,
      notes = "off-target: no group-unit effect planted here")
  }
  for (group_ratio in c(1.2, 1.5)) {
    cat(sprintf("[c] power group_ratio=%.1f (set=1), %d reps\n", group_ratio, nrep(150)))
    b <- efrm_batch(1, group_ratio, nrep(150), seed0 = 2.3e6 + group_ratio * 1e5)
    ok_a <- !is.na(b$p_alpha); ok_g <- !is.na(b$p_phi)
    rows[[length(rows) + 1]] <- sv_row(STUDY,
      sprintf("efrm power group_ratio=%.1f set=1 N=200/grp", group_ratio),
      "power (group units phi omnibus, target)", n_reps = sum(ok_g),
      power = mean(b$p_phi[ok_g] < ALPHA), effect = group_ratio,
      n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv)
    rows[[length(rows) + 1]] <- sv_row(STUDY,
      sprintf("efrm power group_ratio=%.1f set=1 N=200/grp", group_ratio),
      "off-target flag rate (set units alpha omnibus)", n_reps = sum(ok_a),
      power = mean(b$p_alpha[ok_a] < ALPHA), effect = group_ratio,
      n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv,
      notes = "off-target: no set-unit effect planted here")
  }
  cat("family (c) elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Family (d): BTL-EFRM unit_omnibus, judge bootstrap
# ---------------------------------------------------------------------------

BOOT_REPS <- 150L   # >= 100 as required; benchmarked ~2.2s/replicate at this size

btl_rep <- function(set_ratio, seed) {
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 6, reps_within = 20, reps_cross = 20,
                         panel_units = c(1, 1), set_units = c(1, set_ratio),
                         set_origins = c(0, 0), seed = seed)
  tr <- attr(d, "truth")
  fit <- tryCatch(btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
                           panels = "panel", object_sets = tr$object_sets,
                           se_method = "judge_bootstrap", boot_reps = BOOT_REPS),
                   error = function(e) NULL)
  if (is.null(fit)) return(list(p_phi = NA, p_alpha = NA, p_kappa = NA,
                                refused = TRUE, nonconv = FALSE))
  if (!isTRUE(fit$converged))
    return(list(p_phi = NA, p_alpha = NA, p_kappa = NA, refused = FALSE, nonconv = TRUE))
  uo <- fit$unit_omnibus
  gp <- function(term) { v <- uo$p[uo$term == term]; if (length(v)) v else NA }
  list(p_phi = gp("panel units (phi)"), p_alpha = gp("set units (alpha)"),
       p_kappa = gp("set origins (kappa)"), refused = FALSE, nonconv = FALSE)
}

btl_batch <- function(set_ratio, n_reps, seed0) {
  p_phi <- p_alpha <- p_kappa <- rep(NA_real_, n_reps)
  n_refused <- 0L; n_nonconv <- 0L
  for (r in seq_len(n_reps)) {
    res <- btl_rep(set_ratio, seed0 + r)
    if (res$refused) { n_refused <- n_refused + 1L; next }
    if (res$nonconv) { n_nonconv <- n_nonconv + 1L; next }
    p_phi[r] <- res$p_phi; p_alpha[r] <- res$p_alpha; p_kappa[r] <- res$p_kappa
  }
  list(p_phi = p_phi, p_alpha = p_alpha, p_kappa = p_kappa, n_attempted = n_reps,
       n_refused = n_refused, n_nonconv = n_nonconv)
}

run_family_d <- function() {
  rows <- list()
  t0 <- Sys.time()

  cat("[d] null set=1 panel=1,", nrep(150), sprintf("reps (judge_bootstrap, boot_reps=%d)\n", BOOT_REPS))
  b <- btl_batch(1, nrep(150), seed0 = 3.1e6)
  for (nm in c("p_phi", "p_alpha", "p_kappa")) {
    v <- b[[nm]]; ok <- !is.na(v)
    label <- c(p_phi = "panel units (phi)", p_alpha = "set units (alpha)",
              p_kappa = "set origins (kappa)")[[nm]]
    rows[[length(rows) + 1]] <- sv_row(STUDY, "btl_efrm null (set=1,panel=1) judge_bootstrap",
      sprintf("type1 (%s omnibus)", label), n_reps = sum(ok), type1 = mean(v[ok] < ALPHA),
      effect = 0, n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv,
      notes = if (nm == "p_phi") sprintf("boot_reps=%d; one fit gives all three nulls", BOOT_REPS) else "")
  }

  for (set_ratio in c(1.3, 1.6)) {
    n_reps <- nrep(80)
    cat(sprintf("[d] power set_ratio=%.1f,", set_ratio), n_reps, "reps\n")
    b <- btl_batch(set_ratio, n_reps, seed0 = 3.2e6 + set_ratio * 1e5)
    ok_a <- !is.na(b$p_alpha); ok_p <- !is.na(b$p_phi)
    rows[[length(rows) + 1]] <- sv_row(STUDY,
      sprintf("btl_efrm power set_ratio=%.1f panel=1", set_ratio),
      "power (set units alpha omnibus, target)", n_reps = sum(ok_a),
      power = mean(b$p_alpha[ok_a] < ALPHA), effect = set_ratio,
      n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv,
      notes = sprintf("boot_reps=%d", BOOT_REPS))
    rows[[length(rows) + 1]] <- sv_row(STUDY,
      sprintf("btl_efrm power set_ratio=%.1f panel=1", set_ratio),
      "off-target flag rate (panel units phi omnibus)", n_reps = sum(ok_p),
      power = mean(b$p_phi[ok_p] < ALPHA), effect = set_ratio,
      n_attempted = b$n_attempted, n_refused = b$n_refused, n_nonconv = b$n_nonconv,
      notes = "off-target: no panel-unit effect planted here")
  }
  cat("family (d) elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

if (chunk == "combine") {
  p1 <- "tools/simval/results/custom-wald-tests-chunk1.csv"
  p2 <- "tools/simval/results/custom-wald-tests-chunk2.csv"
  r1 <- read.csv(p1, stringsAsFactors = FALSE)
  r2 <- read.csv(p2, stringsAsFactors = FALSE)
  sv_write(rbind(r1, r2), STUDY)
} else {
  rows_a <- NULL; rows_bcd <- NULL
  if (chunk %in% c("1", "all")) rows_a <- run_family_a()
  if (chunk %in% c("2", "all")) {
    rb <- run_family_b(); rc <- run_family_c(); rd <- run_family_d()
    rows_bcd <- rbind(rb, rc, rd)
  }
  if (chunk == "1") sv_write(rows_a, paste0(STUDY, "-chunk1"))
  if (chunk == "2") sv_write(rows_bcd, paste0(STUDY, "-chunk2"))
  if (chunk == "all") sv_write(rbind(rows_a, rows_bcd), STUDY)
}

cat("\nTOTAL elapsed:", as.numeric(Sys.time() - t_study, units = "secs"), "s\n")
