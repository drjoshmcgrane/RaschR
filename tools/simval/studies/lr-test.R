#!/usr/bin/env Rscript
# Simulation validation: lr_test() -- the composite-likelihood RSM-vs-PCM
# test with Kent (1982) / Godambe adjustment. See R/compare.R for the
# implementation and its documentation.
#
# (a) SIZE: empirical rejection rate of p_adj (and, for comparison, the raw
#     unadjusted composite statistic's p) at alpha = .05 under data
#     simulated from a TRUE RSM. Principal cell 500 persons x 8 items x 3
#     categories at >= 1000 replicates; a secondary crossing of
#     n_persons {300, 800} x n_items {6, 12} x n_categories {3, 4} (8 cells)
#     at 400 replicates each.
# (b) POWER: TRUE PCM data whose item threshold SHAPES diverge by design --
#     odd items get one pattern, even items its mirror image, with the
#     divergence set by an effect size `eff` (mild vs strong) -- something a
#     true RSM (one shared pattern) could never produce. Hand-built via the
#     package's own .sim_thresholds()/.sim_item()/.sim_theta() primitives
#     (used internally by simulate_rasch(model = "PCM")) so the truth is
#     deterministic given the scenario's arguments and held fixed across
#     replicates, as the methodology requires; only persons/responses are
#     redrawn each replicate.
# (c) Refusals (errors) and non-convergence (PCM or its RSM refit) are
#     counted per cell, not silently dropped.
#
# Usage:
#   Rscript tools/simval/studies/lr-test.R                  # full run -> results/lr-test.csv
#   Rscript tools/simval/studies/lr-test.R size_principal    # partial chunk
#   Rscript tools/simval/studies/lr-test.R size_secondary    # partial chunk
#   Rscript tools/simval/studies/lr-test.R power              # partial chunk
#   Rscript tools/simval/studies/lr-test.R combine             # merge the 3 partials
#
# During development the three chunks were run as (at most 2 concurrent)
# background Rscript processes and merged with `combine`; running the
# script with no arguments reproduces the identical full CSV in one
# (slower, single-process) pass -- this is the form to trust for
# reproducibility.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

args <- commandArgs(trailingOnly = TRUE)
chunk <- if (length(args)) args[1] else "all"
stopifnot(chunk %in% c("all", "size_principal", "size_secondary", "power", "combine"))

alpha <- 0.05

# ---------------------------------------------------------------------------
# (a) SIZE: one replicate under a TRUE RSM. simulate_rasch(model = "RSM")
# builds item locations from `difficulty` and the single shared threshold
# pattern from `threshold_spread` deterministically (no RNG touches the
# structural parameters for RSM); only the person draw and the item
# responses use the RNG stream seeded per replicate. So looping seed = seed0
# + r holds the RSM truth fixed and redraws persons/responses only.
# ---------------------------------------------------------------------------
size_once <- function(n_persons, n_items, n_categories, seed) {
  d <- simulate_rasch(n_persons = n_persons, n_items = n_items, model = "RSM",
                       n_categories = n_categories, difficulty = c(-2, 2),
                       threshold_spread = 1.2, seed = seed)
  fit <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
  if (inherits(fit, "error")) return(list(status = "refusal"))
  if (!isTRUE(fit$est$converged)) return(list(status = "nonconv"))
  lr <- tryCatch(lr_test(fit), error = function(e) e)
  if (inherits(lr, "error")) return(list(status = "refusal"))
  if (!isTRUE(lr$fit_rsm$est$converged)) return(list(status = "nonconv"))
  if (!is.finite(lr$p_adj)) return(list(status = "no_adjustment"))
  list(status = "ok", p_adj = lr$p_adj, p_raw = lr$p)
}

run_size_cell <- function(label, n_persons, n_items, n_categories, nrep, seed0) {
  p_adj <- p_raw <- numeric(0)
  n_refusal <- n_nonconv <- n_noadj <- 0L
  for (r in seq_len(nrep)) {
    out <- size_once(n_persons, n_items, n_categories, seed0 + r)
    if (identical(out$status, "refusal")) { n_refusal <- n_refusal + 1L; next }
    if (identical(out$status, "nonconv")) { n_nonconv <- n_nonconv + 1L; next }
    if (identical(out$status, "no_adjustment")) { n_noadj <- n_noadj + 1L; next }
    p_adj <- c(p_adj, out$p_adj); p_raw <- c(p_raw, out$p_raw)
  }
  n_ok <- length(p_adj)
  cat(sprintf("[size %-22s] np=%d ni=%d ncat=%d ok=%d/%d refusal=%d nonconv=%d noadj=%d type1_adj=%.4f type1_raw=%.4f\n",
              label, n_persons, n_items, n_categories, n_ok, nrep,
              n_refusal, n_nonconv, n_noadj,
              if (n_ok) mean(p_adj < alpha) else NA_real_,
              if (n_ok) mean(p_raw < alpha) else NA_real_))
  rbind(
    sv_row("lr-test", label, "type1_p_adj", n_ok,
           type1 = if (n_ok) mean(p_adj < alpha) else NA_real_,
           n_attempted = nrep, n_refused = n_refusal + n_noadj,
           n_nonconv = n_nonconv,
           notes = sprintf(
             "np=%d ni=%d ncat=%d; no_adjustment(df<=0)=%d/%d; Kent-adjusted p_adj, nominal alpha=.05",
             n_persons, n_items, n_categories, n_noadj, nrep)),
    sv_row("lr-test", label, "type1_p_raw", n_ok,
           type1 = if (n_ok) mean(p_raw < alpha) else NA_real_,
           n_attempted = nrep, n_refused = n_refusal + n_noadj,
           n_nonconv = n_nonconv,
           notes = "raw unadjusted composite chi-square p (conventional display); documents the anticonservatism the Godambe adjustment corrects")
  )
}

# ---------------------------------------------------------------------------
# (b) POWER: TRUE PCM data with hand-built threshold-SHAPE heterogeneity.
# step0 is the common (RSM) rating pattern; odd items get step0 + u, even
# items step0 - u, where u scales with `eff`. eff = 0 collapses to the RSM
# null (patterns identical); eff > 0 is a genuine per-item shape divergence
# no shared rating pattern can reproduce, growing with eff.
# ---------------------------------------------------------------------------
build_pcm_truth <- function(n_items, n_categories, spread, eff,
                             delta_range = c(-2, 2)) {
  m <- n_categories - 1L
  delta <- stats::setNames(seq(delta_range[1], delta_range[2], length.out = n_items),
                            sprintf("I%02d", seq_len(n_items)))
  step0 <- seq(-1, 1, length.out = m) * spread
  u <- seq(-1, 1, length.out = m) * eff
  tau <- lapply(seq_len(n_items), function(i) {
    stepi <- if (i %% 2 == 1) step0 + u else step0 - u
    rasch:::.sim_thresholds(delta[i], m, spread, FALSE, pattern = stepi)
  })
  list(delta = delta, tau = tau, m = m, n_items = n_items)
}

gen_pcm_data <- function(truth, n_persons, seed) {
  set.seed(seed)
  theta <- rasch:::.sim_theta(n_persons, 0, 1, "normal")
  X <- matrix(NA_integer_, n_persons, truth$n_items)
  for (i in seq_len(truth$n_items)) X[, i] <- rasch:::.sim_item(theta, truth$tau[[i]])
  colnames(X) <- sprintf("I%02d", seq_len(truth$n_items))
  as.data.frame(X)
}

power_once <- function(truth, n_persons, seed) {
  d <- gen_pcm_data(truth, n_persons, seed)
  fit <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
  if (inherits(fit, "error")) return(list(status = "refusal"))
  if (!isTRUE(fit$est$converged)) return(list(status = "nonconv"))
  lr <- tryCatch(lr_test(fit), error = function(e) e)
  if (inherits(lr, "error")) return(list(status = "refusal"))
  if (!isTRUE(lr$fit_rsm$est$converged)) return(list(status = "nonconv"))
  if (!is.finite(lr$p_adj)) return(list(status = "no_adjustment"))
  list(status = "ok", p_adj = lr$p_adj)
}

run_power_cell <- function(label, eff, n_persons, n_items, n_categories,
                            spread, nrep, seed0) {
  truth <- build_pcm_truth(n_items, n_categories, spread, eff)
  p_adj <- numeric(0)
  n_refusal <- n_nonconv <- n_noadj <- 0L
  for (r in seq_len(nrep)) {
    out <- power_once(truth, n_persons, seed0 + r)
    if (identical(out$status, "refusal")) { n_refusal <- n_refusal + 1L; next }
    if (identical(out$status, "nonconv")) { n_nonconv <- n_nonconv + 1L; next }
    if (identical(out$status, "no_adjustment")) { n_noadj <- n_noadj + 1L; next }
    p_adj <- c(p_adj, out$p_adj)
  }
  n_ok <- length(p_adj)
  cat(sprintf("[power %-10s] eff=%.2f ok=%d/%d refusal=%d nonconv=%d noadj=%d power=%.4f\n",
              label, eff, n_ok, nrep, n_refusal, n_nonconv, n_noadj,
              if (n_ok) mean(p_adj < alpha) else NA_real_))
  sv_row("lr-test", label, "power_p_adj", n_ok, effect = eff,
         power = if (n_ok) mean(p_adj < alpha) else NA_real_,
         n_attempted = nrep, n_refused = n_refusal + n_noadj,
         n_nonconv = n_nonconv,
         notes = sprintf(
           "np=%d ni=%d ncat=%d spread=%.2f; odd/even items' threshold SHAPE diverges by +-%.2f logits (hand-built, not simulate_rasch's random PCM patterns); no_adjustment=%d/%d",
           n_persons, n_items, n_categories, spread, eff, n_noadj, nrep))
}

# ---------------------------------------------------------------------------
# dispatch
# ---------------------------------------------------------------------------
t0 <- Sys.time()

if (chunk %in% c("all", "size_principal")) {
  principal_rows <- run_size_cell("size_principal_500x8x3", 500, 8, 3, 2000, 100000)
  if (chunk == "size_principal") sv_write(principal_rows, "lr-test-size_principal")
}

if (chunk %in% c("all", "size_secondary")) {
  secondary_grid <- expand.grid(n_persons = c(300, 800), n_items = c(6, 12),
                                 n_categories = c(3, 4))
  secondary_rows <- do.call(rbind, lapply(seq_len(nrow(secondary_grid)), function(k) {
    g <- secondary_grid[k, ]
    lbl <- sprintf("size_sec_%dx%dx%d", g$n_persons, g$n_items, g$n_categories)
    run_size_cell(lbl, g$n_persons, g$n_items, g$n_categories, 400, 200000 + k * 10000)
  }))
  if (chunk == "size_secondary") sv_write(secondary_rows, "lr-test-size_secondary")
}

if (chunk %in% c("all", "power")) {
  power_rows <- rbind(
    run_power_cell("power_mild",   eff = 0.10, n_persons = 500, n_items = 8,
                    n_categories = 4, spread = 1.2, nrep = 500, seed0 = 900000),
    run_power_cell("power_strong", eff = 0.18, n_persons = 500, n_items = 8,
                    n_categories = 4, spread = 1.2, nrep = 500, seed0 = 910000)
  )
  if (chunk == "power") sv_write(power_rows, "lr-test-power")
}

if (chunk == "all") {
  rows <- rbind(principal_rows, secondary_rows, power_rows)
  sv_write(rows, "lr-test")
}

if (chunk == "combine") {
  # The partial CSVs may have been written by concurrent Rscript processes
  # that sourced different revisions of this shared harness.R (its sv_row()
  # schema has been evolving across the round-2 battery while this study
  # ran); align on the union of columns rather than assume identical
  # schemas, so a harness.R column added/removed between two chunks' runs
  # can't crash the merge.
  parts <- c("lr-test-size_principal", "lr-test-size_secondary", "lr-test-power")
  paths <- file.path("tools/simval/results", paste0(parts, ".csv"))
  miss <- !file.exists(paths)
  if (any(miss)) stop("missing partial result(s): ", paste(paths[miss], collapse = ", "))
  dfs <- lapply(paths, utils::read.csv, stringsAsFactors = FALSE)
  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    for (nm in setdiff(all_cols, names(d))) d[[nm]] <- NA
    d[all_cols]
  })
  rows <- do.call(rbind, dfs)
  sv_write(rows, "lr-test")
}

cat("TOTAL TIME:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
