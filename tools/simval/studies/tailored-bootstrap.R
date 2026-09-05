# Simulation validation: calibration of tailored_analysis() bootstrap
# inference (R/tailored.R) -- the sign-count bootstrap p, Holm adjustment
# across items, p-value resolution floor, and automatic anchor selection.
#
# Run from the package root:
#   Rscript tools/simval/studies/tailored-bootstrap.R          # everything, serially (slow)
#   Rscript tools/simval/studies/tailored-bootstrap.R 1        # chunk 1: family (a), FWER principal
#   Rscript tools/simval/studies/tailored-bootstrap.R 2        # chunk 2: family (b), power grid
#   Rscript tools/simval/studies/tailored-bootstrap.R combine  # stitch the two chunk CSVs
#
# Chunks 1 and 2 were run as two concurrent Rscript processes (the package's
# "at most 2 concurrent" budget) and then stitched with `combine`.
#
# ---------------------------------------------------------------------------
# BENCHMARK (this machine, single core, N=300 8-item dichotomous design --
# the "typical design" the task nominates):
#   tailored_analysis(se_method="bootstrap") at the DEFAULT boot_reps=999:
#       151.0s  (boot_reps_used = 997/999)
#   At 200 replicates that is 30,200s (~8.4h) for the FWER principal claim
#   alone -- utterly infeasible against a ~70-minute total budget. Dropping
#   to boot_reps=399 (the task's own suggested fallback; Holm floor for
#   m=8 items is 2*8/400 = 0.040 < 0.05, so detection stays possible):
#       70.8s/replicate (N=300, I=8, boot_reps=399, boot_reps_used=399)
#   The principal cell therefore uses 55 predeclared replicates -- well under
#   the aspirational ">= 300" for the principal claim. Reducing
#   N to 150 saves only ~16% (59.7s/replicate) -- the cost is dominated by
#   the number of bootstrap iterations x item-pair count, not person count,
#   so shrinking N is not a useful lever and was NOT used (it would also
#   confound the across-design comparisons). ALL scenarios therefore keep
#   N=300 exactly as benchmarked.
#
#   CRITICAL, task-relevant discovery: the task's own boot_reps=399
#   fallback (sized for the m=8 principal design) VIOLATES ITS OWN safety
#   criterion once test length becomes m=15 (one of the two "test length"
#   levels required for the power grid): the Holm floor is
#   2*15/400 = 0.075 > 0.05, so at boot_reps=399 NO 15-item design could
#   EVER flag a single item significant, regardless of effect size -- the
#   package's own floor warning fires and says so verbatim (confirmed by
#   benchmark: "with 396 usable replicates and 15 items the smallest
#   achievable Holm-adjusted p is 0.076 (> 0.05) ... use boot_reps >= 599").
#   The 15-item power cells therefore use boot_reps=601 (floor
#   2*15/602 = 0.0498 < 0.05), benchmarked at 231.4s/replicate (N=300,
#   I=15, boot_reps=399->601 rescaled from a direct 153.6s/396-used
#   benchmark at boot_reps=399) -- confirmed directly at boot_reps=601,
#   N=150: 186.5s (vs the 300-persons number extrapolated at ~231s;
#   N-scaling is mild, consistent with the I=8 finding above).
#
#   CONSEQUENCE FOR SCOPE: the full requested 2(guess) x 2(chance) x
#   2(length) x 2(targeting) = 16-cell power grid, at the boot_reps each
#   cell's item count requires to keep detection even POSSIBLE, costs
#   8 cells x 70.8s (I=8) + 8 cells x 231.4s (I=15) per replicate = 2418s
#   PER REPLICATE across the full grid -- ~150 replicates/cell (the task's
#   aspirational range) would need ~100 CPU-hours. That is not recoverable
#   inside a ~70-minute budget by any reallocation. The power grid is
#   therefore SCOPED DOWN, deliberately and transparently:
#     - the full 2x2x2 (guess x chance x targeting) grid at I=8 keeps ALL
#       8 cells, at boot_reps=399, with 4 predeclared replicates per cell;
#     - I=15 is run as a smaller CONFIRMATORY check (guess rate only,
#       chance=0.25 and on-target held at their I=8 defaults) rather than
#       the full 8-cell crossing, at boot_reps=601, so it still says
#       something about whether detection survives a longer test with a
#       correctly-sized bootstrap, without pretending to full-grid
#       precision it cannot afford.
#   Replicate counts per cell are consequently small (I=8: 4/cell = 32
#   total; I=15: 6/cell = 12 total) and the resulting rates carry large
#   Monte Carlo SEs, reported honestly per the harness convention. Main-
#   effect marginals (pooling replicates across the factor levels held
#   fixed) are reported alongside the raw per-cell rates to extract more
#   precision on the primary guessing-rate contrast from the same draws.
#
# METHODOLOGY: true generating item difficulties and (where planted) the
# guessing floor are fixed within a scenario; only person theta and
# responses are redrawn each replicate (via simulate_rasch(..., seed=)).
# A replicate that errors during simulation/fit/tailoring is a REFUSAL; one
# whose initial fit does not converge is a NONCONV; both are excluded from
# n_reps and reported via n_attempted/n_refused/n_nonconv. Bootstrap-
# internal non-convergence and other errors are tracked separately. For 30
# or more requested draws, tailored_analysis() refuses inference unless at
# least 90% (and at least 30) are usable.
#
# "2 low items" (the task's planted-guessing items) are read as the 2
# HARDEST items (highest difficulty; last two of the increasing-difficulty
# vector) -- this matches the package's own worked example and the round-1
# battery (tools/simval/round1/tailored/tailored_battery.R), and matches
# the mechanism: tailoring's P(success) < chance cutoff concentrates on
# hard items, so guessing signal is most detectable there.
#
# "chance parameter" {0.20, 0.25} is the ASSUMED per-item chance floor
# passed to tailored_analysis(chance=), independent of the TRUE generating
# guessing floor ("guessing rate" {0.15, 0.30}) planted via
# simulate_rasch(guessing=). This lets the grid show whether detection is
# sensitive to a chance assumption that does not exactly match the truth
# (e.g. assuming 4-option items when the true floor behaves like 5-option).
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

STUDY <- "tailored-bootstrap"
ALPHA <- 0.05
t_study <- Sys.time()

.sv_row <- function(..., n_withheld = 0L, n_metric_unavailable = 0L) {
  sv_row(..., n_withheld = n_withheld,
         n_metric_unavailable = n_metric_unavailable)
}

args <- commandArgs(trailingOnly = TRUE)
chunk <- if (length(args) >= 1) args[1] else "all"

DIFF8  <- seq(-2, 2.5, length.out = 8)
DIFF15 <- seq(-2.5, 3, length.out = 15)

# --- shared: one attempt (simulate -> fit -> tailored bootstrap) -----------
# Returns a list with $status in {"ok", "refused", "nonconv", "error"} and, when ok,
# the comparison table, boot_reps_used, responses removed, anchor count and
# whether the p-floor warning fired (it should not, given the boot_reps
# chosen below).
run_one <- function(N, diff, guess_vec, theta_mean, theta_sd, chance,
                    boot_reps, seed) {
  dat <- simulate_rasch(N, length(diff), difficulty = diff,
                        guessing = guess_vec, theta_mean = theta_mean,
                        theta_sd = theta_sd, seed = seed)
  fit0 <- tryCatch(rasch(dat, id = "id"), error = function(e) e)
  if (inherits(fit0, "error"))
    return(list(status = if (inherits(fit0, "rasch_refusal")) "refused" else "error",
                stage = "fit", msg = conditionMessage(fit0)))
  if (!isTRUE(fit0$est$converged))
    return(list(status = "nonconv"))
  floor_warned <- FALSE
  ta <- withCallingHandlers(
    tryCatch(tailored_analysis(fit0, chance = chance, se_method = "bootstrap",
                               boot_reps = boot_reps), error = function(e) e),
    warning = function(w) {
      if (grepl("smallest achievable", conditionMessage(w))) floor_warned <<- TRUE
      invokeRestart("muffleWarning")
    })
  if (inherits(ta, "error")) {
    structured <- inherits(ta, "rasch_fit_bootstrap_refusal")
    return(list(status = if (inherits(ta, "rasch_refusal")) "refused" else "error",
                stage = "tailored", msg = conditionMessage(ta),
                boot_reps = if (structured) ta$B else NULL,
                boot_used = if (structured) ta$B_used else NULL,
                boot_nonconv = if (structured) ta$B_nonconverged else NULL,
                boot_errors = if (structured) ta$B_errors else NULL))
  }
  list(status = "ok", table = ta$table, boot_used = ta$boot_reps_used,
       boot_reps = boot_reps, floor_warned = floor_warned,
       boot_nonconv = ta$boot_reps_nonconverged,
       boot_errors = ta$boot_reps_errors,
       boot_minimum = ta$boot_minimum_usable,
       n_removed = ta$n_removed, n_anchor = length(ta$anchor_items))
}

# ---------------------------------------------------------------------------
# Family (a): FWER under NO guessing -- the PRINCIPAL claim
# ---------------------------------------------------------------------------
run_family_a <- function() {
  t0 <- Sys.time()
  MAX_REPS <- 55L
  n_attempted <- 0L; n_refused <- 0L; n_nonconv <- 0L; n_error <- 0L
  tabs <- list(); boot_used_v <- integer(0); boot_reqd_v <- integer(0)
  boot_nonconv_v <- integer(0); boot_errors_v <- integer(0)
  floor_warned_any <- logical(0)
  first_ok_seed <- NA_integer_
  cat("[a] FWER principal: N=300, I=8, chance=0.25, no guessing, boot_reps=399\n")
  repeat {
    if (n_attempted >= MAX_REPS) break
    n_attempted <- n_attempted + 1L
    seed <- 90000L + n_attempted
    out <- run_one(300, DIFF8, rep(0, 8), theta_mean = 0, theta_sd = 1,
                   chance = 0.25, boot_reps = 399, seed = seed)
    if (!is.null(out$boot_reps)) {
      boot_used_v <- c(boot_used_v, out$boot_used)
      boot_reqd_v <- c(boot_reqd_v, out$boot_reps)
      boot_nonconv_v <- c(boot_nonconv_v, out$boot_nonconv)
      boot_errors_v <- c(boot_errors_v, out$boot_errors)
    }
    if (out$status == "refused") { n_refused <- n_refused + 1L; next }
    if (out$status == "nonconv") { n_nonconv <- n_nonconv + 1L; next }
    if (out$status == "error") { n_error <- n_error + 1L; next }
    if (is.na(first_ok_seed)) first_ok_seed <- seed
    tabs[[length(tabs) + 1L]] <- out$table
    floor_warned_any <- c(floor_warned_any, out$floor_warned)
  }
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("[a] done: %d attempted, %d ok, %d refused, %d nonconv, %d error, %.1fs (%.1fs/rep)\n",
              n_attempted, length(tabs), n_refused, n_nonconv, n_error, elapsed,
              elapsed / max(1, n_attempted)))

  rows <- list()
  n_reps <- length(tabs)
  # guard n_reps == 0 (every attempt refused/nonconverged): mean() of an
  # empty vector is NaN, not NA -- keep the CSV clean if that ever happens
  per_rep_type1 <- if (n_reps) vapply(tabs, function(t) mean(t$significant %in% TRUE), 0) else numeric(0)
  any_sig <- if (n_reps) vapply(tabs, function(t) any(t$significant %in% TRUE), logical(1)) else logical(0)
  degen_loss <- (boot_reqd_v - boot_used_v) / boot_reqd_v   # fraction of inner draws lost
  degen_any <- boot_used_v < boot_reqd_v
  inner <- list(n_boot_attempted = sum(boot_reqd_v),
                n_boot_used = sum(boot_used_v),
                n_boot_nonconv = sum(boot_nonconv_v),
                n_boot_errors = sum(boot_errors_v))

  rows[[1]] <- .sv_row(STUDY, "FWER null: N=300 I=8 chance=.25 boot_reps=399 no guessing",
    "familywise (any item significant per replicate, Holm-adjusted alpha=.05)",
    n_reps = n_reps, familywise = if (n_reps) mean(any_sig) else NA_real_, effect = 0,
    n_attempted = n_attempted, n_refused = n_refused,
    n_nonconv = n_nonconv, n_error = n_error,
    n_boot_attempted = inner$n_boot_attempted,
    n_boot_used = inner$n_boot_used, n_boot_nonconv = inner$n_boot_nonconv,
    n_boot_errors = inner$n_boot_errors,
    notes = "PRINCIPAL claim; 55 predeclared replicates at N=300, I=8 and boot_reps=399; see the script header for the boot_reps=999 benchmark.")
  rows[[2]] <- .sv_row(STUDY, "FWER null: N=300 I=8 chance=.25 boot_reps=399 no guessing",
    "type1_item (per item, Holm-adjusted, pooled)",
    n_reps = n_reps, type1 = if (n_reps) mean(per_rep_type1) else NA_real_,
    mc_override = list(type1 = if (n_reps > 1) stats::sd(per_rep_type1) / sqrt(n_reps) else NA_real_),
    effect = 0, n_attempted = n_attempted, n_refused = n_refused,
    n_nonconv = n_nonconv, n_error = n_error,
    n_boot_attempted = inner$n_boot_attempted, n_boot_used = inner$n_boot_used,
    n_boot_nonconv = inner$n_boot_nonconv, n_boot_errors = inner$n_boot_errors,
    notes = "pooled over 8 items/replicate; MC SE is cluster-robust (sd of per-replicate proportions / sqrt(n_reps)), not the plug-in binomial SE, since items within a replicate share one fit.")
  rows[[3]] <- .sv_row(STUDY, "FWER null: N=300 I=8 chance=.25 boot_reps=399 no guessing",
    "bootstrap inner-draw degeneracy (fraction of B inner resamples that failed)",
    n_reps = length(degen_loss),
    bias = if (length(degen_loss)) mean(degen_loss) else NA_real_,
    emp_sd = if (length(degen_loss) > 1) stats::sd(degen_loss) else NA_real_,
    n_attempted = length(degen_loss), n_refused = 0L, n_nonconv = 0L,
    n_error = 0L, n_boot_attempted = inner$n_boot_attempted,
    n_boot_used = inner$n_boot_used, n_boot_nonconv = inner$n_boot_nonconv,
    n_boot_errors = inner$n_boot_errors,
    notes = if (length(degen_loss)) sprintf(
      paste0("mean fraction lost = %.4f; %d/%d bootstrap runs lost >=1 ",
             "inner draw; denominator is every run that reached the inner ",
             "bootstrap, including a run later refused by the usable-draw rule"),
      mean(degen_loss), sum(degen_any), length(degen_loss)) else
        "no outer replicate reached the bootstrap")
  rows[[4]] <- .sv_row(STUDY, "FWER null: N=300 I=8 chance=.25 boot_reps=399 no guessing",
    "p-floor warning fired during principal loop (should never, floor=0.040<0.05)",
    n_reps = n_reps, type1 = if (n_reps) mean(floor_warned_any) else NA_real_,
    mc_override = list(type1 = 0),
    n_attempted = n_attempted, n_refused = n_refused, n_nonconv = n_nonconv,
    n_error = n_error, n_boot_attempted = inner$n_boot_attempted,
    n_boot_used = inner$n_boot_used, n_boot_nonconv = inner$n_boot_nonconv,
    n_boot_errors = inner$n_boot_errors,
    notes = "deterministic given boot_reps=399, m=8: floor 2*8/400=0.040<0.05, so the floor warning is not expected to fire in ANY replicate.")

  # deterministic floor-warning demonstration: same data, boot_reps=50
  # (floor 2*8/51=0.314>0.05, must fire) vs the boot_reps=399 runs above
  # (must not fire) -- reuses the first successfully-fit replicate's seed
  if (!is.na(first_ok_seed)) {
    demo <- run_one(300, DIFF8, rep(0, 8), theta_mean = 0, theta_sd = 1,
                    chance = 0.25, boot_reps = 50, seed = first_ok_seed)
    fired <- if (demo$status == "ok") demo$floor_warned else NA
    rows[[5]] <- .sv_row(STUDY, "floor-warning behavioural check: same draw, boot_reps=50",
      "p-floor warning fires when boot_reps too small for m=8 (floor 2*8/51=0.314>0.05)",
      n_reps = 1L, type1 = if (isTRUE(fired)) 1 else 0, mc_override = list(type1 = 0),
      n_attempted = 1L, n_refused = if (demo$status == "refused") 1L else 0L,
      n_nonconv = if (demo$status == "nonconv") 1L else 0L,
      n_error = if (demo$status == "error") 1L else 0L,
      n_boot_attempted = if (is.null(demo$boot_reps)) 0L else demo$boot_reps,
      n_boot_used = if (is.null(demo$boot_used)) 0L else demo$boot_used,
      n_boot_nonconv = if (is.null(demo$boot_nonconv)) 0L else demo$boot_nonconv,
      n_boot_errors = if (is.null(demo$boot_errors)) 0L else demo$boot_errors,
      notes = sprintf("deterministic behavioural check (n=1 by construction, not a rate): fired=%s (expect TRUE)",
                      if (is.na(fired)) "NA (replicate unavailable)" else fired))
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Family (b): POWER grid -- guessing planted on the 2 hardest items
# ---------------------------------------------------------------------------
guessed_idx  <- function(n_items) c(n_items - 1L, n_items)
guessed_nm   <- function(n_items) sprintf("I%02d", guessed_idx(n_items))

run_cell <- function(n_items, diff, guess_rate, chance, theta_mean, boot_reps,
                     n_reps, seed0) {
  guess_vec <- rep(0, n_items); guess_vec[guessed_idx(n_items)] <- guess_rate
  gnm <- guessed_nm(n_items)
  n_attempted <- 0L; n_refused <- 0L; n_nonconv <- 0L; n_error <- 0L
  det <- numeric(0); ff <- numeric(0)   # per-replicate proportions
  boot_attempted <- integer(0); boot_used <- integer(0)
  boot_nonconv <- integer(0); boot_errors <- integer(0)
  for (r in seq_len(n_reps)) {
    n_attempted <- n_attempted + 1L
    out <- run_one(300, diff, guess_vec, theta_mean = theta_mean, theta_sd = 1,
                   chance = chance, boot_reps = boot_reps, seed = seed0 + r)
    if (!is.null(out$boot_reps)) {
      boot_attempted <- c(boot_attempted, out$boot_reps)
      boot_used <- c(boot_used, out$boot_used)
      boot_nonconv <- c(boot_nonconv, out$boot_nonconv)
      boot_errors <- c(boot_errors, out$boot_errors)
    }
    if (out$status == "refused") { n_refused <- n_refused + 1L; next }
    if (out$status == "nonconv") { n_nonconv <- n_nonconv + 1L; next }
    if (out$status == "error") { n_error <- n_error + 1L; next }
    tab <- out$table
    sig <- tab$significant %in% TRUE
    is_g <- tab$item %in% gnm
    det <- c(det, mean(sig[is_g]))
    ff  <- c(ff, mean(sig[!is_g]))
  }
  list(det = det, ff = ff, n_attempted = n_attempted, n_refused = n_refused,
       n_nonconv = n_nonconv, n_error = n_error,
       n_boot_attempted = sum(boot_attempted), n_boot_used = sum(boot_used),
       n_boot_nonconv = sum(boot_nonconv), n_boot_errors = sum(boot_errors),
       n_items = n_items, guess_rate = guess_rate,
       chance = chance, theta_mean = theta_mean, boot_reps = boot_reps)
}

cell_row <- function(cell, label) {
  n_reps <- length(cell$det)
  list(
    .sv_row(STUDY, label, "detection rate (planted items, Holm-adjusted, pooled over the 2 guessed items)",
      n_reps = n_reps, power = if (n_reps) mean(cell$det) else NA_real_,
      mc_override = list(power = if (n_reps > 1) stats::sd(cell$det) / sqrt(n_reps) else NA_real_),
      effect = cell$guess_rate, n_attempted = cell$n_attempted,
      n_refused = cell$n_refused, n_nonconv = cell$n_nonconv,
      n_error = cell$n_error, n_boot_attempted = cell$n_boot_attempted,
      n_boot_used = cell$n_boot_used,
      n_boot_nonconv = cell$n_boot_nonconv,
      n_boot_errors = cell$n_boot_errors,
      notes = sprintf("I=%d guess=%.2f chance=%.2f theta_mean=%.0f boot_reps=%d",
                      cell$n_items, cell$guess_rate, cell$chance, cell$theta_mean, cell$boot_reps)),
    .sv_row(STUDY, label, "false-flag rate (clean items, Holm-adjusted, pooled)",
      n_reps = n_reps, type1 = if (n_reps) mean(cell$ff) else NA_real_,
      mc_override = list(type1 = if (n_reps > 1) stats::sd(cell$ff) / sqrt(n_reps) else NA_real_),
      effect = 0, n_attempted = cell$n_attempted,
      n_refused = cell$n_refused, n_nonconv = cell$n_nonconv,
      n_error = cell$n_error, n_boot_attempted = cell$n_boot_attempted,
      n_boot_used = cell$n_boot_used,
      n_boot_nonconv = cell$n_boot_nonconv,
      n_boot_errors = cell$n_boot_errors,
      notes = sprintf("I=%d guess=%.2f chance=%.2f theta_mean=%.0f boot_reps=%d",
                      cell$n_items, cell$guess_rate, cell$chance, cell$theta_mean, cell$boot_reps))
  )
}

run_family_b <- function() {
  t0 <- Sys.time()
  rows <- list()

  # --- I=8 full 2x2x2 grid: guess x chance x targeting, boot_reps=399 ------
  REPS_8 <- 4L
  cells8 <- list()
  seed_base <- 100000L
  cat(sprintf("[b] I=8 grid: 8 cells x %d reps, boot_reps=399\n", REPS_8))
  for (guess_rate in c(0.15, 0.30)) for (chance in c(0.20, 0.25))
    for (theta_mean in c(0, 1)) {
      key <- sprintf("g%.2f_c%.2f_t%.0f", guess_rate, chance, theta_mean)
      seed_base <- seed_base + 1000L
      cell <- run_cell(8, DIFF8, guess_rate, chance, theta_mean, 399L,
                       REPS_8, seed_base)
      cells8[[key]] <- cell
      label <- sprintf("I=8 power: guess=%.2f chance=%.2f target=%s",
                       guess_rate, chance, if (theta_mean == 0) "on" else "+1logit")
      rows <- c(rows, cell_row(cell, label))
    }
  cat(sprintf("[b] I=8 grid done: %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))

  # main-effect marginals (pool raw per-replicate proportions across the
  # other two factors, for more precision on the primary contrast)
  pool_by <- function() {
    for (g in c(0.15, 0.30)) {
      keys <- names(cells8)[vapply(cells8, function(c) c$guess_rate == g, TRUE)]
      det <- unlist(lapply(cells8[keys], `[[`, "det"))
      ff  <- unlist(lapply(cells8[keys], `[[`, "ff"))
      na  <- sum(vapply(cells8[keys], `[[`, 0L, "n_attempted"))
      nr  <- sum(vapply(cells8[keys], `[[`, 0L, "n_refused"))
      nn  <- sum(vapply(cells8[keys], `[[`, 0L, "n_nonconv"))
      ne  <- sum(vapply(cells8[keys], `[[`, 0L, "n_error"))
      nba <- sum(vapply(cells8[keys], `[[`, 0L, "n_boot_attempted"))
      nbu <- sum(vapply(cells8[keys], `[[`, 0L, "n_boot_used"))
      nbn <- sum(vapply(cells8[keys], `[[`, 0L, "n_boot_nonconv"))
      nbe <- sum(vapply(cells8[keys], `[[`, 0L, "n_boot_errors"))
      n_reps <- length(det)
      rows[[length(rows) + 1L]] <<- .sv_row(STUDY,
        sprintf("I=8 power MARGINAL: guess=%.2f (pooled over chance, targeting)", g),
        "detection rate (planted items, pooled marginal)",
        n_reps = n_reps, power = if (n_reps) mean(det) else NA_real_,
        mc_override = list(power = if (n_reps > 1) stats::sd(det) / sqrt(n_reps) else NA_real_),
        effect = g, n_attempted = na, n_refused = nr, n_nonconv = nn,
        n_error = ne, n_boot_attempted = nba, n_boot_used = nbu,
        n_boot_nonconv = nbn, n_boot_errors = nbe,
        notes = "main effect of guessing rate, pooling the 4 chance x targeting cells at this rate")
      rows[[length(rows) + 1L]] <<- .sv_row(STUDY,
        sprintf("I=8 power MARGINAL: guess=%.2f (pooled over chance, targeting)", g),
        "false-flag rate (clean items, pooled marginal)",
        n_reps = n_reps, type1 = if (n_reps) mean(ff) else NA_real_,
        mc_override = list(type1 = if (n_reps > 1) stats::sd(ff) / sqrt(n_reps) else NA_real_),
        effect = 0, n_attempted = na, n_refused = nr, n_nonconv = nn,
        n_error = ne, n_boot_attempted = nba, n_boot_used = nbu,
        n_boot_nonconv = nbn, n_boot_errors = nbe,
        notes = "main effect of guessing rate, pooling the 4 chance x targeting cells at this rate")
    }
  }
  pool_by()

  # --- I=15 confirmatory check: guess rate only, boot_reps=601 -------------
  # chance=0.25 and on-target held at the I=8 default cell; NOT a full
  # 8-cell crossing (see header: infeasible within budget). boot_reps=601
  # is the smallest value keeping the Holm floor (2*15/(B+1)) below 0.05
  # for m=15 items -- boot_reps=399 would make detection IMPOSSIBLE here.
  REPS_15 <- 6L
  cat(sprintf("[b] I=15 confirmatory: 2 cells x %d reps, boot_reps=601\n", REPS_15))
  seed_base <- 200000L
  for (guess_rate in c(0.15, 0.30)) {
    seed_base <- seed_base + 1000L
    cell <- run_cell(15, DIFF15, guess_rate, 0.25, 0, 601L, REPS_15, seed_base)
    label <- sprintf("I=15 power (confirmatory): guess=%.2f chance=0.25 target=on", guess_rate)
    rows <- c(rows, cell_row(cell, label))
  }
  cat(sprintf("[b] I=15 confirmatory done, total family (b): %.1fs\n",
              as.numeric(Sys.time() - t0, units = "secs")))
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
if (chunk == "combine") {
  p1 <- "tools/simval/results/tailored-bootstrap-chunk1.csv"
  p2 <- "tools/simval/results/tailored-bootstrap-chunk2.csv"
  r1 <- read.csv(p1, stringsAsFactors = FALSE)
  r2 <- read.csv(p2, stringsAsFactors = FALSE)
  sv_write(sv_bind_rows(r1, r2), STUDY)
} else {
  rows_a <- NULL; rows_b <- NULL
  if (chunk %in% c("1", "all")) rows_a <- run_family_a()
  if (chunk %in% c("2", "all")) rows_b <- run_family_b()
  if (chunk == "1") sv_write(rows_a, paste0(STUDY, "-chunk1"))
  if (chunk == "2") sv_write(rows_b, paste0(STUDY, "-chunk2"))
  if (chunk == "all") sv_write(rbind(rows_a, rows_b), STUDY)
}

cat("\nTOTAL elapsed:", as.numeric(Sys.time() - t_study, units = "secs"), "s\n")
