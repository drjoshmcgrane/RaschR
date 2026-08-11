# STUDY: structural-missingness
#
# The round-1 battery exercised blanket 25% MCAR. This study exercises the
# STRUCTURAL missingness cases the package explicitly supports and
# documents (identification/connectedness guards in R/estimation.R;
# ignorability discussion in R/rasch.R; the DIF within-cell centring note
# in R/dif.R; residual_pca's disjoint-column refusal in R/dimensionality.R):
#
#   (a) WEAK LINKS: two 10-item forms sharing 3 common items, at several
#       total sample sizes -- does rasch() proceed or refuse, and when it
#       proceeds, item recovery/SE calibration for common vs form-unique
#       items, benchmarked against a complete-data reference design.
#   (b) UNEQUAL EXPOSURE: one item answered by only 10% of persons -- bias/
#       RMSE/SE-ratio/coverage for that item specifically, against a fully
#       observed contrast item from the same fits.
#   (c) MISSINGNESS DEPENDING ON OBSERVED GROUP (PRINCIPAL): group B
#       structurally skips the 3 hardest items (booklet determined
#       entirely by the observed group factor -- textbook MAR-given-group).
#       Tests whether dif_anova() false-flags the SHARED items merely
#       because group B's booklet is shorter, and whether group-B person
#       measures are biased by the shorter test.
#   (d) SPARSE CATEGORIES x MISSINGNESS: PCM, 4 categories with a
#       deliberately rare top category, plus 30% MCAR -- threshold recovery
#       and the weak-category guard's (.pcml_weak_thresholds, R/estimation.R)
#       firing rate, against a no-missingness reference and a
#       moderate-category contrast.
#   (e) DIMENSIONALITY UNDER STRUCTURAL MISSINGNESS: residual_pca /
#       dimensionality_test false-positive behaviour with NO planted second
#       dimension, under two structural designs -- the literal two-form
#       disjoint booklet of (a) (e1), and a 4-form overlapping spiral
#       booklet designed so the diagnostic can actually run (e2).
#
# Every scenario reports refusal/non-convergence as first-class results.
#
# ---------------------------------------------------------------------
# BENCHMARKED replicate costs (this machine, single core; see the
# scratch benchmarking that preceded this script):
#   (a) 17-item two-form booklet fit, N=20..600:       0.05-0.40s
#   (b) 15-item fit, one item at 10% exposure, N=500:  ~0.15-0.20s
#   (c) 15-item booklet + dif_anova(), N=600:          ~0.37s/rep
#   (d) PCM, 10 items, 4 categories, N=500, 30% MCAR:  ~2.7-2.9s/rep
#       (dominates the budget; run as its own concurrent chunk)
#   (e1) 17-item disjoint booklet + dimensionality_test, N=300: ~0.40s
#   (e2) 4-form 16-item overlapping-spiral booklet, N=500:      ~0.18s
#
# An early attempt at (e2) used a genuinely per-person RANDOM item subset
# (~unique missingness pattern per person): cost exploded to ~9s/rep
# because person_estimates() root-finds once per DISTINCT missing-data
# pattern, and unique-per-person patterns defeat that sharing. The
# overlapping-spiral design (a handful of distinct forms) avoids this and
# is also the realistic structural case (planned/matrix-sample booklets
# reuse a small number of forms) -- documented here as a discovered
# performance characteristic, not a package defect.
#
# The PRINCIPAL claim is (c)'s shared-item DIF false-flag rate: cost
# ~0.37s/rep is comfortably under the 1.5s/rep threshold, so it runs at
# 1500 replicates (>= the required 1000).
#
# Chunk 1 (this process alone) = scenario (d), the expensive one (observed
# 20.9 CPU-minutes). Chunk 2 = scenarios (a), (b), (c), (e) (observed 14.1
# CPU-minutes). Run concurrently (as committed) the wall clock is ~21
# minutes, comfortably inside the ~70-minute budget; running with no
# argument reproduces the same rows serially in one process (~35
# CPU-minutes).
#
# Run from the package root:
#   Rscript tools/simval/studies/structural-missingness.R          # everything, serially
#   Rscript tools/simval/studies/structural-missingness.R 1        # chunk 1: (d)
#   Rscript tools/simval/studies/structural-missingness.R 2        # chunk 2: (a)(b)(c)(e)

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

STUDY <- "structural-missingness"
args  <- commandArgs(trailingOnly = TRUE)
chunk <- if (length(args) >= 1) suppressWarnings(as.integer(args[1])) else NA

t_start <- Sys.time()
rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- sv_row(STUDY, ...)

logit <- function(p) log(p / (1 - p))

# =======================================================================
# (a) WEAK LINKS: two 10-item forms sharing 3 common items {C1,C2,C3}.
# =======================================================================
run_weak_links <- function(n_reps, seed_base) {
  common <- setNames(c(-1, 0, 1), c("C1", "C2", "C3"))
  uniqA  <- setNames(seq(-3, 3, length.out = 7), paste0("A", 1:7))
  uniqB  <- setNames(seq(-3.3, 3.3, length.out = 7), paste0("B", 1:7))
  dA <- c(common, uniqA); dB <- c(common, uniqB)
  allitems <- c(names(common), names(uniqA), names(uniqB))
  truth <- c(common, uniqA, uniqB)
  roles <- list(common = names(common), uniqueA = names(uniqA),
                uniqueB = names(uniqB))

  sim_booklet <- function(Ntot) {
    N <- Ntot %/% 2L
    thetaA <- rnorm(N); thetaB <- rnorm(N)
    XA <- matrix(rbinom(N * 10, 1, plogis(outer(thetaA, dA, "-"))), N, 10,
                 dimnames = list(NULL, names(dA)))
    XB <- matrix(rbinom(N * 10, 1, plogis(outer(thetaB, dB, "-"))), N, 10,
                 dimnames = list(NULL, names(dB)))
    XAf <- matrix(NA_real_, N, length(allitems), dimnames = list(NULL, allitems))
    XBf <- matrix(NA_real_, N, length(allitems), dimnames = list(NULL, allitems))
    XAf[, names(dA)] <- XA; XBf[, names(dB)] <- XB
    rbind(XAf, XBf)
  }
  sim_complete <- function(Ntot) {
    theta <- rnorm(Ntot)
    X <- matrix(rbinom(Ntot * length(truth), 1,
                       plogis(outer(theta, truth, "-"))), Ntot, length(truth),
                dimnames = list(NULL, allitems))
    X
  }

  for (Ntot in c(20L, 40L, 100L, 600L)) {
    est_b <- se_b <- matrix(NA_real_, n_reps, length(allitems),
                            dimnames = list(NULL, allitems))
    est_r <- se_r <- est_b
    n_ref_b <- 0L; n_ref_r <- 0L
    for (i in seq_len(n_reps)) {
      set.seed(seed_base + i)
      Xb <- sim_booklet(Ntot)
      Xr <- sim_complete(Ntot)
      fb <- tryCatch(suppressWarnings(rasch(Xb)), error = function(e) e)
      fr <- tryCatch(suppressWarnings(rasch(Xr)), error = function(e) e)
      if (inherits(fb, "error")) n_ref_b <- n_ref_b + 1L else {
        idx <- match(fb$items$item, allitems)
        est_b[i, idx] <- fb$items$location; se_b[i, idx] <- fb$items$se
      }
      if (inherits(fr, "error")) n_ref_r <- n_ref_r + 1L else {
        idx <- match(fr$items$item, allitems)
        est_r[i, idx] <- fr$items$location; se_r[i, idx] <- fr$items$se
      }
    }
    scen <- sprintf("weak_links N=%d", Ntot)
    add(scen, "booklet_fit_refusal", n_reps - n_ref_b,
        n_attempted = n_reps, n_refused = n_ref_b,
        notes = "connectivity is structural (3 common items always shared): near-zero refusal is the expected finding, not a bug")
    add(scen, "reference_fit_refusal", n_reps - n_ref_r,
        n_attempted = n_reps, n_refused = n_ref_r,
        notes = "complete-data reference (every person answers all 17 items)")
    for (role in names(roles)) {
      its <- roles[[role]]
      tv <- truth[its]
      for (design in c("booklet", "reference")) {
        E <- if (design == "booklet") est_b[, its, drop = FALSE] else est_r[, its, drop = FALSE]
        S <- if (design == "booklet") se_b[, its, drop = FALSE] else se_r[, its, drop = FALSE]
        d <- sweep(E, 2, tv, "-")
        dv <- as.vector(d); sv <- as.vector(S)
        add(scen, sprintf("item_%s_%s", role, design), sum(is.finite(dv)),
            bias = mean(dv, na.rm = TRUE),
            emp_sd = stats::sd(dv, na.rm = TRUE),
            mean_se = mean(sv, na.rm = TRUE),
            coverage95 = sv_coverage(dv, sv, 0),
            notes = sprintf("%s: %d item(s), pooled over items x replicates; %s is answered by %s",
                            design, length(its),
                            paste(its, collapse = ","),
                            if (design == "booklet")
                              (if (role == "common") "BOTH forms (full N)" else "one form only (half N)")
                            else "everyone (full N, all items)"))
      }
    }
  }
  invisible(NULL)
}

# =======================================================================
# (b) UNEQUAL EXPOSURE: item L1 answered by only 10% of persons.
# =======================================================================
run_unequal_exposure <- function(n_reps, seed_base) {
  I <- 14L
  uniq <- setNames(seq(-2.6, 2.6, length.out = I), paste0("U", 1:I))
  L1 <- c(L1 = 0)
  truth <- c(L1, uniq)
  contrast_item <- names(uniq)[which.min(abs(uniq))]   # fully-observed item nearest L1's difficulty
  N <- 500L
  scen <- "unequal_exposure N=500"
  est <- se <- matrix(NA_real_, n_reps, length(truth),
                      dimnames = list(NULL, names(truth)))
  n_refused <- 0L; n_L1_dropped <- 0L
  for (i in seq_len(n_reps)) {
    set.seed(seed_base + i)
    theta <- rnorm(N)
    X <- matrix(rbinom(N * length(truth), 1,
                       plogis(outer(theta, truth, "-"))), N, length(truth),
                dimnames = list(NULL, names(truth)))
    exposed <- sample.int(N, size = round(0.10 * N))
    X[setdiff(seq_len(N), exposed), "L1"] <- NA
    f <- tryCatch(suppressWarnings(rasch(X)), error = function(e) e)
    if (inherits(f, "error")) { n_refused <- n_refused + 1L; next }
    if (!"L1" %in% f$items$item) n_L1_dropped <- n_L1_dropped + 1L
    idx <- match(f$items$item, names(truth))
    est[i, idx] <- f$items$location; se[i, idx] <- f$items$se
  }
  add(scen, "fit_refusal", n_reps - n_refused, n_attempted = n_reps,
      n_refused = n_refused)
  add(scen, "L1_dropped_as_constant", n_reps, n_attempted = n_reps,
      n_refused = n_L1_dropped,
      notes = "L1 answered by only 50/500 persons; occasionally constant among its few respondents and dropped by .prepare_X -- excluded from the L1 recovery row below, not silently averaged in")
  for (it in c("L1", contrast_item)) {
    d <- est[, it] - truth[it]; s <- se[, it]
    add(scen, paste0("item_", it), sum(is.finite(d)),
        bias = mean(d, na.rm = TRUE), emp_sd = stats::sd(d, na.rm = TRUE),
        mean_se = mean(s, na.rm = TRUE),
        coverage95 = sv_coverage(d, s, 0),
        notes = if (it == "L1") "10% exposure (~50/500 persons); location 0"
                else sprintf("fully-observed contrast item nearest L1's true difficulty (%s)", it))
  }
  invisible(NULL)
}

# =======================================================================
# (c) PRINCIPAL: group B structurally skips the 3 hardest items.
# =======================================================================
run_group_mar <- function(n_reps, seed_base) {
  I <- 15L
  # mean-zero item bank: rasch() recentres the ESTIMATED item locations to
  # mean 0 every replicate (see tools/simval/studies/wle-coverage.R's
  # identification-indeterminacy note); an off-centre TRUE bank would
  # mechanically bias the person-theta comparison against true theta by
  # exactly the true bank's mean offset, confounding this scenario's actual
  # question (does the shorter group-B booklet bias person measures) with a
  # pure origin artifact. An earlier version of this scenario used
  # seq(-2, 2.4, ...) (true mean +0.2) and both groups showed bias ~= -0.20,
  # i.e. exactly -offset, confirming the mechanism before this fix.
  delta <- seq(-2.2, 2.2, length.out = I); names(delta) <- sprintf("I%02d", 1:I)
  shared <- names(delta)[1:12]; hard3 <- names(delta)[13:15]
  N <- 300L                               # per group
  scen <- "group_MAR N=600"

  flag_mat <- matrix(NA, n_reps, length(shared), dimnames = list(NULL, shared))
  hard_tested <- rep(NA, n_reps); hard_flagged <- rep(NA, n_reps)
  thetaA_est <- thetaA_se <- thetaA_true <- numeric(0)
  thetaB_est <- thetaB_se <- thetaB_true <- numeric(0)
  n_refused <- 0L

  for (i in seq_len(n_reps)) {
    set.seed(seed_base + i)
    thA <- rnorm(N); thB <- rnorm(N)
    XA <- matrix(rbinom(N * I, 1, plogis(outer(thA, delta, "-"))), N, I,
                dimnames = list(NULL, names(delta)))
    XB <- matrix(rbinom(N * I, 1, plogis(outer(thB, delta, "-"))), N, I,
                dimnames = list(NULL, names(delta)))
    XB[, hard3] <- NA
    X <- rbind(XA, XB)
    grp <- rep(c("A", "B"), each = N)
    f <- tryCatch(suppressWarnings(rasch(X, factors = grp)), error = function(e) e)
    if (inherits(f, "error") || !isTRUE(f$est$converged)) { n_refused <- n_refused + 1L; next }
    dd <- tryCatch(dif_anova(f), error = function(e) e)
    if (inherits(dd, "error")) { n_refused <- n_refused + 1L; next }
    s <- dd$summary
    flag_mat[i, ] <- s$uniform_DIF[match(shared, s$item)]
    tested3 <- hard3 %in% s$item
    hard_tested[i] <- any(tested3)
    if (any(tested3)) hard_flagged[i] <- any(s$uniform_DIF[s$item %in% hard3])

    pA <- f$person[f$person$grp == "A", ]
    pB <- f$person[f$person$grp == "B", ]
    thetaA_est <- c(thetaA_est, pA$theta); thetaA_se <- c(thetaA_se, pA$se)
    thetaA_true <- c(thetaA_true, thA)
    thetaB_est <- c(thetaB_est, pB$theta); thetaB_se <- c(thetaB_se, pB$se)
    thetaB_true <- c(thetaB_true, thB)
  }
  n_ok <- n_reps - n_refused
  add(scen, "fit_refusal", n_ok, n_attempted = n_reps, n_refused = n_refused)

  # PRINCIPAL: shared-item DIF Type I, pooled cluster-robust MC SE (12
  # correlated item-flags per replicate; the harness's own prescription
  # for pooled-over-units rates).
  per_rep_rate <- rowMeans(flag_mat, na.rm = TRUE)
  ok_rep <- is.finite(per_rep_rate)
  overall_rate <- mean(flag_mat, na.rm = TRUE)
  cluster_se <- stats::sd(per_rep_rate[ok_rep]) / sqrt(sum(ok_rep))
  add(scen, "dif_type1_shared_items", sum(ok_rep), type1 = overall_rate,
      mc_override = list(type1 = cluster_se),
      notes = "PRINCIPAL claim: uniform_DIF flag rate (BH-adjusted alpha=0.05, dif_anova default) on the 12 items BOTH groups answer, under the null (item locations identical for both groups; group B's booklet omits only the 3 hardest items). mc SE is the cluster-robust sd(per-replicate rate)/sqrt(n_reps), not the naive binomial plug-in, because the 12 item flags per replicate are correlated.")
  fam <- apply(flag_mat, 1, function(v) if (all(is.na(v))) NA else any(v, na.rm = TRUE))
  add(scen, "dif_familywise_shared_items", sum(!is.na(fam)),
      familywise = mean(fam, na.rm = TRUE),
      notes = "any of the 12 shared items flagged, per replicate")

  add(scen, "dif_skipped_items_tested", n_ok,
      refusal_rate = 1 - mean(hard_tested, na.rm = TRUE), n_attempted = n_ok,
      n_refused = sum(!hard_tested, na.rm = TRUE),
      notes = "rate the 3 group-B-skipped items even ENTER dif_anova()'s summary table (group B has zero data for them, so complete.cases() drops every group-B row and the item's group factor has a single remaining level -- the guarantee under test is that this is refused/untested, never silently flagged)")
  n_tested_hard <- sum(hard_tested, na.rm = TRUE)
  add(scen, "dif_skipped_items_flagged_when_tested", n_tested_hard,
      type1 = if (n_tested_hard > 0) mean(hard_flagged, na.rm = TRUE) else NA_real_,
      notes = if (n_tested_hard == 0)
        "never tested in any replicate (see dif_skipped_items_tested); no false-flag opportunity arose"
      else sprintf("tested in %d/%d replicates; flag rate among those", n_tested_hard, n_ok))

  dA <- thetaA_est - thetaA_true; sA <- thetaA_se
  add(scen, "person_theta_bias_groupA", sum(is.finite(dA)),
      bias = mean(dA, na.rm = TRUE), emp_sd = stats::sd(dA, na.rm = TRUE),
      mean_se = mean(sA, na.rm = TRUE), coverage95 = sv_coverage(dA, sA, 0),
      notes = "group A: full 15-item test, WLE vs true simulated theta, pooled over persons x replicates")
  dB <- thetaB_est - thetaB_true; sB <- thetaB_se
  add(scen, "person_theta_bias_groupB", sum(is.finite(dB)),
      bias = mean(dB, na.rm = TRUE), emp_sd = stats::sd(dB, na.rm = TRUE),
      mean_se = mean(sB, na.rm = TRUE), coverage95 = sv_coverage(dB, sB, 0),
      notes = "group B: shortened 12-item test (3 hardest items structurally missing), WLE vs true simulated theta -- same theta distribution as group A, so any bias here is attributable to the shorter test, not a true group difference")
  invisible(NULL)
}

# =======================================================================
# (d) SPARSE CATEGORIES x MISSINGNESS: PCM, 4 categories, rare top
#     category, 30% MCAR (+ a no-missingness reference and a
#     moderate-category contrast).
# =======================================================================
run_sparse_categories <- function(n_reps_main, n_reps_ref, n_reps_contrast, seed_base) {
  I <- 10L; m <- 3L
  delta <- seq(-1.2, 1.2, length.out = I); names(delta) <- sprintf("I%02d", 1:I)
  step_rare <- c(-1.8, 0.2, 4.2)       # pattern C: marginal category props ~ 7/33/55/6%
  step_moderate <- c(-1.2, 0, 1.2)     # even spacing: no rare category
  taus_rare <- lapply(seq_len(I), function(i)
    .sim_thresholds(delta[i], m, 1, FALSE, pattern = step_rare))
  taus_mod <- lapply(seq_len(I), function(i)
    .sim_thresholds(delta[i], m, 1, FALSE, pattern = step_moderate))

  run_cell <- function(scen, taus, n_reps, miss_prop, seed0) {
    N <- 500L
    est <- se <- weak <- array(NA_real_, c(n_reps, I, m),
                               dimnames = list(NULL, names(delta), 1:m))
    n_refused <- 0L; n_nonconv <- 0L; any_weak <- rep(NA, n_reps)
    for (i in seq_len(n_reps)) {
      set.seed(seed0 + i)
      theta <- rnorm(N)
      X <- sapply(seq_len(I), function(j) .sim_item(theta, taus[[j]]))
      colnames(X) <- names(delta)
      if (miss_prop > 0) X[matrix(runif(N * I) < miss_prop, N, I)] <- NA
      f <- tryCatch(suppressWarnings(rasch(X, model = "PCM")), error = function(e) e)
      if (inherits(f, "error")) { n_refused <- n_refused + 1L; next }
      if (!isTRUE(f$est$converged)) { n_nonconv <- n_nonconv + 1L; next }
      th <- f$thresholds
      for (r in seq_len(nrow(th))) {
        est[i, th$item[r], th$k[r]] <- th$tau[r]
        se[i, th$item[r], th$k[r]] <- th$se[r]
        weak[i, th$item[r], th$k[r]] <- th$weak[r]
      }
      any_weak[i] <- any(th$weak)
    }
    n_ok <- n_reps - n_refused - n_nonconv
    add(scen, "fit_refusal", n_ok, n_attempted = n_reps, n_refused = n_refused)
    add(scen, "fit_nonconvergence", n_ok, n_attempted = n_reps, n_nonconv = n_nonconv)

    truth_arr <- array(NA_real_, c(I, m))
    for (j in seq_len(I)) truth_arr[j, ] <- taus[[j]]
    for (k in seq_len(m)) {
      Ek <- est[, , k]; d <- sweep(Ek, 2, truth_arr[, k], "-")
      Sk <- se[, , k]
      dv <- as.vector(d); sv <- as.vector(Sk)
      ok <- is.finite(dv) & is.finite(sv)
      add(scen, sprintf("threshold_k%d", k), sum(is.finite(dv)),
          bias = mean(dv, na.rm = TRUE), emp_sd = stats::sd(dv, na.rm = TRUE),
          mean_se = mean(sv, na.rm = TRUE),
          coverage95 = if (any(ok)) mean(abs(dv[ok]) <= 1.96 * sv[ok]) else NA_real_,
          notes = sprintf("pooled over %d items x replicates; NA se cells (weak-guard fires) excluded from mean_se/coverage but counted in n_reps via finite point estimates; POOLED emp_sd/mean_se MIXES heterogeneous item precisions (Jensen: sqrt(mean sd^2)/mean se > 1 even when every item calibrates) and, post-guard, selection -- read the per-item row below for calibration; %s",
                          I, if (k == m) "k=3 is the deliberately rare top-category threshold" else "not the rare threshold"))
      per_item <- vapply(seq_len(I), function(j) {
        okjk <- is.finite(est[, j, k]) & is.finite(se[, j, k])
        if (sum(okjk) < 30) return(NA_real_)
        stats::sd(est[okjk, j, k]) / mean(se[okjk, j, k])
      }, 0)
      fin <- is.finite(per_item)
      add(scen, sprintf("threshold_k%d_per_item_ratio", k), sum(fin),
          notes = sprintf(
            "per-item empirical-SD/mean-SE among REPORTED (non-withheld) cells: median %.2f, max %.2f (item %s); items with <30 reported cells (guard-withheld) excluded: %s",
            stats::median(per_item[fin]), if (any(fin)) max(per_item[fin]) else NA,
            if (any(fin)) names(delta)[which.max(replace(per_item, !fin, -Inf))] else "none",
            if (all(fin)) "none" else paste(names(delta)[!fin], collapse = "/")))
      Wk <- weak[, , k]
      per_rep <- rowMeans(Wk, na.rm = TRUE)
      okr <- is.finite(per_rep)
      add(scen, sprintf("weak_guard_k%d", k), sum(okr),
          type1 = mean(Wk, na.rm = TRUE),
          mc_override = list(type1 = stats::sd(per_rep[okr]) / sqrt(sum(okr))),
          notes = "guard-firing rate (SE reported NA), pooled over items x replicates, cluster-robust mc SE")
    }
    add(scen, "weak_guard_any_threshold", sum(!is.na(any_weak)),
        type1 = mean(any_weak, na.rm = TRUE),
        notes = "rate at least one threshold anywhere in the 10-item test triggers the weak-category guard (among fitted, converged replicates)")
    invisible(NULL)
  }

  run_cell("pcm_rare_top+30MCAR", taus_rare, n_reps_main, 0.30, seed_base)
  run_cell("pcm_rare_top+noMiss", taus_rare, n_reps_ref, 0.00, seed_base + 500000)
  run_cell("pcm_moderate_top+30MCAR", taus_mod, n_reps_contrast, 0.30, seed_base + 1000000)
  invisible(NULL)
}

# =======================================================================
# (e) DIMENSIONALITY UNDER STRUCTURAL MISSINGNESS.
# =======================================================================
run_dimensionality <- function(n_reps_e1, n_reps_e2, seed_base) {
  # e1: literal two-form disjoint booklet of (a), no planted 2nd dimension.
  common <- setNames(c(-1, 0, 1), c("C1", "C2", "C3"))
  uniqA  <- setNames(seq(-3, 3, length.out = 7), paste0("A", 1:7))
  uniqB  <- setNames(seq(-3.3, 3.3, length.out = 7), paste0("B", 1:7))
  dA <- c(common, uniqA); dB <- c(common, uniqB)
  allitems <- c(names(common), names(uniqA), names(uniqB))
  scen1 <- "dim_disjoint_booklet N=300"
  n_fit_refused <- 0L; n_dimtest_na <- 0L; n_ok <- 0L
  for (i in seq_len(n_reps_e1)) {
    set.seed(seed_base + i)
    N <- 150L
    thetaA <- rnorm(N); thetaB <- rnorm(N)
    XA <- matrix(rbinom(N * 10, 1, plogis(outer(thetaA, dA, "-"))), N, 10,
                dimnames = list(NULL, names(dA)))
    XB <- matrix(rbinom(N * 10, 1, plogis(outer(thetaB, dB, "-"))), N, 10,
                dimnames = list(NULL, names(dB)))
    XAf <- matrix(NA_real_, N, length(allitems), dimnames = list(NULL, allitems))
    XBf <- matrix(NA_real_, N, length(allitems), dimnames = list(NULL, allitems))
    XAf[, names(dA)] <- XA; XBf[, names(dB)] <- XB
    X <- rbind(XAf, XBf)
    f <- tryCatch(suppressWarnings(rasch(X)), error = function(e) e)
    if (inherits(f, "error")) { n_fit_refused <- n_fit_refused + 1L; next }
    n_ok <- n_ok + 1L
    dt <- dimensionality_test(f)
    if (is.na(dt$multidimensional)) n_dimtest_na <- n_dimtest_na + 1L
  }
  add(scen1, "fit_refusal", n_ok, n_attempted = n_reps_e1, n_refused = n_fit_refused)
  add(scen1, "dimtest_refusal", n_ok, refusal_rate = n_dimtest_na / max(n_ok, 1),
      n_attempted = n_ok, n_refused = n_dimtest_na,
      notes = "residual_pca refuses (multidimensional=NA, 'no respondents in common') whenever two form-unique items are never co-answered by the same person -- expected near-total refusal for a strictly two-form disjoint booklet with NO overlap beyond the 3 common items; this IS the package's documented guard, not a defect")

  # e2: 4-form overlapping spiral booklet (window covers every item pair),
  # no planted 2nd dimension -- the diagnostic CAN run here.
  n_items <- 16L; n_forms <- 4L; window <- 12L
  delta2 <- seq(-2, 2, length.out = n_items); names(delta2) <- sprintf("I%02d", 1:n_items)
  step <- n_items / n_forms
  forms <- lapply(seq_len(n_forms), function(fm) {
    start <- round((fm - 1) * step)
    ((start + 0:(window - 1)) %% n_items) + 1
  })
  scen2 <- "dim_spiral_booklet N=500"
  n_fit_refused2 <- 0L; n_dimtest_na2 <- 0L; n_ok2 <- 0L
  flagged <- logical(0); eigs <- numeric(0)
  for (i in seq_len(n_reps_e2)) {
    set.seed(seed_base + 2000000L + i)
    N <- 500L
    theta <- rnorm(N)
    form_id <- rep(seq_len(n_forms), length.out = N)
    X <- matrix(NA_real_, N, n_items, dimnames = list(NULL, names(delta2)))
    for (p in seq_len(N)) {
      its <- forms[[form_id[p]]]
      X[p, its] <- rbinom(length(its), 1, plogis(theta[p] - delta2[its]))
    }
    f <- tryCatch(suppressWarnings(rasch(X)), error = function(e) e)
    if (inherits(f, "error")) { n_fit_refused2 <- n_fit_refused2 + 1L; next }
    n_ok2 <- n_ok2 + 1L
    dt <- dimensionality_test(f)
    if (is.na(dt$multidimensional)) { n_dimtest_na2 <- n_dimtest_na2 + 1L; next }
    flagged <- c(flagged, isTRUE(dt$multidimensional))
    eigs <- c(eigs, dt$first_eigen)
  }
  add(scen2, "fit_refusal", n_ok2, n_attempted = n_reps_e2, n_refused = n_fit_refused2)
  add(scen2, "dimtest_refusal", n_ok2, refusal_rate = n_dimtest_na2 / max(n_ok2, 1),
      n_attempted = n_ok2, n_refused = n_dimtest_na2,
      notes = "every item pair is co-answered within at least one of the 4 forms by construction, so this should be near zero")
  add(scen2, "dimtest_type1_flag", length(flagged), type1 = mean(flagged),
      notes = sprintf("dimensionality_test's multidimensional flag rate under the null (no 2nd dimension planted; alpha=0.05 nominal); mean first_eigen=%.3f, sd=%.3f over %d estimable replicates -- if the spiral pattern masqueraded as a dimension this rate would sit well above 0.05",
                      mean(eigs), stats::sd(eigs), length(eigs)))
  invisible(NULL)
}

# =======================================================================
# Run
# =======================================================================
if (is.na(chunk) || chunk == 2L) {
  run_weak_links(n_reps = 300L, seed_base = 1e6)
  run_unequal_exposure(n_reps = 800L, seed_base = 2e6)
  run_group_mar(n_reps = 1500L, seed_base = 3e6)
  run_dimensionality(n_reps_e1 = 300L, n_reps_e2 = 400L, seed_base = 5e6)
}
if (is.na(chunk) || chunk == 1L) {
  run_sparse_categories(n_reps_main = 300L, n_reps_ref = 150L,
                        n_reps_contrast = 150L, seed_base = 4e6)
}

out <- do.call(rbind, rows)
suffix <- if (!is.na(chunk)) paste0("-chunk", chunk) else ""
sv_write(out, paste0(STUDY, suffix))
cat(sprintf("total wall time: %.1fs\n", as.numeric(Sys.time() - t_start, units = "secs")))
