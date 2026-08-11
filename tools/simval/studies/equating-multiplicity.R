# STUDY: equating-multiplicity
#
# Ordinary equating (equate_tests(), R/equating.R) and its paired-comparison
# analogue (btl_equate(), R/btl-equating.R). Two 12-item forms share a set
# of common ("anchor") items with common TRUE locations; the rest of each
# form is unique to that form. equate_tests() estimates a precision-weighted
# scale shift over the common items and tests each one against the shifted
# identity line, BH-adjusted at a hard-coded alpha = 0.05 (the adjustment
# method and alpha are not parameters of equate_tests()). btl_equate() is
# the paired-comparison counterpart: Holm-adjusted at alpha = 0.05 by
# default (both ARE parameters there).
#
# Because BH controls the FDR and, under the GLOBAL null (every common item
# truly unmoved), FDR = P(any rejection) (V = R when every null is true),
# the replicate-level any-false-flag ("familywise") rate at alpha is
# expected to sit at or below 0.05 even though BH is not a familywise
# procedure in general. That is the principal calibration claim.
#
# Design
# ------
# (a) NULL, anchor count k in {3, 5, 10}: two independent 500-person
#     samples per replicate, one per form, common anchors' true locations
#     held fixed and IDENTICAL across forms (no drift anywhere). Records
#     the replicate-level familywise any-flag rate and the per-anchor flag
#     rate at k = 3, 5, 10; k = 5 is the principal cell (most replicates).
# (b) DRIFT POWER, k = 5: exactly one anchor (A3) has a true location that
#     differs by delta in {0.4, 0.8} logits between the two forms; the
#     other four anchors are unmoved. Detection rate on A3 and false-flag
#     rate pooled over the four clean anchors.
# (c) CONTAMINATION, k = 5: two of five anchors (A2, A4) drift +0.6 logits.
#     Are both found? At least one? False-flag rate on the three clean
#     anchors? Bias of the estimated shift away from its true value of 0
#     (the shift is only ever a MINORITY-drift device; the print-time
#     documentation of equate_tests()/btl_equate() already warns that a
#     shift estimated from a badly contaminated anchor set is pulled
#     toward the movers -- this scenario quantifies that for a textbook
#     "2 of 5" contamination).
# (d) SHIFT INFERENCE, comparing the three reference modes equate_tests()
#     supports (read from .equate_loc_cov()/.equate_bank_cov() in
#     R/equating.R): (i) two separately fitted rasch() calibrations
#     (independent = TRUE), which draws each side's covariance from its own
#     cov_tau; (ii) a FIXED bank (se = 0, the true population anchor
#     locations, no sampling noise contributed); (iii) a bank carrying its
#     own joint location covariance in attr(bank, "cov_location") (built
#     here from an actual second calibration's .equate_loc_cov(), so mode
#     (iii) is a bank-shaped restatement of the same information mode (i)
#     uses live -- the two should closely agree, which is itself a useful
#     cross-check of the bank-covariance code path). Because equate_tests()
#     does not return a shift SE directly, this study reconstructs it with
#     the package's OWN unexported .equate_loc_cov() and the same
#     precision weights the drift tests use (see equate_shift_se() below);
#     nothing here is an independent re-derivation of the variance formula.
#     Evaluated under the (k = 5) null and under the (k = 5, delta = 0.8)
#     single-anchor drift of (b), reusing those replicates rather than
#     refitting.
# (e) BTL ANALOGUE, k = 5 anchors + 5 unique objects per side: null
#     familywise/per-item rates, and one drift-power cell (delta = 0.8 on
#     one anchor), using btl_equate()'s default Holm/alpha = 0.05.
#
# Method
# ------
# True generating locations are fixed once per scenario (deterministic,
# no RNG); only person samples (and hence responses) are redrawn each
# replicate. A replicate that fails to converge, or whose equating comes
# back non-inferential (fewer than 3 usable common items -- should never
# happen here since k >= 3 by construction, but checked), is excluded from
# the substantive rates and counted in nonconv_rate / refusal_rate instead
# of being silently dropped.
#
# Run from the package root:
#   Rscript tools/simval/studies/equating-multiplicity.R          # everything, sequential
#   Rscript tools/simval/studies/equating-multiplicity.R A        # chunk A (heavier)
#   Rscript tools/simval/studies/equating-multiplicity.R B        # chunk B (lighter)
# Chunks A and B write suffixed CSVs so they can run concurrently; running
# with no argument runs every scenario in one process and writes the single
# canonical CSV (bench: ~0.13-0.20s/replicate, so the full unchunked run is
# comfortably inside the ~70-minute wall-clock budget on its own).

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

STUDY <- "equating-multiplicity"
N_PERSONS <- 500L                 # per form, per side

args  <- commandArgs(trailingOnly = TRUE)
chunk <- if (length(args) >= 1) toupper(args[1]) else NA_character_

t_start <- Sys.time()
rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- sv_row(STUDY, ...)

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

# Independent draw of one dichotomous Rasch form: N persons ~ N(0,1),
# items named after `delta_named`, whose names ARE the item names.
sim_form <- function(delta_named, N) {
  I <- length(delta_named)
  theta <- rnorm(N, 0, 1)
  X <- matrix(rbinom(N * I, 1, plogis(outer(theta, delta_named, "-"))), N, I)
  colnames(X) <- names(delta_named)
  X
}

# Shift SE reconstructed from the SAME weights and joint covariance
# equate_tests() itself assembles for the drift-test denominators
# (Sg = .equate_loc_cov(fit) + .equate_loc_cov(reference); u = w / sum(w);
# Var(c0) = u' Sg u). `ref` may be a rasch fit or a bank data frame, exactly
# as accepted by equate_tests().
equate_shift_se <- function(fit1, ref, eq) {
  tab <- eq$table
  v <- tab$se_1^2 + tab$se_2^2
  usable <- is.finite(tab$difference) & is.finite(v)
  if (sum(usable) < 2) return(NA_real_)
  w <- 1 / pmax(v[usable], 1e-10)
  u <- w / sum(w)
  S1 <- .equate_loc_cov(fit1, tab$item)
  S2 <- .equate_loc_cov(ref, tab$item)
  Sg <- (S1 + S2)[usable, usable, drop = FALSE]
  sqrt(max(drop(t(u) %*% Sg %*% u), 0))
}

# Fixed (se = 0) bank at given locations.
fixed_bank <- function(items, locs) data.frame(item = items, location = unname(locs), se = 0)

# Covariance-carrying bank built from an actual second calibration `fit2`,
# restricted to `items` (a bank-shaped restatement of fit2's own
# location covariance, via the package's own .equate_loc_cov()).
cov_bank <- function(fit2, items) {
  idx <- match(items, fit2$items$item)
  bank <- data.frame(item = items, location = fit2$items$location[idx])
  attr(bank, "cov_location") <- .equate_loc_cov(fit2, items)
  bank
}

# ---------------------------------------------------------------------
# (a) NULL familywise / per-item calibration, k in {3, 5, 10}; k = 5 is
#     principal. Also drives the null half of (d).
# ---------------------------------------------------------------------
run_null <- function(k, n_reps, seed_base, do_shift_inference = FALSE) {
  anchors_true <- setNames(seq(-1.2, 1.2, length.out = k), paste0("A", seq_len(k)))
  n_uniq <- 12L - k
  uniq1 <- setNames(seq(-2, 2, length.out = n_uniq), paste0("U1_", seq_len(n_uniq)))
  uniq2 <- setNames(seq(-2.3, 2.3, length.out = n_uniq), paste0("U2_", seq_len(n_uniq)))
  d1 <- c(anchors_true, uniq1); d2 <- c(anchors_true, uniq2)
  anames <- names(anchors_true)

  any_flag <- logical(0); item_flags <- list()
  shift <- list(fitfit = numeric(0), bankfixed = numeric(0), bankcov = numeric(0))
  se_l  <- list(fitfit = numeric(0), bankfixed = numeric(0), bankcov = numeric(0))
  n_nonconv <- 0L; n_refused <- 0L; n_attempt <- 0L

  for (i in seq_len(n_reps)) {
    n_attempt <- n_attempt + 1L
    ok <- tryCatch({
      set.seed(seed_base + i)
      X1 <- sim_form(d1, N_PERSONS); X2 <- sim_form(d2, N_PERSONS)
      f1 <- suppressWarnings(rasch(X1)); f2 <- suppressWarnings(rasch(X2))
      if (!isTRUE(f1$est$converged) || !isTRUE(f2$est$converged)) {
        n_nonconv <- n_nonconv + 1L
        FALSE
      } else {
        eq <- equate_tests(f1, f2, independent = TRUE)
        if (!isTRUE(eq$inferential)) { n_refused <- n_refused + 1L; FALSE } else {
          fl <- eq$table$drift[match(anames, eq$table$item)]
          if (anyNA(fl)) { n_refused <- n_refused + 1L; FALSE } else {
            any_flag[length(any_flag) + 1] <- any(fl)
            item_flags[[length(item_flags) + 1]] <- fl
            if (do_shift_inference) {
              shift$fitfit[length(shift$fitfit) + 1] <- eq$shift
              se_l$fitfit[length(se_l$fitfit) + 1] <- equate_shift_se(f1, f2, eq)

              bf <- fixed_bank(anames, anchors_true)
              eqf <- equate_tests(f1, bf, independent = TRUE)
              shift$bankfixed[length(shift$bankfixed) + 1] <- eqf$shift
              se_l$bankfixed[length(se_l$bankfixed) + 1] <- equate_shift_se(f1, bf, eqf)

              bc <- cov_bank(f2, anames)
              eqc <- equate_tests(f1, bc, independent = TRUE)
              shift$bankcov[length(shift$bankcov) + 1] <- eqc$shift
              se_l$bankcov[length(se_l$bankcov) + 1] <- equate_shift_se(f1, bc, eqc)
            }
            TRUE
          }
        }
      }
    }, error = function(e) { n_refused <<- n_refused + 1L; FALSE })
  }

  fam <- mean(any_flag)
  # per-replicate proportion of flagged anchors: anchors within one
  # replicate share a fit, so the MC SE must come from the spread of
  # per-replicate proportions, not the plug-in binomial over all cells
  itm_rep <- vapply(item_flags, mean, 0)
  itm <- mean(itm_rep)
  scen <- sprintf("null k=%d", k)
  add(scen, "familywise_flag", length(any_flag), familywise = fam,
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = sprintf("BH-adjusted alpha=0.05 (hard-coded in equate_tests); N=%d/form; anchors truth identical both forms", N_PERSONS))
  add(scen, "per_item_flag", length(any_flag), type1 = itm,
      mc_override = list(type1 = stats::sd(itm_rep) / sqrt(length(itm_rep))),
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = sprintf("pooled over %d anchors x %d replicates; MC SE from per-replicate proportions", k, length(any_flag)))

  if (do_shift_inference) {
    for (mode in c("fitfit", "bankfixed", "bankcov")) {
      s <- shift[[mode]]; se <- se_l[[mode]]
      add("shift inference: null (k=5)", paste0("shift_", mode), length(s),
          bias = mean(s), emp_sd = sd(s), mean_se = mean(se, na.rm = TRUE),
          coverage95 = sv_coverage(s, se, 0),
          notes = switch(mode,
            fitfit = "reference = a second live rasch() fit, independent=TRUE",
            bankfixed = "reference = fixed bank (se=0) at the TRUE anchor locations",
            bankcov = "reference = bank carrying attr(,'cov_location') from the second fit's own .equate_loc_cov()"))
    }
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------
# (b)+(d) DRIFT POWER on one anchor (A3 of 5), k = 5; delta in {0.4, 0.8}.
#     delta = 0.8 also drives the drift half of (d).
# ---------------------------------------------------------------------
run_drift <- function(delta, n_reps, seed_base, do_shift_inference = FALSE) {
  k <- 5L
  anchors_true <- setNames(seq(-1.2, 1.2, length.out = k), paste0("A", seq_len(k)))
  n_uniq <- 12L - k
  uniq1 <- setNames(seq(-2, 2, length.out = n_uniq), paste0("U1_", seq_len(n_uniq)))
  uniq2 <- setNames(seq(-2.3, 2.3, length.out = n_uniq), paste0("U2_", seq_len(n_uniq)))
  drift_item <- "A3"
  anchors_form1 <- anchors_true
  anchors_form1[drift_item] <- anchors_form1[drift_item] + delta   # form 1 drifts; form 2/bank stay at truth
  d1 <- c(anchors_form1, uniq1); d2 <- c(anchors_true, uniq2)
  anames <- names(anchors_true); clean <- setdiff(anames, drift_item)

  detect <- logical(0); false_flag <- list()
  shift <- list(fitfit = numeric(0), bankfixed = numeric(0), bankcov = numeric(0))
  se_l  <- list(fitfit = numeric(0), bankfixed = numeric(0), bankcov = numeric(0))
  n_nonconv <- 0L; n_refused <- 0L; n_attempt <- 0L

  for (i in seq_len(n_reps)) {
    n_attempt <- n_attempt + 1L
    ok <- tryCatch({
      set.seed(seed_base + i)
      X1 <- sim_form(d1, N_PERSONS); X2 <- sim_form(d2, N_PERSONS)
      f1 <- suppressWarnings(rasch(X1)); f2 <- suppressWarnings(rasch(X2))
      if (!isTRUE(f1$est$converged) || !isTRUE(f2$est$converged)) {
        n_nonconv <- n_nonconv + 1L; FALSE
      } else {
        eq <- equate_tests(f1, f2, independent = TRUE)
        if (!isTRUE(eq$inferential)) { n_refused <- n_refused + 1L; FALSE } else {
          fl <- eq$table$drift[match(anames, eq$table$item)]
          if (anyNA(fl)) { n_refused <- n_refused + 1L; FALSE } else {
            names(fl) <- anames
            detect[length(detect) + 1] <- isTRUE(fl[[drift_item]])
            false_flag[[length(false_flag) + 1]] <- fl[clean]
            if (do_shift_inference) {
              shift$fitfit[length(shift$fitfit) + 1] <- eq$shift
              se_l$fitfit[length(se_l$fitfit) + 1] <- equate_shift_se(f1, f2, eq)

              bf <- fixed_bank(anames, anchors_true)
              eqf <- equate_tests(f1, bf, independent = TRUE)
              shift$bankfixed[length(shift$bankfixed) + 1] <- eqf$shift
              se_l$bankfixed[length(se_l$bankfixed) + 1] <- equate_shift_se(f1, bf, eqf)

              bc <- cov_bank(f2, anames)
              eqc <- equate_tests(f1, bc, independent = TRUE)
              shift$bankcov[length(shift$bankcov) + 1] <- eqc$shift
              se_l$bankcov[length(se_l$bankcov) + 1] <- equate_shift_se(f1, bc, eqc)
            }
            TRUE
          }
        }
      }
    }, error = function(e) { n_refused <<- n_refused + 1L; FALSE })
  }

  scen <- sprintf("drift k=5 delta=%.1f", delta)
  add(scen, "detect_drifted_anchor", length(detect), power = mean(detect), effect = delta,
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = sprintf("anchor A3 drifts +%.1f logits in form 1 only; N=%d/form", delta, N_PERSONS))
  ff_rep <- vapply(false_flag, mean, 0)
  add(scen, "false_flag_clean_anchors", length(detect), type1 = mean(ff_rep), effect = delta,
      mc_override = list(type1 = stats::sd(ff_rep) / sqrt(length(ff_rep))),
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = "pooled over the 4 undrifted anchors; MC SE from per-replicate proportions")

  if (do_shift_inference) {
    for (mode in c("fitfit", "bankfixed", "bankcov")) {
      s <- shift[[mode]]; se <- se_l[[mode]]
      add(sprintf("shift inference: drift (k=5, delta=%.1f)", delta), paste0("shift_", mode),
          length(s), bias = mean(s), emp_sd = sd(s), mean_se = mean(se, na.rm = TRUE),
          coverage95 = sv_coverage(s, se, 0), effect = delta,
          notes = "truth = 0: form 2 / bank never drift, only form 1's A3 does; a single-anchor drift is expected to pull the precision-weighted shift off 0 -- this is documented package behaviour, not a defect")
    }
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------
# (c) CONTAMINATION: 2 of 5 anchors (A2, A4) drift +0.6, k = 5
# ---------------------------------------------------------------------
run_contamination <- function(n_reps, seed_base) {
  k <- 5L; delta <- 0.6
  anchors_true <- setNames(seq(-1.2, 1.2, length.out = k), paste0("A", seq_len(k)))
  n_uniq <- 12L - k
  uniq1 <- setNames(seq(-2, 2, length.out = n_uniq), paste0("U1_", seq_len(n_uniq)))
  uniq2 <- setNames(seq(-2.3, 2.3, length.out = n_uniq), paste0("U2_", seq_len(n_uniq)))
  drift_items <- c("A2", "A4")
  anchors_form1 <- anchors_true
  anchors_form1[drift_items] <- anchors_form1[drift_items] + delta
  d1 <- c(anchors_form1, uniq1); d2 <- c(anchors_true, uniq2)
  anames <- names(anchors_true); clean <- setdiff(anames, drift_items)

  both <- logical(0); at_least_one <- logical(0); false_flag <- list(); shift <- numeric(0)
  n_nonconv <- 0L; n_refused <- 0L; n_attempt <- 0L

  for (i in seq_len(n_reps)) {
    n_attempt <- n_attempt + 1L
    tryCatch({
      set.seed(seed_base + i)
      X1 <- sim_form(d1, N_PERSONS); X2 <- sim_form(d2, N_PERSONS)
      f1 <- suppressWarnings(rasch(X1)); f2 <- suppressWarnings(rasch(X2))
      if (!isTRUE(f1$est$converged) || !isTRUE(f2$est$converged)) {
        n_nonconv <- n_nonconv + 1L
      } else {
        eq <- equate_tests(f1, f2, independent = TRUE)
        if (!isTRUE(eq$inferential)) { n_refused <- n_refused + 1L } else {
          fl <- eq$table$drift[match(anames, eq$table$item)]
          if (anyNA(fl)) { n_refused <- n_refused + 1L } else {
            names(fl) <- anames
            both[length(both) + 1] <- all(isTRUE(fl[[drift_items[1]]]), isTRUE(fl[[drift_items[2]]]))
            at_least_one[length(at_least_one) + 1] <- any(isTRUE(fl[[drift_items[1]]]), isTRUE(fl[[drift_items[2]]]))
            false_flag[[length(false_flag) + 1]] <- fl[clean]
            shift[length(shift) + 1] <- eq$shift
          }
        }
      }
    }, error = function(e) { n_refused <<- n_refused + 1L })
  }

  scen <- "contamination k=5 (2 of 5 drift +0.6)"
  add(scen, "detect_both_drifted", length(both), power = mean(both),
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = sprintf("A2 and A4 drift +%.1f logits in form 1 only; N=%d/form", delta, N_PERSONS))
  add(scen, "detect_at_least_one_drifted", length(at_least_one), power = mean(at_least_one),
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv)
  ff_rep <- vapply(false_flag, mean, 0)
  add(scen, "false_flag_clean_anchors", length(both), type1 = mean(ff_rep),
      mc_override = list(type1 = stats::sd(ff_rep) / sqrt(length(ff_rep))),
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = "pooled over the 3 undrifted anchors; MC SE from per-replicate proportions")
  add(scen, "shift_bias", length(shift), bias = mean(shift), emp_sd = sd(shift),
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = "truth = 0; the precision-weighted shift is expected to be pulled toward the 2 contaminating anchors (40% of the anchor set) -- documented package behaviour (see btl_equate()'s Details), not a defect")
  invisible(NULL)
}

# ---------------------------------------------------------------------
# (e) BTL analogue: null familywise/per-item + one drift-power cell
# ---------------------------------------------------------------------
sim_btl_form <- function(beta_named, reps_per_pair = 25L) {
  objs <- names(beta_named)
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = reps_per_pair), b = rep(pr[, 2], each = reps_per_pair),
                  stringsAsFactors = FALSE)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta_named[d$a] - beta_named[d$b]), d$a, d$b)
  d
}

run_btl_null <- function(n_reps, seed_base) {
  k <- 5L; n_uniq <- 5L
  anchors_true <- setNames(seq(-1.2, 1.2, length.out = k), paste0("A", seq_len(k)))
  uniq1 <- setNames(seq(-2, 2, length.out = n_uniq), paste0("U1_", seq_len(n_uniq)))
  uniq2 <- setNames(seq(-2.3, 2.3, length.out = n_uniq), paste0("U2_", seq_len(n_uniq)))
  b1 <- c(anchors_true, uniq1); b2 <- c(anchors_true, uniq2)
  anames <- names(anchors_true)

  any_flag <- logical(0); item_flags <- list()
  n_nonconv <- 0L; n_refused <- 0L; n_attempt <- 0L
  for (i in seq_len(n_reps)) {
    n_attempt <- n_attempt + 1L
    tryCatch({
      set.seed(seed_base + i)
      d1 <- sim_btl_form(b1); d2 <- sim_btl_form(b2)
      f1 <- suppressWarnings(btl(d1, "a", "b", "win"))
      f2 <- suppressWarnings(btl(d2, "a", "b", "win"))
      if (!isTRUE(f1$converged) || !isTRUE(f2$converged)) {
        n_nonconv <- n_nonconv + 1L
      } else {
        eq <- btl_equate(f1, f2, independent = TRUE)
        if (!isTRUE(eq$inferential)) { n_refused <- n_refused + 1L } else {
          fl <- eq$table$drifting[match(anames, eq$table$object)]
          if (anyNA(fl)) { n_refused <- n_refused + 1L } else {
            any_flag[length(any_flag) + 1] <- any(fl)
            item_flags[[length(item_flags) + 1]] <- fl
          }
        }
      }
    }, error = function(e) { n_refused <<- n_refused + 1L })
  }
  add("btl null k=5", "familywise_flag", length(any_flag), familywise = mean(any_flag),
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = "btl_equate() default: Holm-adjusted alpha=0.05; 5 anchors + 5 unique objects/side, reps_per_pair=25")
  itm_rep <- vapply(item_flags, mean, 0)
  add("btl null k=5", "per_item_flag", length(any_flag), type1 = mean(itm_rep),
      mc_override = list(type1 = stats::sd(itm_rep) / sqrt(length(itm_rep))),
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = "MC SE from per-replicate proportions")
  invisible(NULL)
}

run_btl_drift <- function(delta, n_reps, seed_base) {
  k <- 5L; n_uniq <- 5L
  anchors_true <- setNames(seq(-1.2, 1.2, length.out = k), paste0("A", seq_len(k)))
  uniq1 <- setNames(seq(-2, 2, length.out = n_uniq), paste0("U1_", seq_len(n_uniq)))
  uniq2 <- setNames(seq(-2.3, 2.3, length.out = n_uniq), paste0("U2_", seq_len(n_uniq)))
  drift_item <- "A3"
  anchors_form1 <- anchors_true; anchors_form1[drift_item] <- anchors_form1[drift_item] + delta
  b1 <- c(anchors_form1, uniq1); b2 <- c(anchors_true, uniq2)
  anames <- names(anchors_true); clean <- setdiff(anames, drift_item)

  detect <- logical(0); false_flag <- list()
  n_nonconv <- 0L; n_refused <- 0L; n_attempt <- 0L
  for (i in seq_len(n_reps)) {
    n_attempt <- n_attempt + 1L
    tryCatch({
      set.seed(seed_base + i)
      d1 <- sim_btl_form(b1); d2 <- sim_btl_form(b2)
      f1 <- suppressWarnings(btl(d1, "a", "b", "win"))
      f2 <- suppressWarnings(btl(d2, "a", "b", "win"))
      if (!isTRUE(f1$converged) || !isTRUE(f2$converged)) {
        n_nonconv <- n_nonconv + 1L
      } else {
        eq <- btl_equate(f1, f2, independent = TRUE)
        if (!isTRUE(eq$inferential)) { n_refused <- n_refused + 1L } else {
          fl <- eq$table$drifting[match(anames, eq$table$object)]
          if (anyNA(fl)) { n_refused <- n_refused + 1L } else {
            names(fl) <- anames
            detect[length(detect) + 1] <- isTRUE(fl[[drift_item]])
            false_flag[[length(false_flag) + 1]] <- fl[clean]
          }
        }
      }
    }, error = function(e) { n_refused <<- n_refused + 1L })
  }
  scen <- sprintf("btl drift k=5 delta=%.1f", delta)
  add(scen, "detect_drifted_anchor", length(detect), power = mean(detect), effect = delta,
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = "object A3 drifts in form 1 only")
  ff_rep <- vapply(false_flag, mean, 0)
  add(scen, "false_flag_clean_anchors", length(detect), type1 = mean(ff_rep), effect = delta,
      mc_override = list(type1 = stats::sd(ff_rep) / sqrt(length(ff_rep))),
      n_attempted = n_attempt, n_refused = n_refused, n_nonconv = n_nonconv,
      notes = "pooled over the 4 undrifted anchor objects; MC SE from per-replicate proportions")
  invisible(NULL)
}

# ---------------------------------------------------------------------
# Scenario schedule (~0.13-0.20s/replicate benchmarked on this machine)
# ---------------------------------------------------------------------
run_A <- function() {
  run_null(k = 5, n_reps = 4000L, seed_base = 100000L, do_shift_inference = TRUE)   # principal
  run_drift(delta = 0.8, n_reps = 1500L, seed_base = 300000L, do_shift_inference = TRUE)
  run_contamination(n_reps = 1500L, seed_base = 500000L)
  run_btl_null(n_reps = 2000L, seed_base = 700000L)
}
run_B <- function() {
  run_null(k = 3, n_reps = 2000L, seed_base = 110000L)
  run_null(k = 10, n_reps = 2000L, seed_base = 120000L)
  run_drift(delta = 0.4, n_reps = 1500L, seed_base = 310000L)
  run_btl_drift(delta = 0.8, n_reps = 1200L, seed_base = 710000L)
}

if (is.na(chunk)) {
  run_A(); run_B()
} else if (chunk == "A") {
  run_A()
} else if (chunk == "B") {
  run_B()
} else stop("unknown chunk: ", chunk, " (use A, B, or no argument)")

res <- do.call(rbind, rows)
suffix <- if (!is.na(chunk)) paste0("-chunk", chunk) else ""
sv_write(res, paste0("equating-multiplicity", suffix))

cat(sprintf("total wall time: %.1fs\n", as.numeric(Sys.time() - t_start, units = "secs")))
