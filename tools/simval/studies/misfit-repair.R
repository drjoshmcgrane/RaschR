# STUDY: misfit-repair
#
# Does the diagnose-and-drop workflow actually restore a planted unit ratio?
#
# alpha-set-misfit.csv and frame-invariance-power.csv each stop half way:
# the first shows that misfit concentrated in one item set distorts the set
# unit, the second shows that frame_invariance() detects items that behave
# differently across person frames. Neither closes the loop. This study
# runs the loop a user would run -- fit, read the diagnostic, drop what it
# flags with drop_items(), refit -- and asks whether the planted ratio comes
# back.
#
# Two arms, because the two units are threatened by different departures
# and read by different diagnostics:
#
#   alpha (item-set units): threatened by misfit CONCENTRATED IN ONE SET,
#     read by the item fit statistics in fit$items. One group, two sets.
#
#   phi (person-frame units): threatened by items that behave DIFFERENTLY
#     ACROSS FRAMES, read by frame_invariance(). Two groups, one set.
#
# Four conditions per damage type, the last three paired on the same
# replicates so the comparison is not confounded by simulation noise:
#
#   clean     no misfit planted            -- what recovery looks like at best
#   damaged   misfit planted, no repair    -- the cost of ignoring it
#   repaired  drop what the diagnostic flags, at two screening thresholds:
#             the multiplicity-adjusted one the test reports, and a lenient
#             unadjusted one, because a screen and a confirmatory test are
#             not the same job
#   oracle    drop the planted items regardless of what is flagged -- the
#             ceiling any diagnostic could reach by dropping
#
# The oracle matters: if repaired falls short of the planted ratio but
# matches the oracle, dropping is the limit, not the diagnostic. If oracle
# recovers and repaired does not, the diagnostic is the limit.
#
# Sensitivity (planted items flagged) and false-positive rate (clean items
# flagged) are recorded per condition, along with how many replicates the
# diagnostic flagged nothing in and how many drops drop_items() refused.
#
# N = 500 per group, planted ratio 1.40, 200 replicates. Serial, ~20 min.
#   Rscript tools/simval/studies/misfit-repair.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "misfit-repair"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

ratio <- 1.40
lt <- log(ratio)
u <- c(ratio^-0.5, ratio^0.5)   # geometrically centred pair
N <- 500L
R <- 200L

# ---------------------------------------------------------------- reporting

report <- function(label, quantity, lr, extra = list()) {
  ok <- is.finite(lr)
  add(label, quantity, sum(ok),
      bias = mean(lr[ok]) - lt, emp_sd = stats::sd(lr[ok]),
      n_attempted = R, n_refused = R - sum(ok),
      notes = c(sprintf("recovered %.3f vs planted %.3f",
                        exp(mean(lr[ok])), ratio),
                unlist(extra)))
  cat(sprintf("[%s] %-46s %.3f  (bias %+.4f, sd %.4f, n=%d)\n",
              format(Sys.time(), "%H:%M"), label,
              exp(mean(lr[ok])), mean(lr[ok]) - lt, stats::sd(lr[ok]), sum(ok)))
}

isTRUE_v <- function(x) !is.na(x) & x

# fit$items names an item by its virtual column, "S1I03:g1", because the
# frame model splits each item by the frame that took it. drop_items() wants
# the source item. Strip the frame suffix, but only where the bare name is
# one the fit holds.
bare_item <- function(f, v) {
  src <- names(f$set_of)
  intersect(unique(ifelse(v %in% src, v, sub(":[^:]*$", "", v))), src)
}

# proportion of the planted items a diagnostic caught, and of the sound ones
# it wrongly caught, averaged over replicates
hit_rates <- function(flagged, planted, all_items) {
  sens <- vapply(flagged, function(f) mean(planted %in% f), 0)
  fp <- vapply(flagged, function(f)
    mean(setdiff(all_items, planted) %in% f), 0)
  c(sens = mean(sens, na.rm = TRUE), fp = mean(fp, na.rm = TRUE))
}

refit_ratio <- function(f, drop, which_unit) {
  if (!length(drop)) f2 <- f
  else {
    f2 <- tryCatch(drop_items(f, drop, boot_reps = 0), error = function(e) NULL)
    if (is.null(f2)) return(NA_real_)
  }
  if (which_unit == "alpha") {
    a <- f2$alpha_table$alpha[match(c("set1", "set2"), f2$alpha_table$set)]
    log(a[2] / a[1])
  } else {
    p <- f2$phi_table$phi[match(c("g1", "g2"), f2$phi_table$group)]
    log(p[2] / p[1])
  }
}

# --------------------------------------------- arm A: item-set units (alpha)
#
# One group of persons takes both sets, so alpha is read person-side from
# the true-score variance ratio over the common persons. Misfit in set 1
# distorts set 1's person estimates and set 2 carries nothing to cancel it.

KA <- 8L
delta_A <- seq(-1.5, 1.5, length.out = KA)
items_A <- c(sprintf("S1I%02d", seq_len(KA)), sprintf("S2I%02d", seq_len(KA)))
isets_A <- list(set1 = items_A[seq_len(KA)], set2 = items_A[KA + seq_len(KA)])

gen_A <- function(seed, disc1, disc2) {
  set.seed(seed)
  th <- stats::rnorm(N, 0, 1.3)
  X <- cbind(
    vapply(seq_len(KA), function(i)
      stats::rbinom(N, 1, stats::plogis(u[1] * disc1[i] * (th - delta_A[i]))),
      numeric(N)),
    vapply(seq_len(KA), function(i)
      stats::rbinom(N, 1, stats::plogis(u[2] * disc2[i] * (th - delta_A[i]))),
      numeric(N)))
  colnames(X) <- items_A
  data.frame(id = sprintf("P%05d", seq_len(N)), X, group = "g1",
             check.names = FALSE)
}

fit_A <- function(d)
  tryCatch(rasch_efrm(d, item_sets = isets_A, groups = "group", id = "id",
                      boot_reps = 0), error = function(e) NULL)

run_arm_A <- function(label, disc1, planted) {
  flat <- rep(1, KA)
  dam <- rep(NA_real_, R); rep_ <- rep(NA_real_, R); orc <- rep(NA_real_, R)
  len <- rep(NA_real_, R)
  flagged <- vector("list", R); flagged2 <- vector("list", R)
  n_empty <- 0L; n_empty2 <- 0L
  for (r in seq_len(R)) {
    d <- gen_A(310000 + r, disc1, flat)
    f <- fit_A(d)
    if (is.null(f)) next
    dam[r] <- refit_ratio(f, character(0), "alpha")
    # strict: the item fit test as the fit reports it, multiplicity adjusted
    fl <- bare_item(f, f$items$item[isTRUE_v(f$items$p_adj < 0.05)])
    flagged[[r]] <- fl
    if (!length(fl)) n_empty <- n_empty + 1L
    rep_[r] <- refit_ratio(f, fl, "alpha")
    # lenient: any item whose infit is two standard units from expectation
    fl2 <- bare_item(f, f$items$item[isTRUE_v(abs(f$items$infit_z) > 2)])
    flagged2[[r]] <- fl2
    if (!length(fl2)) n_empty2 <- n_empty2 + 1L
    len[r] <- refit_ratio(f, fl2, "alpha")
    orc[r] <- refit_ratio(f, planted, "alpha")
  }
  nn <- function(x) x[!vapply(x, is.null, TRUE)]
  hr <- hit_rates(nn(flagged), planted, items_A)
  hr2 <- hit_rates(nn(flagged2), planted, items_A)
  cap <- function(h, ne) list(
    sprintf("flagged %.0f%% of planted, %.0f%% of sound items",
            100 * h["sens"], 100 * h["fp"]),
    sprintf("%d/%d replicates flagged nothing", ne, R))
  report(paste0("alpha | ", label, " | damaged"), "log alpha ratio", dam)
  report(paste0("alpha | ", label, " | repaired, item fit p_adj < .05"),
         "log alpha ratio", rep_, cap(hr, n_empty))
  report(paste0("alpha | ", label, " | repaired, |infit z| > 2"),
         "log alpha ratio", len, cap(hr2, n_empty2))
  report(paste0("alpha | ", label, " | oracle drop"), "log alpha ratio", orc)
}

# ------------------------------------------ arm B: person-frame units (phi)
#
# Two groups take the same single set, so phi is identified within frame.
# An item that behaves differently across frames breaks the invariance the
# frame model assumes, and frame_invariance() is the test for it.

KB <- 10L
delta_B <- seq(-1.6, 1.6, length.out = KB)
items_B <- sprintf("I%02d", seq_len(KB))
isets_B <- list(set1 = items_B)

gen_B <- function(seed, dif2, disc2) {
  set.seed(seed)
  mk <- function(g, shift, disc) {
    th <- stats::rnorm(N, 0, 1.3)
    X <- vapply(seq_len(KB), function(i)
      stats::rbinom(N, 1, stats::plogis(
        u[g] * disc[i] * (th - (delta_B[i] + shift[i])))), numeric(N))
    colnames(X) <- items_B
    X
  }
  X <- rbind(mk(1, rep(0, KB), rep(1, KB)), mk(2, dif2, disc2))
  data.frame(id = sprintf("P%05d", seq_len(2L * N)), X,
             group = rep(c("g1", "g2"), each = N), check.names = FALSE)
}

fit_B <- function(d)
  tryCatch(rasch_efrm(d, item_sets = isets_B, groups = "group", id = "id",
                      boot_reps = 0), error = function(e) NULL)

run_arm_B <- function(label, dif2, disc2, planted) {
  dam <- rep(NA_real_, R); rep_ <- rep(NA_real_, R); orc <- rep(NA_real_, R)
  len <- rep(NA_real_, R)
  flagged <- vector("list", R); flagged2 <- vector("list", R)
  n_empty <- 0L; n_empty2 <- 0L
  for (r in seq_len(R)) {
    d <- gen_B(320000 + r, dif2, disc2)
    f <- fit_B(d)
    if (is.null(f)) next
    dam[r] <- refit_ratio(f, character(0), "phi")
    inv <- tryCatch(frame_invariance(f), error = function(e) NULL)
    # strict: flagged as the test reports it, multiplicity adjusted
    fl <- if (is.null(inv)) character(0) else unique(c(
      inv$locations$item[isTRUE_v(inv$locations$flagged)],
      inv$discrimination$item[isTRUE_v(inv$discrimination$flagged)]))
    fl <- bare_item(f, fl)
    flagged[[r]] <- fl
    if (!length(fl)) n_empty <- n_empty + 1L
    rep_[r] <- refit_ratio(f, fl, "phi")
    # lenient: unadjusted probabilities, which is what screening rather than
    # confirming calls for
    fl2 <- if (is.null(inv)) character(0) else unique(c(
      inv$locations$item[isTRUE_v(inv$locations$p < 0.05)],
      inv$discrimination$item[isTRUE_v(inv$discrimination$p < 0.05)]))
    fl2 <- bare_item(f, fl2)
    flagged2[[r]] <- fl2
    if (!length(fl2)) n_empty2 <- n_empty2 + 1L
    len[r] <- refit_ratio(f, fl2, "phi")
    orc[r] <- refit_ratio(f, planted, "phi")
  }
  nn <- function(x) x[!vapply(x, is.null, TRUE)]
  hr <- hit_rates(nn(flagged), planted, items_B)
  hr2 <- hit_rates(nn(flagged2), planted, items_B)
  cap <- function(h, ne) list(
    sprintf("flagged %.0f%% of planted, %.0f%% of sound items",
            100 * h["sens"], 100 * h["fp"]),
    sprintf("%d/%d replicates flagged nothing", ne, R))
  report(paste0("phi | ", label, " | damaged"), "log phi ratio", dam)
  report(paste0("phi | ", label, " | repaired, frame_invariance p_adj < .05"),
         "log phi ratio", rep_, cap(hr, n_empty))
  report(paste0("phi | ", label, " | repaired, unadjusted p < .05"),
         "log phi ratio", len, cap(hr2, n_empty2))
  report(paste0("phi | ", label, " | oracle drop"), "log phi ratio", orc)
}


# ------------------------------------------------------------------- cells

flatA <- rep(1, KA)
flatB <- rep(1, KB)

# clean references, one fit each
{
  lr <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    f <- fit_A(gen_A(310000 + r, flatA, flatA))
    if (!is.null(f)) lr[r] <- refit_ratio(f, character(0), "alpha")
  }
  report("alpha | clean", "log alpha ratio", lr)
}
{
  lr <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    f <- fit_B(gen_B(320000 + r, rep(0, KB), flatB))
    if (!is.null(f)) lr[r] <- refit_ratio(f, character(0), "phi")
  }
  report("phi | clean", "log phi ratio", lr)
}

# alpha: misfit concentrated in set 1
run_arm_A("2 over-discriminating in set 1",
          replace(flatA, c(3, 6), 2.0), items_A[c(3, 6)])
run_arm_A("2 under-discriminating in set 1",
          replace(flatA, c(3, 6), 0.5), items_A[c(3, 6)])

# phi: items behaving differently across frames
run_arm_B("2 items with DIF across frames (0.8 logits)",
          replace(rep(0, KB), c(3, 7), 0.8), flatB, items_B[c(3, 7)])
run_arm_B("2 items discriminating differently across frames (1.8x)",
          rep(0, KB), replace(flatB, c(3, 7), 1.8), items_B[c(3, 7)])
run_arm_B("both departures on 4 items",
          replace(rep(0, KB), c(3, 7), 0.8), replace(flatB, c(2, 9), 1.8),
          items_B[c(2, 3, 7, 9)])

sv_write(do.call(rbind, rows), "misfit-repair")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
