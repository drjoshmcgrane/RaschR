# STUDY: alpha-correction
#
# Validates the truncated-score-moment correction of the EFRM item-set
# unit (alpha) linking estimator, which replaces the naive
# var(u_hat) - mean(se^2) true-score variance. Diagnosis (2026-08-17):
# at 8 dichotomous items/set the naive correction under-recovers the
# true-score variance by >50% per set (reported SE^2 overstates the
# actual error variance; WLE shrinkage makes cov(u, error) < 0), and
# the imperfect cancellation between sets biased the log unit ratio by
# +0.046 (externally confirmed against an unbiased TAM 2PL slope-group
# anchor in cross-package-validation.csv). The corrected estimator was
# unbiased in the diagnosis harness (-0.0015 +/- 0.0131).
# Cells: null size + SE calibration, alternative bias + power, the old
# estimator's worst case (5 items/set), the full two-group frame grid,
# and the TAM-anchored bias comparison rerun.
# Serial. Rscript tools/simval/studies/alpha-correction.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "alpha-correction"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M")),
                          sprintf(...), "\n")

run_cell <- function(scen, n_reps, seed0, ratio, ips = 8, ng = 1,
                     gur = 1, npg = 500, boot = 300) {
  la2 <- se2 <- pom <- rep(NA_real_, n_reps); n_ref <- 0L
  for (r in seq_len(n_reps)) {
    s <- simulate_efrm(n_per_group = npg, items_per_set = ips, n_sets = 2,
                       n_groups = ng, set_unit_ratio = ratio,
                       group_unit_ratio = gur, seed = seed0 + r)
    tr <- attr(s, "truth")
    f <- tryCatch(rasch_efrm(s, item_sets = tr$item_sets, groups = "group",
                             id = "id", boot_reps = boot),
                  error = function(e) NULL)
    if (is.null(f)) { n_ref <- n_ref + 1L; next }
    at <- f$alpha_table
    la2[r] <- log(at$alpha[at$set == "set2"])
    se2[r] <- at$se_log_alpha[at$set == "set2"]
    om <- f$efrm_vs_rasch$unit_omnibus
    pom[r] <- om$p[grepl("alpha", om$term)][1]
  }
  ok <- is.finite(la2) & is.finite(se2)
  truth_la2 <- log(ratio) / 2                   # geometric-mean-centred
  list(bias = mean(la2[ok]) - truth_la2, emp_sd = sd(la2[ok]),
       mean_se = mean(se2[ok]),
       size = mean(pom[ok] < 0.05),
       cover = mean(abs(la2[ok] - truth_la2) <= 1.96 * se2[ok]),
       n_ok = sum(ok), n_ref = n_ref, n_reps = n_reps, scen = scen)
}
put <- function(z, what) {
  add(z$scen, sprintf("log alpha[set2] %s", what), z$n_ok,
      bias = z$bias, emp_sd = z$emp_sd, mean_se = z$mean_se,
      coverage95 = z$cover,
      type1 = if (what == "null") z$size else NA_real_,
      power = if (what != "null") z$size else NA_real_,
      n_attempted = z$n_reps, n_refused = z$n_ref,
      notes = sprintf("alpha omnibus rejection %.4f; corrected linking estimator", z$size))
  tick("%s: bias %+.4f se_ratio %.3f size/power %.3f cover %.3f",
       z$scen, z$bias, z$emp_sd / z$mean_se, z$size, z$cover)
}

put(run_cell("null, 8 items/set, 1 group, N=500", 400L, 20e3, 1), "null")
put(run_cell("ratio 1.4, 8 items/set, 1 group, N=500", 200L, 21e3, 1.4), "alt")
put(run_cell("ratio 1.4, 5 items/set (old worst case), N=500", 150L, 22e3, 1.4,
             ips = 5), "alt")
put(run_cell("null, 8 items/set, 2 groups (phi 1.15), N=500/group", 200L,
             23e3, 1, ng = 2, gur = 1.15), "null")

## TAM-anchored bias comparison, post-correction (cross-package cell C rerun)
R <- 25L; lr <- matrix(NA_real_, R, 2); n_refC <- 0L
for (r in seq_len(R)) {
  s <- simulate_efrm(n_per_group = 1000, items_per_set = 8, n_sets = 2,
                     n_groups = 1, set_unit_ratio = 1.4, seed = 900 + r)
  tr <- attr(s, "truth"); its <- unlist(tr$item_sets, use.names = FALSE)
  fe <- tryCatch(rasch_efrm(s, item_sets = tr$item_sets, groups = "group",
                            id = "id", boot_reps = 50), error = function(e) NULL)
  if (is.null(fe)) { n_refC <- n_refC + 1L; next }
  al <- fe$alpha_table$alpha[match(c("set1","set2"), fe$alpha_table$set)]
  mt <- TAM::tam.mml.2pl(resp = as.matrix(s[, its]), irtmodel = "2PL.groups",
                         est.slopegroups = rep(1:2, each = 8),
                         control = list(progress = FALSE))
  a <- mt$B[its, 2, 1]
  lr[r, ] <- c(log(al[2]/al[1]), log(mean(a[9:16])/mean(a[1:8])))
}
okC <- stats::complete.cases(lr); lt <- log(1.4)
add("TAM-anchored log unit ratio, post-correction (truth 0.3365)",
    "bias of corrected rasch_efrm vs TAM", sum(okC),
    bias = mean(lr[okC,1]) - lt, emp_sd = sd(lr[okC,1]),
    n_attempted = R, n_refused = n_refC,
    notes = sprintf("TAM anchor bias %+.4f (sd %.4f); pre-correction efrm was +0.0462 on these seeds (cross-package-validation.csv); mean |log-ratio diff| %.4f",
      mean(lr[okC,2]) - lt, sd(lr[okC,2]), mean(abs(lr[okC,1] - lr[okC,2]))))
tick("TAM cell: efrm %+.4f, TAM %+.4f", mean(lr[okC,1]) - lt, mean(lr[okC,2]) - lt)

sv_write(do.call(rbind, rows), "alpha-correction")
cat(sprintf("TOTAL elapsed: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
