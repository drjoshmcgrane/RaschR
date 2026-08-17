# STUDY: alpha-n-sweep
#
# The corrected EFRM set-unit estimator across linking sample sizes
# N = 250..10,000 at the shipped configuration (hybrid, boot_reps = 300;
# 8 dichotomous items/set, ratio 1.4). Tracks: net bias (the small-sample
# Jensen term cancelling the fixed-design floor at practical N, the floor
# emerging as N grows), empirical SD vs mean reported SE (calibration
# across the range), and coverage (onset of floor-driven erosion).
# Serial. Rscript tools/simval/studies/alpha-n-sweep.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "alpha-n-sweep"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M")),
                          sprintf(...), "\n")

lt <- log(1.4) / 2
for (n in c(250, 500, 1000, 2500, 5000, 10000)) {
  R <- 60L; la <- se <- rep(NA_real_, R); n_ref <- 0L
  for (r in seq_len(R)) {
    s <- simulate_efrm(n_per_group = n, items_per_set = 8, n_sets = 2,
                       n_groups = 1, set_unit_ratio = 1.4,
                       seed = round(56e3 + n / 10) + r)
    tr <- attr(s, "truth")
    f <- tryCatch(rasch_efrm(s, item_sets = tr$item_sets, groups = "group",
                             id = "id", boot_reps = 300),
                  error = function(e) NULL)
    if (is.null(f)) { n_ref <- n_ref + 1L; next }
    at <- f$alpha_table
    la[r] <- log(at$alpha[at$set == "set2"])
    se[r] <- at$se_log_alpha[at$set == "set2"]
  }
  ok <- is.finite(la) & is.finite(se)
  add(sprintf("N = %d, 8 items/set, ratio 1.4, boot_reps 300", n),
      "log alpha[set2] bias / SE calibration / coverage", sum(ok),
      bias = mean(la[ok]) - lt, emp_sd = sd(la[ok]), mean_se = mean(se[ok]),
      coverage95 = mean(abs(la[ok] - lt) <= 1.96 * se[ok]),
      effect = n, n_attempted = R, n_refused = n_ref)
  tick("N %5d: bias %+.4f  emp_sd %.4f  mean_se %.4f  se_ratio %.3f  cover %.3f",
       n, mean(la[ok]) - lt, sd(la[ok]), mean(se[ok]),
       sd(la[ok]) / mean(se[ok]), mean(abs(la[ok] - lt) <= 1.96 * se[ok]))
}

sv_write(do.call(rbind, rows), "alpha-n-sweep")
cat(sprintf("TOTAL elapsed: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
