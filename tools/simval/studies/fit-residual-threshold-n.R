# STUDY: fit-residual-threshold-n
#
# Why a fixed cut on a fit statistic does not survive a real sample size.
#
# fit-residual-screens.csv recommended screening item sets on
# |fit_resid| > 2, measured at 500 persons. Run against the Rosenberg
# Self-Esteem data the wording case study uses -- 6,000 respondents, ten
# items -- that screen selects seven of the ten items, |infit z| > 2
# selects eight, and the chi-square test selects all ten, so drop_items()
# refuses every one of them for emptying a set. The recommendation did not
# survive contact with the data it was written for.
#
# The obvious explanation is that real items never fit exactly while
# simulated ones do. That explanation is wrong, and this study shows why:
# the cut degrades with N even when the sound items are generated from the
# model exactly. Two of eight items in a set discriminating twice as
# steeply forces the estimated model to a compromise, under which the
# remaining items genuinely depart; that departure is fixed in size, so
# only its detectability grows. A fixed threshold on any fit statistic is
# therefore a statement about power, not about magnitude, and at survey
# sample sizes it flags everything that departs at all.
#
# Two generating conditions at 500, 2,000 and 6,000 persons:
#
#   exact    sound items generated from the model exactly
#   jittered every item carries its own small slope departure
#            (log-normal, sd 0.15), which is what real data looks like
#
# Reported per cell: how many sound items clear the cut, how many planted
# ones do, and how often the two planted items are the two largest
# departures -- the ranking, which is what a screen can actually rely on.
#
# 100 replicates. Serial, ~20 min.
#   Rscript tools/simval/studies/fit-residual-threshold-n.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "fit-residual-threshold-n"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

ratio <- 1.40
u <- c(ratio^-0.5, ratio^0.5)
K <- 8L
delta <- seq(-1.5, 1.5, length.out = K)
items <- c(sprintf("S1I%02d", seq_len(K)), sprintf("S2I%02d", seq_len(K)))
isets <- list(set1 = items[seq_len(K)], set2 = items[K + seq_len(K)])
planted <- items[c(3, 6)]
n_sound <- 2L * K - 2L
R <- 100L

run_cell <- function(N, jitter_sd) {
  lab <- sprintf("%s, N = %d",
                 if (jitter_sd > 0) "every item departs a little" else
                   "sound items fit exactly", N)
  sd_flag <- pl_flag <- rank_ok <- refused <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    set.seed(500000 + r)
    j1 <- exp(stats::rnorm(K, 0, jitter_sd))
    j2 <- exp(stats::rnorm(K, 0, jitter_sd))
    d1 <- replace(rep(1, K), c(3, 6), 2.0) * j1
    th <- stats::rnorm(N, 0, 1.3)
    X <- cbind(
      vapply(seq_len(K), function(i)
        stats::rbinom(N, 1, stats::plogis(u[1] * d1[i] * (th - delta[i]))),
        numeric(N)),
      vapply(seq_len(K), function(i)
        stats::rbinom(N, 1, stats::plogis(u[2] * j2[i] * (th - delta[i]))),
        numeric(N)))
    colnames(X) <- items
    d <- data.frame(id = sprintf("P%06d", seq_len(N)), X, group = "g1",
                    check.names = FALSE)
    f <- tryCatch(rasch_efrm(d, item_sets = isets, groups = "group",
                             id = "id", boot_reps = 0), error = function(e) NULL)
    if (is.null(f)) next
    it <- f$items
    nm <- sub(":[^:]*$", "", it$item)
    fr <- abs(it$fit_resid)
    pl <- nm %in% planted
    sd_flag[r] <- sum(fr[!pl] > 2, na.rm = TRUE)
    pl_flag[r] <- sum(fr[pl] > 2, na.rm = TRUE)
    o <- order(fr, decreasing = TRUE, na.last = TRUE)
    rank_ok[r] <- mean(pl[utils::head(o, 2L)])
    # what the cut actually does when carried through to a repair
    sel <- unique(nm[!is.na(fr) & fr > 2])
    refused[r] <- if (!length(sel)) 0 else
      is.null(tryCatch(drop_items(f, sel, boot_reps = 0),
                       error = function(e) NULL))
  }
  ok <- is.finite(sd_flag)
  add(lab, "sound items clearing |fit_resid| > 2", sum(ok),
      effect = mean(sd_flag[ok]),
      n_attempted = R, n_refused = R - sum(ok),
      notes = c(sprintf("%.1f of %d sound items flagged", mean(sd_flag[ok]),
                        n_sound),
                sprintf("%.1f of 2 planted items flagged", mean(pl_flag[ok])),
                sprintf("planted are the two largest departures in %.0f%% of replicates",
                        100 * mean(rank_ok[ok])),
                sprintf("the resulting drop is refused in %.0f%% of replicates",
                        100 * mean(refused[ok]))))
  cat(sprintf("[%s] %-42s sound %4.1f/%d  planted %3.1f/2  ranked %3.0f%%  refused %3.0f%%\n",
              format(Sys.time(), "%H:%M"), lab, mean(sd_flag[ok]), n_sound,
              mean(pl_flag[ok]), 100 * mean(rank_ok[ok]),
              100 * mean(refused[ok])))
}

for (jit in c(0, 0.15))
  for (N in c(500L, 2000L, 6000L)) run_cell(N, jit)

sv_write(do.call(rbind, rows), "fit-residual-threshold-n")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
