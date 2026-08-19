# STUDY: fit-residual-screens
#
# Which standardised fit statistic should drive a screen?
#
# frame_invariance() compares items across frames on infit_z, and
# misfit-repair.csv screens item sets on infit_z or the chi-square test.
# Both are the Wilson-Hilferty cube-root standardisation of a mean square.
# The package also computes the log-transformed fit residual,
# f (log y2 - log f) / sqrt(v), which is the form the Rasch measurement
# literature normally reports, and nothing here has ever compared them.
#
# The comparison matters because a screen is judged on ranking, not on any
# single value: when two of ten items discriminate more steeply, the Rasch
# compromise makes the remaining eight look flatter, so every item departs
# and what a screen needs is for the planted ones to depart most. A single
# dataset suggested infit_z fails exactly there, ranking two sound items
# above both planted ones while the fit residual ranked the planted ones
# first and second.
#
# Part A: the across-frames test in frame_invariance(). Type I error under
# a true null (frames differing by a 1.40 unit ratio and nothing else), and
# power against two items discriminating 1.8 times as steeply in frame 2,
# at 500, 1,000 and 2,000 persons per frame. Both statistics are put
# through the same rule: the difference of the two frames' values treated
# as having variance 2, two-sided at 5 per cent, unadjusted.
#
# Part B: the item-set screen in the repair loop. Same design as
# misfit-repair.csv arm A, adding |fit_resid| > 2 beside the chi-square
# test and |infit_z| > 2, and reporting the recovered unit ratio so the
# screens are judged on the quantity a user cares about rather than on
# detection alone.
#
# 200 replicates. Serial, ~35 min.
#   Rscript tools/simval/studies/fit-residual-screens.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "fit-residual-screens"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

ratio <- 1.40
lt <- log(ratio)
u <- c(ratio^-0.5, ratio^0.5)
R <- 200L

# ================================================== part A: across frames

KB <- 10L
delta_B <- seq(-1.6, 1.6, length.out = KB)
items_B <- sprintf("I%02d", seq_len(KB))

# each frame is calibrated on its own persons, so its fit statistics are
# computed against its own model and carry no unit to cancel
frame_stats <- function(seed, disc2, N) {
  set.seed(seed)
  one <- function(g, disc) {
    th <- stats::rnorm(N, 0, 1.3)
    X <- vapply(seq_len(KB), function(i)
      stats::rbinom(N, 1, stats::plogis(u[g] * disc[i] * (th - delta_B[i]))),
      numeric(N))
    colnames(X) <- items_B
    d <- data.frame(id = sprintf("P%05d", seq_len(N)), X, check.names = FALSE)
    f <- tryCatch(rasch(d, id = "id"), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    f$items[match(items_B, f$items$item), c("fit_resid", "infit_z")]
  }
  a <- one(1, rep(1, KB)); b <- one(2, disc2)
  if (is.null(a) || is.null(b)) return(NULL)
  list(fit_resid = (a$fit_resid - b$fit_resid) / sqrt(2),
       infit_z = (a$infit_z - b$infit_z) / sqrt(2))
}

run_part_A <- function(label, disc2, planted, N) {
  keep <- c("fit_resid", "infit_z")
  hit <- setNames(rep(0, 2), keep); fp <- setNames(rep(0, 2), keep)
  top <- setNames(rep(0, 2), keep)
  n_ok <- 0L
  pl <- items_B %in% planted
  for (r in seq_len(R)) {
    z <- frame_stats(410000 + r, disc2, N)
    if (is.null(z)) next
    n_ok <- n_ok + 1L
    for (s in keep) {
      v <- z[[s]]
      flag <- is.finite(v) & 2 * stats::pnorm(-abs(v)) < 0.05
      if (any(pl)) hit[s] <- hit[s] + mean(flag[pl])
      fp[s] <- fp[s] + mean(flag[!pl])
      # ranking is what a screen actually depends on: are the planted
      # items the largest departures, whatever the threshold
      if (any(pl)) {
        o <- order(abs(v), decreasing = TRUE, na.last = TRUE)
        top[s] <- top[s] + mean(pl[utils::head(o, sum(pl))])
      }
    }
  }
  for (s in keep) {
    is_null <- !any(pl)
    # sv_row builds its Monte Carlo table by name, so the rates must be
    # bare scalars rather than elements carrying the statistic's name
    hs <- unname(hit[s] / n_ok); fs <- unname(fp[s] / n_ok)
    ts <- unname(top[s] / n_ok)
    add(paste0("across frames | ", label, " | ", s),
        if (is_null) "type I error" else "power", n_ok,
        power = if (is_null) NA_real_ else hs,
        type1 = if (is_null) fs else NA_real_,
        n_attempted = R, n_refused = R - n_ok,
        notes = if (is_null)
          sprintf("flagged %.1f%% of sound items", 100 * fs)
        else c(sprintf("flagged %.0f%% of planted, %.1f%% of sound items",
                       100 * hs, 100 * fs),
               sprintf("planted items are the largest departures %.0f%% of the time",
                       100 * ts)))
    cat(sprintf("[%s] %-52s %-10s flag %5.1f%%  false %4.1f%%  top %5.1f%%\n",
                format(Sys.time(), "%H:%M"), label, s,
                100 * hs, 100 * fs, 100 * ts))
  }
}

flatB <- rep(1, KB)
steep <- replace(flatB, c(3, 7), 1.8)

for (N in c(500L, 1000L, 2000L)) {
  run_part_A(sprintf("null, unit ratio only, N = %d", N), flatB,
             character(0), N)
  run_part_A(sprintf("2 items 1.8x steeper in frame 2, N = %d", N), steep,
             items_B[c(3, 7)], N)
}

# ================================================ part B: item-set screen

KA <- 8L
N <- 500L
delta_A <- seq(-1.5, 1.5, length.out = KA)
items_A <- c(sprintf("S1I%02d", seq_len(KA)), sprintf("S2I%02d", seq_len(KA)))
isets_A <- list(set1 = items_A[seq_len(KA)], set2 = items_A[KA + seq_len(KA)])

bare_item <- function(f, v) {
  src <- names(f$set_of)
  intersect(unique(ifelse(v %in% src, v, sub(":[^:]*$", "", v))), src)
}
isTRUE_v <- function(x) !is.na(x) & x

gen_A <- function(seed, disc1) {
  set.seed(seed)
  th <- stats::rnorm(N, 0, 1.3)
  X <- cbind(
    vapply(seq_len(KA), function(i)
      stats::rbinom(N, 1, stats::plogis(u[1] * disc1[i] * (th - delta_A[i]))),
      numeric(N)),
    vapply(seq_len(KA), function(i)
      stats::rbinom(N, 1, stats::plogis(u[2] * (th - delta_A[i]))), numeric(N)))
  colnames(X) <- items_A
  data.frame(id = sprintf("P%05d", seq_len(N)), X, group = "g1",
             check.names = FALSE)
}

alpha_ratio <- function(f2) {
  a <- f2$alpha_table$alpha[match(c("set1", "set2"), f2$alpha_table$set)]
  log(a[2] / a[1])
}

run_part_B <- function(label, disc1, planted) {
  screens <- c("chi-square p_adj < .05", "|infit z| > 2", "|fit resid| > 2")
  lr <- lapply(screens, function(x) rep(NA_real_, R))
  names(lr) <- screens
  fl_all <- lapply(screens, function(x) vector("list", R))
  names(fl_all) <- screens
  none <- rep(NA_real_, R); orc <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    f <- tryCatch(rasch_efrm(gen_A(310000 + r, disc1), item_sets = isets_A,
                             groups = "group", id = "id", boot_reps = 0),
                  error = function(e) NULL)
    if (is.null(f)) next
    none[r] <- alpha_ratio(f)
    it <- f$items
    picks <- list(
      bare_item(f, it$item[isTRUE_v(it$p_adj < 0.05)]),
      bare_item(f, it$item[isTRUE_v(abs(it$infit_z) > 2)]),
      bare_item(f, it$item[isTRUE_v(abs(it$fit_resid) > 2)]))
    for (j in seq_along(screens)) {
      fl_all[[j]][[r]] <- picks[[j]]
      f2 <- if (!length(picks[[j]])) f else
        tryCatch(drop_items(f, picks[[j]], boot_reps = 0),
                 error = function(e) NULL)
      if (!is.null(f2)) lr[[j]][r] <- alpha_ratio(f2)
    }
    f3 <- tryCatch(drop_items(f, planted, boot_reps = 0), error = function(e) NULL)
    if (!is.null(f3)) orc[r] <- alpha_ratio(f3)
  }
  rep_row <- function(v, nm, fls) {
    ok <- is.finite(v)
    fls <- fls[!vapply(fls, is.null, TRUE)]
    sens <- mean(vapply(fls, function(x) mean(planted %in% x), 0))
    fpr <- mean(vapply(fls, function(x)
      mean(setdiff(items_A, planted) %in% x), 0))
    add(paste0("item set | ", label, " | ", nm), "log alpha ratio", sum(ok),
        bias = mean(v[ok]) - lt, emp_sd = stats::sd(v[ok]),
        n_attempted = R, n_refused = R - sum(ok),
        notes = c(sprintf("recovered %.3f vs planted %.3f",
                          exp(mean(v[ok])), ratio),
                  sprintf("flagged %.0f%% of planted, %.0f%% of sound items",
                          100 * sens, 100 * fpr)))
    cat(sprintf("[%s] %-30s %-24s %.3f  (n=%d, found %.0f%%, false %.0f%%)\n",
                format(Sys.time(), "%H:%M"), label, nm,
                exp(mean(v[ok])), sum(ok), 100 * sens, 100 * fpr))
  }
  add(paste0("item set | ", label, " | no repair"), "log alpha ratio",
      sum(is.finite(none)), bias = mean(none[is.finite(none)]) - lt,
      emp_sd = stats::sd(none[is.finite(none)]), n_attempted = R,
      notes = sprintf("recovered %.3f vs planted %.3f",
                      exp(mean(none[is.finite(none)])), ratio))
  cat(sprintf("[%s] %-30s %-24s %.3f\n", format(Sys.time(), "%H:%M"),
              label, "no repair", exp(mean(none[is.finite(none)]))))
  for (j in seq_along(screens)) rep_row(lr[[j]], screens[j], fl_all[[j]])
  add(paste0("item set | ", label, " | oracle drop"), "log alpha ratio",
      sum(is.finite(orc)), bias = mean(orc[is.finite(orc)]) - lt,
      emp_sd = stats::sd(orc[is.finite(orc)]), n_attempted = R,
      notes = sprintf("recovered %.3f vs planted %.3f",
                      exp(mean(orc[is.finite(orc)])), ratio))
  cat(sprintf("[%s] %-30s %-24s %.3f\n", format(Sys.time(), "%H:%M"),
              label, "oracle drop", exp(mean(orc[is.finite(orc)]))))
}

flatA <- rep(1, KA)
run_part_B("2 over-discriminating", replace(flatA, c(3, 6), 2.0),
           items_A[c(3, 6)])
run_part_B("2 under-discriminating", replace(flatA, c(3, 6), 0.5),
           items_A[c(3, 6)])

sv_write(do.call(rbind, rows), "fit-residual-screens")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
