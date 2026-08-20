# STUDY: frame-invariance-power
#
# The power of the two comparisons frame_invariance() reports.
#
# The help page for frame_invariance() states that its two channels are far
# from equally sensitive -- a location difference found about 97 per cent of
# the time at 500 persons per frame where a difference in discrimination is
# found about 16 per cent, the two converging only near 2,000 -- and those
# figures had no committed artefact behind them. Every other claim in this
# battery does. This study supplies the missing provenance, and will correct
# the help page if it disagrees.
#
# The two departures are matched on the damage they do rather than on their
# nominal size, because that is the comparison a reader needs: each moves
# the planted group-unit ratio by roughly six or seven per cent. Two items
# shifted a logit in frame 2 against two items discriminating half again as
# steeply there.
#
# Cells: a null (frames differing by a unit ratio and nothing else) for the
# type I error of each channel, then each departure, at 500, 1,000 and 2,000
# persons per frame. Flags are read at the function's default, Holm-adjusted
# within each table.
#
# 8 items, one set, two frames, planted group-unit ratio 1.40, 200
# replicates. Serial, ~40 min.
#   Rscript tools/simval/studies/frame-invariance-power.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "frame-invariance-power"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

ratio <- 1.40
u <- c(ratio^-0.5, ratio^0.5)
K <- 8L
delta <- seq(-1.5, 1.5, length.out = K)
items <- sprintf("I%02d", seq_len(K))
HIT <- c(3L, 6L)                     # the two items carrying the departure
R <- 200L

gen <- function(seed, shift, disc, N) {
  set.seed(seed)
  one <- function(g, sh, ds) {
    th <- stats::rnorm(N, 0, 1.3)
    X <- vapply(seq_len(K), function(i)
      stats::rbinom(N, 1, stats::plogis(
        u[g] * ds[i] * (th - delta[i] - sh[i]))), numeric(N))
    colnames(X) <- items
    X
  }
  X <- rbind(one(1, rep(0, K), rep(1, K)), one(2, shift, disc))
  data.frame(id = sprintf("P%05d", seq_len(2L * N)), X,
             group = rep(c("g1", "g2"), each = N), check.names = FALSE)
}

run_cell <- function(label, shift, disc, N, is_null) {
  loc_hit <- dsc_hit <- loc_fp <- dsc_fp <- rep(NA_real_, R)
  n_ref <- 0L
  planted <- items[HIT]
  sound <- setdiff(items, planted)
  for (r in seq_len(R)) {
    f <- tryCatch(rasch_efrm(gen(800000 + r, shift, disc, N),
                             items = items, item_sets = list(set1 = items),
                             groups = "group", id = "id", boot_reps = 0),
                  error = function(e) NULL)
    if (is.null(f)) { n_ref <- n_ref + 1L; next }
    inv <- tryCatch(frame_invariance(f), error = function(e) NULL)
    if (is.null(inv)) { n_ref <- n_ref + 1L; next }
    fl <- function(d) d$item[!is.na(d$flagged) & d$flagged]
    L <- fl(inv$locations); D <- fl(inv$discrimination)
    loc_hit[r] <- mean(planted %in% L); dsc_hit[r] <- mean(planted %in% D)
    loc_fp[r] <- mean(sound %in% L);    dsc_fp[r] <- mean(sound %in% D)
  }
  ok <- is.finite(loc_hit)
  emit <- function(channel, hit, fp) {
    add(sprintf("%s | N = %d | %s", label, N, channel),
        if (is_null) "type I error" else "power", sum(ok),
        power = if (is_null) NA_real_ else mean(hit[ok]),
        type1 = if (is_null) mean(fp[ok]) else NA_real_,
        mc_override = list(
          power = if (is_null) NA_real_ else stats::sd(hit[ok]) / sqrt(sum(ok)),
          type1 = if (is_null) stats::sd(fp[ok]) / sqrt(sum(ok)) else NA_real_),
        n_attempted = R, n_refused = n_ref,
        notes = sprintf("flagged %.0f%% of planted, %.1f%% of sound items",
                        100 * mean(hit[ok]), 100 * mean(fp[ok])))
  }
  emit("location", loc_hit, loc_fp)
  emit("discrimination", dsc_hit, dsc_fp)
  cat(sprintf("[%s] %-38s N=%-5d location %5.1f%% (false %4.1f%%)  discrimination %5.1f%% (false %4.1f%%)\n",
              format(Sys.time(), "%H:%M"), label, N,
              100 * mean(loc_hit[ok]), 100 * mean(loc_fp[ok]),
              100 * mean(dsc_hit[ok]), 100 * mean(dsc_fp[ok])))
}

flat <- rep(1, K); none <- rep(0, K)
shift_1 <- replace(none, HIT, 1.0)     # a logit, in frame 2
disc_15 <- replace(flat, HIT, 1.5)     # half again as steep, in frame 2

for (N in c(500L, 1000L, 2000L)) {
  run_cell("null: unit ratio only", none, flat, N, TRUE)
  run_cell("2 items shifted a logit", shift_1, flat, N, FALSE)
  run_cell("2 items 1.5x as steep", none, disc_15, N, FALSE)
}

sv_write(do.call(rbind, rows), "frame-invariance-power")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
