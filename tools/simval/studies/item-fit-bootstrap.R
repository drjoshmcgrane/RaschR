# Calibration and power of the parametric bootstrap null for the item fit
# statistics (fit_bootstrap), against the asymptotic references rasch()
# reports. Every one of those statistics is computed at estimated person
# locations and referred to a distribution derived as though they were known,
# and each is miscalibrated by an amount that moves with the sample: the
# class-interval chi-square grows against a fixed reference, while the fit
# residual and the Wilson-Hilferty standardisations are referred to a normal
# whose spread they do not have. This study measures the size of each at five
# sample sizes, checks that the bootstrap restores nominal error, and
# separates what the correction can do from what the statistics can.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP  <- as.integer(Sys.getenv("SV_REPS", "100"))
B     <- as.integer(Sys.getenv("SV_B", "99"))
CORES <- as.integer(Sys.getenv("SV_CORES", "8"))
L     <- 8L
TARGET <- 4L        # difficulty -0.29 with the spread below

# Placement decides what a power scenario measures. With difficulty spread
# over (-2, 2) the first item is the easiest on the test, and a departure
# planted there sits where almost nobody's responses can discriminate it: it
# reads as no power when the truth is no leverage. Departures go on a central
# item.
one <- function(r, n, slope, tag) {
  disc <- rep(1, L); disc[TARGET] <- slope
  d <- simulate_rasch(n, L, difficulty = c(-2, 2), discrimination = disc,
                      seed = 1e6 + 1e4 * tag + r)
  f <- tryCatch(rasch(d, id = "id"), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  bs <- tryCatch(suppressWarnings(
    fit_bootstrap(f, B = B, seed = 2e6 + 1e4 * tag + r, workers = 1L)),
    error = function(e) NULL)
  if (is.null(bs)) return(NULL)
  data.frame(rep = r, n = n, effect = slope,
             item = seq_len(L), planted = seq_len(L) == TARGET,
             # asymptotic readings, each at its nominal 5%
             a_chisq = f$items$p < .05,
             a_fitres = abs(f$items$fit_resid) > 1.96,
             a_fitres_conv = abs(f$items$fit_resid) > 2.5,
             a_infit_z = abs(f$items$infit_z) > 1.96,
             a_outfit_z = abs(f$items$outfit_z) > 1.96,
             # bootstrap readings of the same statistics
             b_chisq = bs$items$chisq_p_boot < .05,
             b_fitres = bs$items$fit_resid_p_boot < .05,
             b_infit_z = bs$items$infit_z_p_boot < .05,
             b_outfit_z = bs$items$outfit_z_p_boot < .05,
             chisq = f$items$chisq, df = f$items$df,
             fit_resid = f$items$fit_resid,
             a_total = f$total_chisq_p < .05,
             b_total = bs$total$chisq_p_boot < .05,
             fr_sd = bs$total$fit_resid_sd,
             b_fr_sd = bs$total$fit_resid_sd_p_boot < .05)
}

run <- function(n, slope, tag) {
  reps <- parallel::mclapply(seq_len(NREP), one, n = n, slope = slope,
                             tag = tag, mc.cores = CORES)
  do.call(rbind, reps[!vapply(reps, is.null, TRUE)])
}

# Rates pool L correlated items per replicate: the Monte Carlo SE is the
# between-replicate SD of the per-replicate proportion, not the binomial
# plug-in on n_reps * L.
cluster_se <- function(d, col) {
  per <- tapply(d[[col]], d$rep, mean)
  stats::sd(per, na.rm = TRUE) / sqrt(sum(!is.na(per)))
}

LABEL <- c(chisq = "item-trait chi-square", fitres = "fit residual",
           infit_z = "infit z", outfit_z = "outfit z")

rows <- list()
for (n in c(250L, 500L, 1000L, 2000L, 4000L)) {
  d <- run(n, 1, tag = n)
  nrep <- length(unique(d$rep))
  scen <- sprintf("null, %d persons x %d items", n, L)
  for (st in names(LABEL)) for (m in c("a", "b")) {
    col <- paste0(m, "_", st)
    rows[[length(rows) + 1L]] <- sv_row(
      "item fit bootstrap", scen,
      sprintf("%s Type I error (%s)", LABEL[[st]],
              if (m == "a") "asymptotic" else "parametric bootstrap"),
      n_reps = nrep, n_attempted = NREP, n_nonconv = NREP - nrep,
      type1 = mean(d[[col]], na.rm = TRUE),
      mc_override = list(type1 = cluster_se(d, col)))
  }
  # the conventional fit-residual cut, reported so its drift is on the record
  rows[[length(rows) + 1L]] <- sv_row(
    "item fit bootstrap", scen,
    "fit residual beyond the conventional 2.5 (asymptotic)",
    n_reps = nrep, type1 = mean(d$a_fitres_conv, na.rm = TRUE),
    mc_override = list(type1 = cluster_se(d, "a_fitres_conv")),
    notes = sprintf("fit residual mean %.3f SD %.3f; mean chisq %.1f on df %d",
                    mean(d$fit_resid, na.rm = TRUE),
                    stats::sd(d$fit_resid, na.rm = TRUE),
                    mean(d$chisq, na.rm = TRUE),
                    stats::median(d$df, na.rm = TRUE)))
  tot <- d[!duplicated(d$rep), ]
  for (m in c("a", "b"))
    rows[[length(rows) + 1L]] <- sv_row(
      "item fit bootstrap", scen,
      sprintf("total item-trait Type I error (%s)",
              if (m == "a") "asymptotic" else "parametric bootstrap"),
      n_reps = nrep, type1 = mean(tot[[paste0(m, "_total")]], na.rm = TRUE))
  rows[[length(rows) + 1L]] <- sv_row(
    "item fit bootstrap", scen,
    "item fit residual SD Type I error (parametric bootstrap)",
    n_reps = nrep, type1 = mean(tot$b_fr_sd, na.rm = TRUE),
    notes = sprintf("observed fit residual SD %.3f (the convention reads it against 1)",
                    mean(tot$fr_sd, na.rm = TRUE)))
}

# Over-discrimination is the departure the class-interval statistic is built
# to see. Under-discrimination is the case it cannot see: a flatter item
# varies less across the intervals, so it carries LESS of the selection bias
# than a fitting item does and its chi-square comes out smaller. That is a
# property of the statistic, not of the correction, and the fit residual
# detects the same departure readily -- which is the argument for calibrating
# both rather than choosing between them.
for (cond in list(list(n = 500L, s = 2.5, tag = 91L),
                  list(n = 2000L, s = 2.5, tag = 92L),
                  list(n = 500L, s = 0.5, tag = 93L),
                  list(n = 2000L, s = 0.5, tag = 94L))) {
  d <- run(cond$n, cond$s, cond$tag)
  nrep <- length(unique(d$rep))
  scen <- sprintf("slope %.1f on item %d, %d persons x %d items",
                  cond$s, TARGET, cond$n, L)
  pl <- d[d$planted, ]; cl <- d[!d$planted, ]
  for (st in names(LABEL)) for (m in c("a", "b")) {
    col <- paste0(m, "_", st)
    lab <- sprintf("%s (%s)", LABEL[[st]],
                   if (m == "a") "asymptotic" else "parametric bootstrap")
    rows[[length(rows) + 1L]] <- sv_row(
      "item fit bootstrap", scen, sprintf("power on the planted item, %s", lab),
      n_reps = nrep, effect = cond$s, power = mean(pl[[col]], na.rm = TRUE))
    rows[[length(rows) + 1L]] <- sv_row(
      "item fit bootstrap", scen, sprintf("Type I error on the clean items, %s", lab),
      n_reps = nrep, effect = cond$s, type1 = mean(cl[[col]], na.rm = TRUE),
      mc_override = list(type1 = cluster_se(cl, col)))
  }
}

out <- do.call(rbind, rows)
sv_write(out, "item-fit-bootstrap")
cat(sprintf("%d rows written\n", nrow(out)))
