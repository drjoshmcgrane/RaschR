# STUDY: efrm-unit-multiplicity-supported
#
# Familywise-error validation for the EFRM unit decisions in a supported
# two-group, two-set design. The separate frame-unit-multiplicity study keeps
# the four-item-set boundary cell and its refusal accounting.
# Run from the package root. Set SV_CORES on Unix-like systems to parallelise.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "efrm-unit-multiplicity-supported"
ALPHA <- 0.05
REPS <- suppressWarnings(as.integer(Sys.getenv("SV_EFRM_REPS", "500")))
BOOT <- suppressWarnings(as.integer(Sys.getenv("SV_EFRM_BOOT", "80")))
CORES <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
if (!is.finite(REPS) || REPS < 1L)
  stop("SV_EFRM_REPS must be a positive integer")
if (!is.finite(BOOT) || BOOT < 30L)
  stop("SV_EFRM_BOOT must be an integer of at least 30")
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L

check_adjustment <- function(tab) {
  if (is.null(tab) || !nrow(tab)) return(FALSE)
  ok <- is.finite(tab$p)
  expected <- rep(NA_real_, nrow(tab))
  expected[ok] <- stats::p.adjust(tab$p[ok], method = "holm")
  isTRUE(all.equal(tab$p_adj, expected, tolerance = 1e-14,
                   check.attributes = FALSE))
}

one <- function(r) {
  blank <- stats::setNames(rep(NA_real_, 4L),
    c("raw_omnibus_any", "adjusted_omnibus_any",
      "raw_followup_any", "adjusted_followup_any"))
  d <- simulate_efrm(200, 6, n_sets = 2, n_groups = 2,
                     set_unit_ratio = 1, group_unit_ratio = 1,
                     seed = 640000L + r)
  err <- NULL
  f <- tryCatch(rasch_efrm(
    d, item_sets = attr(d, "truth")$item_sets, groups = "group",
    se_method = "hybrid", boot_reps = BOOT, workers = 1L,
    seed = 740000L + r),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    })
  if (is.null(f))
    return(list(value = blank, refused = 1L,
                nonconv = 0L, reason = err))
  if (!isTRUE(f$est$converged))
    return(list(value = blank, refused = 0L,
                nonconv = 1L, reason = "calibration did not converge"))
  omnibus <- f$efrm_vs_rasch$unit_omnibus
  followups <- f$efrm_vs_rasch$unit_tests
  if (!check_adjustment(omnibus) || !check_adjustment(followups))
    stop("reported adjusted probabilities do not match Holm adjustment")
  list(value = c(
    raw_omnibus_any = any(omnibus$p < ALPHA, na.rm = TRUE),
    adjusted_omnibus_any = any(omnibus$p_adj < ALPHA, na.rm = TRUE),
    raw_followup_any = any(followups$p < ALPHA, na.rm = TRUE),
    adjusted_followup_any = any(followups$p_adj < ALPHA, na.rm = TRUE)),
    refused = 0L, nonconv = 0L, reason = "")
}

cat(sprintf("Supported EFRM null family: %d replicates, B=%d, on %d core(s)\n",
            REPS, BOOT, CORES))
z <- if (CORES > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(seq_len(REPS), one, mc.cores = CORES,
                     mc.set.seed = FALSE)
} else lapply(seq_len(REPS), one)
bad <- vapply(z, inherits, logical(1), "try-error")
if (any(bad)) stop("worker failure: ", as.character(z[[which(bad)[1L]]]))
values <- do.call(rbind, lapply(z, function(x) x$value))
refused <- sum(vapply(z, function(x) x$refused, integer(1)))
nonconv <- sum(vapply(z, function(x) x$nonconv, integer(1)))
ok <- is.finite(values[, "adjusted_omnibus_any"]) &
  is.finite(values[, "adjusted_followup_any"])
reason <- vapply(z, function(x) x$reason, character(1))
reason <- sort(table(reason[nzchar(reason)]), decreasing = TRUE)
reason_note <- if (length(reason)) {
  paste(paste0(names(reason), " (", as.integer(reason), ")"),
        collapse = "; ")
} else "none"
note <- sprintf(paste(
  "hybrid covariance B=%d; Holm is applied across two omnibus tests and",
  "separately across all individual unit contrasts; refusal/non-convergence",
  "reasons: %s"), BOOT, reason_note)
make <- function(field, quantity) sv_row(
  STUDY,
  "EFRM: 2 groups x 2 sets, 200 persons/group, 6 items/set",
  quantity, n_reps = sum(ok), familywise = mean(values[ok, field]),
  n_attempted = REPS, n_refused = refused, n_nonconv = nonconv,
  notes = note)
rows <- rbind(
  make("raw_omnibus_any", "probability of any raw omnibus p below .05"),
  make("adjusted_omnibus_any", "Holm familywise error: omnibus decisions"),
  make("raw_followup_any", "probability of any raw follow-up p below .05"),
  make("adjusted_followup_any",
       "Holm familywise error: individual unit follow-ups"))
sv_write(rows, STUDY)
