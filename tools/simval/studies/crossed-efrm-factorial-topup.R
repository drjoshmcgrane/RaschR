# STUDY: crossed-efrm-factorial-topup
#
# Fresh-seed adjudication of the 6.7% crossed-EFRM factorial familywise rate
# seen in audit-adjusted-dependence.csv. Records term-level tail behaviour as
# well as the complete Holm family. Run from the package root.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "crossed-efrm-factorial-topup"
ALPHA <- 0.05
REPS <- suppressWarnings(as.integer(Sys.getenv("SV_EFRM_REPS", "2000")))
CORES <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
if (!is.finite(REPS) || REPS < 1L) stop("SV_EFRM_REPS must be positive")
if (!is.finite(CORES) || CORES < 1L) stop("SV_CORES must be positive")

gen_null <- function(seed) {
  set.seed(seed)
  N <- 400L
  region <- factor(rep(c("North", "South"), each = N / 2L))
  cohort <- factor(rep(rep(c("A", "B"), each = N / 4L), 2L))
  theta <- stats::rnorm(N, sd = 1.3)
  delta <- seq(-1.8, 1.8, length.out = 12L)
  X <- vapply(delta, function(b)
    stats::rbinom(N, 1, stats::plogis(theta - b)), integer(N))
  colnames(X) <- sprintf("I%02d", seq_len(ncol(X)))
  data.frame(X, region, cohort, check.names = FALSE)
}

one <- function(r) {
  blank <- c(region_p = NA, cohort_p = NA, interaction_p = NA,
             raw_any = NA, holm_any = NA, region_w = NA, cohort_w = NA,
             interaction_w = NA, refused = 0, nonconv = 0)
  d <- gen_null(7100000L + r)
  f <- tryCatch(rasch_efrm(
    d, item_sets = list(scale = sprintf("I%02d", 1:12)),
    groups = c("region", "cohort"), boot_reps = 0),
    error = function(e) NULL)
  if (is.null(f)) { blank["refused"] <- 1; return(blank) }
  if (!isTRUE(f$est$converged)) {
    blank["nonconv"] <- 1
    return(blank)
  }
  tt <- f$phi_factorial_tests
  expected <- stats::p.adjust(tt$p, method = "holm")
  if (!isTRUE(all.equal(tt$p_adj, expected, tolerance = 1e-14,
                        check.attributes = FALSE)))
    stop("reported p_adj does not match Holm adjustment")
  ip <- match(c("region", "cohort", "region:cohort"), tt$term)
  c(region_p = tt$p[ip[1]], cohort_p = tt$p[ip[2]],
    interaction_p = tt$p[ip[3]], raw_any = any(tt$p < ALPHA),
    holm_any = any(tt$p_adj < ALPHA), region_w = tt$wald[ip[1]],
    cohort_w = tt$wald[ip[2]], interaction_w = tt$wald[ip[3]],
    refused = 0, nonconv = 0)
}

cat(sprintf("Crossed-EFRM fresh-seed top-up: %d replicates on %d core(s)\n",
            REPS, CORES))
z <- if (CORES > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(seq_len(REPS), one, mc.cores = CORES,
                     mc.set.seed = FALSE)
} else lapply(seq_len(REPS), one)
bad <- vapply(z, inherits, logical(1), "try-error")
if (any(bad)) stop("worker failure: ", as.character(z[[which(bad)[1L]]]))
z <- do.call(rbind, z)
ok <- is.finite(z[, "holm_any"])
note <- paste(
  "fresh seeds; 400 persons, 100 in each balanced crossed cell, 12",
  "dichotomous items in one set; analytic pairwise Godambe covariance")
rows <- list(
  sv_row(STUDY, "crossed EFRM null", "probability of any raw p below .05",
         sum(ok), familywise = mean(z[ok, "raw_any"]),
         n_attempted = REPS, n_refused = sum(z[, "refused"]),
         n_nonconv = sum(z[, "nonconv"]), notes = note),
  sv_row(STUDY, "crossed EFRM null", "Holm factorial familywise error",
         sum(ok), familywise = mean(z[ok, "holm_any"]),
         n_attempted = REPS, n_refused = sum(z[, "refused"]),
         n_nonconv = sum(z[, "nonconv"]), notes = note))
for (term in c("region", "cohort", "interaction")) {
  pv <- z[ok, paste0(term, "_p")]
  w <- z[ok, paste0(term, "_w")]
  rows[[length(rows) + 1L]] <- sv_row(
    STUDY, paste("crossed EFRM null", term), "marginal Type I error",
    sum(ok), type1 = mean(pv < ALPHA), n_attempted = REPS,
    n_refused = sum(z[, "refused"]), n_nonconv = sum(z[, "nonconv"]),
    notes = paste(note, sprintf("mean Wald %.3f; empirical 95th percentile %.3f",
                               mean(w), unname(stats::quantile(w, .95)))))
}
sv_write(do.call(rbind, rows), STUDY)
