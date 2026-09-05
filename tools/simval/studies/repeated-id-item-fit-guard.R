# Validation of the repeated-ID item-fit guard. Exact row duplication should
# leave the calibration and person-clustered covariance unchanged, while the
# row-based fit references and their independent-row bootstrap are withheld.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

N_REPS <- as.integer(Sys.getenv("SV_REPS", "100"))
if (!is.finite(N_REPS) || N_REPS < 1L)
  stop("SV_REPS must be a positive whole number")

one <- function(r) {
  set.seed(910000L + r)
  N <- 240L; L <- 8L
  theta <- stats::rnorm(N)
  delta <- seq(-1.4, 1.4, length.out = L)
  X <- matrix(stats::rbinom(N * L, 1L,
    stats::plogis(outer(theta, delta, "-"))), N, L,
    dimnames = list(NULL, paste0("I", seq_len(L))))
  id <- sprintf("P%03d", seq_len(N))
  a <- tryCatch(rasch(X, id = id, n_groups = 4L), error = identity)
  b <- tryCatch(rasch(rbind(X, X), id = rep(id, 2L), n_groups = 4L),
                error = identity)
  if (inherits(a, "error") || inherits(b, "error") ||
      !isTRUE(a$est$converged) || !isTRUE(b$est$converged))
    return(list(status = "failed"))
  refused <- tryCatch({
    fit_bootstrap(b, B = 1L, workers = 1L, seed = 1L)
    FALSE
  }, rasch_refusal = function(e) TRUE, error = function(e) FALSE)
  list(
    status = "ok",
    location_diff = max(abs(a$thresholds$tau - b$thresholds$tau)),
    covariance_diff = max(abs(a$est$cov_tau - b$est$cov_tau)),
    probabilities_withheld = all(is.na(b$items$p)) &&
      all(is.na(b$items$p_adj)) && all(is.na(b$items$p_anova)) &&
      is.na(b$total_chisq_p),
    descriptive_retained = any(is.finite(b$items$chisq)) &&
      any(is.finite(b$items$fit_resid)),
    bootstrap_refused = refused)
}

z <- lapply(seq_len(N_REPS), one)
ok <- vapply(z, function(x) identical(x$status, "ok"), logical(1))
good <- z[ok]
n_ok <- length(good)
n_failed <- N_REPS - n_ok
if (!n_ok) stop("no repeated-ID validation replicate was usable")

v <- function(name) vapply(good, `[[`, logical(1), name)
d <- function(name) vapply(good, `[[`, numeric(1), name)

rows <- rbind(
  sv_row("repeated-id-item-fit-guard", "exact row duplication",
    "maximum absolute threshold difference", n_ok,
    effect = max(d("location_diff")), n_attempted = N_REPS,
    n_error = n_failed,
    notes = "threshold estimates from the original and duplicated response matrices"),
  sv_row("repeated-id-item-fit-guard", "exact row duplication",
    "maximum absolute clustered-covariance difference", n_ok,
    effect = max(d("covariance_diff")), n_attempted = N_REPS,
    n_error = n_failed,
    notes = "person-clustered threshold covariance from the original and duplicated response matrices"),
  sv_row("repeated-id-item-fit-guard", "exact row duplication",
    "item-fit probabilities withheld", n_ok,
    n_attempted = N_REPS, n_error = n_failed,
    n_withheld = sum(v("probabilities_withheld")),
    notes = "item-trait, ANOVA and total probabilities all unavailable"),
  sv_row("repeated-id-item-fit-guard", "exact row duplication",
    "descriptive fit statistics retained", n_ok,
    effect = mean(v("descriptive_retained")), n_attempted = N_REPS,
    n_error = n_failed,
    notes = "finite item chi-square and fit-residual summaries remain available"),
  sv_row("repeated-id-item-fit-guard", "exact row duplication",
    "independent-row fit bootstrap refused", n_ok,
    n_attempted = N_REPS, n_refused = sum(v("bootstrap_refused")),
    n_error = n_failed,
    notes = "refusal occurs before any bootstrap replicate is generated")
)

sv_write(rows, "repeated-id-item-fit-guard")
