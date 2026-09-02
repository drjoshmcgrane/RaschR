#!/usr/bin/env Rscript
# Validation of the score-conditional parallel-analysis reference used by
# plot_scree(). Null cells should centre the observed/reference first-
# eigenvalue ratio on one; planted second dimensions should raise it.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "20"))
INNER <- as.integer(Sys.getenv("SV_INNER", "15"))
CORES <- as.integer(Sys.getenv("SV_CORES", "4"))
if (!is.finite(NREP) || NREP < 2L) stop("SV_REPS must be at least 2")
if (!is.finite(INNER) || INNER < 2L) stop("SV_INNER must be at least 2")
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L

scenarios <- expand.grid(
  design = c("dichotomous", "PCM", "RSM", "booklet", "explanatory"),
  departure = c("null", "second dimension"),
  stringsAsFactors = FALSE)

make_fit <- function(design, departure, seed) {
  sim_model <- if (design %in% c("PCM", "RSM")) design else "dichotomous"
  fit_model <- if (sim_model == "dichotomous") "PCM" else sim_model
  second <- if (departure == "second dimension")
    list(items = sprintf("I%02d", 6:10), rho = 0.2) else NULL
  d <- simulate_rasch(
    n_persons = 400, n_items = 10, model = sim_model, n_categories = 4,
    difficulty = c(-1.8, 1.8), second_dim = second, seed = seed)
  items <- sprintf("I%02d", 1:10)
  if (design == "booklet") {
    form <- (seq_len(nrow(d)) - 1L) %% 3L
    d[form == 1L, items[8:10]] <- NA
    d[form == 2L, items[1:3]] <- NA
  }
  if (design == "explanatory") {
    predictors <- data.frame(item = items,
                             feature = seq(-1.8, 1.8, length.out = 10))
    rasch_explanatory(d, predictors, ~ feature, id = "id", items = items)
  } else {
    rasch(d, model = fit_model, id = "id", items = items)
  }
}

one <- function(r, design, departure, tag) {
  fit <- tryCatch(make_fit(design, departure, 710000L + tag * 1000L + r),
                  error = identity)
  if (inherits(fit, "condition"))
    return(c(ratio = NA_real_, refused = 1, nonconv = 0, error = 0))
  if (!isTRUE(fit$est$converged))
    return(c(ratio = NA_real_, refused = 0, nonconv = 1, error = 0))
  observed <- tryCatch(residual_pca(fit, 1)$first_eigen, error = identity)
  reference <- tryCatch(.scree_reference(fit, 1, INNER)[1L], error = identity)
  if (inherits(observed, "condition") || inherits(reference, "condition") ||
      !is.finite(observed) || !is.finite(reference) || reference <= 0)
    return(c(ratio = NA_real_, refused = 0, nonconv = 0, error = 1))
  c(ratio = observed / reference, refused = 0, nonconv = 0, error = 0)
}

run_cell <- function(design, departure, tag) {
  f <- function(r) one(r, design, departure, tag)
  z <- if (CORES > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(NREP), f, mc.cores = CORES,
                       mc.set.seed = FALSE) else lapply(seq_len(NREP), f)
  do.call(rbind, z)
}

rows <- vector("list", nrow(scenarios))
for (s in seq_len(nrow(scenarios))) {
  design <- scenarios$design[s]
  departure <- scenarios$departure[s]
  z <- run_cell(design, departure, s)
  ratio <- z[, "ratio"]
  ok <- is.finite(ratio)
  rows[[s]] <- sv_row(
    "scree conditional reference",
    paste(design, departure, sep = ", "),
    "observed/reference first eigenvalue ratio",
    n_reps = sum(ok), n_attempted = NREP,
    n_refused = sum(z[, "refused"]), n_nonconv = sum(z[, "nonconv"]),
    n_error = sum(z[, "error"]),
    bias = mean(ratio[ok]) - 1, emp_sd = stats::sd(ratio[ok]),
    notes = sprintf(
      "mean ratio %.4f; median %.4f; P(observed > mean reference) %.3f; %d conditional draws per outer replicate",
      mean(ratio[ok]), stats::median(ratio[ok]), mean(ratio[ok] > 1), INNER))
}

out <- do.call(rbind, rows)
sv_write(out, "scree-conditional-reference")
print(out[, c("scenario", "n_reps", "bias", "emp_sd", "notes")],
      row.names = FALSE)
