#!/usr/bin/env Rscript
# Validation of the score-conditional parallel-analysis reference used by
# plot_scree(). The fitted-model maximum-statistic reference should control the
# familywise null rejection rate across the displayed components, while a
# planted second dimension should clear it with useful power.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "200"))
INNER <- as.integer(Sys.getenv("SV_INNER", "50"))
CORES <- as.integer(Sys.getenv("SV_CORES", "4"))
N_COMPONENTS <- as.integer(Sys.getenv("SV_COMPONENTS", "10"))
if (!is.finite(NREP) || NREP < 2L) stop("SV_REPS must be at least 2")
if (!is.finite(INNER) || INNER < 20L) stop("SV_INNER must be at least 20")
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L
if (!is.finite(N_COMPONENTS) || N_COMPONENTS < 2L)
  stop("SV_COMPONENTS must be at least 2")

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
    return(c(ratio_mean = NA_real_, ratio_critical = NA_real_, p_upper = NA_real_,
             p_adjusted = NA_real_, rejected_pc1 = NA_real_,
             rejected_family = NA_real_, refused = 1, nonconv = 0, error = 0))
  if (!isTRUE(fit$est$converged))
    return(c(ratio_mean = NA_real_, ratio_critical = NA_real_, p_upper = NA_real_,
             p_adjusted = NA_real_, rejected_pc1 = NA_real_,
             rejected_family = NA_real_, refused = 0, nonconv = 1, error = 0))
  observed <- tryCatch(
    residual_pca(fit, N_COMPONENTS)$eigen_table$eigenvalue,
    error = identity)
  set.seed(810000L + tag * 1000L + r)
  reference <- tryCatch(.scree_reference(fit, N_COMPONENTS, INNER),
                        error = identity)
  if (inherits(observed, "condition") || inherits(reference, "condition") ||
      length(observed) != N_COMPONENTS || any(!is.finite(observed)) ||
      !is.finite(reference[1L]) ||
      reference[1L] <= 0)
    return(c(ratio_mean = NA_real_, ratio_critical = NA_real_, p_upper = NA_real_,
             p_adjusted = NA_real_, rejected_pc1 = NA_real_,
             rejected_family = NA_real_, refused = 0, nonconv = 0, error = 1))
  inference <- tryCatch(
    .sim_upper_family(observed, attr(reference, "draws"), 0.05),
    error = identity)
  if (inherits(inference, "condition"))
    return(c(ratio_mean = NA_real_, ratio_critical = NA_real_, p_upper = NA_real_,
             p_adjusted = NA_real_, rejected_pc1 = NA_real_,
             rejected_family = NA_real_, refused = 0, nonconv = 0, error = 1))
  c(ratio_mean = observed[1L] / attr(reference, "mean")[1L],
    ratio_critical = observed[1L] / inference$critical[1L],
    p_upper = inference$p[1L], p_adjusted = inference$p_adjusted[1L],
    rejected_pc1 = as.numeric(inference$significant[1L]),
    rejected_family = as.numeric(any(inference$significant)),
    refused = 0, nonconv = 0, error = 0)
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
  ratio_mean <- z[, "ratio_mean"]
  ratio_critical <- z[, "ratio_critical"]
  p_upper <- z[, "p_upper"]
  p_adjusted <- z[, "p_adjusted"]
  rejected_pc1 <- z[, "rejected_pc1"]
  rejected_family <- z[, "rejected_family"]
  ok <- is.finite(ratio_mean) & is.finite(ratio_critical) &
    is.finite(p_upper) & is.finite(p_adjusted)
  rate_pc1 <- mean(rejected_pc1[ok])
  rate_family <- mean(rejected_family[ok])
  rows[[s]] <- sv_row(
    "scree conditional reference",
    paste(design, departure, sep = ", "),
    "simulated upper-tail dimensionality decision",
    n_reps = sum(ok), n_attempted = NREP,
    n_refused = sum(z[, "refused"]), n_nonconv = sum(z[, "nonconv"]),
    n_error = sum(z[, "error"]),
    bias = mean(ratio_mean[ok]) - 1, emp_sd = stats::sd(ratio_mean[ok]),
    type1 = if (departure == "null") rate_pc1 else NA_real_,
    familywise = if (departure == "null") rate_family else NA_real_,
    power = if (departure != "null") rate_family else NA_real_,
    notes = sprintf(
      paste("PC1 observed/null-mean ratio %.4f; PC1 observed/null-critical",
            "ratio %.4f; PC1 rejection %.3f; family rejection %.3f;",
            "%d components and %d conditional draws per outer replicate;",
            "largest attainable level %.4f"),
      mean(ratio_mean[ok]), mean(ratio_critical[ok]), rate_pc1, rate_family,
      N_COMPONENTS, INNER, floor(0.05 * (INNER + 1)) / (INNER + 1)))
}

out <- do.call(rbind, rows)
sv_write(out, "scree-conditional-reference")
print(out[, c("scenario", "n_reps", "bias", "emp_sd", "type1",
              "familywise", "power", "notes")],
      row.names = FALSE)
