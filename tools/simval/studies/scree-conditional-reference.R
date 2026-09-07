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
ONLY_DESIGNS <- trimws(strsplit(Sys.getenv("SV_DESIGNS", ""), ",",
                                fixed = TRUE)[[1L]])
ONLY_DESIGNS <- ONLY_DESIGNS[nzchar(ONLY_DESIGNS)]
OUTPUT <- Sys.getenv("SV_OUTPUT", "scree-conditional-reference")
if (!is.finite(NREP) || NREP < 2L) stop("SV_REPS must be at least 2")
if (!is.finite(INNER) || INNER < 20L) stop("SV_INNER must be at least 20")
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L
if (!is.finite(N_COMPONENTS) || N_COMPONENTS < 2L)
  stop("SV_COMPONENTS must be at least 2")
if (!grepl("^[A-Za-z0-9._-]+$", OUTPUT))
  stop("SV_OUTPUT must be a plain result-file stem")

scenarios <- expand.grid(
  design = c("dichotomous", "PCM", "RSM", "booklet", "explanatory",
             "moderately sparse PCM", "sparse PCM"),
  departure = c("null", "second dimension"),
  stringsAsFactors = FALSE)
scenarios <- scenarios[!scenarios$design %in%
                         c("moderately sparse PCM", "sparse PCM") |
                         scenarios$departure == "null", , drop = FALSE]
if (length(ONLY_DESIGNS)) {
  unknown <- setdiff(ONLY_DESIGNS, unique(scenarios$design))
  if (length(unknown))
    stop("unknown SV_DESIGNS value(s): ", paste(unknown, collapse = ", "))
  scenarios <- scenarios[scenarios$design %in% ONLY_DESIGNS, , drop = FALSE]
}

make_fit <- function(design, departure, seed) {
  sim_model <- if (design == "RSM") "RSM" else
    if (design %in% c("PCM", "moderately sparse PCM", "sparse PCM"))
      "PCM" else "dichotomous"
  fit_model <- if (sim_model == "dichotomous") "PCM" else sim_model
  second <- if (departure == "second dimension")
    list(items = sprintf("I%02d", 6:10), rho = 0.2) else NULL
  sparse <- design %in% c("moderately sparse PCM", "sparse PCM")
  severe_sparse <- design == "sparse PCM"
  d <- simulate_rasch(
    n_persons = if (severe_sparse) 35 else if (sparse) 300 else 400,
    n_items = 10,
    model = sim_model, n_categories = 4,
    difficulty = if (severe_sparse) c(-3, 3) else if (sparse)
      c(-1.5, 1.5) else c(-1.8, 1.8),
    second_dim = second, seed = seed)
  items <- sprintf("I%02d", 1:10)
  if (sparse) {
    set.seed(seed + 1L)
    xd <- as.matrix(d[items])
    missing_rate <- if (severe_sparse) .10 else .05
    xd[matrix(stats::runif(length(xd)) < missing_rate,
              nrow(xd), ncol(xd))] <- NA
    d[items] <- xd
  }
  if (design == "booklet") {
    form <- (seq_len(nrow(d)) - 1L) %% 3L
    d[form == 1L, items[8:10]] <- NA
    d[form == 2L, items[1:3]] <- NA
  }
  fit <- if (design == "explanatory") {
    predictors <- data.frame(item = items,
                             feature = seq(-1.8, 1.8, length.out = 10))
    rasch_explanatory(d, predictors, ~ feature, id = "id", items = items)
  } else {
    rasch(d, model = fit_model, id = "id", items = items)
  }
  if (!identical(fit$model, fit_model))
    stop("scenario model mismatch: requested ", fit_model,
         " but fitted ", fit$model)
  fit
}

empty_result <- function(refused = 0, nonconv = 0, error = 0,
                         reference = NULL) {
  boot_value <- function(name) {
    if (is.null(reference)) return(NA_real_)
    if (inherits(reference, "rasch_fit_bootstrap_refusal")) {
      field <- switch(name, inner_requested = "B", inner_used = "B_used",
                      inner_nonconv = "B_nonconverged",
                      inner_error = "B_errors")
      return(as.numeric(reference[[field]] %||% NA_real_))
    }
    as.numeric(attr(reference, switch(
      name, inner_requested = "n_requested", inner_used = "n_used",
      inner_nonconv = "n_nonconverged", inner_error = "n_errors")) %||%
        NA_real_)
  }
  c(ratio_mean = NA_real_, ratio_critical = NA_real_, p_upper = NA_real_,
    p_adjusted = NA_real_, rejected_pc1 = NA_real_,
    rejected_family = NA_real_,
    inner_requested = boot_value("inner_requested"),
    inner_used = boot_value("inner_used"),
    inner_nonconv = boot_value("inner_nonconv"),
    inner_error = boot_value("inner_error"),
    refused = refused, nonconv = nonconv, error = error)
}

one <- function(r, design, departure, tag) {
  fit <- tryCatch(make_fit(design, departure, 710000L + tag * 1000L + r),
                  error = identity)
  if (inherits(fit, "condition")) {
    if (inherits(fit, "rasch_refusal")) return(empty_result(refused = 1))
    return(empty_result(error = 1))
  }
  if (!isTRUE(fit$est$converged))
    return(empty_result(nonconv = 1))
  observed <- tryCatch(
    residual_pca(fit, N_COMPONENTS)$eigen_table$eigenvalue,
    error = identity)
  reference <- tryCatch(.scree_reference(
    fit, N_COMPONENTS, INNER, seed = 810000L + tag * 1000L + r),
                        error = identity)
  if (inherits(observed, "condition")) return(empty_result(error = 1))
  if (inherits(reference, "condition")) {
    if (inherits(reference, "rasch_fit_bootstrap_refusal"))
      return(empty_result(refused = 1, reference = reference))
    return(empty_result(error = 1))
  }
  if (length(observed) != N_COMPONENTS || any(!is.finite(observed)) ||
      length(reference) < 1L || !is.finite(reference[1L]) ||
      reference[1L] <= 0)
    return(empty_result(error = 1, reference = reference))
  inference <- tryCatch(
    .sim_upper_family(observed, attr(reference, "draws"), 0.05),
    error = identity)
  if (inherits(inference, "condition"))
    return(empty_result(error = 1, reference = reference))
  c(ratio_mean = observed[1L] / attr(reference, "mean")[1L],
    ratio_critical = observed[1L] / inference$critical[1L],
    p_upper = inference$p[1L], p_adjusted = inference$p_adjusted[1L],
    rejected_pc1 = as.numeric(inference$significant[1L]),
    rejected_family = as.numeric(any(inference$significant)),
    inner_requested = attr(reference, "n_requested"),
    inner_used = attr(reference, "n_used"),
    inner_nonconv = attr(reference, "n_nonconverged"),
    inner_error = attr(reference, "n_errors"),
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
mean_or_na <- function(x) if (length(x)) mean(x) else NA_real_
sd_or_na <- function(x) if (length(x) > 1L) stats::sd(x) else NA_real_
sum_or_na <- function(x) if (any(is.finite(x))) sum(x, na.rm = TRUE) else NA_real_
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
  n_ok <- sum(ok)
  rate_pc1 <- mean_or_na(rejected_pc1[ok])
  rate_family <- mean_or_na(rejected_family[ok])
  boot_requested <- sum_or_na(z[, "inner_requested"])
  boot_used <- sum_or_na(z[, "inner_used"])
  boot_nonconv <- sum_or_na(z[, "inner_nonconv"])
  boot_error <- sum_or_na(z[, "inner_error"])
  guarded_sparse <- design == "sparse PCM" && departure == "null"
  rows[[s]] <- sv_row(
    "scree conditional reference",
    paste(design, departure, sep = ", "),
    if (guarded_sparse) "conditional-reference support guard" else
      "simulated upper-tail dimensionality decision",
    n_reps = n_ok, n_attempted = NREP,
    n_refused = sum(z[, "refused"]), n_nonconv = sum(z[, "nonconv"]),
    n_error = sum(z[, "error"]),
    n_boot_attempted = boot_requested, n_boot_used = boot_used,
    n_boot_nonconv = boot_nonconv, n_boot_errors = boot_error,
    bias = if (n_ok) mean(ratio_mean[ok]) - 1 else NA_real_,
    emp_sd = sd_or_na(ratio_mean[ok]),
    type1 = if (departure == "null") rate_pc1 else NA_real_,
    familywise = if (departure == "null") rate_family else NA_real_,
    power = if (departure != "null") rate_family else NA_real_,
    notes = sprintf(
      paste("generating/fitted model %s; PC1 observed/null-mean ratio %s;",
            "PC1 observed/null-critical ratio %s; PC1 rejection %s;",
            "family rejection %s; %d components; %s/%s conditional draws",
            "used in total, with %s non-converged and %s other failures%s"),
      if (design == "RSM") "RSM" else if (design == "explanatory")
        "explanatory PCM" else if (design %in%
          c("PCM", "moderately sparse PCM", "sparse PCM"))
          "PCM" else "dichotomous PCM",
      if (n_ok) sprintf("%.4f", mean(ratio_mean[ok])) else "unavailable",
      if (n_ok) sprintf("%.4f", mean(ratio_critical[ok])) else "unavailable",
      if (is.finite(rate_pc1)) sprintf("%.3f", rate_pc1) else "unavailable",
      if (is.finite(rate_family)) sprintf("%.3f", rate_family) else "unavailable",
      N_COMPONENTS,
      if (is.finite(boot_used)) format(boot_used, scientific = FALSE) else "unknown",
      if (is.finite(boot_requested)) format(boot_requested, scientific = FALSE) else "unknown",
      if (is.finite(boot_nonconv)) format(boot_nonconv, scientific = FALSE) else "unknown",
      if (is.finite(boot_error)) format(boot_error, scientific = FALSE) else "unknown",
      if (guarded_sparse) paste0(
        "; this deliberately severe sparse-category design tests safe refusal, ",
        "not null rejection calibration") else ""))
}

out <- do.call(rbind, rows)
sv_write(out, OUTPUT)
print(out[, c("scenario", "n_reps", "bias", "emp_sd", "type1",
              "familywise", "power", "notes")],
      row.names = FALSE)
