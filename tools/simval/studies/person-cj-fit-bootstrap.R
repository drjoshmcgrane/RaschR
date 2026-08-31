# Null calibration and planted-misfit detection for the person and
# paired-comparison extensions of fit_bootstrap(). Run from the package root.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "200"))
B <- as.integer(Sys.getenv("SV_B", "399"))
CORES <- as.integer(Sys.getenv("SV_CORES", "8"))

person_failure <- function(r, model, careless, status, message,
                           B_used = NA_integer_, B_nonconverged = NA_integer_,
                           B_errors = NA_integer_) data.frame(
  rep = r, model = model, careless = careless, status = status,
  message = message, marginal_clean = NA_real_, familywise_clean = NA,
  power_raw = NA_real_, power_adjusted = NA_real_, any_adjusted = NA,
  B_used = B_used, B_nonconverged = B_nonconverged, B_errors = B_errors)

person_one <- function(r, model = "dichotomous", careless = 0) {
  d <- simulate_rasch(240, 12, model = model,
    n_categories = if (model == "dichotomous") 2 else 4,
    careless = careless, seed = 110000L + r + 1000L * (model != "dichotomous") +
      10000L * (careless > 0))
  truth <- attr(d, "truth")
  f <- tryCatch(rasch(d, id = "id", model = if (model == "dichotomous")
    "PCM" else model), error = function(e) e)
  if (inherits(f, "error"))
    return(person_failure(r, model, careless, "refused", conditionMessage(f)))
  if (!isTRUE(f$est$converged))
    return(person_failure(r, model, careless, "nonconv",
                          "the observed-data fit did not converge"))
  z <- tryCatch(suppressWarnings(fit_bootstrap(
    f, B = B, workers = 1L, seed = 210000L + r)), error = function(e) e)
  if (inherits(z, "error"))
    return(person_failure(r, model, careless,
      if (inherits(z, "rasch_refusal")) "refused" else "error",
      conditionMessage(z), z$B_used %||% NA_integer_,
      z$B_nonconverged %||% NA_integer_, z$B_errors %||% NA_integer_))
  if (is.null(z$persons))
    return(person_failure(r, model, careless, "error",
                          "the conditional bootstrap omitted person results"))
  bad <- seq_len(nrow(d)) %in% (truth$careless_idx %||% integer(0))
  p <- z$persons$fit_resid_p_boot
  pa <- z$persons$fit_resid_p_boot_adj
  if (!any(is.finite(pa)))
    return(person_failure(
      r, model, careless, "refused",
      "the person fit-residual adjusted family was unavailable",
      z$B_used, z$B_nonconverged, z$B_errors))
  data.frame(rep = r, model = model, careless = careless,
    status = "analysed", message = "",
    marginal_clean = mean(p[!bad] < .05, na.rm = TRUE),
    familywise_clean = any(pa[!bad] < .05, na.rm = TRUE),
    power_raw = if (any(bad)) mean(p[bad] < .05, na.rm = TRUE) else NA_real_,
    power_adjusted = if (any(bad)) mean(pa[bad] < .05, na.rm = TRUE) else NA_real_,
    any_adjusted = any(pa < .05, na.rm = TRUE),
    B_used = z$B_used, B_nonconverged = z$B_nonconverged,
    B_errors = z$B_errors)
}

cj_failure <- function(r, model, erratic, status, message,
                       B_used = NA_integer_, B_nonconverged = NA_integer_,
                       B_errors = NA_integer_) data.frame(
  rep = r, model = model, erratic = erratic, status = status,
  message = message, total = NA, pair_familywise = NA,
  object_familywise = NA, judge_familywise_clean = NA,
  judge_marginal_clean = NA_real_, judge_power_raw = NA_real_,
  judge_power_adjusted = NA_real_, B_used = B_used,
  B_nonconverged = B_nonconverged, B_errors = B_errors)

cj_one <- function(r, model = "dichotomous", erratic = 0) {
  d <- simulate_btl(6, 20, reps_per_pair = 10, model = model,
    n_categories = if (model == "dichotomous") 2 else 4,
    erratic_judges = erratic,
    seed = 310000L + r + 1000L * (model != "dichotomous") +
      10000L * (erratic > 0))
  truth <- attr(d, "truth")
  f <- tryCatch(if (model == "dichotomous")
    btl(d, "object_a", "object_b", winner = "winner", judge = "judge") else
    btl(d, "object_a", "object_b", response = "response", judge = "judge"),
    error = function(e) e)
  if (inherits(f, "error"))
    return(cj_failure(r, model, erratic, "refused", conditionMessage(f)))
  if (!isTRUE(f$converged))
    return(cj_failure(r, model, erratic, "nonconv",
                      "the observed-data fit did not converge"))
  z <- tryCatch(suppressWarnings(fit_bootstrap(
    f, B = B, workers = 1L, seed = 410000L + r)), error = function(e) e)
  if (inherits(z, "error"))
    return(cj_failure(r, model, erratic,
      if (inherits(z, "rasch_refusal")) "refused" else "error",
      conditionMessage(z), z$B_used %||% NA_integer_,
      z$B_nonconverged %||% NA_integer_, z$B_errors %||% NA_integer_))
  bad <- z$judges$judge %in% (truth$erratic %||% character(0))
  jp <- z$judges$fit_resid_p_boot
  ja <- z$judges$fit_resid_p_boot_adj
  adjusted <- list(z$pairs$chisq_p_boot_adj,
                   z$objects$fit_resid_p_boot_adj, ja)
  if (any(!vapply(adjusted, function(p) any(is.finite(p)), logical(1))))
    return(cj_failure(
      r, model, erratic, "refused",
      "at least one paired-comparison adjusted family was unavailable",
      z$B_used, z$B_nonconverged, z$B_errors))
  data.frame(rep = r, model = model, erratic = erratic,
    status = "analysed", message = "",
    total = z$total$chisq_p_boot < .05,
    pair_familywise = any(z$pairs$chisq_p_boot_adj < .05, na.rm = TRUE),
    object_familywise = any(z$objects$fit_resid_p_boot_adj < .05, na.rm = TRUE),
    judge_familywise_clean = any(ja[!bad] < .05, na.rm = TRUE),
    judge_marginal_clean = mean(jp[!bad] < .05, na.rm = TRUE),
    judge_power_raw = if (any(bad)) mean(jp[bad] < .05, na.rm = TRUE) else NA_real_,
    judge_power_adjusted = if (any(bad)) mean(ja[bad] < .05, na.rm = TRUE) else NA_real_,
    B_used = z$B_used, B_nonconverged = z$B_nonconverged,
    B_errors = z$B_errors)
}

run_outer <- function(fun, ...) {
  z <- parallel::mclapply(seq_len(NREP), fun, ...,
                          mc.cores = min(CORES, NREP))
  d <- do.call(rbind, z)
  if (any(d$status == "error"))
    stop("unexpected simulation error: ",
         paste(unique(d$message[d$status == "error"]), collapse = "; "))
  if (any(d$status != "analysed"))
    print(unique(d[d$status != "analysed", c("status", "message")]))
  d
}

analysed <- function(d) d[d$status == "analysed", , drop = FALSE]
accounting <- function(d) list(
  attempted = nrow(d), refused = sum(d$status == "refused"),
  nonconv = sum(d$status == "nonconv"), error = sum(d$status == "error"))

cluster_mcse <- function(d, name) {
  z <- d[[name]]
  stats::sd(z, na.rm = TRUE) / sqrt(sum(is.finite(z)))
}

inner_counts <- function(d) {
  used <- sum(d$B_used, na.rm = TRUE)
  nonconv <- sum(d$B_nonconverged, na.rm = TRUE)
  errors <- sum(d$B_errors, na.rm = TRUE)
  list(attempted = used + nonconv + errors, used = used,
       nonconv = nonconv, errors = errors)
}

inner_accounting <- function(d) {
  z <- inner_counts(d)
  sprintf("%d of %d bootstrap refits did not converge; %d otherwise failed",
          z$nonconv, z$attempted, z$errors)
}

rows <- list()
for (model in c("dichotomous", "PCM")) {
  attempts <- run_outer(person_one, model = model, careless = 0)
  ac <- accounting(attempts); ic <- inner_counts(attempts)
  d <- analysed(attempts)
  nr <- nrow(d)
  for (metric in c("marginal_clean", "familywise_clean"))
    rows[[length(rows) + 1L]] <- sv_row(
      "person fit bootstrap", sprintf("%s, 240 persons x 12 items", model),
      if (metric == "marginal_clean")
        "person fit-residual marginal Type I error" else
        "person fit-residual joint-adjusted familywise error",
      n_reps = nr, n_attempted = ac$attempted,
      n_refused = ac$refused, n_nonconv = ac$nonconv,
      n_error = ac$error, n_boot_attempted = ic$attempted,
      n_boot_used = ic$used, n_boot_nonconv = ic$nonconv,
      n_boot_errors = ic$errors,
      type1 = if (metric == "marginal_clean") mean(d[[metric]]) else NA_real_,
      familywise = if (metric == "familywise_clean") mean(d[[metric]]) else NA_real_,
      mc_override = if (metric == "marginal_clean")
        list(type1 = cluster_mcse(d, metric)) else list(),
      notes = paste(sprintf("B = %d; fit residual, score-conditional maxT", B),
                    inner_accounting(attempts), sep = "; "))
}

attempts <- run_outer(person_one, model = "dichotomous", careless = .10)
ac <- accounting(attempts); ic <- inner_counts(attempts)
d <- analysed(attempts)
nr <- nrow(d)
for (metric in c("familywise_clean", "power_raw", "power_adjusted"))
  rows[[length(rows) + 1L]] <- sv_row(
    "person fit bootstrap", "10% random responders, 240 persons x 12 items",
    switch(metric, familywise_clean = "clean-person familywise error",
      power_raw = "careless-person marginal detection",
      power_adjusted = "careless-person joint-adjusted detection"),
    n_reps = nr, n_attempted = ac$attempted,
    n_refused = ac$refused, n_nonconv = ac$nonconv,
    n_error = ac$error, n_boot_attempted = ic$attempted,
    n_boot_used = ic$used, n_boot_nonconv = ic$nonconv,
    n_boot_errors = ic$errors,
    familywise = if (metric == "familywise_clean") mean(d[[metric]]) else NA_real_,
    power = if (metric != "familywise_clean") mean(d[[metric]]) else NA_real_,
    mc_override = if (metric != "familywise_clean")
      list(power = cluster_mcse(d, metric)) else list(),
    notes = paste(sprintf("B = %d", B), inner_accounting(attempts), sep = "; "))

for (model in c("dichotomous", "polytomous")) {
  attempts <- run_outer(cj_one, model = model, erratic = 0)
  ac <- accounting(attempts); ic <- inner_counts(attempts)
  d <- analysed(attempts)
  nr <- nrow(d)
  for (metric in c("total", "pair_familywise", "object_familywise",
                   "judge_familywise_clean", "judge_marginal_clean"))
    rows[[length(rows) + 1L]] <- sv_row(
      "paired-comparison fit bootstrap",
      sprintf("%s, 6 objects x 20 judges", model),
      switch(metric, total = "total pairwise chi-square Type I error",
        pair_familywise = "pair chi-square joint-adjusted familywise error",
        object_familywise = "object fit-residual joint-adjusted familywise error",
        judge_familywise_clean = "judge fit-residual joint-adjusted familywise error",
        judge_marginal_clean = "judge fit-residual marginal Type I error"),
      n_reps = nr, n_attempted = ac$attempted,
      n_refused = ac$refused, n_nonconv = ac$nonconv,
      n_error = ac$error, n_boot_attempted = ic$attempted,
      n_boot_used = ic$used, n_boot_nonconv = ic$nonconv,
      n_boot_errors = ic$errors,
      type1 = if (metric %in% c("total", "judge_marginal_clean"))
        mean(d[[metric]]) else NA_real_,
      familywise = if (!metric %in% c("total", "judge_marginal_clean"))
        mean(d[[metric]]) else NA_real_,
      mc_override = if (metric == "judge_marginal_clean")
        list(type1 = cluster_mcse(d, metric)) else list(),
      notes = paste(sprintf("B = %d", B), inner_accounting(attempts), sep = "; "))
}

attempts <- run_outer(cj_one, model = "dichotomous", erratic = .20)
ac <- accounting(attempts); ic <- inner_counts(attempts)
d <- analysed(attempts)
nr <- nrow(d)
for (metric in c("judge_familywise_clean", "judge_power_raw",
                 "judge_power_adjusted"))
  rows[[length(rows) + 1L]] <- sv_row(
    "paired-comparison fit bootstrap", "20% erratic judges, 6 objects x 20 judges",
    switch(metric, judge_familywise_clean = "clean-judge familywise error",
      judge_power_raw = "erratic-judge marginal detection",
      judge_power_adjusted = "erratic-judge joint-adjusted detection"),
    n_reps = nr, n_attempted = ac$attempted,
    n_refused = ac$refused, n_nonconv = ac$nonconv,
    n_error = ac$error, n_boot_attempted = ic$attempted,
    n_boot_used = ic$used, n_boot_nonconv = ic$nonconv,
    n_boot_errors = ic$errors,
    familywise = if (metric == "judge_familywise_clean") mean(d[[metric]]) else NA_real_,
    power = if (metric != "judge_familywise_clean") mean(d[[metric]]) else NA_real_,
    mc_override = if (metric != "judge_familywise_clean")
      list(power = cluster_mcse(d, metric)) else list(),
    notes = paste(sprintf("B = %d", B), inner_accounting(attempts), sep = "; "))

sv_write(do.call(rbind, rows), if (NREP < 30L)
  "person-cj-fit-bootstrap-smoke" else "person-cj-fit-bootstrap")
