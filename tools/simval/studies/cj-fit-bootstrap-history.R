# Null calibration of the paired-comparison fitted-design bootstrap when the
# fitted model contains exposure and carry-over history terms. Run from the
# package root.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "200"))
B <- as.integer(Sys.getenv("SV_B", "199"))
CORES <- as.integer(Sys.getenv("SV_CORES", "4"))

failure <- function(model, status, message, B_used = NA_integer_,
                    B_nonconverged = NA_integer_, B_errors = NA_integer_) {
  data.frame(
    model = model, status = status, message = message,
    total = NA, pair_familywise = NA, object_familywise = NA,
    judge_familywise = NA, judge_marginal = NA_real_,
    B_used = B_used, B_nonconverged = B_nonconverged,
    B_errors = B_errors, stringsAsFactors = FALSE)
}

one <- function(r, model) {
  polytomous <- identical(model, "polytomous")
  d <- simulate_btl(
    n_objects = 6, n_judges = 30, reps_per_pair = 24,
    model = model, n_categories = if (polytomous) 4 else 2,
    dependence = list(exposure = 0.5, carry_over = 0.4),
    seed = 910000L + 10000L * polytomous + r)
  f <- tryCatch(
    if (polytomous)
      btl(d, "object_a", "object_b", response = "response",
          judge = "judge", order = "order")
    else
      btl(d, "object_a", "object_b", winner = "winner",
          judge = "judge", order = "order"),
    error = function(e) e)
  if (inherits(f, "error"))
    return(failure(model,
      if (inherits(f, "rasch_refusal")) "refused" else "error",
      conditionMessage(f)))
  if (!isTRUE(f$converged))
    return(failure(model, "nonconv",
                   "the observed-data fit did not converge"))

  z <- tryCatch(suppressWarnings(fit_bootstrap(
    f, B = B, workers = 1L, seed = 920000L + 10000L * polytomous + r)),
    error = function(e) e)
  if (inherits(z, "error"))
    return(failure(
      model, if (inherits(z, "rasch_refusal")) "refused" else "error",
      conditionMessage(z), z$B_used %||% NA_integer_,
      z$B_nonconverged %||% NA_integer_, z$B_errors %||% NA_integer_))

  adjusted <- list(z$pairs$chisq_p_boot_adj,
                   z$objects$fit_resid_p_boot_adj,
                   z$judges$fit_resid_p_boot_adj)
  if (any(!vapply(adjusted, function(p) any(is.finite(p)), logical(1))))
    return(failure(
      model, "refused",
      "at least one paired-comparison adjusted family was unavailable",
      z$B_used, z$B_nonconverged, z$B_errors))

  data.frame(
    model = model, status = "analysed", message = "",
    total = z$total$chisq_p_boot < .05,
    pair_familywise = any(z$pairs$chisq_p_boot_adj < .05, na.rm = TRUE),
    object_familywise = any(z$objects$fit_resid_p_boot_adj < .05,
                            na.rm = TRUE),
    judge_familywise = any(z$judges$fit_resid_p_boot_adj < .05,
                           na.rm = TRUE),
    judge_marginal = mean(z$judges$fit_resid_p_boot < .05, na.rm = TRUE),
    B_used = z$B_used, B_nonconverged = z$B_nonconverged,
    B_errors = z$B_errors, stringsAsFactors = FALSE)
}

accounting <- function(d) list(
  attempted = nrow(d), refused = sum(d$status == "refused"),
  nonconv = sum(d$status == "nonconv"), error = sum(d$status == "error"),
  boot_used = sum(d$B_used, na.rm = TRUE),
  boot_nonconv = sum(d$B_nonconverged, na.rm = TRUE),
  boot_errors = sum(d$B_errors, na.rm = TRUE))

rows <- list()
for (model in c("dichotomous", "polytomous")) {
  attempts <- parallel::mclapply(
    seq_len(NREP), one, model = model, mc.cores = min(CORES, NREP))
  attempts <- do.call(rbind, attempts)
  if (any(attempts$status == "error"))
    stop("unexpected simulation error: ", paste(unique(
      attempts$message[attempts$status == "error"]), collapse = "; "))
  a <- accounting(attempts)
  d <- attempts[attempts$status == "analysed", , drop = FALSE]
  nr <- nrow(d)
  if (!nr) stop("no analysed replicates for ", model)
  boot_attempted <- a$boot_used + a$boot_nonconv + a$boot_errors
  note <- sprintf(paste(
    "B = %d; exposure 0.5 and carry-over 0.4 fitted jointly;",
    "%d of %d bootstrap refits did not converge; %d otherwise failed"),
    B, a$boot_nonconv, boot_attempted, a$boot_errors)
  for (metric in c("total", "pair_familywise", "object_familywise",
                   "judge_familywise", "judge_marginal")) {
    marginal <- identical(metric, "judge_marginal")
    rate <- mean(d[[metric]])
    rows[[length(rows) + 1L]] <- sv_row(
      "paired-comparison history fit bootstrap",
      sprintf("%s, 6 objects x 30 judges", model),
      switch(metric,
        total = "total pairwise chi-square Type I error",
        pair_familywise = "pair chi-square joint-adjusted familywise error",
        object_familywise = "object fit-residual joint-adjusted familywise error",
        judge_familywise = "judge fit-residual joint-adjusted familywise error",
        judge_marginal = "judge fit-residual marginal Type I error"),
      n_reps = nr, n_attempted = a$attempted,
      n_refused = a$refused, n_nonconv = a$nonconv, n_error = a$error,
      n_boot_attempted = boot_attempted, n_boot_used = a$boot_used,
      n_boot_nonconv = a$boot_nonconv, n_boot_errors = a$boot_errors,
      type1 = if (metric %in% c("total", "judge_marginal")) rate else NA_real_,
      familywise = if (!metric %in% c("total", "judge_marginal")) rate else NA_real_,
      mc_override = if (marginal)
        list(type1 = stats::sd(d[[metric]]) / sqrt(nr)) else list(),
      notes = note)
  }
}

sv_write(do.call(rbind, rows), if (NREP < 50L)
  "cj-fit-bootstrap-history-smoke" else "cj-fit-bootstrap-history")
