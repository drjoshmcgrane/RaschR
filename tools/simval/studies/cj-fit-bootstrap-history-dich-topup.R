# Fresh-seed top-up for the dichotomous history-dependent judge-family result
# in cj-fit-bootstrap-history.R. Seeds 201:1000 do not overlap its first 200.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "800"))
B <- as.integer(Sys.getenv("SV_B", "199"))
CORES <- as.integer(Sys.getenv("SV_CORES", "8"))
OFFSET <- as.integer(Sys.getenv("SV_OFFSET", "200"))

blank <- function(status, message, B_used = NA_integer_,
                  B_nonconverged = NA_integer_, B_errors = NA_integer_) {
  data.frame(status = status, message = message, total = NA,
    pair_familywise = NA, object_familywise = NA, judge_familywise = NA,
    judge_marginal = NA_real_, B_used = B_used,
    B_nonconverged = B_nonconverged, B_errors = B_errors,
    stringsAsFactors = FALSE)
}

one <- function(k) {
  r <- OFFSET + k
  d <- simulate_btl(
    n_objects = 6, n_judges = 30, reps_per_pair = 24,
    dependence = list(exposure = 0.5, carry_over = 0.4),
    seed = 910000L + r)
  f <- tryCatch(btl(
    d, "object_a", "object_b", winner = "winner", judge = "judge",
    order = "order"), error = function(e) e)
  if (inherits(f, "error"))
    return(blank(if (inherits(f, "rasch_refusal")) "refused" else "error",
                 conditionMessage(f)))
  if (!isTRUE(f$converged))
    return(blank("nonconv", "the observed-data fit did not converge"))
  z <- tryCatch(suppressWarnings(fit_bootstrap(
    f, B = B, workers = 1L, seed = 920000L + r)),
    error = function(e) e)
  if (inherits(z, "error"))
    return(blank(if (inherits(z, "rasch_refusal")) "refused" else "error",
      conditionMessage(z), z$B_used %||% NA_integer_,
      z$B_nonconverged %||% NA_integer_, z$B_errors %||% NA_integer_))
  adjusted <- list(z$pairs$chisq_p_boot_adj,
                   z$objects$fit_resid_p_boot_adj,
                   z$judges$fit_resid_p_boot_adj)
  if (any(!vapply(adjusted, function(p) any(is.finite(p)), logical(1))))
    return(blank("refused",
      "at least one paired-comparison adjusted family was unavailable",
      z$B_used, z$B_nonconverged, z$B_errors))
  data.frame(
    status = "analysed", message = "",
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

attempts <- parallel::mclapply(seq_len(NREP), one,
                               mc.cores = min(CORES, NREP))
attempts <- do.call(rbind, attempts)
if (any(attempts$status == "error"))
  stop("unexpected simulation error: ", paste(unique(
    attempts$message[attempts$status == "error"]), collapse = "; "))
d <- attempts[attempts$status == "analysed", , drop = FALSE]
if (!nrow(d)) stop("no analysed replicates")
a <- list(
  attempted = nrow(attempts), refused = sum(attempts$status == "refused"),
  nonconv = sum(attempts$status == "nonconv"),
  error = sum(attempts$status == "error"),
  boot_used = sum(attempts$B_used, na.rm = TRUE),
  boot_nonconv = sum(attempts$B_nonconverged, na.rm = TRUE),
  boot_errors = sum(attempts$B_errors, na.rm = TRUE))
boot_attempted <- a$boot_used + a$boot_nonconv + a$boot_errors
note <- sprintf(paste(
  "fresh seeds %d:%d; B = %d; exposure 0.5 and carry-over 0.4 fitted jointly;",
  "%d of %d bootstrap refits did not converge; %d otherwise failed"),
  OFFSET + 1L, OFFSET + NREP, B, a$boot_nonconv, boot_attempted,
  a$boot_errors)

rows <- lapply(c("total", "pair_familywise", "object_familywise",
                 "judge_familywise", "judge_marginal"), function(metric) {
  rate <- mean(d[[metric]])
  sv_row(
    "paired-comparison history fit bootstrap top-up",
    "dichotomous, 6 objects x 30 judges",
    switch(metric,
      total = "total pairwise chi-square Type I error",
      pair_familywise = "pair chi-square joint-adjusted familywise error",
      object_familywise = "object fit-residual joint-adjusted familywise error",
      judge_familywise = "judge fit-residual joint-adjusted familywise error",
      judge_marginal = "judge fit-residual marginal Type I error"),
    n_reps = nrow(d), n_attempted = a$attempted,
    n_refused = a$refused, n_nonconv = a$nonconv, n_error = a$error,
    n_boot_attempted = boot_attempted, n_boot_used = a$boot_used,
    n_boot_nonconv = a$boot_nonconv, n_boot_errors = a$boot_errors,
    type1 = if (metric %in% c("total", "judge_marginal")) rate else NA_real_,
    familywise = if (!metric %in% c("total", "judge_marginal")) rate else NA_real_,
    mc_override = if (metric == "judge_marginal")
      list(type1 = stats::sd(d[[metric]]) / sqrt(nrow(d))) else list(),
    notes = note)
})

sv_write(do.call(rbind, rows), if (NREP < 50L)
  "cj-fit-bootstrap-history-dich-topup-smoke" else
  "cj-fit-bootstrap-history-dich-topup")
