# Top-up for fitted-design goodness-of-fit calibration in four-category
# comparative judgement. Run from the package root.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "200"))
B <- as.integer(Sys.getenv("SV_B", "199"))
CORES <- as.integer(Sys.getenv("SV_CORES", "4"))

one <- function(r) {
  d <- simulate_btl(6, 20, reps_per_pair = 10, model = "polytomous",
                    n_categories = 4, seed = 710000L + r)
  f <- tryCatch(btl(d, "object_a", "object_b", response = "response",
                    judge = "judge"), error = function(e) e)
  blank <- function(status, message, B_used = NA_integer_,
                    B_nonconverged = NA_integer_, B_errors = NA_integer_) data.frame(
    status = status, message = message, total = NA,
    pair_familywise = NA, object_familywise = NA,
    judge_familywise = NA, judge_marginal = NA_real_, B_used = B_used,
    B_nonconverged = B_nonconverged, B_errors = B_errors)
  if (inherits(f, "error")) return(blank("refused", conditionMessage(f)))
  if (!isTRUE(f$converged))
    return(blank("nonconv", "the observed-data fit did not converge"))
  z <- tryCatch(suppressWarnings(fit_bootstrap(
    f, B = B, workers = 1L, seed = 810000L + r)), error = function(e) e)
  if (inherits(z, "error"))
    return(blank(if (inherits(z, "rasch_refusal")) "refused" else "error",
                 conditionMessage(z), z$B_used %||% NA_integer_,
                 z$B_nonconverged %||% NA_integer_,
                 z$B_errors %||% NA_integer_))
  adjusted <- list(z$pairs$chisq_p_boot_adj,
                   z$objects$fit_resid_p_boot_adj,
                   z$judges$fit_resid_p_boot_adj)
  if (any(!vapply(adjusted, function(p) any(is.finite(p)), logical(1))))
    return(blank("refused",
      "at least one paired-comparison adjusted family was unavailable",
      z$B_used, z$B_nonconverged, z$B_errors))
  data.frame(
    status = "analysed", message = "", total = z$total$chisq_p_boot < .05,
    pair_familywise = any(z$pairs$chisq_p_boot_adj < .05, na.rm = TRUE),
    object_familywise = any(z$objects$fit_resid_p_boot_adj < .05, na.rm = TRUE),
    judge_familywise = any(z$judges$fit_resid_p_boot_adj < .05, na.rm = TRUE),
    judge_marginal = mean(z$judges$fit_resid_p_boot < .05, na.rm = TRUE),
    B_used = z$B_used, B_nonconverged = z$B_nonconverged,
    B_errors = z$B_errors)
}

z <- parallel::mclapply(seq_len(NREP), one,
                        mc.cores = min(CORES, NREP))
d <- do.call(rbind, z)
if (any(d$status == "error"))
  stop("unexpected simulation error: ",
       paste(unique(d$message[d$status == "error"]), collapse = "; "))
if (any(d$status != "analysed"))
  print(unique(d[d$status != "analysed", c("status", "message")]))
a <- d[d$status == "analysed", , drop = FALSE]
nr <- nrow(a)
acct <- list(attempted = nrow(d), refused = sum(d$status == "refused"),
             nonconv = sum(d$status == "nonconv"))
inner_note <- sprintf("%d of %d bootstrap refits did not converge; %d otherwise failed",
  sum(d$B_nonconverged, na.rm = TRUE),
  sum(d$B_used + d$B_nonconverged + d$B_errors, na.rm = TRUE),
  sum(d$B_errors, na.rm = TRUE))

make_row <- function(metric) sv_row(
  "paired-comparison fit bootstrap",
  "polytomous, 6 objects x 20 judges, four categories",
  switch(metric,
    total = "total pairwise chi-square Type I error",
    pair_familywise = "pair chi-square joint-adjusted familywise error",
    object_familywise = "object fit-residual joint-adjusted familywise error",
    judge_familywise = "judge fit-residual joint-adjusted familywise error",
    judge_marginal = "judge fit-residual marginal Type I error"),
  n_reps = nr, n_attempted = acct$attempted, n_refused = acct$refused,
  n_nonconv = acct$nonconv,
  type1 = if (metric %in% c("total", "judge_marginal"))
    mean(a[[metric]]) else NA_real_,
  familywise = if (!metric %in% c("total", "judge_marginal"))
    mean(a[[metric]]) else NA_real_,
  mc_override = if (metric == "judge_marginal")
    list(type1 = stats::sd(a[[metric]]) / sqrt(nr)) else list(),
  notes = paste(sprintf("B = %d; fitted four-category scale retained", B),
                inner_note, sep = "; "))

rows <- do.call(rbind, lapply(
  c("total", "pair_familywise", "object_familywise",
    "judge_familywise", "judge_marginal"), make_row))
sv_write(rows, if (NREP < 30L)
  "cj-fit-bootstrap-polytomous-topup-smoke" else
    "cj-fit-bootstrap-polytomous-topup")
