# Supported-design top-up for score-conditional person-fit calibration in a
# four-category partial credit model. Run from the package root.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "200"))
B <- as.integer(Sys.getenv("SV_B", "199"))
CORES <- as.integer(Sys.getenv("SV_CORES", "8"))

one <- function(r) {
  d <- simulate_rasch(500, 12, model = "PCM", n_categories = 4,
                      seed = 510000L + r)
  f <- tryCatch(rasch(d, id = "id", model = "PCM"), error = function(e) e)
  if (inherits(f, "error"))
    return(data.frame(status = "refused", message = conditionMessage(f),
                      marginal = NA_real_, familywise = NA, B_used = NA_integer_,
                      B_nonconverged = NA_integer_, B_errors = NA_integer_))
  if (!isTRUE(f$est$converged))
    return(data.frame(status = "nonconv",
                      message = "the observed-data fit did not converge",
                      marginal = NA_real_, familywise = NA, B_used = NA_integer_,
                      B_nonconverged = NA_integer_, B_errors = NA_integer_))
  z <- tryCatch(suppressWarnings(fit_bootstrap(
    f, B = B, workers = 1L, seed = 610000L + r)), error = function(e) e)
  if (inherits(z, "error"))
    return(data.frame(
      status = if (inherits(z, "rasch_refusal")) "refused" else "error",
      message = conditionMessage(z), marginal = NA_real_, familywise = NA,
      B_used = z$B_used %||% NA_integer_,
      B_nonconverged = z$B_nonconverged %||% NA_integer_,
      B_errors = z$B_errors %||% NA_integer_))
  p <- z$persons$fit_resid_p_boot
  pa <- z$persons$fit_resid_p_boot_adj
  if (!any(is.finite(pa)))
    return(data.frame(
      status = "refused",
      message = "the person fit-residual adjusted family was unavailable",
      marginal = NA_real_, familywise = NA, B_used = z$B_used,
      B_nonconverged = z$B_nonconverged, B_errors = z$B_errors))
  data.frame(status = "analysed", message = "",
             marginal = mean(p < .05, na.rm = TRUE),
             familywise = any(pa < .05, na.rm = TRUE), B_used = z$B_used,
             B_nonconverged = z$B_nonconverged, B_errors = z$B_errors)
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

rows <- rbind(
  sv_row("person fit bootstrap", "PCM, 500 persons x 12 four-category items",
    "person fit-residual marginal Type I error", n_reps = nr,
    n_attempted = acct$attempted, n_refused = acct$refused,
    n_nonconv = acct$nonconv, type1 = mean(a$marginal),
    mc_override = list(type1 = stats::sd(a$marginal) / sqrt(nr)),
    notes = paste(sprintf("B = %d; score-conditional null", B), inner_note,
                  sep = "; ")),
  sv_row("person fit bootstrap", "PCM, 500 persons x 12 four-category items",
    "person fit-residual joint-adjusted familywise error", n_reps = nr,
    n_attempted = acct$attempted, n_refused = acct$refused,
    n_nonconv = acct$nonconv, familywise = mean(a$familywise),
    notes = paste(sprintf("B = %d; maximum-statistic adjustment", B),
                  inner_note, sep = "; ")))

sv_write(rows, if (NREP < 30L)
  "person-fit-bootstrap-pcm-topup-smoke" else "person-fit-bootstrap-pcm-topup")
