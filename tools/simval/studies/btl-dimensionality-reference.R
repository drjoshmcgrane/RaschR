#!/usr/bin/env Rscript
# Validation of the finite-simulation upper-tail decision used by
# btl_dimensionality(). The ordinary and Extended Frames paths share the same
# one-statistic rule, so each null design is checked at the minimum and a more
# stable reference size. The planted second-attribute cell checks that replacing
# an interpolated percentile with the finite probability does not remove useful
# power.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "1000"))
CORES <- as.integer(Sys.getenv("SV_CORES", "4"))
REFS <- as.numeric(strsplit(Sys.getenv("SV_INNER", "20,200"), ",",
                            fixed = TRUE)[[1L]])
if (!is.finite(NREP) || NREP < 2L) stop("SV_REPS must be at least 2")
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L
if (!length(REFS) || any(!is.finite(REFS) | REFS < 20L | REFS != floor(REFS)))
  stop("SV_INNER must contain comma-separated whole numbers of at least 20")
REFS <- as.integer(REFS)

scenarios <- rbind(
  expand.grid(model = "BTL", departure = c("null", "second attribute"),
              reference_reps = REFS, stringsAsFactors = FALSE),
  expand.grid(model = "BTL-EFRM", departure = "null",
              reference_reps = REFS, stringsAsFactors = FALSE))

one <- function(r, model, departure, reference_reps, tag) {
  seed <- 930000L + tag * 10000L + r
  fit <- tryCatch({
    if (model == "BTL") {
      d <- simulate_btl(
        n_objects = 8, n_judges = 24, reps_per_pair = 70,
        object_sd = if (departure == "second attribute") 2.4 else 1,
        second_attribute = if (departure == "second attribute")
          list(rho = 0.3) else NULL,
        seed = seed)
      btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
    } else {
      d <- simulate_btl_efrm(
        n_objects_per_set = 6, n_sets = 2, n_panels = 2,
        n_judges_per_panel = 10, reps_within = 20, reps_cross = 20,
        seed = seed)
      btl_efrm(d, "object_a", "object_b", winner = "winner",
               judge = "judge", panels = "panel",
               object_sets = attr(d, "truth")$object_sets,
               se_method = "conditional")
    }
  }, error = identity)
  if (inherits(fit, "condition"))
    return(c(reject = NA_real_, p = NA_real_, refused = 1,
             nonconv = 0, error = 0))
  converged <- fit$converged
  if (!isTRUE(converged))
    return(c(reject = NA_real_, p = NA_real_, refused = 0,
             nonconv = 1, error = 0))
  set.seed(1430000L + tag * 10000L + r)
  ans <- tryCatch(btl_dimensionality(fit, reps = reference_reps),
                  error = identity)
  if (inherits(ans, "condition") || is.na(ans$leading_structured) ||
      !is.finite(ans$reference$p_adj))
    return(c(reject = NA_real_, p = NA_real_, refused = 0,
             nonconv = 0, error = 1))
  c(reject = as.numeric(ans$leading_structured), p = ans$reference$p_adj,
    refused = 0, nonconv = 0, error = 0)
}

run_cell <- function(model, departure, reference_reps, tag) {
  f <- function(r) one(r, model, departure, reference_reps, tag)
  z <- if (CORES > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(NREP), f, mc.cores = CORES,
                       mc.set.seed = FALSE) else lapply(seq_len(NREP), f)
  do.call(rbind, z)
}

rows <- vector("list", nrow(scenarios))
for (s in seq_len(nrow(scenarios))) {
  model <- scenarios$model[s]
  departure <- scenarios$departure[s]
  reference_reps <- scenarios$reference_reps[s]
  z <- run_cell(model, departure, reference_reps, s)
  ok <- is.finite(z[, "reject"]) & is.finite(z[, "p"])
  rate <- mean(z[ok, "reject"])
  rows[[s]] <- sv_row(
    "BTL dimensionality reference",
    sprintf("%s, %s, B = %d", model, departure, reference_reps),
    "finite simulated upper-tail decision",
    n_reps = sum(ok), n_attempted = NREP,
    n_refused = sum(z[, "refused"]), n_nonconv = sum(z[, "nonconv"]),
    n_error = sum(z[, "error"]),
    type1 = if (departure == "null") rate else NA_real_,
    power = if (departure != "null") rate else NA_real_,
    notes = sprintf(
      "rejection %.4f using (1 + exceedances) / (B + 1); mean p %.4f",
      rate, mean(z[ok, "p"])))
}

out <- do.call(rbind, rows)
sv_write(out, "btl-dimensionality-reference")
print(out[, c("scenario", "n_reps", "type1", "power", "notes")],
      row.names = FALSE)
