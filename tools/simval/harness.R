# Shared reporting helpers for the simulation validation studies.
# Source from the package root:  source("tools/simval/harness.R")
# Studies load the in-tree package with  pkgload::load_all(".").
#
# CONVENTIONS
# - se_ratio is empirical SD / mean reported SE throughout the battery
#   (> 1 means the reported SE understates the sampling variability).
# - Replicate accounting distinguishes n_attempted (replicates started),
#   n_refused (identification/guard refusals), n_nonconv (fitted but not
#   converged), n_error (other failures), and n_reps = the ANALYSED replicates every conditional
#   performance quantity (bias, SD, coverage, rates) is computed on.
#   Rates carry Monte Carlo standard errors on their own denominators.

# Monte Carlo standard error of an estimated proportion.
mc_se_prop <- function(p, n) sqrt(pmax(p * (1 - p), 0) / n)

# Provenance, captured once when the harness is sourced: the package git
# SHA (marked +dirty when R code differs from HEAD), the run date, the
# running study script (from --file=, overridable via
# options(simval.script = "tools/...")), and that script's md5 hash --
# the package SHA does not identify an uncommitted study script, the
# content hash does.
.sv_prov <- local({
  sha <- tryCatch(suppressWarnings(
    system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)),
    error = function(e) character(0))
  dirty <- tryCatch(length(suppressWarnings(
    system("git status --porcelain -- R DESCRIPTION", intern = TRUE,
           ignore.stderr = TRUE))) > 0, error = function(e) NA)
  farg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- getOption("simval.script",
    if (length(farg) == 1) {
      f <- sub("^--file=", "", farg[1])
      tryCatch({
        rel <- sub(paste0("^", normalizePath(getwd()), "/"), "",
                   normalizePath(f))
        rel
      }, error = function(e) f)
    } else NA_character_)
  hash <- if (!is.na(script) && file.exists(script))
    unname(tools::md5sum(script)) else NA_character_
  rfiles <- sort(list.files("R", pattern = "[.]R$", full.names = TRUE))
  rtree <- if (length(rfiles)) {
    tf <- tempfile()
    writeLines(paste(unname(tools::md5sum(rfiles)), collapse = ""), tf)
    h <- substr(unname(tools::md5sum(tf)), 1, 12)
    unlink(tf)
    h
  } else NA_character_
  list(sha = if (length(sha) == 1)
    paste0(sha, if (isTRUE(dirty)) "+dirty" else "") else NA_character_,
    date = format(Sys.Date()), script = script, hash = hash, rtree = rtree)
})

# One result row per scenario x quantity. Everything optional except the
# identifiers; unused fields stay NA so the CSVs align across studies.
sv_row <- function(study, scenario, quantity, n_reps,
                   bias = NA_real_, emp_sd = NA_real_, mean_se = NA_real_,
                   se_ratio = if (is.finite(emp_sd) && is.finite(mean_se) &&
                                  mean_se > 0) emp_sd / mean_se else NA_real_,
                   coverage95 = NA_real_, type1 = NA_real_,
                   familywise = NA_real_, power = NA_real_,
                   mc_override = list(),
                   effect = NA_real_,
                   n_attempted = NA_integer_, n_refused = NA_integer_,
                   n_nonconv = NA_integer_, n_error = NA_integer_,
                   n_boot_attempted = NA_integer_, n_boot_used = NA_integer_,
                   n_boot_nonconv = NA_integer_, n_boot_errors = NA_integer_,
                   refusal_rate = if (is.finite(n_attempted) &&
                                      n_attempted > 0)
                     n_refused / n_attempted else NA_real_,
                   nonconv_rate = if (is.finite(n_attempted) &&
                                      n_attempted > 0)
                     n_nonconv / n_attempted else NA_real_,
                   error_rate = if (is.finite(n_attempted) &&
                                    n_attempted > 0)
                     n_error / n_attempted else NA_real_,
                   notes = "") {
  rate <- c(type1 = type1, familywise = familywise, power = power,
            coverage95 = coverage95)
  mc <- vapply(names(rate), function(q) {
    if (!is.null(mc_override[[q]])) return(mc_override[[q]])
    if (is.finite(rate[[q]])) mc_se_prop(rate[[q]], n_reps) else NA_real_
  }, 0)
  # mc_override: for rates POOLED over multiple correlated units per
  # replicate (e.g. per item x contrast cells), pass the cluster-robust
  # Monte Carlo SE  sd(per-replicate proportions)/sqrt(n_reps)  instead of
  # the plug-in binomial formula, which understates it.
  mc_ref <- if (is.finite(refusal_rate) && is.finite(n_attempted))
    mc_se_prop(refusal_rate, n_attempted) else NA_real_
  mc_nc <- if (is.finite(nonconv_rate) && is.finite(n_attempted))
    mc_se_prop(nonconv_rate, n_attempted) else NA_real_
  mc_err <- if (is.finite(error_rate) && is.finite(n_attempted))
    mc_se_prop(error_rate, n_attempted) else NA_real_
  data.frame(study = study, scenario = scenario, quantity = quantity,
             n_reps = n_reps, n_attempted = n_attempted,
             n_refused = n_refused, n_nonconv = n_nonconv,
             n_error = n_error,
             n_boot_attempted = n_boot_attempted,
             n_boot_used = n_boot_used,
             n_boot_nonconv = n_boot_nonconv,
             n_boot_errors = n_boot_errors,
             effect = effect, bias = bias, emp_sd = emp_sd,
             mean_se = mean_se, se_ratio = se_ratio,
             coverage95 = coverage95, mc_se_coverage = mc[["coverage95"]],
             type1 = type1, mc_se_type1 = mc[["type1"]],
             familywise = familywise, mc_se_familywise = mc[["familywise"]],
             power = power, mc_se_power = mc[["power"]],
             refusal_rate = refusal_rate, mc_se_refusal = mc_ref,
             nonconv_rate = nonconv_rate, mc_se_nonconv = mc_nc,
             error_rate = error_rate, mc_se_error = mc_err,
             script = .sv_prov$script, script_md5 = .sv_prov$hash,
             package_sha = .sv_prov$sha, r_tree_md5 = .sv_prov$rtree,
             executed = .sv_prov$date,
             notes = notes, stringsAsFactors = FALSE)
}

# Append rows to a study's CSV under tools/simval/results/. Before
# writing, re-hash the study script: if it changed on disk since the run
# started, the rows describe the LAUNCHED version (whose md5 they carry),
# not the current file -- warn loudly so the mismatch cannot pass silently.
sv_write <- function(rows, name) {
  if (!is.na(.sv_prov$script) && file.exists(.sv_prov$script)) {
    now <- unname(tools::md5sum(.sv_prov$script))
    if (!identical(now, .sv_prov$hash))
      warning(sprintf(paste(
        "%s changed on disk since this run started (md5 %s at launch, %s",
        "now); the CSV's script_md5 identifies the launched version"),
        .sv_prov$script, .sv_prov$hash, now), call. = FALSE)
  }
  dir.create("tools/simval/results", recursive = TRUE, showWarnings = FALSE)
  path <- file.path("tools/simval/results", paste0(name, ".csv"))
  utils::write.csv(rows, path, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", path, nrow(rows)))
  invisible(path)
}

# Empirical 95% CI coverage from per-replicate estimates and SEs.
sv_coverage <- function(est, se, truth, z = 1.96) {
  ok <- is.finite(est) & is.finite(se) & se > 0
  mean(abs(est[ok] - truth) <= z * se[ok])
}
