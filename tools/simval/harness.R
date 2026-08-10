# Shared reporting helpers for the simulation validation studies.
# Source from the package root:  source("tools/simval/harness.R")
# Studies load the in-tree package with  pkgload::load_all(".").

# Monte Carlo standard error of an estimated proportion.
mc_se_prop <- function(p, n) sqrt(pmax(p * (1 - p), 0) / n)

# One result row per scenario x quantity. Everything optional except the
# identifiers; unused fields stay NA so the CSVs align across studies.
sv_row <- function(study, scenario, quantity, n_reps,
                   bias = NA_real_, emp_sd = NA_real_, mean_se = NA_real_,
                   se_ratio = if (is.finite(emp_sd) && is.finite(mean_se) &&
                                  mean_se > 0) emp_sd / mean_se else NA_real_,
                   coverage95 = NA_real_, type1 = NA_real_,
                   familywise = NA_real_, power = NA_real_,
                   effect = NA_real_, refusal_rate = NA_real_,
                   nonconv_rate = NA_real_, notes = "") {
  rate <- c(type1 = type1, familywise = familywise, power = power,
            coverage95 = coverage95)
  mc <- vapply(rate, function(p)
    if (is.finite(p)) mc_se_prop(p, n_reps) else NA_real_, 0)
  data.frame(study = study, scenario = scenario, quantity = quantity,
             n_reps = n_reps, effect = effect, bias = bias, emp_sd = emp_sd,
             mean_se = mean_se, se_ratio = se_ratio,
             coverage95 = coverage95, mc_se_coverage = mc[["coverage95"]],
             type1 = type1, mc_se_type1 = mc[["type1"]],
             familywise = familywise, mc_se_familywise = mc[["familywise"]],
             power = power, mc_se_power = mc[["power"]],
             refusal_rate = refusal_rate, nonconv_rate = nonconv_rate,
             notes = notes, stringsAsFactors = FALSE)
}

# Append rows to a study's CSV under tools/simval/results/.
sv_write <- function(rows, name) {
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
