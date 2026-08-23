#!/usr/bin/env Rscript
# Numerical parity of the compiled EFRM NPML kernel and its R reference.
# Run from the package root. Timing is deliberately not treated as a
# statistical result because it depends on compiler flags and hardware.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

d <- rasch:::.app_example_data("efrm")
items <- setdiff(names(d), c("person_id", "year_group"))
set_map <- setNames(sub("[_. -]*[0-9]+$", "", items), items)

fit_one <- function(cpp, reps, seed) {
  old <- options(rasch.efrm_cpp = cpp)
  on.exit(options(old), add = TRUE)
  set.seed(seed)
  rasch_efrm(
    d, item_sets = set_map, groups = "year_group", id = "person_id",
    items = items, maxit = 60, tol = 1e-8, boot_reps = reps
  )
}

reference <- fit_one(FALSE, 30L, 4401L)
compiled <- fit_one(TRUE, 30L, 4401L)

differences <- c(
  alpha = max(abs(reference$alpha_table$alpha - compiled$alpha_table$alpha)),
  se_log_alpha = max(abs(reference$alpha_table$se_log_alpha -
                           compiled$alpha_table$se_log_alpha)),
  origin = max(abs(reference$set_table$mu - compiled$set_table$mu)),
  edge_loglik = max(abs(reference$linking$alpha_edges$loglik -
                           compiled$linking$alpha_edges$loglik)),
  thresholds = max(abs(reference$thresholds_arbitrary$delta -
                          compiled$thresholds_arbitrary$delta))
)

rows <- do.call(rbind, lapply(names(differences), function(quantity)
  sv_row(
    "efrm-cpp-parity", "demo_hybrid_30", quantity, 1L,
    bias = unname(differences[[quantity]]),
    n_attempted = 1L, n_refused = 0L, n_nonconv = 0L,
    notes = "absolute difference: compiled kernel versus R reference; identical seed"
  )))

# Convergence is a discrete parity check, retained as a numeric zero/one row.
same_convergence <- identical(reference$linking$alpha_edges$converged,
                              compiled$linking$alpha_edges$converged)
rows <- rbind(rows, sv_row(
  "efrm-cpp-parity", "demo_hybrid_30", "convergence_mismatch", 1L,
  bias = as.numeric(!same_convergence), n_attempted = 1L,
  n_refused = 0L, n_nonconv = 0L,
  notes = "zero means every set-pair convergence flag agreed"
))

sv_write(rows, "efrm-cpp-parity")
