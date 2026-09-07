# Randomised algorithm-conformance checks, not tests of nominal error rates.
# Compare the structural wrapper with an explicit refit that may rescore data.
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "structural-score-preservation"
run_cell <- function(model, boundary) {
  ans <- lapply(seq_len(20L), function(r) {
    d <- simulate_rasch(300, 4, model = "PCM", n_categories = 3,
      n_groups = 2, seed = 984000L + r)
    if (boundary)
      d[d$I01 == 2, c("I02", "I03")] <- 2
    tryCatch(suppressWarnings({
      fit <- if (model == "PCM") rasch(d, id = "id", factors = "group") else
        rasch_efrm(d, item_sets = list(s = sprintf("I%02d", 1:4)),
          groups = "group", id = "id", boot_reps = 0)
      keep <- c("I01", "I02", "I03")
      source <- if (model == "PCM") fit$X[, keep] else
        .efrm_source_matrix(fit, keep)
      direct <- if (model == "PCM") rasch(source) else
        rasch_efrm(data.frame(source, group = d$group),
          item_sets = list(s = keep), groups = "group", boot_reps = 0)
      actual <- if (model == "PCM") direct$X else .efrm_source_matrix(direct, keep)
      changed <- !isTRUE(all.equal(unname(actual), unname(source),
                                   check.attributes = FALSE, tolerance = 0))
      result <- tryCatch(drop_items(fit, "I04", boot_reps = 0), error = identity)
      guarded <- inherits(result, "error") &&
        grepl("cannot preserve the fitted score structure", conditionMessage(result))
      if (inherits(result, "error") && !guarded) stop(result)
      c(changed = changed, guarded = guarded)
    }), error = function(e) c(changed = NA, guarded = NA))
  })
  x <- do.call(rbind, ans)
  ok <- complete.cases(x)
  bad <- sum(x[ok, "changed"] != x[ok, "guarded"])
  stopifnot(bad == 0L, any(ok))
  row <- sv_row(STUDY, paste(model, if (boundary) "boundary categories" else "control"),
    "agreement with explicit refit scoring", sum(ok),
    n_attempted = nrow(x), n_refused = 0L,
    n_error = sum(!ok),
    notes = "20 randomised input designs; compared observed scores, not estimated parameters; guard refusals are intended outcomes")
  row$n_changed <- sum(x[ok, "changed"])
  row$n_guard_refused <- sum(x[ok, "guarded"])
  row$n_preserved <- sum(!x[ok, "changed"])
  row$n_guard_disagreements <- bad
  row
}
out <- do.call(rbind, lapply(c("PCM", "EFRM"), function(model)
  rbind(run_cell(model, FALSE), run_cell(model, TRUE))))
sv_write(out, STUDY)
