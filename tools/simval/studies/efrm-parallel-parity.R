# EFRM parallel-bootstrap parity and timing.
#
# The in-tree package is installed into an isolated temporary library because
# socket workers must load the same compiled namespace as the coordinating R
# process. Random samples are generated before dispatch, so each worker count
# receives exactly the same bootstrap jobs.

options(simval.script = "tools/simval/studies/efrm-parallel-parity.R")
source("tools/simval/harness.R")
if (!requireNamespace("callr", quietly = TRUE))
  stop("this study needs the suggested package 'callr'")

lib <- tempfile("rasch-parallel-lib-")
dir.create(lib)
on.exit(unlink(lib, recursive = TRUE, force = TRUE), add = TRUE)
install_log <- tempfile("rasch-parallel-install-", fileext = ".log")
on.exit(unlink(install_log), add = TRUE)
status <- system2(file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--preclean", "--clean",
    paste0("--library=", lib), normalizePath(".")),
  stdout = install_log, stderr = install_log)
if (!identical(status, 0L))
  stop("isolated package installation failed:\n",
       paste(readLines(install_log, warn = FALSE), collapse = "\n"))

result <- callr::r(function(lib) {
  .libPaths(c(lib, .libPaths()))
  library(rasch)

  numeric_table_difference <- function(a, b, fields) {
    z <- unlist(lapply(fields, function(nm) {
      xa <- a[[nm]]; xb <- b[[nm]]
      cols <- intersect(names(xa)[vapply(xa, is.numeric, logical(1))],
                        names(xb)[vapply(xb, is.numeric, logical(1))])
      unlist(lapply(cols, function(k) abs(xa[[k]] - xb[[k]])))
    }))
    max(z[is.finite(z)], 0)
  }
  fit_difference <- function(a, b) {
    max(numeric_table_difference(a, b,
          c("alpha_table", "set_table", "phi_table", "frames",
            "thresholds_arbitrary", "item_arbitrary")),
        abs(a$linking$alpha_edges$loglik - b$linking$alpha_edges$loglik),
        na.rm = TRUE)
  }
  converged_same <- function(a, b)
    identical(a$linking$alpha_edges$converged,
              b$linking$alpha_edges$converged)

  demo <- rasch:::.app_example_data("efrm")
  items <- grep("^(Number|Algebra|Space)_", names(demo), value = TRUE)
  sets <- setNames(sub("_.*", "", items), items)
  hybrid <- lapply(list(1L, 2L, NULL), function(w) {
    args <- list(data = demo, item_sets = sets, groups = "year_group",
                 id = "person_id", boot_reps = 300, seed = 7123)
    if (!is.null(w)) args$workers <- w
    tm <- system.time(fit <- do.call(rasch_efrm, args))[["elapsed"]]
    list(fit = fit, elapsed = unname(tm), workers = fit$workers,
         default = is.null(w))
  })

  dat <- simulate_efrm(n_per_group = 150, items_per_set = 6, n_sets = 2,
                       n_groups = 2, seed = 119)
  truth <- attr(dat, "truth")
  full <- lapply(c(1L, 2L), function(w) {
    tm <- system.time(fit <- rasch_efrm(
      dat, item_sets = truth$item_sets, groups = "group", id = "id",
      se_method = "bootstrap", boot_reps = 30, workers = w,
      seed = 8124))[["elapsed"]]
    list(fit = fit, elapsed = unname(tm), workers = w)
  })

  hybrid_diff <- vapply(hybrid[-1L], function(x)
    fit_difference(hybrid[[1L]]$fit, x$fit), 0)
  hybrid_conv <- vapply(hybrid[-1L], function(x)
    converged_same(hybrid[[1L]]$fit, x$fit), TRUE)
  full_diff <- fit_difference(full[[1L]]$fit, full[[2L]]$fit)
  full_conv <- converged_same(full[[1L]]$fit, full[[2L]]$fit)
  if (any(hybrid_diff != 0) || !all(hybrid_conv) || full_diff != 0 ||
      !full_conv || !identical(hybrid[[3L]]$workers, 4L) ||
      any(vapply(full, function(x)
        !identical(x$fit$boot_reps_used, 30L), TRUE)))
    stop("serial and parallel EFRM bootstrap results did not agree exactly")

  list(
    hybrid = lapply(seq_along(hybrid[-1L]), function(i) {
      x <- hybrid[-1L][[i]]
      list(
      workers = x$workers, default = x$default, elapsed = x$elapsed,
      serial_elapsed = hybrid[[1L]]$elapsed,
      max_difference = hybrid_diff[i],
      convergence_mismatch = !hybrid_conv[i])
    }),
    full = list(workers = 2L, elapsed = full[[2L]]$elapsed,
      serial_elapsed = full[[1L]]$elapsed,
      max_difference = full_diff,
      convergence_mismatch = !full_conv)
  )
}, args = list(lib = lib), libpath = .libPaths(), timeout = Inf)

rows <- do.call(rbind, c(lapply(result$hybrid, function(x)
  sv_row("efrm-parallel-parity",
    sprintf("hybrid_300_%sworkers_%d",
            if (isTRUE(x$default)) "default_" else "", x$workers),
    "reported_estimates", n_reps = 300, n_attempted = 300,
    effect = x$serial_elapsed / x$elapsed,
    bias = x$max_difference,
    notes = sprintf(paste("maximum absolute serial/parallel difference;",
      "convergence mismatch %s; elapsed %.3fs serial and %.3fs parallel;",
      "timing is machine-specific"), x$convergence_mismatch,
      x$serial_elapsed, x$elapsed))),
  list(sv_row("efrm-parallel-parity", "full_30_workers_2",
    "reported_estimates", n_reps = 30, n_attempted = 30,
    effect = result$full$serial_elapsed / result$full$elapsed,
    bias = result$full$max_difference,
    notes = sprintf(paste("maximum absolute serial/parallel difference;",
      "convergence mismatch %s; elapsed %.3fs serial and %.3fs parallel;",
      "timing is machine-specific"), result$full$convergence_mismatch,
      result$full$serial_elapsed, result$full$elapsed)))))

sv_write(rows, "efrm-parallel-parity")
