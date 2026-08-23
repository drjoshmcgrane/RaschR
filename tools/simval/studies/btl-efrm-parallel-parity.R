# BTL-EFRM judge-bootstrap parity and timing.
#
# The in-tree package is installed into an isolated temporary library because
# socket workers must load the same namespace as the coordinating R process.
# Judge resamples are generated before dispatch, so each worker count receives
# exactly the same bootstrap samples.

options(simval.script = "tools/simval/studies/btl-efrm-parallel-parity.R")
source("tools/simval/harness.R")
if (!requireNamespace("callr", quietly = TRUE))
  stop("this study needs the suggested package 'callr'")

lib <- tempfile("rasch-btlef-parallel-lib-")
dir.create(lib)
on.exit(unlink(lib, recursive = TRUE, force = TRUE), add = TRUE)
install_log <- tempfile("rasch-btlef-parallel-install-", fileext = ".log")
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
  options(rasch.max_workers = 4L)

  d <- simulate_btl_efrm(
    n_objects_per_set = 8, n_sets = 3, n_panels = 2,
    n_judges_per_panel = 10, reps_within = 20, reps_cross = 20,
    set_units = c(1, 1.2, 0.85), set_origins = c(0, 0.5, -0.3),
    seed = 7311)
  object_sets <- attr(d, "truth")$object_sets
  run <- function(workers) {
    elapsed <- system.time(fit <- btl_efrm(
      d, "object_a", "object_b", "winner", "judge", "panel", object_sets,
      se_method = "judge_bootstrap", boot_reps = 200,
      workers = workers, seed = 7312))[["elapsed"]]
    list(fit = fit, elapsed = unname(elapsed))
  }
  serial <- run(1L)
  parallel <- run(4L)

  a <- serial$fit
  b <- parallel$fit
  a$workers <- b$workers <- NULL
  same <- identical(a, b)
  numeric_difference <- function(x, y) {
    nx <- unlist(x, recursive = TRUE, use.names = TRUE)
    ny <- unlist(y, recursive = TRUE, use.names = TRUE)
    keep <- suppressWarnings(is.finite(as.numeric(nx)) &
                             is.finite(as.numeric(ny)))
    if (!any(keep)) return(0)
    max(abs(as.numeric(nx[keep]) - as.numeric(ny[keep])))
  }
  difference <- numeric_difference(a, b)
  if (!same || difference != 0 ||
      !identical(serial$fit$boot_reps_used, 200L) ||
      !identical(parallel$fit$boot_reps_used, 200L) ||
      !identical(parallel$fit$workers, 4L))
    stop("serial and parallel BTL-EFRM results did not agree exactly")

  list(serial_elapsed = serial$elapsed,
       parallel_elapsed = parallel$elapsed,
       max_difference = difference,
       exact = same,
       workers = parallel$fit$workers)
}, args = list(lib = lib), libpath = .libPaths(), timeout = Inf)

rows <- sv_row(
  "btl-efrm-parallel-parity", "judge_200_default_workers_4",
  "reported_estimates", n_reps = 200, n_attempted = 200,
  effect = result$serial_elapsed / result$parallel_elapsed,
  bias = result$max_difference,
  notes = sprintf(paste(
    "exact full-object agreement after excluding the recorded worker count:",
    "%s; elapsed %.3fs serial and %.3fs with %d workers;",
    "timing is machine-specific"), result$exact, result$serial_elapsed,
    result$parallel_elapsed, result$workers))

sv_write(rows, "btl-efrm-parallel-parity")
