# Item fit bootstrap: worker parity and timing.
#
# The in-tree package is installed into an isolated temporary library because
# socket workers must load the same compiled namespace as the coordinating R
# process -- which is why this cannot be answered from the test suite, where
# rasch.max_workers is pinned to 1 and the parallel path never runs. Person
# locations and the per-replicate seeds are drawn before dispatch, so each
# worker count must receive exactly the same bootstrap jobs and return
# exactly the same statistics.

options(simval.script = "tools/simval/studies/item-fit-parallel-parity.R")
source("tools/simval/harness.R")
if (!requireNamespace("callr", quietly = TRUE))
  stop("this study needs the suggested package 'callr'")

lib <- tempfile("rasch-fitboot-lib-")
dir.create(lib)
on.exit(unlink(lib, recursive = TRUE, force = TRUE), add = TRUE)
install_log <- tempfile("rasch-fitboot-install-", fileext = ".log")
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

  designs <- list(
    "500 persons x 10 dichotomous items" =
      list(n = 500L, L = 10L, k = 2L, model = "dichotomous", B = 99L),
    "4000 persons x 30 four-category items" =
      list(n = 4000L, L = 30L, k = 4L, model = "PCM", B = 25L))
  rows <- list()
  for (nm in names(designs)) {
    g <- designs[[nm]]
    d <- simulate_rasch(g$n, g$L, model = g$model, n_categories = g$k, seed = 4)
    f <- rasch(d, id = "id",
               model = if (g$model == "dichotomous") "PCM" else g$model)
    ref <- NULL
    for (w in c(1L, 2L, 4L, 8L)) {
      t0 <- proc.time()[["elapsed"]]
      bs <- suppressWarnings(fit_bootstrap(f, B = g$B, seed = 7, workers = w))
      el <- proc.time()[["elapsed"]] - t0
      if (is.null(ref)) ref <- bs
      diff <- max(abs(unlist(bs$items[vapply(bs$items, is.numeric, TRUE)]) -
                      unlist(ref$items[vapply(ref$items, is.numeric, TRUE)])),
                  0, na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        design = nm, B = g$B, workers = w, seconds = el,
        max_abs_difference = diff, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}, args = list(lib = lib))

rows <- list()
for (nm in unique(result$design)) {
  z <- result[result$design == nm, ]
  serial <- z$seconds[z$workers == 1L]
  for (i in seq_len(nrow(z)))
    rows[[length(rows) + 1L]] <- sv_row(
      "item fit bootstrap parity", nm,
      sprintf("%d worker%s", z$workers[i], if (z$workers[i] == 1L) "" else "s"),
      n_reps = z$B[i], bias = z$max_abs_difference[i],
      notes = sprintf("%.1f s, %.2fx the serial run; largest absolute difference from the serial result %.3g",
                      z$seconds[i], serial / z$seconds[i], z$max_abs_difference[i]))
}
out <- do.call(rbind, rows)
sv_write(out, "item-fit-parallel-parity")
print(result, row.names = FALSE, digits = 3)
