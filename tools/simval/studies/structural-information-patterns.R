# Information conformance under partial frame/facet administrations.
# Fixed calibration parameters; this does not estimate coverage or Type I error.
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "structural-information-patterns"
m <- rep(c(1L, 2L, 3L, 1L), 2)
tau <- lapply(seq_along(m), function(i)
  seq(-0.7, 0.7, length.out = m[i]) + (i - 4.5) / 5)
grid <- c(-1, 0, 1)
h <- 1e-4

rows <- lapply(c("efrm", "mfrm"), function(kind) {
  disc <- if (kind == "efrm") rep(c(0.7, 1 / 0.7), each = 4) else rep(1, 8)
  log_partition <- function(theta, i) {
    eta <- disc[i] * ((0:m[i]) * theta - c(0, cumsum(tau[[i]])))
    max(eta) + log(sum(exp(eta - max(eta))))
  }
  curvature <- vapply(seq_along(m), function(i) vapply(grid, function(th)
    (log_partition(th + h, i) - 2 * log_partition(th, i) +
       log_partition(th - h, i)) / h^2, 0), numeric(length(grid)))
  errors <- vapply(seq_len(50L), function(r) {
    set.seed(917400L + r)
    observed <- matrix(runif(40L * 8L) > 0.3, 40L, 8L)
    observed[1, ] <- FALSE
    X <- matrix(0L, 40L, 8L); X[!observed] <- NA
    vm <- data.frame(vkey = paste0("v", 1:8), item = paste0("I", 1:8),
      group = "G", set = rep(c("A", "B"), each = 4),
      rater = rep(c("R1", "R2"), each = 4))
    if (kind == "mfrm") vm$item <- rep(paste0("I", 1:4), 2)
    fit <- structure(list(est = list(converged = TRUE), X = X,
      alpha_table = data.frame(set = "reference", alpha = 1),
      linking = list(alpha_edges = data.frame(converged = TRUE)),
      tau_list = tau, m = m, disc = disc, virtual_map = vm,
      facet_spec = "rater", items = data.frame(item = vm$vkey)),
      class = c(paste0("rasch_", kind), "rasch"))
    blocks <- .design_blocks(fit)
    key <- function(ii) paste(ii, collapse = ",")
    expected_keys <- unique(apply(observed[rowSums(observed) > 0, , drop = FALSE],
                                  1, function(z) key(which(z))))
    stopifnot(setequal(vapply(blocks, key, ""), expected_keys),
              length(blocks) == length(expected_keys))
    info <- test_information(fit, grid)
    expected <- unlist(lapply(blocks, function(ii)
      rowSums(curvature[, ii, drop = FALSE])), use.names = FALSE)
    max(abs(info$info - expected))
  }, 0)
  stopifnot(all(is.finite(errors)), max(errors) < 1e-5)
  row <- sv_row(STUDY, kind, "exact patterns and numerical likelihood curvature",
    length(errors), n_attempted = length(errors), n_refused = 0L,
    n_nonconv = 0L, n_error = 0L,
    notes = paste("fixed parameters; eight response cells with mixed category",
      "counts; unequal set units in EFRM; 40 missingness patterns per design;",
      "three theta points; conformance only, no error-rate calibration"))
  row$max_absolute_curvature_error <- max(errors)
  row
})
sv_write(do.call(rbind, rows), STUDY)
