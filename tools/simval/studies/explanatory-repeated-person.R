#!/usr/bin/env Rscript
# Repeated-person coefficient inference for explanatory Rasch models.
#
# The item design and true coefficient are fixed within each cell. Repeated
# response rows share a person location, while responses are redrawn on every
# row. The study checks the linearised delete-one-person covariance and its
# finite-cluster t reference. Run from the package root:
#
#   Rscript tools/simval/studies/explanatory-repeated-person.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

STUDY <- "explanatory-repeated-person"

predictors <- data.frame(
  item = sprintf("I%02d", seq_len(12L)),
  feature = seq(-1, 1, length.out = 12L),
  stringsAsFactors = FALSE)

one_rep <- function(seed, n_person, repeat_counts, beta) {
  set.seed(seed)
  if (length(repeat_counts) == 1L)
    repeat_counts <- rep(repeat_counts, n_person)
  if (length(repeat_counts) != n_person || any(repeat_counts < 1L))
    stop("the repeat design must give every person at least one row")
  id0 <- sprintf("P%03d", seq_len(n_person))
  id <- rep(id0, repeat_counts)
  theta0 <- rnorm(n_person)
  theta <- rep(theta0, repeat_counts)
  delta <- beta * predictors$feature
  pr <- plogis(outer(theta, delta, "-"))
  X <- matrix(rbinom(length(pr), 1L, pr), nrow = length(theta),
              ncol = nrow(predictors))
  colnames(X) <- predictors$item

  fit <- tryCatch(
    rasch_explanatory(X, predictors, ~ feature, id = id),
    error = function(e) e)
  if (inherits(fit, "error"))
    return(list(status = "error", message = conditionMessage(fit)))
  if (!isTRUE(fit$est$converged) ||
      !isTRUE(fit$reference_fit$est$converged))
    return(list(status = "nonconv"))
  if (!isTRUE(fit$est$cluster_inference))
    return(list(status = "withheld"))
  cf <- fit$est$coefficients["feature", , drop = FALSE]
  if (!all(is.finite(unlist(cf[c("estimate", "se", "df", "p")],
                              use.names = FALSE))) || cf$se <= 0)
    return(list(status = "unavailable"))
  list(status = "analysed", estimate = cf$estimate, se = cf$se,
       df = cf$df, reject = cf$p < 0.05)
}

run_cell <- function(label, n_person, repeats, beta, n_reps, seed0) {
  ans <- lapply(seq_len(n_reps), function(r)
    one_rep(seed0 + r, n_person, repeats, beta))
  status <- vapply(ans, `[[`, character(1L), "status")
  ok <- status == "analysed"
  est <- vapply(ans[ok], `[[`, numeric(1L), "estimate")
  se <- vapply(ans[ok], `[[`, numeric(1L), "se")
  df <- vapply(ans[ok], `[[`, numeric(1L), "df")
  reject <- vapply(ans[ok], `[[`, logical(1L), "reject")
  coverage <- abs(est - beta) <= stats::qt(0.975, df) * se
  account <- list(
    n_attempted = n_reps,
    n_refused = 0L,
    n_nonconv = sum(status == "nonconv"),
    n_error = sum(status == "error"),
    n_withheld = sum(status == "withheld"),
    n_metric_unavailable = sum(status == "unavailable"))
  notes <- paste(
    "12 fixed dichotomous items; one centred continuous item predictor;",
    "response rows for each person share theta; linearised delete-one-person",
    "covariance and t reference with contributing-person df")
  common <- c(list(
    study = STUDY, scenario = label, n_reps = length(est),
    effect = beta, bias = mean(est - beta), emp_sd = stats::sd(est),
    mean_se = mean(se), coverage95 = mean(coverage), notes = notes), account)
  if (beta == 0) {
    do.call(sv_row, c(common, list(
      quantity = "coefficient null rejection", type1 = mean(reject))))
  } else {
    do.call(sv_row, c(common, list(
      quantity = "coefficient power", power = mean(reject))))
  }
}

cells <- list(
  list(label = "20 persons; 2 rows/person; null", n_person = 20L,
       repeats = 2L, beta = 0, n_reps = 1000L, seed0 = 9310000L),
  list(label = "40 persons; 2 rows/person; null", n_person = 40L,
       repeats = 2L, beta = 0, n_reps = 500L, seed0 = 9320000L),
  list(label = "80 persons; 4 rows/person; null", n_person = 80L,
       repeats = 4L, beta = 0, n_reps = 500L, seed0 = 9330000L),
  list(label = "20 persons; alternating 1/4 rows; null", n_person = 20L,
       repeats = rep(c(1L, 4L), 10L), beta = 0, n_reps = 2000L,
       seed0 = 9340000L),
  list(label = "40 persons; 2 rows/person; beta=0.30", n_person = 40L,
       repeats = 2L, beta = 0.30, n_reps = 500L, seed0 = 9350000L))

workers <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (length(workers) != 1L || !is.finite(workers) || workers < 1L)
  workers <- 4L
rows <- parallel::mclapply(cells, function(z)
  do.call(run_cell, z), mc.cores = min(4L, length(cells), workers))
rows <- do.call(rbind, rows)
sv_write(rows, STUDY)
print(rows[, c("scenario", "n_reps", "bias", "se_ratio", "coverage95",
               "type1", "power", "n_nonconv", "n_error", "n_withheld")],
      row.names = FALSE)
