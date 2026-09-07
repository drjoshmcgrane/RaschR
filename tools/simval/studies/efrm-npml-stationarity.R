#!/usr/bin/env Rscript
# Numerical stationarity of the EFRM NPML set link. The production estimate is
# compared with further alternating mass and link updates at the same grid.
# This checks whether the final nuisance-mass refresh can leave the reported
# transformation materially short of the joint coordinate optimum.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

STUDY <- "efrm-npml-stationarity"
NREP <- as.integer(Sys.getenv("SV_REPS", "100"))
if (!is.finite(NREP) || NREP < 2L) stop("SV_REPS must be at least 2")

# Instrument a local copy of the internal fitter so this study can continue
# from the exact nuisance masses used by the production result. Nothing in the
# package namespace is changed.
npml_debug <- .efrm_npml_pair
fn_body <- body(npml_debug)
return_expression <- fn_body[[length(fn_body)]]
fn_body[[length(fn_body)]] <- bquote({
  .result <- .(return_expression)
  .result$debug_state <- list(
    logw = logw, grid = grid, La = La, oa = oa, ob = ob,
    score_a = score_a, score_b = score_b, count = count,
    mix_idx = mix_idx, tau_b = tau_v[cb], disc_b = db)
  .result
})
body(npml_debug) <- fn_body

generate <- function(seed, ratio, population, n = 400L, ips = 6L) {
  set.seed(seed)
  theta <- if (population == "normal") rnorm(n, 0, 1.2) else
    sample(c(-1.8, 1.8), n, replace = TRUE) + rnorm(n, 0, 0.5)
  difficulty <- seq(-1.5, 1.5, length.out = ips)
  Xa <- vapply(difficulty, function(d)
    rbinom(n, 1L, plogis(theta - d)), numeric(n))
  Xb <- vapply(difficulty, function(d)
    rbinom(n, 1L, plogis(ratio * theta - d)), numeric(n))
  list(X = cbind(Xa, Xb), difficulty = rep(difficulty, 2L), ips = ips)
}

polish <- function(dat, ratio) {
  ips <- dat$ips
  Xm <- dat$X
  vmap <- data.frame(set = rep(c("a", "b"), each = ips), group = "g")
  tau <- lapply(dat$difficulty, function(x) x)
  disc <- rep(1, 2L * ips)
  fitted <- npml_debug(
    Xm, vmap, tau, disc, c("a", "b"), 1L, 2L, seq_len(nrow(Xm)),
    init_log_ratio = log(ratio), init_offset = 0,
    min_link_persons = 30L, grid_n = 61L)
  if (is.null(fitted) || !isTRUE(fitted$converged)) return(NULL)

  state <- fitted$debug_state
  grid <- state$grid
  La <- state$La
  ob <- state$ob
  score_b <- state$score_b
  count <- state$count
  mix_idx <- state$mix_idx
  par <- c(fitted$log_ratio, fitted$offset)
  initial <- par
  ll_initial <- fitted$loglik
  logw <- state$logw
  objective <- function(z, current_logw) {
    Lb <- .efrm_npml_likelihood_r(
      exp(z[1L]) * grid + z[2L], ob, score_b,
      state$tau_b, state$disc_b)
    A <- La + Lb + current_logw[mix_idx, , drop = FALSE]
    mx <- apply(A, 1L, max)
    -sum(count * (mx + log(rowSums(exp(A - mx)))))
  }
  final_ll <- ll_initial
  for (iteration in seq_len(10L)) {
    Lb <- .efrm_npml_likelihood_r(
      exp(par[1L]) * grid + par[2L], ob, score_b,
      state$tau_b, state$disc_b)
    mass <- .efrm_npml_fit_weights_r(
      La + Lb, logw, mix_idx, count, maxit = 5000L,
      tol = 1e-7)
    if (!isTRUE(mass$converged)) return(NULL)
    logw <- mass$logw
    op <- stats::optim(
      par, objective, current_logw = logw, method = "L-BFGS-B",
      lower = c(log(.1), -20), upper = c(log(10), 20),
      control = list(maxit = 100L, factr = 1e5))
    if (op$convergence != 0L || any(!is.finite(op$par))) return(NULL)
    change <- max(abs(op$par - par))
    par <- op$par
    final_ll <- -op$value
    if (change < 1e-6) break
  }
  c(log_ratio_change = par[1L] - initial[1L],
    offset_change = par[2L] - initial[2L],
    loglik_gain = final_ll - ll_initial)
}

run_cell <- function(population, ratio, seed0) {
  draws <- lapply(seq_len(NREP), function(r) {
    dat <- generate(seed0 + r, ratio, population)
    polish(dat, ratio)
  })
  ok <- !vapply(draws, is.null, logical(1L))
  values <- if (any(ok)) do.call(rbind, draws[ok]) else
    matrix(numeric(0), 0L, 3L, dimnames = list(NULL,
      c("log_ratio_change", "offset_change", "loglik_gain")))
  do.call(rbind, lapply(colnames(values), function(quantity) {
    x <- values[, quantity]
    sv_row(
      STUDY, sprintf("%s persons; ratio %.1f; N=400", population, ratio),
      quantity, n_reps = length(x), n_attempted = NREP,
      n_refused = NREP - length(x), effect = 0,
      bias = if (length(x)) mean(x) else NA_real_,
      emp_sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
      notes = if (length(x)) sprintf(
        "median %.6g; maximum absolute %.6g",
        stats::median(x), max(abs(x))) else "no usable link")
  }))
}

rows <- rbind(
  run_cell("normal", 1.0, 9510000L),
  run_cell("normal", 1.4, 9520000L),
  run_cell("bimodal", 1.4, 9530000L))
sv_write(rows, STUDY)
print(rows[, c("scenario", "quantity", "n_reps", "bias", "emp_sd",
               "notes")], row.names = FALSE)
