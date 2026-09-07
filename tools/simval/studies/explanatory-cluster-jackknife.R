#!/usr/bin/env Rscript
# Small-cluster uncertainty for an explanatory Rasch coefficient when response
# occasions are unevenly distributed over persons. The production linearised
# delete-one-person covariance is compared with the former cluster sandwich,
# its CR1 rescaling, and a full delete-one-person refit jackknife.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "500"))
if (!is.finite(NREP) || NREP < 2L) stop("SV_REPS must be at least 2")

predictors <- data.frame(
  item = sprintf("I%02d", seq_len(12L)),
  feature = seq(-1, 1, length.out = 12L),
  stringsAsFactors = FALSE)

# First-order delete-cluster covariance.  For cluster g, let A_g be its
# contribution to the sensitivity matrix and s_g its score at the full-data
# estimate.  Solving (A - A_g) d_g = s_g gives the linearised change produced
# by deleting that cluster.  The ordinary jackknife covariance of the d_g is
# therefore a cheap CR3-type approximation to the refit jackknife.
linearised_cluster_se <- function(X, design, full, cluster) {
  B <- design$B
  tau <- full$thr$tau
  thr <- full$thr
  m <- full$m
  pairs <- .pair_counts(X, m)
  N <- nrow(X)
  P <- ncol(B)
  group_all <- .pcml_cluster_index(cluster, N)
  G_all <- max(group_all)
  S <- matrix(0, N, P)
  A_g <- array(0, c(G_all, P, P))
  pair_load <- integer(N)
  cum <- lapply(seq_along(m), function(i)
    cumsum(tau[thr$item == i]))
  ids <- lapply(seq_along(m), function(i) thr$id[thr$item == i])

  for (pc in pairs) {
    i <- pc$i; j <- pc$j; mi <- m[i]; mj <- m[j]
    Li <- c(0, cum[[i]]); Lj <- c(0, cum[[j]])
    idx <- c(ids[[i]], ids[[j]])
    B_pair <- B[idx, , drop = FALSE]
    score_cell <- matrix(0, (mi + 1L) * (mj + 1L), P)
    info_total <- vector("list", mi + mj + 1L)
    for (r in seq_len(mi + mj - 1L)) {
      ks <- max(0L, r - mj):min(mi, r)
      if (length(ks) < 2L) next
      lp <- -(Li[ks + 1L] + Lj[r - ks + 1L])
      lp <- lp - max(lp); prob <- exp(lp) / sum(exp(lp))
      U <- cbind(-outer(ks, seq_len(mi), ">="),
                 -outer(r - ks, seq_len(mj), ">="))
      storage.mode(U) <- "double"
      Ub <- U %*% B_pair
      mean_u <- drop(crossprod(Ub, prob))
      centred <- sweep(Ub, 2L, mean_u)
      score_cell[ks * (mj + 1L) + (r - ks) + 1L, ] <- centred
      info_total[[r + 1L]] <- crossprod(centred, prob * centred)
    }
    both <- which(!is.na(X[, i]) & !is.na(X[, j]))
    if (!length(both)) next
    total <- X[both, i] + X[both, j]
    informative <- pmin(mi, total) > pmax(0, total - mj)
    if (any(informative)) {
      rows <- both[informative]
      pair_load[rows] <- pair_load[rows] + 1L
      for (h in which(informative)) {
        g <- group_all[both[h]]
        A_g[g, , ] <- A_g[g, , ] + info_total[[total[h] + 1L]]
      }
    }
    cell <- X[both, i] * (mj + 1L) + X[both, j] + 1L
    S[both, ] <- S[both, , drop = FALSE] + score_cell[cell, , drop = FALSE]
  }

  contributes <- pair_load > 0L
  used <- unique(group_all[contributes])
  group <- match(group_all[contributes], used)
  S_g <- rowsum(S[contributes, , drop = FALSE], group, reorder = FALSE)
  A_g <- A_g[used, , , drop = FALSE]
  A <- apply(A_g, c(2L, 3L), sum)
  if (!isTRUE(all.equal(A, -full$H_beta, tolerance = 1e-7,
                        check.attributes = FALSE)))
    stop("cluster sensitivity decomposition did not reproduce the fit: ",
         paste(signif(c(candidate = A, fitted = -full$H_beta), 8),
               collapse = " / "))
  delta <- matrix(NA_real_, length(used), P)
  for (g in seq_along(used)) {
    Ag <- matrix(A_g[g, , , drop = FALSE], P, P)
    delta[g, ] <- tryCatch(solve(A - Ag, S_g[g, ]),
                            error = function(e) rep(NA_real_, P))
  }
  if (any(!is.finite(delta))) return(rep(NA_real_, P))
  delta <- sweep(delta, 2L, colMeans(delta))
  Ainv <- solve(A)
  cov_hc0 <- Ainv %*% crossprod(S_g) %*% Ainv
  cov_cr1 <- length(used) / (length(used) - 1L) * cov_hc0
  cov_linearised <- (length(used) - 1) / length(used) * crossprod(delta)
  list(
    hc0 = stats::setNames(sqrt(diag(cov_hc0)), colnames(B)),
    cr1 = stats::setNames(sqrt(diag(cov_cr1)), colnames(B)),
    linearised = stats::setNames(sqrt(diag(cov_linearised)), colnames(B)))
}

one <- function(seed, beta = 0) {
  set.seed(seed)
  repeats <- rep(c(1L, 4L), 10L)
  id0 <- sprintf("P%03d", seq_len(20L))
  id <- rep(id0, repeats)
  theta <- rep(rnorm(20L), repeats)
  delta <- beta * predictors$feature
  pr <- plogis(outer(theta, delta, "-"))
  X <- matrix(rbinom(length(pr), 1L, pr), nrow = length(theta),
              ncol = nrow(predictors),
              dimnames = list(NULL, predictors$item))
  design <- tryCatch(.explanatory_metadata(predictors, ~ feature, X, "item"),
                     error = identity)
  if (inherits(design, "condition")) return(NULL)
  full <- tryCatch(.pcml_design(
    X, design$B, parameter_names = colnames(design$B), cluster = id),
    error = identity)
  if (inherits(full, "condition") || !isTRUE(full$converged) ||
      !isTRUE(full$cluster_inference)) return(NULL)
  estimate <- unname(full$beta["feature"])
  candidates <- linearised_cluster_se(X, design, full, id)
  se_hc0 <- unname(candidates$hc0["feature"])
  se_cr1 <- unname(candidates$cr1["feature"])
  se_linearised <- sqrt(full$cov_beta["feature", "feature"])
  if (!is.finite(estimate) || !is.finite(se_cr1) || se_cr1 <= 0) return(NULL)

  loo <- vapply(id0, function(g) {
    keep <- id != g
    z <- tryCatch(.pcml_design(
      X[keep, , drop = FALSE], design$B,
      parameter_names = colnames(design$B), cluster = id[keep]),
      error = identity)
    if (inherits(z, "condition") || !isTRUE(z$converged)) NA_real_ else
      unname(z$beta["feature"])
  }, 0)
  se_jk <- if (all(is.finite(loo))) {
    sqrt((length(loo) - 1) / length(loo) * sum((loo - mean(loo))^2))
  } else NA_real_
  c(estimate = estimate, se_hc0 = se_hc0,
    se_cr1 = se_cr1,
    se_linearised = se_linearised,
    se_jackknife = se_jk)
}

run_cell <- function(beta, seed0) {
  z <- lapply(seq_len(NREP), function(r) one(seed0 + r, beta))
  ok <- !vapply(z, is.null, logical(1))
  M <- if (any(ok)) do.call(rbind, z[ok]) else
    matrix(numeric(0), 0L, 5L, dimnames = list(NULL,
      c("estimate", "se_hc0", "se_cr1", "se_linearised",
        "se_jackknife")))
  crit <- stats::qt(.975, 19L)
  methods <- c(`former HC0` = "se_hc0", `CR1 candidate` = "se_cr1",
               `production linearised delete-one-person` = "se_linearised",
               `delete-one-person jackknife` = "se_jackknife")
  do.call(rbind, lapply(names(methods), function(method) {
    se <- M[, methods[[method]]]
    use <- is.finite(M[, "estimate"]) & is.finite(se) & se > 0
    reject <- abs((M[use, "estimate"] - 0) / se) > crit
    sv_row(
      "explanatory-cluster-jackknife",
      if (beta == 0) "20 persons; alternating 1/4 rows; null" else
        "20 persons; alternating 1/4 rows; beta=0.30",
      paste(method, if (beta == 0) "coefficient Type I" else
        "coefficient power"),
      n_reps = sum(use), n_attempted = NREP,
      n_refused = NREP - sum(use), effect = beta,
      bias = if (any(use)) mean(M[use, "estimate"] - beta) else NA_real_,
      emp_sd = if (sum(use) > 1L) stats::sd(M[use, "estimate"]) else NA_real_,
      mean_se = if (any(use)) mean(se[use]) else NA_real_,
      coverage95 = if (any(use))
        mean(abs(M[use, "estimate"] - beta) <= crit * se[use]) else NA_real_,
      type1 = if (beta == 0 && any(use)) mean(reject) else NA_real_,
      power = if (beta != 0 && any(use))
        mean(abs(M[use, "estimate"] / se[use]) > crit) else NA_real_,
      notes = paste(
        "12 fixed dichotomous items; person-cluster t reference on 19 df;",
        "candidate comparison under highly unequal occasion counts"))
  }))
}

rows <- rbind(run_cell(0, 9410000L), run_cell(.30, 9420000L))
sv_write(rows, "explanatory-cluster-jackknife")
print(rows[, c("scenario", "quantity", "n_reps", "bias", "se_ratio",
               "coverage95", "type1", "power")], row.names = FALSE)
