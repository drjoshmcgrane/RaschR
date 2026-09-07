# Conformance of the public conditional DIF bootstrap to an independently
# orchestrated refit loop. This is not an operating-characteristic study: the
# null-size and power studies are dif-conditional-bootstrap*.R. It verifies
# that the exported function preserves scores and missingness, repeats the
# declared design, aligns the complete item-by-term family, and calculates the
# documented minimum-p probabilities.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

dif_source_md5 <- unname(tools::md5sum("R/dif-bootstrap.R"))

B <- as.integer(Sys.getenv("SV_BOOT", "49"))
workers <- max(1L, as.integer(Sys.getenv("SV_CORES", "1")))

simulate_item <- function(theta, tau) vapply(theta, function(th)
  sample.int(length(tau) + 1L, 1L,
             prob = item_moments(th, tau)$P) - 1L, integer(1))

make_cases <- function() {
  set.seed(9011)
  N <- 240L; L <- 6L
  group <- factor(rep(c("A", "B"), each = N / 2L))
  region <- factor(ifelse(runif(N) < ifelse(group == "A", .35, .65),
                          "North", "South"))
  theta <- rnorm(N) + ifelse(group == "B", .35, 0)
  d <- seq(-1.2, 1.2, length.out = L)
  X <- matrix(rbinom(N * L, 1, plogis(outer(theta, d, "-"))), N, L)
  colnames(X) <- paste0("D", seq_len(L))
  X[seq(3, N, by = 11), 2] <- NA
  f1 <- rasch(X, factors = data.frame(group, region), n_groups = 3)
  a1 <- dif_anova(f1, effects = "factorial", n_groups = 3)

  set.seed(9012)
  people <- 110L; id <- rep(sprintf("P%03d", seq_len(people)), 2L)
  occasion <- factor(rep(c("pre", "post"), each = people))
  group2 <- factor(rep(rep(c("A", "B"), each = people / 2L), 2L))
  th2 <- rep(rnorm(people), 2L) + ifelse(occasion == "post", .2, 0)
  X2 <- matrix(rbinom(length(th2) * L, 1,
                      plogis(outer(th2, d, "-"))), length(th2), L)
  colnames(X2) <- paste0("W", seq_len(L))
  f2 <- rasch(X2, id = id, factors = data.frame(group = group2, occasion),
              n_groups = 3)
  a2 <- dif_anova(f2, effects = "factorial", within = "occasion",
                  n_groups = 3)

  set.seed(9013)
  N3 <- 260L; group3 <- factor(rep(c("A", "B"), each = N3 / 2L))
  th3 <- rnorm(N3)
  tau <- lapply(seq_len(5L), function(j)
    seq(-1, 1, length.out = 3L) + seq(-.6, .6, length.out = 5L)[j])
  X3 <- vapply(tau, function(tt) simulate_item(th3, tt), integer(N3))
  colnames(X3) <- paste0("P", seq_len(ncol(X3)))
  f3 <- rasch(X3, model = "PCM", factors = data.frame(group = group3),
              n_groups = 3)
  a3 <- dif_anova(f3, n_groups = 3)

  set.seed(9014)
  N4 <- 260L; group4 <- factor(rep(c("A", "B"), each = N4 / 2L))
  th4 <- rnorm(N4)
  common <- c(-.9, 0, .9)
  loc <- seq(-.7, .7, length.out = 5L)
  X4 <- vapply(loc, function(x) simulate_item(th4, common + x), integer(N4))
  colnames(X4) <- paste0("R", seq_len(ncol(X4)))
  f4 <- rasch(X4, model = "RSM", factors = data.frame(group = group4),
              n_groups = 3)
  a4 <- dif_anova(f4, n_groups = 3)

  list(multifactor = list(fit = f1, dif = a1),
       repeated = list(fit = f2, dif = a2),
       PCM = list(fit = f3, dif = a3),
       RSM = list(fit = f4, dif = a4))
}

manual_run <- function(fit, dif, seed) {
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  design <- dif$bootstrap_design
  family <- !dif$term_ids %in% c("Residuals", "ci")
  obs_key <- .factor_keys(data.frame(
    item = as.character(dif$terms$item[family]),
    term = dif$term_ids[family], stringsAsFactors = FALSE))
  F <- P <- matrix(NA_real_, B, length(obs_key))
  score_ok <- missing_ok <- logical(B)
  spec <- fit$refit_spec %||% list()
  for (b in seq_len(B)) {
    set.seed(seeds[b])
    xb <- .fit_gen_conditional(fit$X, fit$tau_list,
      if (anyNA(fit$X)) is.na(fit$X) else NULL)
    score_ok[b] <- identical(as.numeric(rowSums(xb, na.rm = TRUE)),
                             as.numeric(rowSums(fit$X, na.rm = TRUE)))
    missing_ok[b] <- identical(is.na(xb), is.na(fit$X))
    bf <- rasch(xb, model = fit$model, id = design$id,
                factors = design$factors,
                n_groups = spec$n_groups %||% fit$n_groups,
                anchors = spec$anchors, maxit = spec$maxit %||% 60L,
                tol = spec$tol %||% 1e-8)
    bd <- dif_anova(bf, factors = design$factors,
                    n_groups = design$n_groups,
                    p_adjust = design$p_adjust, alpha = design$alpha,
                    effects = design$effects, id = design$id,
                    within = design$within,
                    pool_facets = design$pool_facets)
    fam_b <- !bd$term_ids %in% c("Residuals", "ci")
    key_b <- .factor_keys(data.frame(
      item = as.character(bd$terms$item[fam_b]),
      term = bd$term_ids[fam_b], stringsAsFactors = FALSE))
    at <- match(obs_key, key_b)
    if (anyNA(at) || length(key_b) != length(obs_key))
      stop("manual refit changed the declared family")
    F[b, ] <- bd$terms$F_value[fam_b][at]
    P[b, ] <- bd$terms$p[fam_b][at]
  }
  list(F = F, P = P, score_ok = score_ok, missing_ok = missing_ok)
}

cases <- make_cases()
rows <- list()
for (j in seq_along(cases)) {
  nm <- names(cases)[j]
  fit <- cases[[j]]$fit; dif <- cases[[j]]$dif
  seed <- 9500L + j
  public <- suppressWarnings(dif_bootstrap(
    fit, dif, B = B, workers = workers, seed = seed))
  manual <- manual_run(fit, dif, seed)
  F <- manual$F; P <- manual$P
  family <- !dif$term_ids %in% c("Residuals", "ci")
  obs <- dif$terms[family, , drop = FALSE]
  p_raw <- vapply(seq_len(ncol(P)), function(k)
    (1 + sum(P[, k] <= obs$p[k])) / (B + 1), numeric(1))
  minp <- apply(P, 1L, min)
  p_adj <- vapply(obs$p, function(p)
    (1 + sum(minp <= p)) / (B + 1), numeric(1))
  legacy_raw <- vapply(seq_len(ncol(F)), function(k)
    (1 + sum(F[, k] >= obs$F_value[k])) / (B + 1), numeric(1))
  legacy_adj <- pmax(legacy_raw, p_adj)

  vals <- c(
    max_abs_F_difference = max(abs(public$replicates$F - F)),
    max_abs_marginal_p_difference = max(abs(
      public$terms$p_boot[family] - p_raw)),
    max_abs_adjusted_p_difference = max(abs(
      public$terms$p_boot_adj[family] - p_adj)),
    legacy_to_reference_marginal_max_change = max(abs(legacy_raw - p_raw)),
    legacy_to_reference_adjusted_max_change = max(abs(legacy_adj - p_adj)),
    legacy_to_reference_flag_changes = sum(
      (legacy_adj < dif$alpha) != (p_adj < dif$alpha)),
    score_preservation_failures = sum(!manual$score_ok),
    missingness_preservation_failures = sum(!manual$missing_ok))
  rows[[j]] <- do.call(rbind, lapply(names(vals), function(q)
    sv_row(study = "dif-bootstrap-public-conformance", scenario = nm,
      quantity = q, n_reps = 1L, n_attempted = 1L,
      n_refused = 0L, n_nonconv = 0L, n_error = 0L,
      n_boot_attempted = B, n_boot_used = public$B_used,
      n_boot_nonconv = public$B_nonconverged,
      n_boot_errors = public$B_errors, bias = unname(vals[q]),
      notes = paste("independent public-versus-manual orchestration;",
                    "exact score, missingness, family and reference-p",
                    "probability check; R/dif-bootstrap.R md5",
                    dif_source_md5))))
}
sv_write(do.call(rbind, rows), "dif-bootstrap-public-conformance")
