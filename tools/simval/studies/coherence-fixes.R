# STUDY: coherence-fixes
#
# Post-fix validation for three model-coherence defects: BTL-EFRM likelihoods
# after panel-unit reconciliation, and ordinary/BTL multifactor DIF magnitudes
# under deliberately correlated, unbalanced factor cells. Run from the package
# root. Set SV_CORES to parallelise on Unix-like systems.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "coherence-fixes"

cores <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
if (!is.finite(cores) || cores < 1L) cores <- 1L
map_reps <- function(n, fun) {
  if (cores > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(n), fun, mc.cores = cores,
                       mc.set.seed = FALSE, mc.preschedule = FALSE)
  else lapply(seq_len(n), fun)
}

btlef_one <- function(r, reps_within = 8L, reps_cross = 4L,
                      seed0 = 610000L) {
  d <- simulate_btl_efrm(
    n_objects_per_set = 5, n_sets = 3, n_panels = 3,
    n_judges_per_panel = 10, reps_within = reps_within,
    reps_cross = reps_cross,
    panel_units = c(0.75, 1, 1.35), set_units = c(1, 1.25, 0.8),
    set_origins = c(0, 0.4, -0.3), seed = seed0 + r)
  tr <- attr(d, "truth")
  fit <- tryCatch(suppressWarnings(btl_efrm(
    d, "object_a", "object_b", "winner", "judge", "panel",
    object_sets = tr$object_sets, se_method = "conditional", min_link = 5)),
    error = function(e) NULL)
  blank <- c(gap = NA, rmse = NA, bias = NA, phi2 = NA, alpha2 = NA,
             alpha3 = NA, kappa2 = NA, kappa3 = NA)
  if (is.null(fit)) return(c(blank, refused = 1, nonconv = 0))
  if (!isTRUE(fit$converged))
    return(c(blank, refused = 0, nonconv = 1))
  y <- fit$comparisons$response; p <- fit$comparisons$expected
  ll <- sum(y * log(p) + (1 - y) * log1p(-p))
  est <- setNames(fit$objects$v, fit$objects$object)[names(tr$v)]
  err <- est - tr$v
  c(gap = fit$loglik - ll, rmse = sqrt(mean(err^2)), bias = mean(err),
    phi2 = log(fit$phi_table$phi[fit$phi_table$panel == "panel2"]),
    alpha2 = log(fit$alpha_table$alpha[fit$alpha_table$set == "set2"]),
    alpha3 = log(fit$alpha_table$alpha[fit$alpha_table$set == "set3"]),
    kappa2 = fit$kappa_table$kappa[fit$kappa_table$set == "set2"],
    kappa3 = fit$kappa_table$kappa[fit$kappa_table$set == "set3"],
    refused = 0, nonconv = 0)
}

rasch_dif_one <- function(r) {
  set.seed(620000L + r)
  cell <- rep(c("AX", "AY", "BX", "BY"), c(600, 200, 200, 600))
  A <- factor(substr(cell, 1, 1)); B <- factor(substr(cell, 2, 2))
  n <- length(cell); theta <- rnorm(n)
  difficulty <- seq(-1.5, 1.5, length.out = 8)
  shift <- (A == "B") * 1 + (B == "Y") * 2
  X <- sapply(seq_along(difficulty), function(j)
    rbinom(n, 1, plogis(theta - difficulty[j] -
                          if (j == 3L) shift else 0)))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- tryCatch(suppressWarnings(rasch(
    data.frame(X, A = A, B = B), factors = c("A", "B"))),
    error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$est$converged))
    return(c(A = NA, Ase = NA, Ap = NA, B = NA, Bse = NA, Bp = NA,
             refused = is.null(fit), nonconv = !is.null(fit)))
  z <- dif_anova(fit)
  pa <- dif_posthoc(fit, "I3", "A", factors = c("A", "B"))$table
  pb <- dif_posthoc(fit, "I3", "B", factors = c("A", "B"))$table
  getp <- function(term) {
    x <- z$summary$p_uniform_adj[z$summary$item == "I3" &
                                   z$summary$term == term]
    if (length(x) == 1L) x else NA_real_
  }
  c(A = pa$estimate, Ase = pa$se, Ap = getp("A"),
    B = pb$estimate, Bse = pb$se, Bp = getp("B"),
    refused = 0, nonconv = 0)
}

btl_dif_one <- function(r) {
  set.seed(630000L + r)
  counts <- c(AX = 24, AY = 8, BX = 8, BY = 24)
  cells <- rep(names(counts), counts)
  judges <- sprintf("J%02d", seq_along(cells))
  A <- setNames(substr(cells, 1, 1), judges)
  B <- setNames(substr(cells, 2, 2), judges)
  objects <- paste0("O", 1:6)
  beta <- setNames(seq(-1, 1, length.out = 6), objects)
  pairs <- t(utils::combn(objects, 2))
  rows <- lapply(judges, function(j) {
    z <- data.frame(a = rep(pairs[, 1], each = 2),
                    b = rep(pairs[, 2], each = 2), judge = j)
    shift <- (A[j] == "B") * 1 + (B[j] == "Y") * 2
    eta <- beta[z$a] - beta[z$b] +
      shift * ((z$a == "O3") - (z$b == "O3"))
    z$winner <- ifelse(runif(nrow(z)) < plogis(eta), z$a, z$b)
    z
  })
  fit <- tryCatch(suppressWarnings(btl(
    do.call(rbind, rows), "a", "b", "winner", judge = "judge")),
    error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged))
    return(c(A = NA, Ase = NA, Ap = NA, B = NA, Bse = NA, Bp = NA,
             refused = is.null(fit), nonconv = !is.null(fit)))
  # A near-one alpha requests both complete-design magnitudes regardless of
  # their screening result; power below is still evaluated at adjusted p < .05.
  z <- tryCatch(btl_dif(fit, list(A = A, B = B), objects = "O3",
                        alpha = 0.999999),
                error = function(e) NULL)
  if (is.null(z))
    return(c(A = NA, Ase = NA, Ap = NA, B = NA, Bse = NA, Bp = NA,
             refused = 1, nonconv = 0))
  get <- function(term, col) {
    x <- z$sizes[z$sizes$term == term, col]
    if (length(x) == 1L) x else NA_real_
  }
  getp <- function(term) {
    x <- z$summary$p_uniform_adj[z$summary$term == term]
    if (length(x) == 1L) x else NA_real_
  }
  c(A = get("A", "difference"), Ase = get("A", "se"), Ap = getp("A"),
    B = get("B", "difference"), Bse = get("B", "se"), Bp = getp("B"),
    refused = 0, nonconv = 0)
}

summarise_dif <- function(z, label, attempted) {
  refused <- sum(z[, "refused"]); nonconv <- sum(z[, "nonconv"])
  truth <- c(A = 1, B = 2)
  out <- lapply(names(truth), function(nm) {
    target <- truth[[nm]]
    se_nm <- paste0(nm, "se"); p_nm <- paste0(nm, "p")
    ok <- is.finite(z[, nm]) & is.finite(z[, se_nm]) & z[, se_nm] > 0
    mag <- sv_row(
      STUDY, label, paste0(nm, " adjusted marginal magnitude"), sum(ok),
      effect = target, bias = mean(z[ok, nm]) - target,
      emp_sd = stats::sd(z[ok, nm]), mean_se = mean(z[ok, se_nm]),
      coverage95 = mean(abs(z[ok, nm] - target) <= 1.96 * z[ok, se_nm]),
      n_attempted = attempted,
      n_refused = refused, n_nonconv = nonconv,
      notes = paste("correlated unbalanced cells AX/AY/BX/BY =",
                    "3/1/1/3; both factors fitted jointly"))
    okp <- is.finite(z[, p_nm])
    det <- sv_row(
      STUDY, label, paste0(nm, " adjusted omnibus detection"), sum(okp),
      effect = target, power = mean(z[okp, p_nm] < 0.05),
      n_attempted = attempted, n_refused = refused, n_nonconv = nonconv,
      notes = "power denominator includes every analysed replicate")
    rbind(mag, det)
  })
  do.call(rbind, out)
}

n_btlef_sparse <- 80L; n_btlef_dense <- 60L; n_btlef_high <- 40L
n_rasch <- 100L; n_btl <- 100L
btlef_sparse <- do.call(rbind, map_reps(n_btlef_sparse, btlef_one))
btlef_dense <- do.call(rbind, map_reps(n_btlef_dense, function(r)
  btlef_one(r, reps_within = 24L, reps_cross = 12L, seed0 = 615000L)))
btlef_high <- do.call(rbind, map_reps(n_btlef_high, function(r)
  btlef_one(r, reps_within = 72L, reps_cross = 36L, seed0 = 617000L)))
rdif <- do.call(rbind, map_reps(n_rasch, rasch_dif_one))
bdif <- do.call(rbind, map_reps(n_btl, btl_dif_one))

phi_truth <- c(0.75, 1, 1.35)
phi_truth <- phi_truth / exp(mean(log(phi_truth)))
btlef_truth <- c(phi2 = log(phi_truth[2]), alpha2 = log(1.25),
                 alpha3 = log(0.8), kappa2 = 0.4, kappa3 = -0.3)
summarise_btlef <- function(btlef, label, attempted) {
  ok_btlef <- is.finite(btlef[, "gap"])
  rbind(
    sv_row(
      STUDY, label, "reported-minus-reconstructed loglik",
      sum(ok_btlef), bias = mean(btlef[ok_btlef, "gap"]),
      emp_sd = stats::sd(btlef[ok_btlef, "gap"]),
      n_attempted = attempted, n_refused = sum(btlef[, "refused"]),
      n_nonconv = sum(btlef[, "nonconv"]),
      notes = "identity check against every stored expected comparison probability"),
    sv_row(
      STUDY, label, "common-scale object recovery", sum(ok_btlef),
      bias = mean(btlef[ok_btlef, "bias"]),
      emp_sd = mean(btlef[ok_btlef, "rmse"]), n_attempted = attempted,
      n_refused = sum(btlef[, "refused"]),
      n_nonconv = sum(btlef[, "nonconv"]),
      notes = "emp_sd field holds mean within-replicate object RMSE"),
    do.call(rbind, lapply(names(btlef_truth), function(nm) {
      target <- btlef_truth[[nm]]
      ok <- is.finite(btlef[, nm])
      sv_row(
        STUDY, label, paste0(nm, " recovery"), sum(ok),
        effect = unname(target), bias = mean(btlef[ok, nm]) - target,
        emp_sd = stats::sd(btlef[ok, nm]), n_attempted = attempted,
        n_refused = sum(btlef[, "refused"]),
        n_nonconv = sum(btlef[, "nonconv"]),
        notes = "point-estimate recovery after refitting all sets at reconciled phi")
    })))
}

rows <- rbind(
  summarise_btlef(btlef_sparse,
                  "three sets, three panels; 8 within/4 cross repeats",
                  n_btlef_sparse),
  summarise_btlef(btlef_dense,
                  "three sets, three panels; 24 within/12 cross repeats",
                  n_btlef_dense),
  summarise_btlef(btlef_high,
                  "three sets, three panels; 72 within/36 cross repeats",
                  n_btlef_high),
  summarise_dif(rdif, "ordinary Rasch multifactor DIF", n_rasch),
  summarise_dif(bdif, "paired-comparison multifactor DIF", n_btl))

sv_write(rows, STUDY)
