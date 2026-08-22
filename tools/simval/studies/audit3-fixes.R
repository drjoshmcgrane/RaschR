# Validation of the third-audit corrections:
#   1. repeated-measures DIF follow-ups retain equal design-cell margins;
#   2. BTL-DIF pairwise inference requires eight effective judges per level;
#   3. the spread diagnostic tests the lower-bound departure rather than
#      treating every point estimate below the bound as dependence.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

rows <- list()

# -------------------------------------------------------------------------
# A. Repeated-measures DIF: equal-cell target under an imbalanced nuisance
#    factor. The true marginal time effect is zero (+1 in a, -1 in b), while
#    the superseded shortcut targets 0.1(+1) + 0.9(-1) = -0.8.
# -------------------------------------------------------------------------
NREP_DIF <- 2000L
n_a <- 20L; n_b <- 180L; n <- n_a + n_b
h <- factor(c(rep("a", n_a), rep("b", n_b)))
id <- rep(sprintf("P%03d", seq_len(n)), 2)
fac <- data.frame(time = factor(rep(c("t1", "t2"), each = n)),
                  h = rep(h, 2))
grp <- .factor_cells(fac, sep = ":")
cellmap <- unique(data.frame(cell = as.character(grp), fac,
                             check.names = FALSE))
cellmap <- cellmap[match(levels(grp), cellmap$cell), , drop = FALSE]
family <- .dif_contrast_family(fac, cellmap, "time")
w <- family$family[["time: t2 - t1"]]

p_equal <- p_legacy <- rep(NA_real_, NREP_DIF)
set.seed(2.6e7)
for (r in seq_len(NREP_DIF)) {
  delta <- ifelse(h == "a", 1, -1) + rnorm(n, 0, 1)
  z <- c(-delta / 2, delta / 2)
  tst <- .dif_paired_cell_contrast(z, fac, grp, id, "time", cellmap, w)
  p_equal[r] <- tst$p
  p_legacy[r] <- stats::t.test(delta)$p.value
}
rows[[length(rows) + 1L]] <- sv_row(
  "audit3-fixes", "DIF time main effect; nuisance cells 10:90",
  "equal-cell repeated-measures Type I", NREP_DIF,
  type1 = mean(p_equal < 0.05),
  n_attempted = NREP_DIF, n_refused = 0L, n_nonconv = 0L,
  notes = "true cell effects +1 and -1; full cell weights; Welch-Satterthwaite reference")
rows[[length(rows) + 1L]] <- sv_row(
  "audit3-fixes", "DIF time main effect; nuisance cells 10:90",
  "superseded person-frequency shortcut rejection", NREP_DIF,
  type1 = mean(p_legacy < 0.05),
  n_attempted = NREP_DIF, n_refused = 0L, n_nonconv = 0L,
  notes = "diagnostic comparison only: this shortcut targeted -0.8 rather than the zero equal-cell estimand")

# The same check for a mixed time-by-group interaction, marginalised over an
# imbalanced between-person nuisance factor.
n_s1 <- 60L; n_s2 <- 540L; n <- n_s1 + n_s2
site <- factor(c(rep("s1", n_s1), rep("s2", n_s2)))
between_group <- factor(unlist(lapply(c(n_s1, n_s2), function(nn)
  rep(c("A", "B"), each = nn / 2))))
id <- rep(sprintf("M%03d", seq_len(n)), 2)
fac <- data.frame(time = factor(rep(c("t1", "t2"), each = n)),
                  group = rep(between_group, 2), site = rep(site, 2))
grp <- .factor_cells(fac, sep = ":")
cellmap <- unique(data.frame(cell = as.character(grp), fac,
                             check.names = FALSE))
cellmap <- cellmap[match(levels(grp), cellmap$cell), , drop = FALSE]
family <- .dif_posthoc_family(fac, cellmap, c("time", "group"), "time")
w <- family$family[[1]]
p_equal <- p_legacy <- rep(NA_real_, NREP_DIF)
set.seed(2.7e7)
for (r in seq_len(NREP_DIF)) {
  noise <- stats::rnorm(n, 0, 1)
  site_effect <- ifelse(site == "s1", 1, -1)
  delta <- noise + ifelse(between_group == "B", site_effect, 0)
  z <- c(-delta / 2, delta / 2)
  tst <- .dif_paired_cell_contrast(z, fac, grp, id, "time", cellmap, w)
  p_equal[r] <- tst$p
  p_legacy[r] <- stats::t.test(delta[between_group == "B"],
                               delta[between_group == "A"])$p.value
}
rows[[length(rows) + 1L]] <- sv_row(
  "audit3-fixes", "DIF time x group; nuisance sites 10:90",
  "equal-cell mixed-interaction Type I", NREP_DIF,
  type1 = mean(p_equal < 0.05),
  n_attempted = NREP_DIF, n_refused = 0L, n_nonconv = 0L,
  notes = "true site-specific interactions +1 and -1; full cell weights; Welch-Satterthwaite reference")
rows[[length(rows) + 1L]] <- sv_row(
  "audit3-fixes", "DIF time x group; nuisance sites 10:90",
  "superseded person-frequency shortcut rejection", NREP_DIF,
  type1 = mean(p_legacy < 0.05),
  n_attempted = NREP_DIF, n_refused = 0L, n_nonconv = 0L,
  notes = "diagnostic comparison only: this shortcut targeted -0.8 rather than the zero equal-cell interaction")

# -------------------------------------------------------------------------
# B. BTL-DIF resolved pairwise contrasts under judge-specific object
#    preferences. Every judge supplies three observations per object pair, so
#    raw and effective judges per factor level coincide.
# -------------------------------------------------------------------------
btl_pair_batch <- function(n_levels, judges_per_level, n_reps, seed0,
                           workload = NULL) {
  K <- 6L
  J <- n_levels * judges_per_level
  objects <- sprintf("O%d", seq_len(K))
  judges <- sprintf("J%d", seq_len(J))
  group <- setNames(rep(sprintf("g%d", seq_len(n_levels)),
                        each = judges_per_level), judges)
  if (is.null(workload)) workload <- rep(1L, J)
  stopifnot(length(workload) == J, all(workload >= 1))
  lev_eff <- vapply(seq_len(n_levels), function(g) {
    ww <- workload[group == sprintf("g%d", g)]
    sum(ww)^2 / sum(ww^2)
  }, 0)
  object_pairs <- t(utils::combn(objects, 2))
  p_pair_rep <- p_global_rep <- vector("list", n_reps)
  n_refused <- n_nonconv <- 0L
  for (r in seq_len(n_reps)) {
    set.seed(seed0 + r)
    judge_effect <- matrix(stats::rnorm(J * K, 0, 0.8), J, K,
                           dimnames = list(judges, objects))
    dat <- do.call(rbind, lapply(seq_along(judges), function(jj) {
      j <- judges[jj]
      d <- data.frame(object_a = rep(object_pairs[, 1],
                                     each = 3 * workload[jj]),
                      object_b = rep(object_pairs[, 2],
                                     each = 3 * workload[jj]),
                      judge = j)
      beta <- as.numeric(scale(seq_len(K)))
      eta <- beta[match(d$object_a, objects)] -
        beta[match(d$object_b, objects)] +
        judge_effect[j, d$object_a] - judge_effect[j, d$object_b]
      d$winner <- ifelse(stats::runif(nrow(d)) < stats::plogis(eta),
                         d$object_a, d$object_b)
      d
    }))
    fit <- tryCatch(
      btl(dat, "object_a", "object_b", winner = "winner", judge = "judge"),
      error = function(e) NULL)
    if (is.null(fit)) {
      n_refused <- n_refused + 1L
      next
    }
    if (!isTRUE(fit$converged)) {
      n_nonconv <- n_nonconv + 1L
      next
    }
    cm <- fit$comparisons
    cell <- unname(group[cm$judge])
    a2 <- ifelse(cm$object_a == "O3",
                 paste0("O3 (", cell, ")"), cm$object_a)
    b2 <- ifelse(cm$object_b == "O3",
                 paste0("O3 (", cell, ")"), cm$object_b)
    resolved <- tryCatch(.btl_graded(
      a2, b2, cm$response, cm$judge, cm$weight, fit$categories,
      60, 1e-8, character(0), thr = fit$thr_structure),
      error = function(e) NULL)
    if (is.null(resolved)) {
      n_refused <- n_refused + 1L
      next
    }
    if (!isTRUE(resolved$converged)) {
      n_nonconv <- n_nonconv + 1L
      next
    }
    idx <- match(paste0("O3 (g", seq_len(n_levels), ")"),
                 resolved$objects$object)
    if (anyNA(idx) || anyNA(resolved$cov_beta[idx, idx])) {
      n_refused <- n_refused + 1L
      next
    }
    loc <- resolved$objects$location[idx]
    V <- resolved$cov_beta[idx, idx, drop = FALSE]
    pr <- t(utils::combn(seq_len(n_levels), 2))
    se <- sqrt(diag(V)[pr[, 1]] + diag(V)[pr[, 2]] -
                 2 * V[cbind(pr[, 1], pr[, 2])])
    df <- lev_eff[pr[, 1]] + lev_eff[pr[, 2]] - 2
    stat <- (loc[pr[, 1]] - loc[pr[, 2]]) / se
    p_pair_rep[[r]] <- 2 * stats::pt(-abs(stat), df = df)
    p_global_rep[[r]] <- 2 * stats::pt(-abs(stat), df = J - 1)
  }
  p_pair_rep <- Filter(Negate(is.null), p_pair_rep)
  p_global_rep <- Filter(Negate(is.null), p_global_rep)
  pair_rep <- vapply(p_pair_rep, function(p) mean(p < 0.05), 0)
  global_rep <- vapply(p_global_rep, function(p) mean(p < 0.05), 0)
  list(type1_pair = mean(pair_rep),
       mcse_pair = stats::sd(pair_rep) / sqrt(length(pair_rep)),
       type1_global = mean(global_rep),
       mcse_global = stats::sd(global_rep) / sqrt(length(global_rep)),
       n = length(pair_rep), n_attempted = n_reps,
       n_refused = n_refused, n_nonconv = n_nonconv,
       effective_judges = lev_eff)
}

for (cell in list(c(2L, 6L, 800L, 3.0e7),
                  c(2L, 8L, 600L, 3.1e7),
                  c(2L, 10L, 600L, 3.2e7),
                  c(5L, 4L, 400L, 3.3e7))) {
  n_levels <- cell[1]; judges_per_level <- cell[2]
  res <- btl_pair_batch(n_levels, judges_per_level, cell[3], cell[4])
  thin <- judges_per_level < 8
  rows[[length(rows) + 1L]] <- sv_row(
    "audit3-fixes",
    sprintf("BTL-DIF %d levels x %d judges", n_levels, judges_per_level),
    if (thin) "underlying raw pairwise Type I in withheld region"
    else "pairwise Type I in inference region",
    res$n, type1 = res$type1_pair,
    mc_override = list(type1 = res$mcse_pair),
    n_attempted = res$n_attempted, n_refused = res$n_refused,
    n_nonconv = res$n_nonconv,
    notes = if (thin)
      "public btl_dif withholds these pairwise SEs and probabilities because each level has fewer than eight effective judges"
    else
      "balanced workload; pair-specific df = effective judges in the two levels minus 2")
  if (thin)
    rows[[length(rows) + 1L]] <- sv_row(
      "audit3-fixes",
      sprintf("BTL-DIF %d levels x %d judges", n_levels, judges_per_level),
      "superseded global-df pairwise Type I", res$n,
      type1 = res$type1_global,
      mc_override = list(type1 = res$mcse_global),
      n_attempted = res$n_attempted, n_refused = res$n_refused,
      n_nonconv = res$n_nonconv,
      notes = "diagnostic comparison only: df used every judge in the fit rather than the judges supporting the two cells")
}

# Mild workload imbalance remains inside the public inference region. One
# judge in each level carries twice the ordinary allocation, leaving 9.3
# effective judges per level.
res <- btl_pair_batch(2L, 10L, 600L, 3.35e7,
                      workload = c(2L, rep(1L, 9L),
                                   2L, rep(1L, 9L)))
rows[[length(rows) + 1L]] <- sv_row(
  "audit3-fixes", "BTL-DIF 2 x 10 judges; one twofold workload per level",
  "pairwise Type I in imbalanced inference region",
  res$n, type1 = res$type1_pair,
  mc_override = list(type1 = res$mcse_pair),
  n_attempted = res$n_attempted, n_refused = res$n_refused,
  n_nonconv = res$n_nonconv,
  notes = sprintf("public inference retained with a caution: %.1f effective judges per level",
                  res$effective_judges[1]))

# Raw counts can also conceal concentrated workloads. Here both levels contain
# ten judges, but one judge carries five times the ordinary allocation,
# reducing that level to 5.8 effective judges. Public inference is withheld.
res <- btl_pair_batch(2L, 10L, 600L, 3.4e7,
                      workload = c(5L, rep(1L, 9L), rep(1L, 10L)))
rows[[length(rows) + 1L]] <- sv_row(
  "audit3-fixes", "BTL-DIF 2 x 10 judges; one fivefold workload",
  "underlying raw pairwise Type I in workload-withheld region",
  res$n, type1 = res$type1_pair,
  mc_override = list(type1 = res$mcse_pair),
  n_attempted = res$n_attempted, n_refused = res$n_refused,
  n_nonconv = res$n_nonconv,
  notes = sprintf("public inference withheld: effective judges %.1f and %.1f",
                  res$effective_judges[1], res$effective_judges[2]))
rows[[length(rows) + 1L]] <- sv_row(
  "audit3-fixes", "BTL-DIF 2 x 10 judges; one fivefold workload",
  "superseded global-df pairwise Type I", res$n,
  type1 = res$type1_global,
  mc_override = list(type1 = res$mcse_global),
  n_attempted = res$n_attempted, n_refused = res$n_refused,
  n_nonconv = res$n_nonconv,
  notes = "diagnostic comparison only: global df ignores workload concentration within the factor level")

# -------------------------------------------------------------------------
# C. Spread lower-bound test. Equal-difficulty independent items sit at the
#    boundary; different difficulties move above it; planted dependence moves
#    below it.
# -------------------------------------------------------------------------
spread_batch <- function(kind = c("boundary", "above", "dependent"),
                         n_reps, seed0) {
  kind <- match.arg(kind)
  below <- reject <- rep(NA, n_reps)
  n_refused <- n_nonconv <- 0L
  difficulty <- if (kind == "above")
    c(-1, 1, -1.5, -0.8, -0.3, 0.3, 0.8, 1.5)
  else c(0, 0, -1.5, -0.8, -0.3, 0.3, 0.8, 1.5)
  for (r in seq_len(n_reps)) {
    dat <- simulate_rasch(
      900, 8, difficulty = difficulty,
      dependence = if (kind == "dependent")
        list(pairs = list(c(1, 2)), strength = 3) else NULL,
      seed = seed0 + r)
    fit <- tryCatch(rasch(dat, id = "id"), error = function(e) NULL)
    if (is.null(fit)) {
      n_refused <- n_refused + 1L
      next
    }
    if (!isTRUE(fit$est$converged)) {
      n_nonconv <- n_nonconv + 1L
      next
    }
    ans <- tryCatch(
      spread_test(combine_items(fit, list(c("I01", "I02")))),
      error = function(e) NULL)
    if (is.null(ans)) {
      n_refused <- n_refused + 1L
      next
    }
    below[r] <- ans$below_bound
    reject[r] <- ans$dependent
  }
  ok <- !is.na(reject)
  list(below = mean(below[ok]), reject = mean(reject[ok]), n = sum(ok),
       n_attempted = n_reps, n_refused = n_refused,
       n_nonconv = n_nonconv)
}

for (cell in list(c("boundary", 1000, 4.1e7),
                  c("above", 400, 4.2e7),
                  c("dependent", 400, 4.3e7))) {
  ans <- spread_batch(cell[1], as.integer(cell[2]), as.numeric(cell[3]))
  rows[[length(rows) + 1L]] <- sv_row(
    "audit3-fixes", sprintf("spread %s", cell[1]),
    "raw point estimate below bound", ans$n, type1 = ans$below,
    n_attempted = ans$n_attempted, n_refused = ans$n_refused,
    n_nonconv = ans$n_nonconv,
    notes = "descriptive below_bound rate; not an inferential verdict")
  rows[[length(rows) + 1L]] <- sv_row(
    "audit3-fixes", sprintf("spread %s", cell[1]),
    if (cell[1] == "dependent") "power of adjusted one-sided test"
    else "Type I of adjusted one-sided test",
    ans$n,
    type1 = if (cell[1] == "dependent") NA_real_ else ans$reject,
    power = if (cell[1] == "dependent") ans$reject else NA_real_,
    n_attempted = ans$n_attempted, n_refused = ans$n_refused,
    n_nonconv = ans$n_nonconv,
    notes = "one eligible superitem, so Holm adjustment leaves the one-sided p unchanged")
}

sv_write(do.call(rbind, rows), "audit3-fixes")
