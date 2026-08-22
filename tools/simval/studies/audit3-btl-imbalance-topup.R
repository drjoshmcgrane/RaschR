# Fresh-seed top-up for the mildly imbalanced BTL-DIF pairwise cell found in
# audit3-fixes.R. Both factor levels have ten judges; one judge in each carries
# twice the ordinary comparison allocation (9.3 effective judges per level).

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

n_reps <- 2000L
K <- 6L
J <- 20L
objects <- sprintf("O%d", seq_len(K))
judges <- sprintf("J%d", seq_len(J))
group <- setNames(rep(c("g1", "g2"), each = 10L), judges)
workload <- c(2L, rep(1L, 9L), 2L, rep(1L, 9L))
object_pairs <- t(utils::combn(objects, 2))
lev_eff <- vapply(c("g1", "g2"), function(g) {
  ww <- workload[group == g]
  sum(ww)^2 / sum(ww^2)
}, 0)
df_pair <- sum(lev_eff) - 2
df_min <- min(lev_eff) - 1

reject <- reject_min <- rep(NA, n_reps)
n_refused <- n_nonconv <- 0L
for (r in seq_len(n_reps)) {
  set.seed(3.36e7 + r)
  judge_effect <- matrix(stats::rnorm(J * K, 0, 0.8), J, K,
                         dimnames = list(judges, objects))
  dat <- do.call(rbind, lapply(seq_along(judges), function(jj) {
    j <- judges[jj]
    d <- data.frame(
      object_a = rep(object_pairs[, 1], each = 3 * workload[jj]),
      object_b = rep(object_pairs[, 2], each = 3 * workload[jj]),
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
  a2 <- ifelse(cm$object_a == "O3", paste0("O3 (", cell, ")"), cm$object_a)
  b2 <- ifelse(cm$object_b == "O3", paste0("O3 (", cell, ")"), cm$object_b)
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
  idx <- match(c("O3 (g1)", "O3 (g2)"), resolved$objects$object)
  if (anyNA(idx) || anyNA(resolved$cov_beta[idx, idx])) {
    n_refused <- n_refused + 1L
    next
  }
  loc <- resolved$objects$location[idx]
  V <- resolved$cov_beta[idx, idx, drop = FALSE]
  se <- sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2])
  stat <- diff(loc) / se
  p <- 2 * stats::pt(-abs(stat), df = df_pair)
  reject[r] <- p < 0.05
  reject_min[r] <- 2 * stats::pt(-abs(stat), df = df_min) < 0.05
}

ok <- !is.na(reject)
rows <- list(sv_row(
  "audit3-btl-imbalance-topup",
  "BTL-DIF 2 x 10 judges; one twofold workload per level",
  "fresh-seed pairwise Type I in imbalanced inference region",
  sum(ok), type1 = mean(reject[ok]),
  n_attempted = n_reps, n_refused = n_refused,
  n_nonconv = n_nonconv,
  notes = sprintf("%.2f effective judges per level; pair df %.2f",
                  lev_eff[1], df_pair)),
  sv_row(
    "audit3-btl-imbalance-topup",
    "BTL-DIF 2 x 10 judges; one twofold workload per level",
    "diagnostic minimum-cell-df Type I",
    sum(ok), type1 = mean(reject_min[ok]),
    n_attempted = n_reps, n_refused = n_refused,
    n_nonconv = n_nonconv,
    notes = sprintf("diagnostic comparison only: df = minimum effective judges minus 1 = %.2f",
                    df_min)))
sv_write(do.call(rbind, rows), "audit3-btl-imbalance-topup")
