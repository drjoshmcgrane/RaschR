suppressWarnings(pkgload::load_all(".", quiet = TRUE))

# ---- custom DIF generator (simulate_mfrm has no `dif=` argument) ----------
# Mirrors simulate_mfrm()'s generation code (R/simulate.R), adding a 2-level
# person group and a location-shift DIF on one named item for group "B".
sim_mfrm_dif <- function(n_persons = 200, n_items = 6, n_raters = 6,
                         n_categories = 4, theta_sd = 1.2, item_sd = 1,
                         rater_severity_sd = 0.6, dif_item = "I3",
                         dif_shift = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- n_persons; I <- n_items; R <- n_raters; m <- n_categories - 1L
  pids <- sprintf("P%03d", seq_len(N)); iids <- sprintf("I%d", seq_len(I))
  rids <- sprintf("R%d", seq_len(R))
  group <- sample(rep(c("A", "B"), length.out = N))
  theta <- .sim_theta(N, 0, theta_sd)
  delta <- setNames(seq(-item_sd, item_sd, length.out = I), iids)
  lambda <- setNames(as.numeric(scale(stats::rnorm(R))) * rater_severity_sd, rids)
  base_tau <- .sim_thresholds(0, m, 1.2)
  grid <- expand.grid(p = seq_len(N), i = seq_len(I), r = seq_len(R))
  score <- integer(nrow(grid))
  for (i in seq_len(I)) for (r in seq_len(R)) {
    rows <- which(grid$i == i & grid$r == r)
    if (iids[i] == dif_item) {
      rowsA <- rows[group[grid$p[rows]] == "A"]
      rowsB <- rows[group[grid$p[rows]] == "B"]
      score[rowsA] <- .sim_item(theta[grid$p[rowsA]], base_tau + delta[i] + lambda[r])
      score[rowsB] <- .sim_item(theta[grid$p[rowsB]], base_tau + delta[i] + dif_shift + lambda[r])
    } else {
      score[rows] <- .sim_item(theta[grid$p[rows]], base_tau + delta[i] + lambda[r])
    }
  }
  d <- data.frame(person = pids[grid$p], item = iids[grid$i], rater = rids[grid$r],
                  score = score, group = group[grid$p], stringsAsFactors = FALSE)
  attr(d, "truth") <- list(theta = theta, difficulty = delta, severity = lambda,
                           base_tau = base_tau, dif_item = dif_item,
                           dif_shift = dif_shift, group = group)
  d
}

# ---- true item thresholds implied by simulate_mfrm()'s generation ---------
# .sim_thresholds(delta = 0, m, spread = 1.2) is DETERMINISTIC (no RNG):
# step <- seq(-1, 1, length.out = m) * 1.2; tau <- step - mean(step) = step
true_base_tau <- function(n_categories) {
  m <- n_categories - 1L
  seq(-1, 1, length.out = m) * 1.2
}

# ---- connected-incomplete judging plan -------------------------------------
# Each person is rated by a random subset of `k` raters (out of n_raters),
# on all items -- a common linked judging-plan pattern. With enough persons
# and overlap this is (almost surely) connected.
make_connected_incomplete <- function(d, rater_col = "rater", person_col = "person",
                                      k = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  raters <- sort(unique(d[[rater_col]]))
  persons <- unique(d[[person_col]])
  assign <- lapply(persons, function(p) sample(raters, k))
  names(assign) <- persons
  keep <- mapply(function(p, r) r %in% assign[[p]], d[[person_col]], d[[rater_col]])
  d[keep, , drop = FALSE]
}

# ---- disconnected judging plan (two non-overlapping blocks) ---------------
make_disconnected <- function(d, rater_col = "rater", person_col = "person") {
  raters <- sort(unique(d[[rater_col]]))
  persons <- sort(unique(d[[person_col]]))
  half_r <- ceiling(length(raters) / 2); half_p <- ceiling(length(persons) / 2)
  blockA_r <- raters[seq_len(half_r)]; blockB_r <- raters[(half_r + 1):length(raters)]
  blockA_p <- persons[seq_len(half_p)]; blockB_p <- persons[(half_p + 1):length(persons)]
  keepA <- d[[person_col]] %in% blockA_p & d[[rater_col]] %in% blockA_r
  keepB <- d[[person_col]] %in% blockB_p & d[[rater_col]] %in% blockB_r
  d[keepA | keepB, , drop = FALSE]
}

mc_se <- function(p, n) sqrt(p * (1 - p) / n)
