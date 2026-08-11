# FIXED TRUTH: simulate_mfrm() redraws the rater severities from every
# seed, which would vary the truth across replicates and fold that
# variation into the replicate distribution. Here the truth -- severities
# (the only random component of simulate_mfrm's truth; its item locations
# and thresholds are deterministic grids), items, thresholds, and the
# planted interaction -- is drawn ONCE per rater count R (7.7e6 + R,
# shared across every N, null, and effect cell), and only
# the persons and their responses are redrawn each replicate, using the
# package's own categorical sampler (.sim_item, exposed by load_all).
mfrm_truth <- function(R, seed_truth) {
  set.seed(seed_truth)
  I <- 6L
  list(lambda = setNames(as.numeric(scale(stats::rnorm(R))) * 0.7,
                         sprintf("R%d", seq_len(R))),
       delta = setNames(seq(-1, 1, length.out = I), sprintf("I%d", seq_len(I))),
       base_tau = c(-1.2, 0, 1.2))          # .sim_thresholds(0, 3, 1.2)
}

mfrm_rep <- function(N, truth, bias = NA, item = "I3", rater = "R3", seed) {
  lambda <- truth$lambda; delta <- truth$delta; base_tau <- truth$base_tau
  I <- length(delta); R <- length(lambda)
  set.seed(seed)
  theta <- stats::rnorm(N, 0, 1.2)
  grid <- expand.grid(p = seq_len(N), i = seq_len(I), r = seq_len(R))
  score <- integer(nrow(grid))
  for (i in seq_len(I)) for (r in seq_len(R)) {
    sel <- grid$i == i & grid$r == r
    tau_ir <- base_tau + delta[i] + lambda[r] +
      if (!is.na(bias) && names(delta)[i] == item && names(lambda)[r] == rater)
        bias else 0
    score[sel] <- .sim_item(theta[grid$p[sel]], tau_ir)
  }
  d <- data.frame(person = sprintf("P%03d", grid$p), item = names(delta)[grid$i],
                  rater = names(lambda)[grid$r], score = score,
                  stringsAsFactors = FALSE)
  fit <- tryCatch(rasch_mfrm(d, person = "person", item = "item", score = "score",
                             facets = "rater", interaction = "rater"),
                   error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA, df = NA, refused = TRUE, nonconv = FALSE))
  if (!isTRUE(fit$est$converged))
    return(list(p = NA, df = NA, refused = FALSE, nonconv = TRUE))
  list(p = fit$interaction_test$p, df = fit$interaction_test$df,
       refused = FALSE, nonconv = FALSE)
}

mfrm_batch <- function(N, R, bias, n_reps, seed0) {
  # ONE truth per rater count R, shared across every N, null, and effect
  # cell (a per-scenario truth would let the planted rater's severity
  # differ between effect sizes, confounding the power comparison);
  # response seeds still differ per cell via seed0
  truth <- mfrm_truth(R, seed_truth = 7.7e6 + R)
  p <- rep(NA_real_, n_reps); n_refused <- 0L; n_nonconv <- 0L; dfq <- NA
  for (r in seq_len(n_reps)) {
    res <- mfrm_rep(N, truth, bias, seed = seed0 + r)
    if (res$refused) { n_refused <- n_refused + 1L; next }
    if (res$nonconv) { n_nonconv <- n_nonconv + 1L; next }
    p[r] <- res$p; dfq <- res$df
  }
  list(p = p[!is.na(p)], n_attempted = n_reps, n_refused = n_refused,
       n_nonconv = n_nonconv, df = dfq)
}
