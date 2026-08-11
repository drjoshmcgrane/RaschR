# Simulation validation: judge-clustered inference for the BTL comparative-
# judgement model (R/btl.R).
#
# btl(..., judge = ) clusters the Godambe sandwich by judge (crossprod of
# per-judge score sums), applies a CR1 finite-cluster correction, and -- per
# the withholding guard (since 1.14.2: at least 10 judges AND at least 8
# EFFECTIVE judges by the inverse Simpson index of comparison shares, with
# more effective judges than parameters) -- WITHHOLDS covariance-based
# inference (object SEs all NA, a note attached) whenever there are fewer
# than 10 judge clusters or no more clusters than fitted parameters, while
# still reporting point estimates and fit summaries. This study asks: (a) is
# the Type I error of pairwise object-contrast tests near .05 once >=10
# judges are present, (b) do clustered 95% CIs on pairwise contrasts cover
# at the nominal rate, (c) how much worse is the naive (unclustered,
# design-based) SE by comparison, and (d) exactly how sharp is the J=9/J=10
# withholding boundary, and how do (a)/(b) compare at J=10 vs J=50.
#
# DESIGN. 8 objects, evenly spaced (object_sd = 1, i.e. beta = scale(1:8));
# 28 pairs, 40 comparisons/pair/replicate (1,120 rows). Two truths per
# replicate, same design and same generating beta/dependence realisation:
#   - FLAT truth (beta = 0 for all 8 objects): every one of the 28 pairwise
#     nulls (object_i - object_j = 0) is TRUE. Feeds Type I error (a)/(c)
#     and familywise error (28 non-independent tests/replicate).
#   - SPACED truth (beta = scale(1:8)): feeds 95% CI coverage (b) pooled
#     over all 28 pairs (clustered vs naive SE), plus power at two
#     contrasts -- ADJACENT (O4 vs O5, |true diff| ~ 0.41 logits) and
#     EXTREME (O1 vs O8, |true diff| ~ 2.86 logits) -- testing H0: diff = 0
#     (false under this truth) with the SAME clustered fit used for
#     coverage, so no extra simulation is needed for power.
# Object-contrast tests use the package's own few-cluster convention (see
# the `dependence` effects table in R/btl.R): clustered tests/CIs use a
# Student-t reference on df = n_judge_clusters - 1; the naive (no `judge =`)
# fit's tests/CIs use the normal reference (df = Inf), matching how the
# package itself treats an unclustered fit as descriptive/independent-rows.
#
# WITHIN-JUDGE DEPENDENCE. Induced by judge-by-object random effects: for
# judge j and object k, eta[j,k] ~ N(0, dep_sd^2), added to both objects'
# logits in every comparison judge j makes involving object k (so a judge's
# repeated comparisons of the same object share a persistent idiosyncratic
# bias -- the textbook source of within-cluster correlation the sandwich is
# meant to absorb). eta is REDRAWN every replicate (judges, like persons in
# the person-facing studies, are resampled each replicate; the fixed
# generating parameters held constant across replicates are the object
# locations beta and the population dependence SD, not one particular
# eta draw). dep_sd in {0 (none), 0.3 (moderate), 0.6 (strong)}.
#
# ALLOCATION. balanced: judge = sample(judges, n, replace = TRUE) (uniform).
# skewed: one judge gets probability 0.5, the rest share the remaining 0.5
# uniformly (~half of all comparisons go to a single judge).
#
# METHODOLOGY. beta (flat and spaced) and the object/pair design are fixed;
# per replicate only eta (judge idiosyncrasies), judge assignment, and
# comparison outcomes are redrawn. A btl() call that errors (stop()) is a
# REFUSAL; one that fits but reports converged = FALSE is a NONCONV; both
# are excluded from the n_reps quantities are computed on and reported via
# n_attempted/n_refused/n_nonconv (harness rates).
#
# BENCHMARK (this machine, single core): one full replicate -- flat
# truth x {clustered, naive} + spaced truth x {clustered, naive}, 4 btl()
# fits, 1,120 rows each -- costs ~0.30s regardless of J/allocation/dep_sd
# (J in {9,50}, both allocations, dep_sd in {0, 0.6} all benchmarked within
# 0.08-0.40s/fit). At 0.30s/replicate the >=1000-replicate rule for the
# principal claim applies with margin to spare.
#
# REPLICATE BUDGET (~70 min wall clock, <=2 concurrent processes):
#   principal  J=10, balanced, dep_sd=0        : 1,200 reps (~360s)
#   boundary   J=9,  balanced, dep_sd=0        :   500 reps (~150s)
#   secondary  4 J x 2 alloc x 3 dep, minus the
#              principal cell already counted   :  23 x 700 reps (~4,830s)
# Total ~5,340 CPU-s (~89 CPU-min); split into two ~45-min concurrent chunks
# (odd/even scenario ids) comfortably inside the 70-minute budget.
#
# Run from the package root:
#   Rscript tools/simval/studies/btl-clustered.R 1        # chunk 1 (odd ids)
#   Rscript tools/simval/studies/btl-clustered.R 2        # chunk 2 (even ids)
#   Rscript tools/simval/studies/btl-clustered.R combine   # stitch chunk CSVs
#   Rscript tools/simval/studies/btl-clustered.R           # everything, serially
#
# The committed results CSV (tools/simval/results/btl-clustered.csv) is the
# FULL design (SV_SCALE unset: principal 1,200 replicates, boundary 500,
# secondary 700 per cell), re-run after the 1.14.2 corrections with
# coverage targeting the population-averaged contrasts. (An earlier
# agent-produced version used SV_SCALE=0.15 after a session time limit;
# it is superseded.) Every reported rate carries its own MC SE computed on
# the actual n_reps used.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

SCALE <- as.numeric(Sys.getenv("SV_SCALE", "1"))   # replicate-count scale, for smoke tests
nrep <- function(n) max(1L, round(n * SCALE))

args <- commandArgs(trailingOnly = TRUE)
chunk <- if (length(args) >= 1) args[1] else "all"

STUDY <- "btl-clustered"

if (identical(chunk, "combine")) {
  f1 <- file.path("tools/simval/results", paste0(STUDY, "_chunk1.csv"))
  f2 <- file.path("tools/simval/results", paste0(STUDY, "_chunk2.csv"))
  r1 <- utils::read.csv(f1, stringsAsFactors = FALSE)
  r2 <- utils::read.csv(f2, stringsAsFactors = FALSE)
  sv_write(rbind(r1, r2), STUDY)
  quit(save = "no")
}

## ---- fixed design -----------------------------------------------------

K <- 8L
objs <- sprintf("O%d", seq_len(K))
beta_spaced <- setNames(as.numeric(scale(seq_len(K))) * 1, objs)
beta_flat <- setNames(rep(0, K), objs)
REPS_PER_PAIR <- 40L

pair_idx <- utils::combn(K, 2)               # 2 x 28, columns i<j
n_pairs <- ncol(pair_idx)
idx_small <- which(objs[pair_idx[1, ]] == "O4" & objs[pair_idx[2, ]] == "O5")
idx_large <- which(objs[pair_idx[1, ]] == "O1" & objs[pair_idx[2, ]] == "O8")
true_diff_flat <- rep(0, n_pairs)
true_diff_spaced <- beta_spaced[pair_idx[1, ]] - beta_spaced[pair_idx[2, ]]

# Population-averaged targets under judge-by-object heterogeneity. With
# eta_jk ~ N(0, dep_sd^2) added to both objects' logits, the marginal
# win probability of each pair is a logistic-normal integral, and the
# marginal (no-random-effects) BTL fit converges to the KL projection of
# those probabilities onto the BTL family -- the population-averaged
# contrasts, attenuated relative to the conditional generating beta
# (about 11% at dep_sd = 0.6 for the extreme pair). Coverage of the
# clustered CIs must target THIS estimand; targeting the conditional
# beta conflates non-collapsibility with covariance error.
marginal_diff_targets <- function(beta, dep_sd) {
  if (dep_sd == 0) return(beta[pair_idx[1, ]] - beta[pair_idx[2, ]])
  pm <- vapply(seq_len(n_pairs), function(e) {
    d <- beta[pair_idx[1, e]] - beta[pair_idx[2, e]]
    stats::integrate(function(u) stats::plogis(d + u) *
                       stats::dnorm(u, 0, sqrt(2) * dep_sd),
                     -Inf, Inf)$value
  }, 0)
  X <- matrix(0, n_pairs, K)
  for (e in seq_len(n_pairs)) {          # pair_idx holds numeric indices
    X[e, pair_idx[1, e]] <- 1
    X[e, pair_idx[2, e]] <- -1
  }
  Cm <- rbind(diag(K - 1L), rep(-1, K - 1L))     # sum-to-zero
  Xr <- X %*% Cm
  b <- rep(0, K - 1L)
  for (it in 1:200) {
    pr <- stats::plogis(drop(Xr %*% b))
    g <- crossprod(Xr, pm - pr)
    H <- crossprod(Xr, Xr * (pr * (1 - pr)))
    step <- solve(H, g)
    b <- b + step
    if (max(abs(step)) < 1e-12) break
  }
  bm <- drop(Cm %*% b)
  bm[pair_idx[1, ]] - bm[pair_idx[2, ]]
}

gen_data <- function(beta, J, allocation, dep_sd) {
  jids <- sprintf("J%d", seq_len(J))
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(object_a = rep(pr[, 1], each = REPS_PER_PAIR),
                  object_b = rep(pr[, 2], each = REPS_PER_PAIR),
                  stringsAsFactors = FALSE)
  n <- nrow(d)
  if (allocation == "balanced") {
    d$judge <- sample(jids, n, replace = TRUE)
  } else {
    prob <- c(0.5, rep(0.5 / (J - 1), J - 1))
    d$judge <- sample(jids, n, replace = TRUE, prob = prob)
  }
  eta <- if (dep_sd > 0)
    matrix(stats::rnorm(J * K, 0, dep_sd), J, K, dimnames = list(jids, objs))
  else matrix(0, J, K, dimnames = list(jids, objs))
  lp <- (beta[d$object_a] + eta[cbind(d$judge, d$object_a)]) -
        (beta[d$object_b] + eta[cbind(d$judge, d$object_b)])
  y <- stats::rbinom(n, 1, stats::plogis(lp))
  d$winner <- ifelse(y == 1, d$object_a, d$object_b)
  d
}

# fit + extract every pairwise contrast's estimate/SE/df in one shot;
# NULL on a refusal (stop()), $converged FALSE on non-convergence, SEs NA
# when clustered inference is withheld (few-cluster guard)
fit_pairs <- function(d, clustered) {
  f <- tryCatch(
    if (clustered) btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
    else btl(d, "object_a", "object_b", winner = "winner"),
    error = function(e) e)
  if (inherits(f, "error")) return(list(refused = TRUE, msg = conditionMessage(f)))
  loc <- f$objects$location
  diffs <- loc[pair_idx[1, ]] - loc[pair_idx[2, ]]
  cb <- f$cov_beta
  var_ij <- diag(cb)[pair_idx[1, ]] + diag(cb)[pair_idx[2, ]] -
    2 * cb[cbind(pair_idx[1, ], pair_idx[2, ])]
  se <- sqrt(pmax(var_ij, 0))
  nc <- if (clustered) f$cl$n_units else Inf
  df <- if (clustered) max(nc - 1, 1) else Inf
  list(refused = FALSE, converged = isTRUE(f$converged),
       inference_available = isTRUE(f$cl$inference_available),
       diffs = diffs, se = se, df = df, nc = nc, notes = f$notes)
}

## ---- one replicate: 4 fits (flat/spaced x clustered/naive) ------------

one_rep <- function(J, allocation, dep_sd) {
  d_flat <- gen_data(beta_flat, J, allocation, dep_sd)
  d_sp   <- gen_data(beta_spaced, J, allocation, dep_sd)
  list(flat_c = fit_pairs(d_flat, TRUE), flat_n = fit_pairs(d_flat, FALSE),
       sp_c   = fit_pairs(d_sp, TRUE),   sp_n   = fit_pairs(d_sp, FALSE))
}

## ---- scenario table (interleaved for a roughly balanced 2-chunk split) --

principal <- data.frame(J = 10L, allocation = "balanced", dep_sd = 0,
                        role = "principal", reps = nrep(1200), stringsAsFactors = FALSE)
boundary <- data.frame(J = 9L, allocation = "balanced", dep_sd = 0,
                       role = "boundary", reps = nrep(500), stringsAsFactors = FALSE)
grid <- expand.grid(J = c(10L, 12L, 20L, 50L), allocation = c("balanced", "skewed"),
                    dep_sd = c(0, 0.3, 0.6), stringsAsFactors = FALSE)
# drop the principal cell (already scheduled with more reps above)
grid <- grid[!(grid$J == 10L & grid$allocation == "balanced" & grid$dep_sd == 0), ]
secondary <- data.frame(grid, role = "secondary", reps = nrep(700), stringsAsFactors = FALSE)

scenarios <- rbind(principal, boundary, secondary)
scenarios$id <- seq_len(nrow(scenarios))
scenarios$label <- sprintf("J%02d_%s_dep%.1f", scenarios$J, scenarios$allocation, scenarios$dep_sd)

if (chunk == "1") scenarios <- scenarios[scenarios$id %% 2 == 1, ]
if (chunk == "2") scenarios <- scenarios[scenarios$id %% 2 == 0, ]

cat(sprintf("chunk=%s: %d scenarios, %d total replicates\n", chunk, nrow(scenarios), sum(scenarios$reps)))

## ---- run ----------------------------------------------------------------

rows <- list()
t_start <- Sys.time()

for (s in seq_len(nrow(scenarios))) {
  sc <- scenarios[s, ]
  seed <- 811000L + sc$id                 # scenario-specific, chunk-independent
  set.seed(seed)
  n_reps <- sc$reps

  # coverage/bias targets for this scenario: the population-averaged
  # contrasts (equal to the generating contrasts when dep_sd = 0)
  target_sp <- marginal_diff_targets(beta_spaced, sc$dep_sd)

  # accumulators
  rej_flat_c <- rej_flat_n <- matrix(NA, n_reps, n_pairs)  # Type I (flat truth)
  any_rej_flat_c <- logical(n_reps)                         # familywise (clustered)
  cov_sp_c <- cov_sp_n <- matrix(NA, n_reps, n_pairs)       # coverage (spaced truth)
  diff_small <- se_small <- diff_large <- se_large <- rep(NA_real_, n_reps)
  pow_small <- pow_large <- rep(NA, n_reps)
  infer_avail_c <- rep(NA, n_reps)                           # withholding guard state

  n_att <- n_reps
  n_ref_flat <- n_ref_sp <- 0L
  n_nc_flat <- n_nc_sp <- 0L

  for (r in seq_len(n_reps)) {
    rep_out <- tryCatch(one_rep(sc$J, sc$allocation, sc$dep_sd),
                        error = function(e) e)
    if (inherits(rep_out, "error")) { n_ref_flat <- n_ref_flat + 1L; n_ref_sp <- n_ref_sp + 1L; next }

    fc <- rep_out$flat_c; fn <- rep_out$flat_n
    sc_ <- rep_out$sp_c;  sn <- rep_out$sp_n

    if (isTRUE(fc$refused) || isTRUE(fn$refused)) n_ref_flat <- n_ref_flat + 1L
    if (isTRUE(sc_$refused) || isTRUE(sn$refused)) n_ref_sp <- n_ref_sp + 1L
    if (!isTRUE(fc$refused) && !isTRUE(fn$refused) && (!fc$converged || !fn$converged))
      n_nc_flat <- n_nc_flat + 1L
    if (!isTRUE(sc_$refused) && !isTRUE(sn$refused) && (!sc_$converged || !sn$converged))
      n_nc_sp <- n_nc_sp + 1L

    if (!isTRUE(fc$refused) && fc$converged && fc$inference_available) {
      crit <- stats::qt(0.975, df = fc$df)
      z <- fc$diffs / fc$se
      rej_flat_c[r, ] <- abs(z) > crit
      any_rej_flat_c[r] <- any(rej_flat_c[r, ], na.rm = TRUE)
    }
    if (!isTRUE(fn$refused) && fn$converged) {
      crit <- stats::qnorm(0.975)
      z <- fn$diffs / fn$se
      rej_flat_n[r, ] <- abs(z) > crit
    }

    infer_avail_c[r] <- if (!isTRUE(sc_$refused)) sc_$inference_available else NA

    if (!isTRUE(sc_$refused) && sc_$converged && sc_$inference_available) {
      crit <- stats::qt(0.975, df = sc_$df)
      cov_sp_c[r, ] <- abs(sc_$diffs - target_sp) <= crit * sc_$se
      diff_small[r] <- sc_$diffs[idx_small]; se_small[r] <- sc_$se[idx_small]
      diff_large[r] <- sc_$diffs[idx_large]; se_large[r] <- sc_$se[idx_large]
      pow_small[r] <- abs(sc_$diffs[idx_small] / sc_$se[idx_small]) > crit
      pow_large[r] <- abs(sc_$diffs[idx_large] / sc_$se[idx_large]) > crit
    }
    if (!isTRUE(sn$refused) && sn$converged) {
      crit <- stats::qnorm(0.975)
      cov_sp_n[r, ] <- abs(sn$diffs - target_sp) <= crit * sn$se
    }
  }

  # per-replicate proportions across the 28 (correlated) pairs -> cluster-
  # robust MC SE on pooled rates (sd of per-rep proportions / sqrt(n_reps))
  pooled_rate <- function(mat) {
    perrep <- rowMeans(mat, na.rm = TRUE)
    ok <- is.finite(perrep)
    p <- mean(perrep[ok]); mcse <- stats::sd(perrep[ok]) / sqrt(sum(ok))
    list(p = p, n = sum(ok), mcse = mcse)
  }

  ty1c <- pooled_rate(rej_flat_c); ty1n <- pooled_rate(rej_flat_n)
  covc <- pooled_rate(cov_sp_c);   covn <- pooled_rate(cov_sp_n)
  fw <- any_rej_flat_c[is.finite(rowMeans(rej_flat_c, na.rm = TRUE))]

  mk <- function(quantity, ...) sv_row(STUDY, sc$label, quantity, n_attempted = n_att, ...)

  rows[[length(rows) + 1]] <- mk("type1_clustered", n_reps = ty1c$n, type1 = ty1c$p,
    mc_override = list(type1 = ty1c$mcse), n_refused = n_ref_flat, n_nonconv = n_nc_flat,
    notes = "pooled over 28 pairs, flat truth (all true diffs = 0), t(df=nJudge-1) test")
  rows[[length(rows) + 1]] <- mk("type1_naive", n_reps = ty1n$n, type1 = ty1n$p,
    mc_override = list(type1 = ty1n$mcse), n_refused = n_ref_flat, n_nonconv = n_nc_flat,
    notes = "pooled over 28 pairs, flat truth, unclustered SE, normal(z) test")
  rows[[length(rows) + 1]] <- mk("familywise_clustered", n_reps = length(fw),
    familywise = mean(fw), mc_override = list(familywise = mc_se_prop(mean(fw), length(fw))),
    n_refused = n_ref_flat, n_nonconv = n_nc_flat,
    notes = "P(>=1 false rejection among 28 non-independent pairwise tests/replicate)")
  rows[[length(rows) + 1]] <- mk("coverage_clustered", n_reps = covc$n, coverage95 = covc$p,
    mc_override = list(coverage95 = covc$mcse), n_refused = n_ref_sp, n_nonconv = n_nc_sp,
    notes = "pooled over 28 pairs, spaced truth, clustered SE, t(df=nJudge-1) CI; target = population-averaged (marginal) contrasts, the estimand of the marginal BTL fit under judge heterogeneity")
  rows[[length(rows) + 1]] <- mk("coverage_naive", n_reps = covn$n, coverage95 = covn$p,
    mc_override = list(coverage95 = covn$mcse), n_refused = n_ref_sp, n_nonconv = n_nc_sp,
    notes = "pooled over 28 pairs, spaced truth, unclustered SE, normal(z) CI; same population-averaged targets")

  ok_s <- is.finite(diff_small); ok_l <- is.finite(diff_large)
  bs <- mean(diff_small[ok_s] - target_sp[idx_small])
  bl <- mean(diff_large[ok_l] - target_sp[idx_large])
  rows[[length(rows) + 1]] <- mk("contrast_small_O4_O5", n_reps = sum(ok_s),
    effect = target_sp[idx_small], bias = bs, emp_sd = stats::sd(diff_small[ok_s]),
    mean_se = mean(se_small[ok_s]), coverage95 = covc$p,
    power = mean(pow_small[ok_s]), mc_override = list(power = mc_se_prop(mean(pow_small[ok_s]), sum(ok_s))),
    n_refused = n_ref_sp, n_nonconv = n_nc_sp,
    notes = "adjacent-object contrast, |true diff| ~0.41 logits; coverage95 is the pooled 28-pair figure")
  rows[[length(rows) + 1]] <- mk("contrast_large_O1_O8", n_reps = sum(ok_l),
    effect = target_sp[idx_large], bias = bl, emp_sd = stats::sd(diff_large[ok_l]),
    mean_se = mean(se_large[ok_l]), coverage95 = covc$p,
    power = mean(pow_large[ok_l]), mc_override = list(power = mc_se_prop(mean(pow_large[ok_l]), sum(ok_l))),
    n_refused = n_ref_sp, n_nonconv = n_nc_sp,
    notes = "extreme-object contrast, |true diff| ~2.86 logits; coverage95 is the pooled 28-pair figure")

  ok_i <- is.finite(infer_avail_c)
  rows[[length(rows) + 1]] <- mk("clustered_inference_available", n_reps = sum(ok_i),
    coverage95 = mean(infer_avail_c[ok_i]),
    mc_override = list(coverage95 = mc_se_prop(mean(infer_avail_c[ok_i]), sum(ok_i))),
    n_refused = n_ref_sp, n_nonconv = n_nc_sp,
    notes = sprintf(paste("fraction of replicates where cluster_inference == TRUE",
                          "(J=%d, balanced-count expectation %s; since 1.14.2 the",
                          "guard also requires >=8 EFFECTIVE judges, so skewed",
                          "allocations withhold regardless of J)"),
                    sc$J, if (sc$J >= 10L) "1 for balanced allocations" else "0"))

  cat(sprintf("[%s] scenario %d/%d (%s, %d reps) done, elapsed %.1fs\n",
              format(Sys.time(), "%H:%M:%S"), s, nrow(scenarios), sc$label, n_reps,
              as.numeric(Sys.time() - t_start, units = "secs")))
}

out <- do.call(rbind, rows)
out_name <- if (chunk %in% c("1", "2")) paste0(STUDY, "_chunk", chunk) else STUDY
sv_write(out, out_name)
