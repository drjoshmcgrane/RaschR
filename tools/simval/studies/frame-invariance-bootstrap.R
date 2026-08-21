# STUDY: frame-invariance-bootstrap
#
# Benchmark the person-within-frame bootstrap used by frame_invariance().
# A principal null cell records standard-error calibration, interval coverage,
# per-item size and combined-family error. Two departure cells record power
# and off-channel error. Set the environment variables only for smoke runs.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "frame-invariance-bootstrap"
rows <- list()
t0 <- Sys.time()

K <- 8L
items <- sprintf("I%02d", seq_len(K))
delta <- seq(-1.5, 1.5, length.out = K)
HIT <- c(3L, 6L)
ratio <- 1.40
unit <- c(ratio^-0.5, ratio^0.5)
R_NULL <- as.integer(Sys.getenv("SV_FRAME_BOOT_NULL_REPS", "300"))
R_POWER <- as.integer(Sys.getenv("SV_FRAME_BOOT_POWER_REPS", "120"))
B <- as.integer(Sys.getenv("SV_FRAME_BOOT_REPS", "200"))
N_CORES <- as.integer(Sys.getenv("SV_CORES", "8"))
if (!is.finite(N_CORES) || N_CORES < 1L) N_CORES <- 1L

run_reps <- function(R, fun) {
  if (.Platform$OS.type == "windows" || N_CORES == 1L)
    return(lapply(seq_len(R), fun))
  parallel::mclapply(seq_len(R), fun, mc.cores = min(N_CORES, R),
                     mc.preschedule = FALSE, mc.set.seed = FALSE)
}

gen <- function(seed, shift, disc, N = 500L) {
  set.seed(seed)
  one <- function(g, sh, ds) {
    theta <- stats::rnorm(N, 0, 1.3)
    X <- vapply(seq_len(K), function(i)
      stats::rbinom(N, 1, stats::plogis(
        unit[g] * ds[i] * (theta - delta[i] - sh[i]))), numeric(N))
    colnames(X) <- items
    X
  }
  X <- rbind(one(1, rep(0, K), rep(1, K)), one(2, shift, disc))
  data.frame(id = sprintf("P%05d", seq_len(2L * N)), X,
             group = rep(c("g1", "g2"), each = N), check.names = FALSE)
}

run_cell <- function(label, shift, disc, target = c("null", "location",
                                                    "discrimination")) {
  target <- match.arg(target)
  R <- if (target == "null") R_NULL else R_POWER
  est_l <- se_l <- est_d <- se_d <- p_l <- pa_l <- p_d <- pa_d <-
    matrix(NA_real_, R, K, dimnames = list(NULL, items))
  one_rep <- function(r) {
    f <- tryCatch(rasch_efrm(
      gen(910000 + r, shift, disc), items = items,
      item_sets = list(set1 = items), groups = "group", id = "id",
      boot_reps = 0), error = function(e) NULL)
    if (is.null(f)) return(list(status = "refused"))
    if (!isTRUE(f$est$converged)) return(list(status = "nonconv"))
    z <- tryCatch(frame_invariance(
      f, se_method = "bootstrap", boot_reps = B, seed = 1200000 + r),
      error = function(e) NULL)
    if (is.null(z)) return(list(status = "refused"))
    if (any(!is.finite(c(
      z$locations$difference, z$locations$se, z$locations$p,
      z$locations$p_adj, log(z$discrimination$disc_ratio),
      z$discrimination$se_log_disc_ratio, z$discrimination$p,
      z$discrimination$p_adj))))
      return(list(status = "refused"))
    list(status = "analysed", locations = z$locations,
         discrimination = z$discrimination)
  }
  rr <- run_reps(R, one_rep)
  status <- vapply(rr, `[[`, "", "status")
  n_ref <- sum(status == "refused")
  n_nonconv <- sum(status == "nonconv")
  for (r in which(status == "analysed")) {
    j <- match(rr[[r]]$locations$item, items)
    est_l[r, j] <- rr[[r]]$locations$difference
    se_l[r, j] <- rr[[r]]$locations$se
    p_l[r, j] <- rr[[r]]$locations$p
    pa_l[r, j] <- rr[[r]]$locations$p_adj
    j <- match(rr[[r]]$discrimination$item, items)
    est_d[r, j] <- log(rr[[r]]$discrimination$disc_ratio)
    se_d[r, j] <- rr[[r]]$discrimination$se_log_disc_ratio
    p_d[r, j] <- rr[[r]]$discrimination$p
    pa_d[r, j] <- rr[[r]]$discrimination$p_adj
  }
  cat(sprintf("[%s] %-28s analysed=%d/%d\n", format(Sys.time(), "%H:%M"),
              label, sum(status == "analysed"), R))
  ok <- stats::complete.cases(est_l, se_l, est_d, se_d, p_l, pa_l, p_d, pa_d)
  common <- list(study = STUDY, scenario = label, n_reps = sum(ok),
                 n_attempted = R, n_refused = n_ref, n_nonconv = n_nonconv)
  put <- function(quantity, ..., mc_override = list(), notes = "")
    rows[[length(rows) + 1L]] <<- do.call(
      sv_row, c(common, list(quantity = quantity,
                            mc_override = mc_override, notes = notes),
                list(...)))
  if (!any(ok)) return(invisible(NULL))

  if (target == "null") {
    for (channel in c("location", "discrimination")) {
      E <- if (channel == "location") est_l else est_d
      S <- if (channel == "location") se_l else se_d
      P <- if (channel == "location") p_l else p_d
      emp <- sqrt(mean(apply(E[ok, , drop = FALSE], 2, stats::var)))
      mse <- mean(S[ok, , drop = FALSE])
      cover <- rowMeans(abs(E[ok, , drop = FALSE]) <=
                          1.96 * S[ok, , drop = FALSE])
      reject <- rowMeans(P[ok, , drop = FALSE] < 0.05)
      put(paste(channel, "SE calibration"), emp_sd = emp, mean_se = mse,
          se_ratio = emp / mse)
      put(paste(channel, "per-item coverage"), coverage95 = mean(cover),
          mc_override = list(coverage95 = stats::sd(cover) / sqrt(sum(ok))))
      put(paste(channel, "raw per-item type I error"), type1 = mean(reject),
          mc_override = list(type1 = stats::sd(reject) / sqrt(sum(ok))))
    }
    fw <- rowSums(cbind(pa_l[ok, , drop = FALSE] < 0.05,
                        pa_d[ok, , drop = FALSE] < 0.05)) > 0
    put("combined Holm familywise type I error", familywise = mean(fw))
  } else {
    planted <- items[HIT]
    sound <- setdiff(items, planted)
    for (channel in c("location", "discrimination")) {
      P <- if (channel == "location") pa_l else pa_d
      M <- P[ok, , drop = FALSE] < 0.05
      if (channel == target) {
        hit <- rowMeans(M[, planted, drop = FALSE])
        put(paste(channel, "Holm per-planted-item power"), power = mean(hit),
            mc_override = list(power = stats::sd(hit) / sqrt(sum(ok))))
        other <- rowMeans(M[, sound, drop = FALSE])
        if (channel == "location")
          put("location Holm per-unplanted-item spillover flag rate",
              power = mean(other),
              mc_override = list(power = stats::sd(other) / sqrt(sum(ok))),
              notes = "Centred frame contrasts give unplanted items non-zero relative differences under concentrated DIF.")
        else
          put("discrimination Holm per-unplanted-item false-positive rate",
              type1 = mean(other),
              mc_override = list(type1 = stats::sd(other) / sqrt(sum(ok))))
      } else {
        per <- rowMeans(M)
        put(paste(channel, "Holm per-item off-channel type I error"),
            type1 = mean(per),
            mc_override = list(type1 = stats::sd(per) / sqrt(sum(ok))))
      }
    }
  }
  invisible(NULL)
}

flat <- rep(1, K)
none <- rep(0, K)
run_cell("null | N = 500 per frame", none, flat, "null")
run_cell("2 items shifted 1 logit | N = 500 per frame",
         replace(none, HIT, 1), flat, "location")
run_cell("2 items 1.5x as steep | N = 500 per frame",
         none, replace(flat, HIT, 1.5), "discrimination")

sv_write(do.call(rbind, rows), "frame-invariance-bootstrap")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
