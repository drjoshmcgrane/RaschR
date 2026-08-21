# STUDY: frame-invariance-power
#
# Frequentist performance of the conditional frame-invariance location test.
# The null cells assess standard-error calibration, coverage, per-item size and
# Holm familywise error. Departure cells report power for the two shifted
# items and the expected relative contrasts among the remaining items.
#
# Principal null cells use 1,000 replicates; power cells use 400. Set
# SV_FRAME_NULL_REPS or SV_FRAME_POWER_REPS only for development smoke runs.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "frame-invariance-power"
rows <- list()
t0 <- Sys.time()

ratio <- 1.40
u <- c(ratio^-0.5, ratio^0.5)
K <- 8L
delta <- seq(-1.5, 1.5, length.out = K)
items <- sprintf("I%02d", seq_len(K))
HIT <- c(3L, 6L)
R_NULL <- as.integer(Sys.getenv("SV_FRAME_NULL_REPS", "1000"))
R_POWER <- as.integer(Sys.getenv("SV_FRAME_POWER_REPS", "400"))
N_CORES <- as.integer(Sys.getenv("SV_CORES", "8"))
if (!is.finite(N_CORES) || N_CORES < 1L) N_CORES <- 1L

run_reps <- function(R, fun) {
  if (.Platform$OS.type == "windows" || N_CORES == 1L)
    return(lapply(seq_len(R), fun))
  parallel::mclapply(seq_len(R), fun, mc.cores = min(N_CORES, R),
                     mc.preschedule = FALSE, mc.set.seed = FALSE)
}

gen <- function(seed, shift, N) {
  set.seed(seed)
  one <- function(g, sh) {
    th <- stats::rnorm(N, 0, 1.3)
    X <- vapply(seq_len(K), function(i)
      stats::rbinom(N, 1, stats::plogis(
        u[g] * (th - delta[i] - sh[i]))), numeric(N))
    colnames(X) <- items
    X
  }
  X <- rbind(one(1, rep(0, K)), one(2, shift))
  data.frame(id = sprintf("P%05d", seq_len(2L * N)), X,
             group = rep(c("g1", "g2"), each = N), check.names = FALSE)
}

run_cell <- function(label, shift, N, target = c("null", "location")) {
  target <- match.arg(target)
  is_null <- target == "null"
  R <- if (is_null) R_NULL else R_POWER
  est <- se <- p_raw <- p_holm <- matrix(
    NA_real_, R, K, dimnames = list(NULL, items))
  planted <- items[HIT]
  sound <- setdiff(items, planted)
  one_rep <- function(r) {
    f <- tryCatch(rasch_efrm(
      gen(800000 + r, shift, N), items = items,
      item_sets = list(set1 = items), groups = "group", id = "id",
      boot_reps = 0), error = function(e) NULL)
    if (is.null(f)) return(list(status = "refused"))
    if (!isTRUE(f$est$converged)) return(list(status = "nonconv"))
    inv <- tryCatch(frame_invariance(f, se_method = "conditional"),
                    error = function(e) NULL)
    if (is.null(inv) || any(!is.finite(c(
      inv$locations$difference, inv$locations$se,
      inv$locations$p, inv$locations$p_adj))))
      return(list(status = "refused"))
    list(status = "analysed", locations = inv$locations)
  }
  rr <- run_reps(R, one_rep)
  status <- vapply(rr, `[[`, "", "status")
  n_ref <- sum(status == "refused")
  n_nonconv <- sum(status == "nonconv")
  for (r in which(status == "analysed")) {
    j <- match(rr[[r]]$locations$item, items)
    est[r, j] <- rr[[r]]$locations$difference
    se[r, j] <- rr[[r]]$locations$se
    p_raw[r, j] <- rr[[r]]$locations$p
    p_holm[r, j] <- rr[[r]]$locations$p_adj
  }
  ok <- stats::complete.cases(est, se, p_raw, p_holm)
  common <- list(
    study = STUDY,
    scenario = sprintf("%s | N = %d", label, N),
    n_reps = sum(ok), n_attempted = R, n_refused = n_ref,
    n_nonconv = n_nonconv)
  put <- function(quantity, ..., mc_override = list(), notes = "")
    rows[[length(rows) + 1L]] <<- do.call(
      sv_row, c(common, list(quantity = quantity,
                            mc_override = mc_override, notes = notes),
                list(...)))
  if (!any(ok)) return(invisible(NULL))

  if (is_null) {
    emp <- sqrt(mean(apply(est[ok, , drop = FALSE], 2, stats::var)))
    mse <- mean(se[ok, , drop = FALSE])
    cover <- rowMeans(abs(est[ok, , drop = FALSE]) <=
                        1.96 * se[ok, , drop = FALSE])
    reject <- rowMeans(p_raw[ok, , drop = FALSE] < 0.05)
    family <- rowSums(p_holm[ok, , drop = FALSE] < 0.05) > 0
    put("location SE calibration", emp_sd = emp, mean_se = mse,
        se_ratio = emp / mse)
    put("location per-item coverage", coverage95 = mean(cover),
        mc_override = list(coverage95 = stats::sd(cover) / sqrt(sum(ok))))
    put("location raw per-item type I error", type1 = mean(reject),
        mc_override = list(type1 = stats::sd(reject) / sqrt(sum(ok))))
    put("location Holm familywise type I error", familywise = mean(family))
  } else {
    for (rule in c("raw", "Holm")) {
      P <- if (rule == "raw") p_raw else p_holm
      M <- P[ok, , drop = FALSE] < 0.05
      hit <- rowMeans(M[, planted, drop = FALSE])
      spill <- rowMeans(M[, sound, drop = FALSE])
      put(paste(rule, "per-planted-item power"), power = mean(hit),
          mc_override = list(power = stats::sd(hit) / sqrt(sum(ok))))
      put(paste(rule, "any-planted-item power"),
          power = mean(rowSums(M[, planted, drop = FALSE]) > 0))
      put(paste(rule, "per-unplanted-item relative-contrast flag rate"),
          power = mean(spill),
          mc_override = list(power = stats::sd(spill) / sqrt(sum(ok))),
          notes = paste("Separate frame origins centre the item differences;",
                        "unplanted items therefore have non-zero relative",
                        "contrasts under concentrated DIF."))
    }
  }
  cat(sprintf("[%s] %-34s N=%-5d analysed=%d/%d\n",
              format(Sys.time(), "%H:%M"), label, N, sum(ok), R))
  invisible(NULL)
}

none <- rep(0, K)
shift_1 <- replace(none, HIT, 1.0)
for (N in c(500L, 1000L, 2000L)) {
  run_cell("null: unit ratio only", none, N, "null")
  run_cell("2 items shifted 1 logit", shift_1, N, "location")
}

sv_write(do.call(rbind, rows), "frame-invariance-power")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
