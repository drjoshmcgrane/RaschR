# Anchored-estimation validation. External anchors are set to their generating
# values, so bias and coverage are conditional on an exact calibration bank.
pkgload::load_all(".", quiet = TRUE)
options(simval.script = "tools/simval/studies/anchored-estimation.R")
source("tools/simval/harness.R")

pcm_reps <- as.integer(Sys.getenv("RASCH_ANCHOR_PCM_REPS", "250"))
btl_reps <- as.integer(Sys.getenv("RASCH_ANCHOR_BTL_REPS", "500"))

draw_item <- function(theta, tau) {
  score <- 0:length(tau)
  vapply(theta, function(th) {
    lp <- score * th - c(0, cumsum(tau))
    p <- exp(lp - max(lp)); p <- p / sum(p)
    sample(score, 1L, prob = p)
  }, 0L)
}

pcm_one <- function(seed, disconnected = FALSE) {
  set.seed(seed)
  m <- c(1L, 2L, 3L, 2L, 1L, 3L, 2L, 3L)
  loc <- seq(-1.5, 1.5, length.out = length(m)) + 0.6
  tau <- Map(function(mi, di)
    di + if (mi == 1L) 0 else seq(-0.8, 0.8, length.out = mi), m, loc)
  if (!disconnected) {
    theta <- rnorm(500, 0.6, 1.2)
    X <- sapply(tau, function(tt) draw_item(theta, tt))
  } else {
    theta <- rnorm(600, 0.6, 1.2)
    X <- matrix(NA_integer_, 600, length(m))
    X[1:300, 1:4] <- sapply(tau[1:4], function(tt)
      draw_item(theta[1:300], tt))
    X[301:600, 5:8] <- sapply(tau[5:8], function(tt)
      draw_item(theta[301:600], tt))
  }
  colnames(X) <- paste0("I", seq_along(m))
  anchors <- if (disconnected)
    data.frame(item = c("I1", "I5"), k = 1,
               tau = c(tau[[1]][1], tau[[5]][1]))
  else data.frame(item = c("I1", "I4", "I7"), k = c(1, NA, 2),
                  tau = c(tau[[1]][1], mean(tau[[4]]), tau[[7]][2]))
  fit <- tryCatch(rasch(X, anchors = anchors, maxit = 100, tol = 1e-9),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(status = "error"))
  if (!isTRUE(fit$est$converged)) return(list(status = "nonconverged"))
  fixed_mean <- if (disconnected) c(1L, 5L) else c(1L, 4L)
  use <- setdiff(seq_along(m), fixed_mean)
  err <- fit$items$location[use] - loc[use]
  se <- fit$items$se[use]
  available <- is.finite(err) & is.finite(se) & se > 0
  anchor_error <- if (disconnected) {
    max(abs(fit$tau_list$I1[1] - tau[[1]][1]),
        abs(fit$tau_list$I5[1] - tau[[5]][1]))
  } else {
    max(abs(fit$tau_list$I1[1] - tau[[1]][1]),
        abs(mean(fit$tau_list$I4) - mean(tau[[4]])),
        abs(fit$tau_list$I7[2] - tau[[7]][2]))
  }
  list(status = "ok", item = use[available], error = err[available],
       se = se[available],
       covered = abs(err[available]) <= 1.96 * se[available],
       metric_missing = any(!available), anchor_error = anchor_error)
}

summarise_pcm <- function(disconnected) {
  n <- pcm_reps
  z <- lapply(seq_len(n), function(r)
    pcm_one(810000L + 10000L * disconnected + r, disconnected))
  status <- vapply(z, `[[`, "", "status")
  good <- status == "ok"
  detail <- do.call(rbind, lapply(which(good), function(i)
    data.frame(replicate = i, item = z[[i]]$item,
               error = z[[i]]$error, se = z[[i]]$se,
               covered = z[[i]]$covered)))
  by_item <- split(detail, detail$item)
  item_bias <- vapply(by_item, function(x) mean(x$error), 0)
  item_sd <- vapply(by_item, function(x) stats::sd(x$error), 0)
  item_se <- vapply(by_item, function(x) mean(x$se), 0)
  cov_rep <- vapply(z[good], function(x)
    if (length(x$covered)) mean(x$covered) else NA_real_, 0)
  cov_rep <- cov_rep[is.finite(cov_rep)]
  n_missing <- sum(vapply(z[good], `[[`, FALSE, "metric_missing"))
  anchor_error <- max(vapply(z[good], `[[`, 0, "anchor_error"))
  scenario <- if (disconnected)
    "PCM, two disconnected four-item blocks, one true anchor per block"
  else "PCM, mixed maximum scores, threshold and location anchors"
  rbind(
    sv_row("anchored-estimation", scenario, "free item-location recovery",
      n_reps = sum(good), bias = max(abs(item_bias)),
      emp_sd = mean(item_sd), mean_se = mean(item_se),
      se_ratio = mean(item_sd / item_se), coverage95 = mean(cov_rep),
      mc_override = list(coverage95 = stats::sd(cov_rep) / sqrt(sum(good))),
      n_attempted = n, n_nonconv = sum(status == "nonconverged"),
      n_error = sum(status == "error"), n_metric_unavailable = n_missing,
      notes = paste("bias is the maximum absolute item bias; the SE ratio is",
                    "the mean itemwise ratio; coverage and standard errors are conditional on exact",
                    "external anchor values; a replicate is marked metric-unavailable",
                    "when any reported free-item SE is withheld")),
    sv_row("anchored-estimation", scenario, "maximum absolute anchor error",
      n_reps = sum(good), bias = anchor_error,
      n_attempted = n, n_nonconv = sum(status == "nonconverged"),
      n_error = sum(status == "error"),
      notes = "fixed thresholds and location means should be exact to numerical tolerance")
  )
}

btl_one <- function(seed) {
  set.seed(seed)
  truth <- setNames(seq(-1.5, 1.5, length.out = 7), paste0("O", 1:7))
  pr <- t(utils::combn(names(truth), 2))
  d <- data.frame(a = rep(pr[, 1], each = 24),
                  b = rep(pr[, 2], each = 24))
  d$winner <- ifelse(runif(nrow(d)) < plogis(truth[d$a] - truth[d$b]),
                     d$a, d$b)
  fit <- tryCatch(btl(d, "a", "b", "winner",
                      anchors = truth[c("O1", "O7")]),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(status = "error"))
  if (!isTRUE(fit$converged)) return(list(status = "nonconverged"))
  use <- match(paste0("O", 2:6), fit$objects$object)
  err <- fit$objects$location[use] - truth[paste0("O", 2:6)]
  se <- fit$objects$se[use]
  list(status = "ok", item = paste0("O", 2:6), error = unname(err), se = se,
       covered = abs(err) <= 1.96 * se,
       anchor_error = max(abs(fit$objects$location[
         match(c("O1", "O7"), fit$objects$object)] - truth[c("O1", "O7")])))
}

z <- lapply(seq_len(btl_reps), function(r) btl_one(920000L + r))
status <- vapply(z, `[[`, "", "status"); good <- status == "ok"
detail <- do.call(rbind, lapply(which(good), function(i)
  data.frame(replicate = i, item = z[[i]]$item,
             error = z[[i]]$error, se = z[[i]]$se,
             covered = z[[i]]$covered)))
by_item <- split(detail, detail$item)
item_bias <- vapply(by_item, function(x) mean(x$error), 0)
item_sd <- vapply(by_item, function(x) stats::sd(x$error), 0)
item_se <- vapply(by_item, function(x) mean(x$se), 0)
cov_rep <- vapply(z[good], function(x) mean(x$covered), 0)
rows_btl <- rbind(
  sv_row("anchored-estimation", "BTL, seven objects, two true fixed anchors",
    "free object-location recovery", n_reps = sum(good),
    bias = max(abs(item_bias)), emp_sd = mean(item_sd), mean_se = mean(item_se),
    se_ratio = mean(item_sd / item_se),
    coverage95 = mean(cov_rep),
    mc_override = list(coverage95 = stats::sd(cov_rep) / sqrt(sum(good))),
    n_attempted = btl_reps, n_nonconv = sum(status == "nonconverged"),
    n_error = sum(status == "error"),
    notes = paste("bias is the maximum absolute object bias; the SE ratio is",
                  "the mean objectwise ratio; coverage and standard errors",
                  "are conditional on exact external anchor values")),
  sv_row("anchored-estimation", "BTL, seven objects, two true fixed anchors",
    "maximum absolute anchor error", n_reps = sum(good),
    bias = max(vapply(z[good], `[[`, 0, "anchor_error")),
    n_attempted = btl_reps, n_nonconv = sum(status == "nonconverged"),
    n_error = sum(status == "error"),
    notes = "fixed object locations should be exact to numerical tolerance")
)

sv_write(sv_bind_rows(summarise_pcm(FALSE), summarise_pcm(TRUE), rows_btl),
         "anchored-estimation")
