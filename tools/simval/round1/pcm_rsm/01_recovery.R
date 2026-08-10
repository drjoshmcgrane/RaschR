suppressWarnings(pkgload::load_all("."), quiet=TRUE))
options(width=140)

## ---------------------------------------------------------------
## helper: run one replicate, PCM or RSM, complete/MCAR/booklet, and
## return pooled threshold + location recovery stats plus one item's
## (location) estimate & se for SE calibration.
## ---------------------------------------------------------------

make_booklet <- function(d) {
  # 10 items: I01-I04 unique to form A, I05-I06 link (everyone), I07-I10 unique to form B
  N <- nrow(d)
  grp <- sample(c("A","B"), N, replace = TRUE)
  d2 <- d
  d2[grp == "A", c("I07","I08","I09","I10")] <- NA
  d2[grp == "B", c("I01","I02","I03","I04")] <- NA
  d2
}

run_one <- function(model, cond, seed, n_persons = 600, n_items = 10) {
  d <- simulate_rasch(n_persons = n_persons, n_items = n_items, model = model,
                       n_categories = 4, difficulty = c(-2, 2),
                       threshold_spread = 1.3, seed = seed)
  tr <- attr(d, "truth")
  if (cond == "mcar") {
    Xcols <- grep("^I", names(d))
    n <- length(Xcols) * nrow(d)
    idx <- sample(seq_len(nrow(d)), size = round(0.25 * nrow(d)))
    for (ic in Xcols) {
      rows <- sample(seq_len(nrow(d)), size = round(0.25 * nrow(d)))
      d[rows, ic] <- NA
    }
  } else if (cond == "booklet") {
    d <- make_booklet(d)
  }
  fit <- tryCatch(rasch(d, model = model), error = function(e) e)
  if (inherits(fit, "error")) return(list(error = conditionMessage(fit)))
  # rare booklet edge case: a subgroup never uses an item's extreme category,
  # so .prepare_X rescores it and the threshold count no longer matches the
  # full-category truth -- not a bug, just incomparable; skip the replicate
  if (nrow(fit$thresholds) != length(unlist(tr$thresholds)) ||
      nrow(fit$items) != length(tr$difficulty))
    return(list(error = "category rescored on a booklet subgroup (skipped, not a fit error)"))
  shift <- mean(unlist(tr$thresholds))
  tau_truth <- unlist(tr$thresholds) - shift
  tau_est <- fit$thresholds$tau
  loc_truth <- tr$difficulty - shift
  loc_est <- fit$items$location
  list(error = NULL,
       cor_tau = suppressWarnings(cor(tau_truth, tau_est)),
       rmse_tau = sqrt(mean((tau_truth - tau_est)^2)),
       cor_loc = suppressWarnings(cor(loc_truth, loc_est)),
       rmse_loc = sqrt(mean((loc_truth - loc_est)^2)),
       # focal item = middle item (I05), stable presence across booklet (link item)
       loc5_est = fit$items$location[fit$items$item == "I05"],
       loc5_se  = fit$items$se[fit$items$item == "I05"],
       loc5_truth = tr$difficulty["I05"] - shift,
       any_weak = any(fit$thresholds$weak, na.rm = TRUE))
}

summarise_batch <- function(model, cond, nrep, seed0) {
  res <- lapply(seq_len(nrep), function(r) run_one(model, cond, seed = seed0 + r))
  errs <- vapply(res, function(x) !is.null(x$error), TRUE)
  ok <- res[!errs]
  n_ok <- length(ok)
  cor_tau <- vapply(ok, `[[`, 0, "cor_tau")
  rmse_tau <- vapply(ok, `[[`, 0, "rmse_tau")
  cor_loc <- vapply(ok, `[[`, 0, "cor_loc")
  rmse_loc <- vapply(ok, `[[`, 0, "rmse_loc")
  loc5_est <- vapply(ok, `[[`, 0, "loc5_est")
  loc5_se <- vapply(ok, `[[`, 0, "loc5_se")
  loc5_truth <- vapply(ok, `[[`, 0, "loc5_truth")[1]
  any_weak <- vapply(ok, `[[`, TRUE, "any_weak")
  se_ratio <- sd(loc5_est) / mean(loc5_se)
  bias5 <- mean(loc5_est) - loc5_truth
  cat(sprintf(
    "[%s | %s] n_ok=%d/%d (errors: %s)\n  cor_tau=%.4f (MCse %.4f) rmse_tau=%.4f\n  cor_loc=%.4f (MCse %.4f) rmse_loc=%.4f\n  item I05: bias=%.4f empSD=%.4f meanSE=%.4f ratio(empSD/meanSE)=%.3f\n  any_weak_flag_rate=%.3f\n\n",
    model, cond, n_ok, nrep, sum(errs),
    mean(cor_tau), sd(cor_tau)/sqrt(n_ok), mean(rmse_tau),
    mean(cor_loc), sd(cor_loc)/sqrt(n_ok), mean(rmse_loc),
    bias5, sd(loc5_est), mean(loc5_se), se_ratio,
    mean(any_weak)))
  invisible(list(model=model, cond=cond, n_ok=n_ok, nrep=nrep,
                  cor_tau=mean(cor_tau), cor_tau_mcse=sd(cor_tau)/sqrt(n_ok),
                  rmse_tau=mean(rmse_tau),
                  cor_loc=mean(cor_loc), cor_loc_mcse=sd(cor_loc)/sqrt(n_ok),
                  rmse_loc=mean(rmse_loc),
                  bias5=bias5, se_ratio=se_ratio, any_weak_rate=mean(any_weak)))
}

t0 <- Sys.time()
out <- list()
out$pcm_complete <- summarise_batch("PCM", "complete", nrep = 100, seed0 = 1000)
out$pcm_mcar     <- summarise_batch("PCM", "mcar",     nrep = 80,  seed0 = 2000)
out$pcm_booklet  <- summarise_batch("PCM", "booklet",  nrep = 60,  seed0 = 3000)
out$rsm_complete <- summarise_batch("RSM", "complete", nrep = 80,  seed0 = 4000)
out$rsm_mcar     <- summarise_batch("RSM", "mcar",     nrep = 60,  seed0 = 5000)
out$rsm_booklet  <- summarise_batch("RSM", "booklet",  nrep = 60,  seed0 = 6000)
t1 <- Sys.time()
cat("TOTAL TIME:", as.numeric(t1-t0, units="secs"), "s\n")

saveRDS(out, "recovery_results.rds")
