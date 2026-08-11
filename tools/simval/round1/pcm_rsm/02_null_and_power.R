suppressWarnings(pkgload::load_all(".", quiet=TRUE))
options(width=140)

make_booklet <- function(d) {
  N <- nrow(d)
  grp <- sample(c("A","B"), N, replace = TRUE)
  d2 <- d
  d2[grp == "A", c("I07","I08","I09","I10")] <- NA
  d2[grp == "B", c("I01","I02","I03","I04")] <- NA
  d2
}

## ---------------------------------------------------------------
## (3) Null calibration: item chi-square (p_adj) and |infit_z|,
## |outfit_z| flag rates under model-true PCM data.
## ---------------------------------------------------------------
null_rates <- function(cond, nrep, seed0, n_persons = 500, n_items = 8) {
  chisq_flags <- c(); infit_flags <- c(); outfit_flags <- c(); errs <- 0
  for (r in seq_len(nrep)) {
    d <- simulate_rasch(n_persons = n_persons, n_items = n_items, model = "PCM",
                         n_categories = 4, difficulty = c(-2, 2),
                         threshold_spread = 1.3, seed = seed0 + r)
    if (cond == "mcar") {
      for (ic in grep("^I", names(d))) {
        rows <- sample(seq_len(nrow(d)), size = round(0.25 * nrow(d)))
        d[rows, ic] <- NA
      }
    } else if (cond == "booklet") d <- make_booklet(d)
    fit <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
    if (inherits(fit, "error")) { errs <- errs + 1; next }
    chisq_flags <- c(chisq_flags, fit$items$p_adj < 0.05)
    infit_flags <- c(infit_flags, abs(fit$items$infit_z) > 1.96)
    outfit_flags <- c(outfit_flags, abs(fit$items$outfit_z) > 1.96)
  }
  cat(sprintf("[null | PCM | %s] reps_ok=%d/%d errors=%d\n  chisq p_adj<.05 rate=%.4f (n=%d, MCse=%.4f)\n  |infit_z|>1.96 rate=%.4f  |outfit_z|>1.96 rate=%.4f\n\n",
              cond, nrep - errs, nrep, errs,
              mean(chisq_flags), length(chisq_flags), sqrt(0.05*0.95/length(chisq_flags)),
              mean(infit_flags), mean(outfit_flags)))
}

## ---------------------------------------------------------------
## (4) Threshold-order diagnostics: null (well-ordered) vs power
## (planted disordered, sparse-middle-category mechanism)
## ---------------------------------------------------------------
threshold_order_rates <- function(planted, nrep, seed0, n_persons = 500, n_items = 8) {
  flag_target <- c(); flag_other <- c(); errs <- 0
  for (r in seq_len(nrep)) {
    dis_arg <- if (planted) "I04" else NULL
    d <- simulate_rasch(n_persons = n_persons, n_items = n_items, model = "PCM",
                         n_categories = 4, difficulty = c(-2, 2),
                         threshold_spread = 1.3, disordered = dis_arg,
                         seed = seed0 + r)
    fit <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
    if (inherits(fit, "error")) { errs <- errs + 1; next }
    ordered_target <- fit$thresholds_diag[["I04"]]$ordered
    flag_target <- c(flag_target, !ordered_target)
    ordered_other <- vapply(fit$thresholds_diag[setdiff(names(fit$thresholds_diag),"I04")],
                             function(x) x$ordered, TRUE)
    flag_other <- c(flag_other, mean(!ordered_other))
  }
  cat(sprintf("[threshold-order | planted=%s] reps_ok=%d/%d errors=%d\n  I04 flagged-disordered rate=%.4f (MCse=%.4f)\n  other items flagged-disordered rate=%.4f\n\n",
              planted, nrep - errs, nrep, errs,
              mean(flag_target), sqrt(mean(flag_target)*(1-mean(flag_target))/length(flag_target)),
              mean(flag_other)))
}

t0 <- Sys.time()
for (cond in c("complete", "mcar", "booklet"))
  null_rates(cond, nrep = 250, seed0 = switch(cond, complete=10000, mcar=20000, booklet=30000))

threshold_order_rates(planted = FALSE, nrep = 150, seed0 = 40000)
threshold_order_rates(planted = TRUE,  nrep = 150, seed0 = 50000)
t1 <- Sys.time()
cat("TOTAL TIME:", as.numeric(t1-t0, units="secs"), "s\n")
