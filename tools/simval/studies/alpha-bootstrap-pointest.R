# STUDY: alpha-bootstrap-pointest
#
# Asked and answered: can a person bootstrap improve the set-unit POINT
# estimate (not its standard error)? Result at 300 persons/group and 8
# items/set, 20 datasets x 30 resamples, against a planted log ratio of
# 0.2624 -- plain (ships): bias -0.0235, sd 0.0886, RMSE 0.0895; bagged:
# -0.0183, 0.0904, 0.0900; bias-corrected: -0.0286, 0.0901, 0.0924.
# Neither variant helps. Bagging trades a little bias for a little
# variance at equal RMSE; bias correction is worse on both counts,
# because it adds resampling noise while removing the small-sample
# Jensen term that partly offsets the fixed-design floor. A bootstrap
# redistributes information already in the sample, and the binding
# constraint here is the thinness of the linking channel itself, not
# statistical inefficiency: longer or polytomous sets, and enough
# linking persons, are the levers (see alpha-n-sweep.csv and the floor
# curve in alpha-correction-limits.csv).
#
# Can a person bootstrap improve the set-unit POINT estimate?
#   plain   theta_hat            (what ships)
#   bagged  mean(theta*)         (bootstrap aggregation)
#   BC      2 theta_hat - mean(theta*)   (bootstrap bias correction)
# Judged on bias and RMSE against the planted log ratio. Design chosen
# where variance is largest and any improvement would matter most.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
D <- 20L; B <- 30L
lt <- log(1.30)
res <- matrix(NA_real_, D, 3, dimnames = list(NULL, c("plain","bagged","BC")))
t0 <- Sys.time()
for (i in seq_len(D)) {
  d <- simulate_efrm(n_per_group = 300, items_per_set = 8,
                     set_unit_ratio = 1.30, group_unit_ratio = 1.10,
                     seed = 94000 + i)
  tr <- attr(d, "truth")
  fit_one <- function(dat) {
    f <- tryCatch(rasch_efrm(dat, item_sets = tr$item_sets, groups = "group",
                             id = "id", boot_reps = 0), error = function(e) NULL)
    if (is.null(f)) return(NA_real_)
    a <- f$alpha_table$alpha
    log(max(a) / min(a))
  }
  th <- fit_one(d)
  if (!is.finite(th)) next
  bs <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    set.seed(500000 + 1000 * i + b)
    idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
    db <- d[idx, ]; db$id <- sprintf("P%05d", seq_len(nrow(db)))
    bs[b] <- fit_one(db)
  }
  mb <- mean(bs, na.rm = TRUE)
  res[i, ] <- c(th, mb, 2 * th - mb)
  cat(sprintf("[%s] dataset %2d: plain %+.4f bagged %+.4f BC %+.4f\n",
      format(Sys.time(), "%H:%M"), i, th - lt, mb - lt, 2 * th - mb - lt))
}
ok <- stats::complete.cases(res)
cat(sprintf("\n%d datasets, B = %d, truth log ratio %.4f\n", sum(ok), B, lt))
for (m in colnames(res))
  cat(sprintf("  %-7s bias %+.4f  sd %.4f  RMSE %.4f\n", m,
      mean(res[ok, m]) - lt, sd(res[ok, m]),
      sqrt(mean((res[ok, m] - lt)^2))))
cat(sprintf("elapsed %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
