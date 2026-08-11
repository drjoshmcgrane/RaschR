a <- commandArgs(TRUE)
offset <- as.integer(a[1]); nreps <- as.integer(a[2]); outfile <- a[3]
suppressWarnings(pkgload::load_all(".", quiet=TRUE))
N <- 800; I <- 10
out <- matrix(NA_real_, nreps, 9, dimnames=list(NULL,
  c("seed","d","se_naive","se_corr","tau_lo","tau_hi","cv_lohi","n_lo1","spread")))
for (r in seq_len(nreps)) {
  sd_ <- offset + r
  s <- simulate_rasch(n_persons = N, n_items = I, seed = sd_)
  f <- rasch(s)
  dm <- tryCatch(dependence_magnitude(f, dependent = "I06", independent = "I05"),
                 error = function(e) NULL)
  if (is.null(dm)) next
  rf <- dm$refit
  thr <- rf$thresholds; item_of <- rf$items$item[thr$item]
  lo_id <- thr$id[item_of == "I06|I05=0" & thr$k == 1]
  hi_id <- thr$id[item_of == "I06|I05=1" & thr$k == 1]
  cv <- rf$est$cov_tau
  se_corr <- sqrt((cv[lo_id,lo_id] + cv[hi_id,hi_id] - 2*cv[lo_id,hi_id]) / 4)
  X <- f$X
  out[r,] <- c(sd_, dm$d, dm$se, se_corr,
               thr$tau[item_of == "I06|I05=0" & thr$k == 1],
               thr$tau[item_of == "I06|I05=1" & thr$k == 1],
               cv[lo_id, hi_id],
               sum(X[,"I05"] == 0 & !is.na(X[,"I06"]), na.rm=TRUE),
               diff(range(f$items$location)))
}
saveRDS(out, outfile)
cat("done", outfile, "\n")
