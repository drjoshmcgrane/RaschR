suppressWarnings(pkgload::load_all("."), quiet=TRUE))
N <- 800; I <- 10; NREPS <- 300
out <- matrix(NA_real_, NREPS, 6,
  dimnames=list(NULL, c("d","se_naive","se_corr","tau_lo","tau_hi","cv_lohi")))
for (r in seq_len(NREPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, seed = 555000 + r)
  f <- rasch(s)
  dm <- dependence_magnitude(f, dependent = "I06", independent = "I05")
  rf <- dm$refit
  thr <- rf$thresholds; item_of <- rf$items$item[thr$item]
  lo_id <- thr$id[item_of == "I06|I05=0" & thr$k == 1]
  hi_id <- thr$id[item_of == "I06|I05=1" & thr$k == 1]
  cv <- rf$est$cov_tau
  se_corr <- sqrt((cv[lo_id,lo_id] + cv[hi_id,hi_id] - 2*cv[lo_id,hi_id]) / 4)
  out[r,] <- c(dm$d, dm$se, se_corr,
               thr$tau[item_of == "I06|I05=0" & thr$k == 1],
               thr$tau[item_of == "I06|I05=1" & thr$k == 1],
               cv[lo_id, hi_id])
}
emp <- sd(out[,"d"])
z_n <- out[,"d"] / out[,"se_naive"]; z_c <- out[,"d"] / out[,"se_corr"]
cat(sprintf("empirical SD(d)        = %.4f\n", emp))
cat(sprintf("mean naive SE          = %.4f  (ratio %.3f)  reject@.05 = %.3f\n",
            mean(out[,"se_naive"]), emp/mean(out[,"se_naive"]),
            mean(2*pnorm(-abs(z_n)) < 0.05)))
cat(sprintf("mean corrected SE      = %.4f  (ratio %.3f)  reject@.05 = %.3f\n",
            mean(out[,"se_corr"]), emp/mean(out[,"se_corr"]),
            mean(2*pnorm(-abs(z_c)) < 0.05)))
cat(sprintf("empirical cor(lo,hi)   = %.3f;  mean sandwich cov(lo,hi) = %.5f (empirical %.5f)\n",
            cor(out[,"tau_lo"], out[,"tau_hi"]), mean(out[,"cv_lohi"]),
            cov(out[,"tau_lo"], out[,"tau_hi"])))
cat(sprintf("MC err on reject rates ~ %.3f\n", sqrt(0.05*0.95/NREPS)))
