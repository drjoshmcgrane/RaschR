suppressWarnings(pkgload::load_all("."), quiet=TRUE))
N <- 800; I <- 10; NREPS <- 400
res <- matrix(NA_real_, NREPS, 3, dimnames=list(NULL, c("d","se","p")))
for (r in seq_len(NREPS)) {
  s <- simulate_rasch(n_persons = N, n_items = I, seed = 77000 + r)
  f <- rasch(s)
  dm <- dependence_magnitude(f, dependent = "I06", independent = "I05")
  res[r,] <- c(dm$d, dm$se, dm$p)
}
cat(sprintf("VERIFY NULL: reject@.05 = %.4f (n=%d, MC err %.4f)\n",
            mean(res[,"p"] < 0.05), NREPS, sqrt(0.05*0.95/NREPS)))
cat(sprintf("mean d = %.4f, empirical SD(d) = %.4f, mean reported SE = %.4f, ratio = %.3f\n",
            mean(res[,"d"]), sd(res[,"d"]), mean(res[,"se"]), sd(res[,"d"])/mean(res[,"se"])))
saveRDS(res, file.path(dirname("tools/simval/round1/verify"), "dm_null_verify.rds"))
