suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/studies/followups/mfrm_fns.R")
for (N in c(50, 200)) {
  b <- mfrm_batch(N, 6, NA, 600, seed0 = 12.5e6 + N)
  r <- mean(b$p < 0.05)
  cat(sprintf("MFRM q=25 recheck N=%d R=6: type1=%.4f (mcse %.4f, n=%d) refused=%d nonconv=%d\n",
      N, r, sqrt(r*(1-r)/length(b$p)), length(b$p), b$n_refused, b$n_nonconv))
}
