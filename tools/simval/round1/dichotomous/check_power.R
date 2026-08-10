suppressWarnings(pkgload::load_all("."), quiet=TRUE))
set.seed(1)
fr <- replicate(40, {
  d <- simulate_rasch(500, 15, difficulty=c(-2,2), discrimination=c(2.8, rep(1,14)))
  it <- rasch(d, id="id")$items
  c(fit_resid=it$fit_resid[1], outfit=it$outfit_ms[1], infit=it$infit_ms[1])
})
print(t(fr))
cat("mean fit_resid:", mean(fr["fit_resid",]), " sd:", sd(fr["fit_resid",]), "\n")
cat("mean outfit:", mean(fr["outfit",]), "\n")

# try a stronger departure and bigger N
fr2 <- replicate(20, {
  d <- simulate_rasch(1000, 15, difficulty=c(-2,2), discrimination=c(4, rep(1,14)))
  it <- rasch(d, id="id")$items
  c(fit_resid=it$fit_resid[1], outfit=it$outfit_ms[1])
})
cat("stronger disc=4, N=1000: mean fit_resid:", mean(fr2["fit_resid",]), "sd:", sd(fr2["fit_resid",]),
    " prop >2.5:", mean(abs(fr2["fit_resid",])>2.5), "\n")
