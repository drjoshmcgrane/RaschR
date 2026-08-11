suppressWarnings(pkgload::load_all(".", quiet=TRUE))
set.seed(2)
for (disc in c(2.8, 4, 6, 8, 12)) {
  fr <- replicate(25, {
    d <- simulate_rasch(500, 15, difficulty=c(-2,2), discrimination=c(disc, rep(1,14)))
    it <- rasch(d, id="id")$items
    c(fit_resid=it$fit_resid[1], p_anova=it$p_anova[1])
  })
  cat(sprintf("disc=%-4s mean fit_resid=%.3f  prop|>2.5|=%.2f  prop p_anova<.05=%.2f\n",
              disc, mean(fr["fit_resid",]), mean(abs(fr["fit_resid",])>2.5),
              mean(fr["p_anova",]<0.05)))
}
