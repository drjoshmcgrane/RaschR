source("helpers.R")
d <- simulate_mfrm(50, 6, 6, n_categories=4, theta_sd=1.2, item_sd=1, rater_severity_sd=0.7, seed=1001)
set.seed(1001+5e5)
d$score[sample.int(nrow(d), round(0.25*nrow(d)))] <- NA
t0<-Sys.time()
mf <- rasch_mfrm(d, person="person", item="item", score="score", facets="rater")
cat("time N=50:", as.numeric(Sys.time()-t0,units="secs"), "\n")
d2 <- simulate_mfrm(50, 6, 6, n_categories=4, theta_sd=1.2, item_sd=1, rater_severity_sd=0.7, seed=1002)
set.seed(1002+5e5)
d2$score[sample.int(nrow(d2), round(0.25*nrow(d2)))] <- NA
t0<-Sys.time()
mf2 <- rasch_mfrm(d2, person="person", item="item", score="score", facets="rater")
cat("time N=50 (2):", as.numeric(Sys.time()-t0,units="secs"), "\n")
