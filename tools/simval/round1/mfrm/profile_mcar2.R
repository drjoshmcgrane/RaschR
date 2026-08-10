source("helpers.R")
d <- simulate_mfrm(200, 6, 6, n_categories=4, theta_sd=1.2, item_sd=1, rater_severity_sd=0.7, seed=1001)
set.seed(1001+5e5)
d$score[sample.int(nrow(d), round(0.25*nrow(d)))] <- NA
Rprof("prof.out")
mf <- rasch_mfrm(d, person="person", item="item", score="score", facets="rater")
Rprof(NULL)
print(summaryRprof("prof.out")$by.self[1:15,])
