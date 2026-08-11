source("tools/simval/round1/mfrm/helpers.R")
d <- simulate_mfrm(150, 6, 8, n_categories=4, theta_sd=1.2, item_sd=1, rater_severity_sd=0.7, halo=0.25, seed=11)
truth <- attr(d,"truth")
cat("halo raters:", truth$halo, "\n")
t0<-Sys.time()
mf <- rasch_mfrm(d, person="person", item="item", score="score", facets="rater", interaction="rater")
cat("fit time:", as.numeric(Sys.time()-t0,units="secs"),"\n")
fe <- mf$facet_effects$rater
fe$halo <- fe$level %in% truth$halo
print(fe[,c("level","severity","fit_resid","fit_resid_pooled","infit_ms","outfit_ms")])
print(fe$halo)
ie <- mf$interaction_effects
ie$halo <- ie$level %in% truth$halo
cat("mean |gamma| halo vs clean:\n")
print(tapply(abs(ie$gamma), ie$halo, mean))
cat("significant interaction rate halo vs clean:\n")
print(tapply(ie$significant, ie$halo, mean))
print(mf$interaction_test)
