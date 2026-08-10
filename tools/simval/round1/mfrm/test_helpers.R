source("helpers.R")

# DIF generator sanity
d <- sim_mfrm_dif(150, 6, 6, dif_item = "I3", dif_shift = 1.0, seed = 1)
mf <- rasch_mfrm(d, person="person", item="item", score="score", facets="rater", factors="group")
da <- dif_anova(mf)
print(da$summary[da$summary$item=="I3", ])
ds <- dif_size(mf, "I3", by="group")
print(ds$pairs)

# connected incomplete
d2 <- simulate_mfrm(150, 6, 6, rater_severity_sd=0.7, seed=2)
dci <- make_connected_incomplete(d2, k=3, seed=1)
cat("rows kept:", nrow(dci), "of", nrow(d2), "\n")
t0<-Sys.time()
mf2 <- rasch_mfrm(dci, person="person", item="item", score="score", facets="rater")
cat("fit time:", as.numeric(Sys.time()-t0,units="secs"),"\n")
print(cor(mf2$facet_effects$rater$severity, attr(d2,"truth")$severity[mf2$facet_effects$rater$level]))

# disconnected
dd <- make_disconnected(d2)
res <- tryCatch(rasch_mfrm(dd, person="person", item="item", score="score", facets="rater"),
                error=function(e) e)
print(class(res))
if(inherits(res,"error")) cat("refused as expected:", conditionMessage(res), "\n")

# true thresholds check
tb <- true_base_tau(4)
print(tb)
truth <- attr(d2,"truth")
true_tau <- outer(tb, truth$difficulty, "+")  # m x I matrix, rows k cols item
print(dim(true_tau))
