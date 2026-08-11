source("tools/simval/round1/mfrm/helpers.R")
for (shift in c(0.4, 0.6, 0.8)) {
  d <- sim_mfrm_dif(200, 6, 6, dif_item = "I3", dif_shift = shift, seed = 1)
  mf <- rasch_mfrm(d, person="person", item="item", score="score", facets="rater", factors="group")
  ds <- dif_size(mf, "I3", by="group")
  cat("shift=", shift, " recovered diff=", ds$pairs$difference, " se=", ds$pairs$se, " weak:", any(ds$levels$weak), "\n")
}
