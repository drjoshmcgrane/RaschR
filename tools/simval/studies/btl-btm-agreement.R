# STUDY: btl-btm-agreement
#
# btl() against sirt::btm on shared dichotomous comparison data. With
# eps = 0 and a fixed home advantage the two maximise the same
# Bradley-Terry likelihood, so agreement should reach solver precision.
# Replicates containing a fully extreme object are recorded separately:
# btl() removes the unidentified object by design, whereas btm keeps it
# and diverges toward the boundary at eps = 0 (its default eps = 0.3
# shrinks it instead), so their solutions differ by policy, not
# estimation.
#   Rscript tools/simval/studies/btl-btm-agreement.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "btl-btm-agreement"
rows <- list()

R <- 25L
md <- rep(NA_real_, R); n_extreme <- 0L
for (r in seq_len(R)) {
  set.seed(710000 + r)
  d <- simulate_btl(n_objects = 8, n_judges = 30, reps_per_pair = 2)
  wins <- tapply(c(d$winner == d$object_a, d$winner == d$object_b),
                 c(d$object_a, d$object_b), mean)
  if (any(wins %in% c(0, 1))) { n_extreme <- n_extreme + 1L; next }
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  b_r <- setNames(f$objects$location, f$objects$object)
  dat <- data.frame(id1 = d$object_a, id2 = d$object_b,
                    result = ifelse(d$winner == d$object_a, 1, 0))
  q <- capture.output(m <- sirt::btm(dat, eps = 0, fix.eta = 0,
                                     maxiter = 5000, conv = 1e-10))
  g <- setNames(m$effects$theta, m$effects$individual)
  common <- intersect(names(b_r), names(g))
  md[r] <- max(abs((b_r[common] - mean(b_r[common])) -
                   (g[common] - mean(g[common]))))
}
ok <- !is.na(md)
rows[[1]] <- sv_row(STUDY, "dichotomous, no extreme objects",
  "max |location difference|", sum(ok),
  bias = mean(md[ok]), emp_sd = max(md[ok]),
  n_attempted = R, n_refused = n_extreme,
  notes = paste("bias holds the mean and emp_sd the worst max-difference;",
                "refusals are replicates with a fully extreme object, where",
                "btl removes the object and btm diverges at eps = 0"))
cat(sprintf("agreement over %d clean replicates: mean %.2e, worst %.2e; %d extreme-object replicates set aside\n",
            sum(ok), mean(md[ok]), max(md[ok]), n_extreme))
sv_write(do.call(rbind, rows), "btl-btm-agreement")
