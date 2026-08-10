suppressWarnings(pkgload::load_all("."), quiet=TRUE))
options(width=140)

## ---------------------------------------------------------------
## (6) RSM vs PCM model comparison: cl_aic prefers the true model.
## ---------------------------------------------------------------
cat("=== compare_fits: cl_aic model selection ===\n")
pick_rsm_on_rsm <- c(); pick_pcm_on_pcm <- c()
nrep <- 15
for (r in seq_len(nrep)) {
  d_rsm <- simulate_rasch(n_persons = 500, n_items = 8, model = "RSM",
                           n_categories = 4, difficulty = c(-2, 2),
                           threshold_spread = 1.2, seed = 90000 + r)
  fp <- rasch(d_rsm, model = "PCM"); fr <- rasch(d_rsm, model = "RSM")
  cmp <- compare_fits(PCM = fp, RSM = fr)
  pick_rsm_on_rsm <- c(pick_rsm_on_rsm, cmp$cl_aic[cmp$label=="RSM"] < cmp$cl_aic[cmp$label=="PCM"])

  d_pcm <- simulate_rasch(n_persons = 500, n_items = 8, model = "PCM",
                           n_categories = 4, difficulty = c(-2, 2),
                           threshold_spread = 1.2, seed = 91000 + r)
  fp2 <- rasch(d_pcm, model = "PCM"); fr2 <- rasch(d_pcm, model = "RSM")
  cmp2 <- compare_fits(PCM = fp2, RSM = fr2)
  pick_pcm_on_pcm <- c(pick_pcm_on_pcm, cmp2$cl_aic[cmp2$label=="PCM"] < cmp2$cl_aic[cmp2$label=="RSM"])
}
cat(sprintf("RSM-true data: cl_aic picks RSM in %d/%d reps\n", sum(pick_rsm_on_rsm), nrep))
cat(sprintf("PCM-true data: cl_aic picks PCM in %d/%d reps\n\n", sum(pick_pcm_on_pcm), nrep))

## ---------------------------------------------------------------
## Identification guard: a DISCONNECTED booklet (no link items) must
## be refused, not silently fit.
## ---------------------------------------------------------------
cat("=== identification guard: disconnected booklet ===\n")
set.seed(99001)
d <- simulate_rasch(n_persons = 400, n_items = 10, model = "PCM",
                     n_categories = 4, difficulty = c(-2, 2),
                     threshold_spread = 1.3, seed = 99001)
N <- nrow(d)
grp <- sample(c("A", "B"), N, replace = TRUE)
d_disc <- d
# NO link items: A only ever sees I01-I05, B only ever sees I06-I10
d_disc[grp == "A", c("I06","I07","I08","I09","I10")] <- NA
d_disc[grp == "B", c("I01","I02","I03","I04","I05")] <- NA
res <- tryCatch(rasch(d_disc, model = "PCM"), error = function(e) e)
cat("class:", class(res)[1], "\n")
if (inherits(res, "error")) cat("refused as expected. message:\n", conditionMessage(res), "\n")
