suppressWarnings(pkgload::load_all(".", quiet=TRUE))
pick_pcm <- function(n_persons, spread, nrep, seed0) {
  hits <- 0; ok <- 0
  for (r in seq_len(nrep)) {
    d <- simulate_rasch(n_persons = n_persons, n_items = 8, model = "PCM",
                         n_categories = 4, difficulty = c(-2, 2),
                         threshold_spread = spread, seed = seed0 + r)
    res <- tryCatch({
      fp <- rasch(d, model = "PCM"); fr <- rasch(d, model = "RSM")
      cmp <- compare_fits(PCM = fp, RSM = fr)
      cmp$cl_aic[cmp$label=="PCM"] < cmp$cl_aic[cmp$label=="RSM"]
    }, error = function(e) NA)
    if (!is.na(res)) { ok <- ok + 1; hits <- hits + res }
  }
  cat(sprintf("n=%d spread=%.1f: PCM picked %d/%d (ok=%d/%d)\n", n_persons, spread, hits, ok, ok, nrep))
}
set.seed(1)
pick_pcm(500, 1.2, 15, 92000)   # replicate the original setting
pick_pcm(1500, 1.2, 15, 93000)  # bigger N, same effect size
pick_pcm(500, 1.8, 15, 94000)   # bigger per-item pattern divergence
pick_pcm(1500, 1.8, 15, 95000)  # both bigger
