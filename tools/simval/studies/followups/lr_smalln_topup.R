# Top-up of the two elevated lr-test secondary cells (N=300, 12 items,
# 3 and 4 categories: 8.0% and 7.55% at 400 replicates) to >= 2,000
# replicates each, per review: a coherent small-N pattern needs the
# replicates to call it. Reuses the study's own generator via source of
# its function block.
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
options(simval.script = "tools/simval/studies/followups/lr_smalln_topup.R")
lines <- readLines("tools/simval/studies/lr-test.R")
a <- grep("^size_once <- function", lines)
b <- grep("^run_size_cell <- function", lines)
eval(parse(text = lines[a:(b - 1)]))          # size_once only
rows <- list()
for (ncat in c(3L, 4L)) {
  NREP <- 2000L
  p_adj <- rep(NA_real_, NREP); n_ref <- n_nc <- n_noadj <- 0L
  for (r in seq_len(NREP)) {
    out <- size_once(300L, 12L, ncat, seed = 41e6 + ncat * 1e5 + r)
    if (identical(out$status, "refusal")) { n_ref <- n_ref + 1L; next }
    if (identical(out$status, "nonconv")) { n_nc <- n_nc + 1L; next }
    if (identical(out$status, "no_adjustment")) { n_noadj <- n_noadj + 1L; next }
    p_adj[r] <- out$p_adj
  }
  ok <- is.finite(p_adj); t1 <- mean(p_adj[ok] < 0.05)
  cat(sprintf("N=300 ni=12 ncat=%d: type1=%.4f (mc %.4f) n=%d ref=%d nc=%d noadj=%d\n",
      ncat, t1, sqrt(t1 * (1 - t1) / sum(ok)), sum(ok), n_ref, n_nc, n_noadj))
  rows[[length(rows) + 1]] <- sv_row("lr-smalln-topup",
    sprintf("size N=300 ni=12 ncat=%d (top-up)", ncat), "type1_p_adj",
    sum(ok), type1 = t1,
    n_attempted = NREP, n_refused = n_ref + n_noadj, n_nonconv = n_nc,
    notes = sprintf("2,000-replicate top-up of the 400-rep secondary cell (was %.1f%%); no_adjustment=%d folded into n_refused",
                    if (ncat == 3) 8.0 else 7.55, n_noadj))
}
sv_write(do.call(rbind, rows), "lr-smalln-topup")
