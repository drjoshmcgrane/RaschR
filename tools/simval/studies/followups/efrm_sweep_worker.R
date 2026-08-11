a <- commandArgs(TRUE)   # name npg K S theta_sd imbalance nreps seed0 outfile
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
name <- a[1]; npg <- as.integer(a[2]); K <- as.integer(a[3]); S <- as.integer(a[4])
tsd <- as.numeric(a[5]); imb <- as.logical(a[6]); nr <- as.integer(a[7])
seed0 <- as.numeric(a[8]); outf <- a[9]
out <- matrix(NA_real_, nr, 4, dimnames = list(NULL, c("la1","se_tab1","p_alpha","p_phi")))
nfail <- 0L
for (r in seq_len(nr)) {
  d <- simulate_efrm(npg, K, n_sets = S, n_groups = 2, set_unit_ratio = 1,
                     group_unit_ratio = 1, theta_sd = tsd, seed = seed0 + r)
  if (imb) {   # keep only a third of group B's persons: 300 vs 100
    gb <- unique(d$group)[2]
    keep_ids <- unique(d$id[d$group == gb])
    drop_ids <- keep_ids[seq_len(floor(2 * length(keep_ids) / 3))]
    d <- d[!(d$id %in% drop_ids), , drop = FALSE]
  }
  fit <- tryCatch(rasch_efrm(d, item_sets = attr(d, "truth")$item_sets,
                             groups = "group", se_method = "hybrid"),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$est$converged)) { nfail <- nfail + 1L; next }
  uo <- fit$efrm_vs_rasch$unit_omnibus
  out[r, ] <- c(log(fit$alpha_table$alpha[1]), fit$alpha_table$se_log_alpha[1],
                uo$p[uo$term == "set units (alpha)"][1],
                uo$p[uo$term == "group units (phi)"][1])
}
saveRDS(list(cell = name, out = out, nfail = nfail), outf)
cat("done", name, "fails", nfail, "\n")
