a <- commandArgs(TRUE); off <- as.integer(a[1]); nr <- as.integer(a[2]); outf <- a[3]
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
out <- matrix(NA_real_, nr, 4, dimnames = list(NULL, c("la1","se_tab1","W","p")))
for (r in seq_len(nr)) {
  d <- simulate_efrm(200, 8, n_sets = 2, n_groups = 2, set_unit_ratio = 1,
                     group_unit_ratio = 1, seed = off + r)
  fit <- tryCatch(rasch_efrm(d, item_sets = attr(d, "truth")$item_sets,
                             groups = "group", se_method = "hybrid"),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$est$converged)) next
  uo <- fit$efrm_vs_rasch$unit_omnibus
  wa <- uo[uo$term == "set units (alpha)", ]
  out[r, ] <- c(log(fit$alpha_table$alpha[1]), fit$alpha_table$se_log_alpha[1],
                wa$wald, wa$p)
}
saveRDS(out, outf)
