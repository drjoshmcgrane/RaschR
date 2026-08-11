suppressWarnings(pkgload::load_all(".", quiet = TRUE))
NR <- 150
out <- matrix(NA_real_, NR, 3, dimnames = list(NULL, c("la1","se1","p_alpha")))
for (r in seq_len(NR)) {
  d <- simulate_efrm(200, 8, n_sets = 2, n_groups = 2, set_unit_ratio = 1,
                     group_unit_ratio = 1, seed = 14.7e6 + r)
  fit <- tryCatch(rasch_efrm(d, item_sets = attr(d, "truth")$item_sets,
                             groups = "group", se_method = "bootstrap",
                             boot_reps = 150),
                  error = function(e) NULL)
  if (is.null(fit)) next
  uo <- fit$efrm_vs_rasch$unit_omnibus
  out[r, ] <- c(log(fit$alpha_table$alpha[1]), fit$alpha_table$se_log_alpha[1],
                uo$p[uo$term == "set units (alpha)"][1])
}
o <- out[complete.cases(out), , drop = FALSE]
r <- mean(o[,"p_alpha"] < .05)
cat(sprintf("EFRM FULL-BOOTSTRAP arm: n=%d SD(la1)=%.4f meanSE=%.4f ratio=%.3f type1=%.4f (mcse %.4f)\n",
    nrow(o), sd(o[,"la1"]), mean(o[,"se1"]), sd(o[,"la1"])/mean(o[,"se1"]),
    r, sqrt(r*(1-r)/nrow(o))))
