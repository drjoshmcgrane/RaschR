suppressWarnings(pkgload::load_all(".", quiet = TRUE))
NREPS <- 400
out <- matrix(NA_real_, NREPS, 5,
  dimnames = list(NULL, c("la1", "la2", "se_tab1", "W", "p")))
for (r in seq_len(NREPS)) {
  d <- simulate_efrm(200, 8, n_sets = 2, n_groups = 2, set_unit_ratio = 1,
                     group_unit_ratio = 1, seed = 4.4e6 + r)
  fit <- tryCatch(rasch_efrm(d, item_sets = attr(d, "truth")$item_sets,
                             groups = "group", se_method = "hybrid"),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$est$converged)) next
  at <- fit$alpha_table
  uo <- fit$efrm_vs_rasch$unit_omnibus
  wa <- uo[uo$term == "set units (alpha)", ]
  out[r, ] <- c(log(at$alpha[1]), log(at$alpha[2]), at$se_log_alpha[1],
                wa$wald, wa$p)
}
out <- out[complete.cases(out), ]
n <- nrow(out)
cat(sprintf("n=%d analysed\n", n))
cat(sprintf("identification: cor(la1, la2) = %.3f, mean(la1+la2) = %.4f\n",
    cor(out[,"la1"], out[,"la2"]), mean(out[,"la1"] + out[,"la2"])))
cat(sprintf("TABLE:   empirical SD(la1) = %.4f vs mean se_log_alpha = %.4f  (ratio %.3f)\n",
    sd(out[,"la1"]), mean(out[,"se_tab1"]), sd(out[,"la1"]) / mean(out[,"se_tab1"])))
d12 <- out[,"la1"] - out[,"la2"]
se_impl <- abs(d12) / sqrt(out[,"W"])          # omnibus-implied SE of the contrast
cat(sprintf("OMNIBUS: empirical SD(la1-la2) = %.4f vs mean implied SE = %.4f  (ratio %.3f)\n",
    sd(d12), mean(se_impl), sd(d12) / mean(se_impl)))
cat(sprintf("omnibus rejection = %.4f; table-z rejection (la1/se_tab1) = %.4f\n",
    mean(out[,"p"] < .05), mean(abs(out[,"la1"] / out[,"se_tab1"]) > 1.96)))
