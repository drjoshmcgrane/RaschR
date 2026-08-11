suppressWarnings(pkgload::load_all(".", quiet = TRUE))
NR <- 400
pk <- rep(NA_real_, NR)
for (r in seq_len(NR)) {
  d <- simulate_btl_efrm(6, 2, set_units = c(1, 1), set_origins = c(0, 0),
                         seed = 13.6e6 + r)
  fit <- tryCatch(btl_efrm(d, "object_a", "object_b", winner = "winner",
                           judge = "judge", panels = "panel",
                           object_sets = attr(d, "truth")$object_sets,
                           se_method = "judge_bootstrap", boot_reps = 150),
                  error = function(e) NULL)
  if (is.null(fit)) next
  uo <- fit$unit_omnibus
  if (!is.null(uo)) {
    pr <- uo$p[grepl("origin|kappa", uo$term, ignore.case = TRUE)]
    if (length(pr)) pk[r] <- pr[1]
  }
}
ok <- is.finite(pk); r <- mean(pk[ok] < .05)
cat(sprintf("BTL-EFRM kappa recheck: type1=%.4f (mcse %.4f, n=%d)\n",
    r, sqrt(r*(1-r)/sum(ok)), sum(ok)))
