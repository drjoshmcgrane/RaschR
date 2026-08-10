suppressWarnings(pkgload::load_all("."), quiet=TRUE))
p_flags <- c(); p_adj_flags <- c(); p_anova_flags <- c()
for (r in seq_len(150)) {
  d <- simulate_rasch(n_persons = 500, n_items = 8, model = "PCM",
                       n_categories = 4, difficulty = c(-2, 2),
                       threshold_spread = 1.3, seed = 15000 + r)
  f <- tryCatch(rasch(d, model = "PCM"), error = function(e) e)
  if (inherits(f, "error")) next
  p_flags <- c(p_flags, f$items$p < 0.05)
  p_adj_flags <- c(p_adj_flags, f$items$p_adj < 0.05)
  p_anova_flags <- c(p_anova_flags, f$items$p_anova < 0.05)
}
cat(sprintf("raw item-trait p<.05 rate=%.4f (n=%d, MCse=%.4f)\n", mean(p_flags), length(p_flags), sqrt(.05*.95/length(p_flags))))
cat(sprintf("p_adj<.05 rate=%.4f\n", mean(p_adj_flags)))
cat(sprintf("class-interval ANOVA p_anova<.05 rate=%.4f\n", mean(p_anova_flags)))
