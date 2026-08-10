suppressWarnings(pkgload::load_all("."), quiet = TRUE))
I <- 15; N <- 500
delta0 <- seq(-2.5, 2.5, length.out = I)
set.seed(51)
n <- 80
raw_flags <- list(); adj_flags <- list()
for (r in seq_len(n)) {
  dA <- simulate_rasch(N, I, difficulty = delta0, seed = 500000 + r)
  dB <- simulate_rasch(N, I, difficulty = delta0, seed = 600000 + r)
  fitA <- rasch(dA, id = "id"); fitB <- rasch(dB, id = "id")
  eq <- equate_tests(fitA, fitB, independent = TRUE)
  raw_flags[[r]] <- eq$table$p < 0.05
  adj_flags[[r]] <- eq$table$p_adj < 0.05
}
rf <- unlist(raw_flags); rf <- rf[!is.na(rf)]
af <- unlist(adj_flags); af <- af[!is.na(af)]
cat(sprintf("raw p<0.05 rate: %.4f (n=%d)\n", mean(rf), length(rf)))
cat(sprintf("BH-adjusted p_adj<0.05 rate: %.4f (n=%d)\n", mean(af), length(af)))
