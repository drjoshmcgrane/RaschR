suppressWarnings(pkgload::load_all(".", quiet=TRUE))

N <- 400; I <- 12
OUTER <- 20; INNER_REPS <- 20

run_cond <- function(gen_fun, label) {
  ratios <- numeric(OUTER); obs1 <- numeric(OUTER); ref1 <- numeric(OUTER)
  above <- logical(OUTER)
  for (r in seq_len(OUTER)) {
    s <- gen_fun(r)
    f <- rasch(s)
    pc <- residual_pca(f)
    # Round 1 assessed centring against the simulated mean. The current plot
    # uses the finite-simulation 5% upper critical value for inference.
    ref <- attr(.scree_reference(f, 5, INNER_REPS), "mean")
    obs1[r] <- pc$first_eigen; ref1[r] <- ref[1]
    ratios[r] <- pc$first_eigen / ref[1]
    above[r] <- pc$first_eigen > ref[1]
  }
  cat(sprintf("[%s] mean obs1=%.3f, mean ref1=%.3f, mean ratio=%.3f, sd ratio=%.3f, P(obs>ref)=%.3f (n=%d, MCerr=%.3f)\n",
              label, mean(obs1), mean(ref1), mean(ratios), sd(ratios), mean(above), OUTER,
              sqrt(mean(above)*(1-mean(above))/OUTER)))
  invisible(list(obs1 = obs1, ref1 = ref1, ratios = ratios, above = above))
}

cat("=== NULL: unidimensional data ===\n")
null_out <- run_cond(function(r) simulate_rasch(n_persons = N, n_items = I, seed = 6000 + r), "null")

cat("\n=== PLANTED: second dimension (rho=0.3) on items 7-12 ===\n")
dim_out <- run_cond(function(r) simulate_rasch(n_persons = N, n_items = I,
                     second_dim = list(items = 7:12, rho = 0.3), seed = 6500 + r), "second_dim rho=0.3")

cat("\n=== PLANTED: second dimension (rho=0.0, orthogonal) on items 7-12 ===\n")
dim_out0 <- run_cond(function(r) simulate_rasch(n_persons = N, n_items = I,
                     second_dim = list(items = 7:12, rho = 0), seed = 7000 + r), "second_dim rho=0.0")
