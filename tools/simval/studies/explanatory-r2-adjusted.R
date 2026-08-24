# STUDY: explanatory-r2-adjusted
#
# The adjusted calibration coefficient explanatory_test() reports beside
# r_squared. The adjustment divides the unexplained proportion by the
# residual share of the degrees of freedom, with the residual dimension
# taken from the rank of the retained explanatory design; it is exact for
# independent homoskedastic least-squares estimates, which calibrations
# are not, so the documentation calls it a descriptive optimism
# adjustment. This study measures how that description behaves: where the
# raw coefficient centres under an uninformative design, where the
# adjusted one centres, and what either costs a true design.
#   Rscript tools/simval/studies/explanatory-r2-adjusted.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "explanatory-r2-adjusted"
rows <- list()
t0 <- Sys.time()

R <- 300L
run_cell <- function(label, K, useless) {
  raw <- adj <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    set.seed(920000 + K * 1000 + r + useless * 500000)
    P <- data.frame(item = sprintf("I%02d", seq_len(K)),
                    steps = rep(1:3, length.out = K),
                    abstract = rep(c(0, 1), length.out = K))
    delta <- if (useless) rnorm(K, 0, 1.0) else
      0.7 * P$steps - 0.5 * P$abstract + rnorm(K, 0, 0.15)
    delta <- delta - mean(delta)
    th <- rnorm(1000, 0, 1.2)
    X <- vapply(seq_len(K), function(i)
      rbinom(1000, 1, plogis(th - delta[i])), numeric(1000))
    colnames(X) <- P$item
    f <- tryCatch(rasch_explanatory(data.frame(id = 1:1000, X),
                                    predictors = P,
                                    formula = ~ steps + abstract, id = "id"),
                  error = function(e) NULL)
    if (is.null(f)) next
    z <- explanatory_test(f)
    raw[r] <- z$r_squared; adj[r] <- z$r_squared_adj
  }
  ok <- is.finite(raw) & is.finite(adj)
  rows[[length(rows) + 1L]] <<- sv_row(STUDY, label, "raw r_squared",
    sum(ok), bias = mean(raw[ok]), emp_sd = sd(raw[ok]),
    n_attempted = R, n_refused = R - sum(ok),
    notes = "bias column holds the mean of the coefficient itself")
  rows[[length(rows) + 1L]] <<- sv_row(STUDY, label, "adjusted r_squared",
    sum(ok), bias = mean(adj[ok]), emp_sd = sd(adj[ok]),
    n_attempted = R, n_refused = R - sum(ok),
    notes = "bias column holds the mean of the coefficient itself")
  cat(sprintf("[%s] %-38s raw %.3f  adjusted %+.3f  (n=%d)\n",
              format(Sys.time(), "%H:%M"), label, mean(raw[ok]),
              mean(adj[ok]), sum(ok)))
}
run_cell("uninformative design, 12 items", 12L, TRUE)
run_cell("uninformative design, 24 items", 24L, TRUE)
run_cell("true design, structural noise 0.15, 12 items", 12L, FALSE)
sv_write(do.call(rbind, rows), "explanatory-r2-adjusted")
cat(sprintf("TOTAL %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
