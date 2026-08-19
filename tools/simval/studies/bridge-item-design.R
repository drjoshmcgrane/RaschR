# STUDY: bridge-item-design
#
# The common-item channel measured in common-item-channel.csv gave each
# frame its own persons. A bridge design for item-set units does not: the
# same person answers the bridge item in both set contexts, so the two
# calibrations share persons and their item-location errors correlate.
# This study measures what that does, and what re-administration risks.
#
# Cells:
#   lambda = 0     shared persons, responses conditionally independent
#                  given theta -- the ideal bridge
#   lambda > 0     shared persons with the second administration partly
#                  following the first (carry-over on the log-odds scale)
#                  -- the realistic bridge, since re-administering an item
#                  invites memory and consistency effects
#
# The estimator is the SD ratio of the two contexts' item calibrations,
# which common-item-channel.csv showed matches conditional ML exactly on
# independent persons. Planted unit ratio 1.30, dichotomous.
#
# What this study cannot address: whether re-administration is
# substantively acceptable in a given instrument. It quantifies the
# statistical cost of the dependence it plants, not the wisdom of the
# design.
# Serial. Rscript tools/simval/studies/bridge-item-design.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "bridge-item-design"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

R <- 40L
ratio <- 1.30
alpha <- c(ratio^-0.5, ratio^0.5)
lt <- log(ratio)

for (K in c(8L, 20L)) {
  delta <- seq(-2, 2, length.out = K)
  for (N in c(250L, 1000L)) {
    for (lambda in c(0, 0.5, 1.0)) {
      lr <- rep(NA_real_, R); n_ref <- 0L
      for (r in seq_len(R)) {
        set.seed(140000 + K * 20000 + N * 3 + round(lambda * 10) + r)
        th <- rnorm(N, 0, 1.3)
        XA <- sapply(delta, function(d)
          rbinom(N, 1, plogis(alpha[1] * (th - d))))
        # second administration of the SAME items to the SAME persons
        XB <- vapply(seq_along(delta), function(i) {
          lp <- alpha[2] * (th - delta[i]) + lambda * (2 * XA[, i] - 1)
          rbinom(N, 1, plogis(lp))
        }, numeric(N))
        colnames(XA) <- colnames(XB) <- sprintf("I%02d", seq_len(K))
        fa <- tryCatch(rasch(XA), error = function(e) NULL)
        fb <- tryCatch(rasch(XB), error = function(e) NULL)
        if (is.null(fa) || is.null(fb)) { n_ref <- n_ref + 1L; next }
        lr[r] <- log(sd(fb$items$location) / sd(fa$items$location))
      }
      ok <- is.finite(lr)
      add(sprintf("%d bridge items, N = %d shared persons, carry-over %.1f",
                  K, N, lambda),
          "log unit ratio, SD-ratio channel", sum(ok),
          bias = mean(lr[ok]) - lt, emp_sd = sd(lr[ok]),
          effect = lambda, n_attempted = R, n_refused = n_ref,
          notes = sprintf("items %d, N %d", K, N))
      cat(sprintf("[%s] K = %2d, N = %4d, carry-over %.1f: bias %+.4f (sd %.4f)\n",
                  format(Sys.time(), "%H:%M"), K, N, lambda,
                  mean(lr[ok]) - lt, sd(lr[ok])))
    }
  }
}

sv_write(do.call(rbind, rows), "bridge-item-design")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
