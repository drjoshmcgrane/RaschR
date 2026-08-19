# STUDY: alpha-set-misfit
#
# What actually threatens an item-set unit.
#
# Item sets partition the items, so an item appears in exactly one frame
# and cannot "behave differently across frames" -- the departure that
# distorts person-group units has no analogue here. The threat to alpha is
# misfit CONCENTRATED IN ONE SET: it distorts that set's person estimates,
# hence that set's true-score variance, and the other set carries nothing
# to cancel it.
#
# Cells: misfit balanced across the two sets (should largely cancel, as in
# misfit-both-channels.csv) against the same amount concentrated in set 1
# (should not). Over- and under-discrimination both, since they distort
# the person estimates in opposite directions.
#
# 8 items per set, N = 500, planted alpha ratio 1.40, estimated by the
# package's corrected linking. Serial.
#   Rscript tools/simval/studies/alpha-set-misfit.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "alpha-set-misfit"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

ratio <- 1.40
alpha <- c(ratio^-0.5, ratio^0.5)
lt <- log(ratio)
K <- 8L
N <- 500L
delta <- seq(-1.5, 1.5, length.out = K)
R <- 100L

run_cell <- function(label, disc1, disc2) {
  lr <- rep(NA_real_, R); n_ref <- 0L
  items <- c(sprintf("S1I%02d", seq_len(K)), sprintf("S2I%02d", seq_len(K)))
  isets <- list(set1 = items[seq_len(K)], set2 = items[K + seq_len(K)])
  for (r in seq_len(R)) {
    set.seed(210000 + nchar(label) * 613 + r)
    th <- rnorm(N, 0, 1.3)
    X <- cbind(
      vapply(seq_len(K), function(i)
        rbinom(N, 1, plogis(alpha[1] * disc1[i] * (th - delta[i]))), numeric(N)),
      vapply(seq_len(K), function(i)
        rbinom(N, 1, plogis(alpha[2] * disc2[i] * (th - delta[i]))), numeric(N)))
    colnames(X) <- items
    d <- data.frame(id = sprintf("P%05d", seq_len(N)), X, group = "g1",
                    check.names = FALSE)
    f <- tryCatch(rasch_efrm(d, item_sets = isets, groups = "group",
                             id = "id", boot_reps = 0), error = function(e) NULL)
    if (is.null(f)) { n_ref <- n_ref + 1L; next }
    a <- f$alpha_table$alpha[match(c("set1", "set2"), f$alpha_table$set)]
    lr[r] <- log(a[2] / a[1])
  }
  ok <- is.finite(lr)
  add(label, "log alpha ratio bias", sum(ok),
      bias = mean(lr[ok]) - lt, emp_sd = sd(lr[ok]),
      n_attempted = R, n_refused = n_ref,
      notes = sprintf("recovered %.3f vs planted %.3f", exp(mean(lr[ok])), ratio))
  cat(sprintf("[%s] %-42s recovered %.3f (bias %+.4f, sd %.4f)\n",
              format(Sys.time(), "%H:%M"), label,
              exp(mean(lr[ok])), mean(lr[ok]) - lt, sd(lr[ok])))
}

flat <- rep(1, K)
two_over  <- replace(rep(1, K), c(3, 6), 2.0)
two_under <- replace(rep(1, K), c(3, 6), 0.5)
four_over <- replace(rep(1, K), c(2, 4, 6, 8), 2.0)

run_cell("clean", flat, flat)
run_cell("balanced: 2 over-discriminating in EACH set", two_over, two_over)
run_cell("concentrated: 2 over-discriminating in set 1", two_over, flat)
run_cell("concentrated: 2 over-discriminating in set 2", flat, two_over)
run_cell("concentrated: 2 under-discriminating in set 1", two_under, flat)
run_cell("concentrated: 4 over-discriminating in set 1", four_over, flat)

sv_write(do.call(rbind, rows), "alpha-set-misfit")
cat(sprintf("TOTAL elapsed: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
