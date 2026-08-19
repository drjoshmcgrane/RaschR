# STUDY: humphry-pgd-replication
#
# Simulate the design of the person-group-discrimination study in Humphry
# (2005, ch. 4) and check whether its estimator recovers the planted unit
# ratio.
#
# The design as reported: WALNA 2003 Reading, 12 common items linking Year
# 5 and Year 7, calibrated separately within each year level, the unit
# ratio read off as the ratio of the standard deviations of the common
# item locations. Reported values: SD ratio 1.22 in the full populations
# and about 1.3 in random samples of roughly 980 per year, giving
# phi_5 = 0.875 and phi_7 = 1.143 under the product constraint -- a ratio
# of 1.306, which is the planted truth here. Common-scale item locations
# are taken from his Table 4.1 Year 5 estimates divided by phi_5.
#
# The feature this design has and earlier studies did not: the two frames
# differ in ability as well as unit, so the common items are differentially
# targeted. That asymmetry is a candidate source of bias in a ratio of
# dispersions, so the ability gap is swept rather than assumed.
# Serial. Rscript tools/simval/studies/humphry-pgd-replication.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "humphry-pgd-replication"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

phi <- c(yr5 = 0.875, yr7 = 1.143)
truth <- unname(phi["yr7"] / phi["yr5"])
lt <- log(truth)
# Table 4.1 Year 5 locations, put on the common scale
yr5_loc <- c(-0.60, -0.28, -0.68, -0.01, -0.12, -0.10,
              0.78,  0.33, -0.09, -1.30,  1.43,  0.65)
delta_common <- yr5_loc / phi["yr5"]
K_common <- length(delta_common)
K_unique <- 18L          # each year also takes its own items (30-item test)
R <- 200L

for (N in c(980L, 5000L)) {
  for (gap in c(0, 0.5, 1.0)) {
    lr <- rep(NA_real_, R); n_ref <- 0L
    for (r in seq_len(R)) {
      set.seed(150000 + N + round(gap * 100) + r)
      loc <- vector("list", 2)
      for (g in 1:2) {
        # ability on the common scale; Year 7 sits `gap` logits higher
        th <- rnorm(N, if (g == 1) 0 else gap, 1.2)
        d_uniq <- seq(-1.5, 1.5, length.out = K_unique) +
          if (g == 1) 0 else gap
        d_all <- c(delta_common, d_uniq)
        X <- vapply(d_all, function(d)
          rbinom(N, 1, plogis(phi[g] * (th - d))), numeric(N))
        colnames(X) <- c(sprintf("C%02d", seq_len(K_common)),
                         sprintf("U%02d", seq_len(K_unique)))
        f <- tryCatch(rasch(X), error = function(e) NULL)
        if (is.null(f)) { loc[[g]] <- NULL; break }
        # the ratio is read off the COMMON items only, as in the source
        loc[[g]] <- f$items$location[seq_len(K_common)]
      }
      if (any(vapply(loc, is.null, TRUE))) { n_ref <- n_ref + 1L; next }
      lr[r] <- log(sd(loc[[2]]) / sd(loc[[1]]))
    }
    ok <- is.finite(lr)
    add(sprintf("12 common items, N = %d per year, ability gap %.1f logits", N, gap),
        "log unit ratio, raw SD ratio of common-item locations", sum(ok),
        bias = mean(lr[ok]) - lt, emp_sd = sd(lr[ok]),
        effect = gap, n_attempted = R, n_refused = n_ref,
        notes = sprintf("planted phi ratio %.3f; recovered mean %.3f",
                        truth, exp(mean(lr[ok]))))
    cat(sprintf("[%s] N = %4d, gap %.1f: bias %+.4f (MC SE %.4f), recovered %.3f vs planted %.3f, sd %.3f\n",
                format(Sys.time(), "%H:%M"), N, gap, mean(lr[ok]) - lt,
                sd(lr[ok]) / sqrt(sum(ok)), exp(mean(lr[ok])), truth, sd(lr[ok])))
  }
}

sv_write(do.call(rbind, rows), "humphry-pgd-replication")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
