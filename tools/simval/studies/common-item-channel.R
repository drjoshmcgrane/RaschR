# STUDY: common-item-channel
#
# When two frames share items, the unit ratio is identified three ways from
# the same data. Which channel should a design aim for?
#
#   cml    the bilinear within-frame conditional fit that rasch_efrm already
#          uses for person-group units: the same items calibrated at two
#          scales, thresholds and units estimated jointly
#   sd     the ratio of the standard deviations of the two frames' item
#          location estimates (Humphry 2005, eq. 2.27)
#   slope  the through-origin least-squares slope of one frame's item
#          locations on the other's
#
# The comparison matters beyond phi. Our item sets partition the items, so
# alpha has no common-item route and must go through the person-side link,
# which carries the error problem the truncated-score-moment correction
# fixes. A bridge design -- some items administered in both set contexts --
# would put alpha on whichever channel wins here, since the identification
# is the same bilinear structure.
#
# Planted group-unit ratio 1.30, dichotomous, independent persons per frame.
# Serial. Rscript tools/simval/studies/common-item-channel.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "common-item-channel"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

R <- 40L
ratio <- 1.30
lt <- log(ratio)

for (K in c(8L, 20L, 40L)) {
  for (N in c(250L, 1000L)) {
    est <- matrix(NA_real_, R, 3,
                  dimnames = list(NULL, c("cml", "sd", "slope")))
    n_ref <- 0L
    for (r in seq_len(R)) {
      d <- simulate_efrm(n_per_group = N, items_per_set = K, n_sets = 1,
                         n_groups = 2, set_unit_ratio = 1,
                         group_unit_ratio = ratio,
                         seed = 130000 + K * 10000 + N + r)
      tr <- attr(d, "truth")
      items <- unlist(tr$item_sets, use.names = FALSE)

      # channel 1: the package's conditional bilinear fit
      f <- tryCatch(rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                               id = "id", boot_reps = 0),
                    error = function(e) NULL)
      if (is.null(f)) { n_ref <- n_ref + 1L; next }
      p <- f$phi_table$phi
      est[r, "cml"] <- max(p) / min(p)

      # channels 2 and 3: separate within-frame calibrations of the same items
      loc <- lapply(c("g1", "g2"), function(g) {
        X <- as.matrix(d[d$group == g, items])
        ff <- tryCatch(rasch(X), error = function(e) NULL)
        if (is.null(ff)) return(NULL)
        ff$items$location
      })
      if (any(vapply(loc, is.null, TRUE))) { n_ref <- n_ref + 1L; next }
      a <- loc[[1]]; b <- loc[[2]]
      est[r, "sd"] <- sd(b) / sd(a)
      est[r, "slope"] <- sum(a * b) / sum(a * a)
    }
    for (m in colnames(est)) {
      ok <- is.finite(est[, m])
      add(sprintf("%d common items, N = %d per frame, unit ratio 1.30", K, N),
          sprintf("log unit ratio, %s channel", m), sum(ok),
          bias = mean(log(est[ok, m])) - lt,
          emp_sd = sd(log(est[ok, m])),
          effect = N, n_attempted = R, n_refused = n_ref,
          notes = sprintf("items %d", K))
    }
    fmt <- function(m) sprintf("%+.4f (sd %.4f)", mean(log(est[, m]), na.rm = TRUE) - lt,
                              sd(log(est[, m]), na.rm = TRUE))
    cat(sprintf("[%s] K = %2d, N = %4d: cml %s | sd %s | slope %s\n",
                format(Sys.time(), "%H:%M"), K, N,
                fmt("cml"), fmt("sd"), fmt("slope")))
  }
}

sv_write(do.call(rbind, rows), "common-item-channel")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
