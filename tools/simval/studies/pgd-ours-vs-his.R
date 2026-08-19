# STUDY: pgd-ours-vs-his
#
# Apples to apples on a common-item linking design. Earlier studies set
# Humphry's item-side estimator against this package's PERSON-side link,
# which is the wrong pairing: given common items and disjoint person
# groups, this package does not use the person side at all. It estimates
# the person-group unit phi by within-frame conditional maximum
# likelihood on the common items -- the same information his SD ratio
# uses, read a different way.
#
# Both estimators are therefore run on his design: the 12 WALNA common
# items of Table 4.1, N = 980 per year level, an ability gap of 0.5
# logits, phi_5 = 0.875 and phi_7 = 1.143 (ratio 1.306).
#
#   his    ratio of the SDs of the item locations, calibrated separately
#          within each year level
#   ours   rasch_efrm phi ratio: one item set, two person groups,
#          conditional ML on the bilinear threshold structure
#
# and against the departures that hurt the SD ratio, to see whether using
# the full threshold pattern rather than its dispersion resists them:
# uniform DIF on two items, and differential discrimination on two items.
# Serial. Rscript tools/simval/studies/pgd-ours-vs-his.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "pgd-ours-vs-his"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

phi <- c(0.875, 1.143)
truth <- phi[2] / phi[1]
lt <- log(truth)
yr5_loc <- c(-0.60, -0.28, -0.68, -0.01, -0.12, -0.10,
              0.78,  0.33, -0.09, -1.30,  1.43,  0.65)
delta <- yr5_loc / phi[1]
K <- length(delta)
N <- 980L
gap <- 0.5
R <- 100L

run_cell <- function(label, kind, size, n_aff) {
  m <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("his", "ours")))
  n_ref <- 0L
  aff <- if (n_aff > 0) seq_len(n_aff) else integer(0)
  sgn <- rep(c(1, -1), length.out = max(n_aff, 1))
  for (r in seq_len(R)) {
    set.seed(200000 + n_aff * 700 + round(size * 100) + r)
    dat <- list(); locs <- list()
    for (g in 1:2) {
      th <- rnorm(N, if (g == 1) 0 else gap, 1.2)
      d <- delta; disc <- rep(1, K)
      if (g == 2 && n_aff > 0) {
        if (kind == "dif") d[aff] <- d[aff] + sgn * size
        if (kind == "disc") disc[aff] <- size
      }
      X <- vapply(seq_len(K), function(i)
        rbinom(N, 1, plogis(phi[g] * disc[i] * (th - d[i]))), numeric(N))
      colnames(X) <- sprintf("C%02d", seq_len(K))
      dat[[g]] <- X
      f <- tryCatch(rasch(X), error = function(e) NULL)
      locs[[g]] <- if (is.null(f)) NULL else f$items$location
    }
    if (any(vapply(locs, is.null, TRUE))) { n_ref <- n_ref + 1L; next }
    m[r, "his"] <- log(sd(locs[[2]]) / sd(locs[[1]]))
    # ours: one item set taken by two person groups, conditional ML
    d_all <- data.frame(id = sprintf("P%05d", seq_len(2L * N)),
                        rbind(dat[[1]], dat[[2]]),
                        group = rep(c("yr5", "yr7"), each = N),
                        check.names = FALSE)
    f <- tryCatch(rasch_efrm(d_all, item_sets = list(common = colnames(dat[[1]])),
                             groups = "group", id = "id", boot_reps = 0),
                  error = function(e) NULL)
    if (is.null(f)) { n_ref <- n_ref + 1L; next }
    p <- f$phi_table$phi
    m[r, "ours"] <- log(p[2] / p[1])
  }
  for (ch in colnames(m)) {
    ok <- is.finite(m[, ch])
    add(label, sprintf("log phi ratio, %s estimator", ch), sum(ok),
        bias = mean(m[ok, ch]) - lt, emp_sd = sd(m[ok, ch]),
        n_attempted = R, n_refused = n_ref,
        notes = sprintf("recovered %.3f vs planted %.3f",
                        exp(mean(m[ok, ch])), truth))
  }
  oh <- is.finite(m[, "his"]); oo <- is.finite(m[, "ours"])
  cat(sprintf("[%s] %-34s his %.3f (bias %+.4f, sd %.4f) | ours %.3f (bias %+.4f, sd %.4f)\n",
      format(Sys.time(), "%H:%M"), label,
      exp(mean(m[oh, "his"])), mean(m[oh, "his"]) - lt, sd(m[oh, "his"]),
      exp(mean(m[oo, "ours"])), mean(m[oo, "ours"]) - lt, sd(m[oo, "ours"])))
}

run_cell("clean", "none", 0, 0L)
run_cell("2 of 12 items, uniform DIF 1.0", "dif", 1.0, 2L)
run_cell("2 of 12 items, discrimination x1.5", "disc", 1.5, 2L)
run_cell("4 of 12 items, discrimination x1.5", "disc", 1.5, 4L)

sv_write(do.call(rbind, rows), "pgd-ours-vs-his")
cat(sprintf("TOTAL elapsed: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
