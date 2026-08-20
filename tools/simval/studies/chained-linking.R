# STUDY: chained-linking
#
# Does a unit ratio recover when the two frames share no items?
#
# A vertical linking design rarely gives every pair of year levels a common
# block. Year 3 and year 5 share one anchor, year 5 and year 7 share a
# different one, and years 3 and 7 share nothing. The model accepts this:
# identification of the person-group units runs over the CONNECTED
# COMPONENTS of the graph whose edges are group pairs sharing at least two
# items within one set, so a chain is enough and only a broken chain is
# refused.
#
# Accepting a design is not the same as recovering from it. The year 3 to
# year 7 comparison is made entirely through year 5, so any error in the two
# estimated links compounds into it, and there is no direct evidence to
# correct it. A single pilot replicate recovered the directly linked ratio
# well and the chained one badly, which is either sampling noise or the
# thing worth knowing about these designs.
#
# Three year levels, units 0.80 / 1.00 / 1.25 before centring, 5 anchor items
# per link, 5 unique items per year, N per year as stated. Reported per cell:
# bias on the log of the DIRECT ratio (year5/year3, one shared anchor) and on
# the log of the INDIRECT ratio (year7/year3, no shared items at all), with
# the mean reported standard error of each so calibration can be judged
# alongside bias.
#
# 200 replicates. Serial, ~20 min.
#   Rscript tools/simval/studies/chained-linking.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "chained-linking"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

phi0 <- c(year3 = 0.80, year5 = 1.00, year7 = 1.25)
phi0 <- phi0 / exp(mean(log(phi0)))
lt_direct <- log(phi0[["year5"]] / phi0[["year3"]])
lt_indirect <- log(phi0[["year7"]] / phi0[["year3"]])
R <- 200L

A <- sprintf("A%02d", 1:5); B <- sprintf("B%02d", 1:5)
U3 <- sprintf("U3_%02d", 1:5); U5 <- sprintf("U5_%02d", 1:4)
U7 <- sprintf("U7_%02d", 1:5)
items <- c(A, B, U3, U5, U7)
takes <- list(year3 = c(U3, A), year5 = c(A, B, U5), year7 = c(B, U7))

gen <- function(seed, N) {
  set.seed(seed)
  delta <- stats::setNames(stats::rnorm(length(items), 0, 1.1), items)
  out <- list()
  for (g in names(takes)) {
    th <- stats::rnorm(N, 0, 1.2)
    M <- matrix(NA_real_, N, length(items), dimnames = list(NULL, items))
    for (it in takes[[g]])
      M[, it] <- stats::rbinom(N, 1, stats::plogis(
        phi0[[g]] * (th - delta[[it]])))
    out[[g]] <- data.frame(M, group = g, check.names = FALSE)
  }
  d <- do.call(rbind, out)
  d$id <- sprintf("P%05d", seq_len(nrow(d)))
  d
}

run_cell <- function(N) {
  lab <- sprintf("chain 3-5-7, N = %d per year", N)
  ld <- li <- sd_ <- si <- rep(NA_real_, R)
  n_ref <- 0L
  for (r in seq_len(R)) {
    f <- tryCatch(rasch_efrm(gen(700000 + r, N), items = items,
                             item_sets = list(set1 = items), groups = "group",
                             id = "id", boot_reps = 0),
                  error = function(e) NULL)
    if (is.null(f)) { n_ref <- n_ref + 1L; next }
    i <- match(names(phi0), f$phi_table$group)
    p <- f$phi_table$phi[i]
    s <- f$phi_table$se_log_phi[i]
    ld[r] <- log(p[2] / p[1]); li[r] <- log(p[3] / p[1])
    # the reported standard error of a ratio of two independent-ish units
    sd_[r] <- sqrt(s[2]^2 + s[1]^2); si[r] <- sqrt(s[3]^2 + s[1]^2)
  }
  ok <- is.finite(ld) & is.finite(li)
  add(paste0(lab, " | direct (year5/year3, shared anchor)"), "log phi ratio",
      sum(ok), bias = mean(ld[ok]) - lt_direct, emp_sd = stats::sd(ld[ok]),
      mean_se = mean(sd_[ok], na.rm = TRUE),
      n_attempted = R, n_refused = n_ref,
      notes = sprintf("recovered %.3f vs planted %.3f", exp(mean(ld[ok])),
                      exp(lt_direct)))
  add(paste0(lab, " | indirect (year7/year3, no shared item)"),
      "log phi ratio", sum(ok), bias = mean(li[ok]) - lt_indirect,
      emp_sd = stats::sd(li[ok]), mean_se = mean(si[ok], na.rm = TRUE),
      n_attempted = R, n_refused = n_ref,
      notes = sprintf("recovered %.3f vs planted %.3f", exp(mean(li[ok])),
                      exp(lt_indirect)))
  cat(sprintf("[%s] N=%-5d direct %.3f (bias %+.4f, sd %.3f, se %.3f) | indirect %.3f (bias %+.4f, sd %.3f, se %.3f)\n",
              format(Sys.time(), "%H:%M"), N, exp(mean(ld[ok])),
              mean(ld[ok]) - lt_direct, stats::sd(ld[ok]),
              mean(sd_[ok], na.rm = TRUE), exp(mean(li[ok])),
              mean(li[ok]) - lt_indirect, stats::sd(li[ok]),
              mean(si[ok], na.rm = TRUE)))
}

for (N in c(300L, 700L, 2000L)) run_cell(N)

sv_write(do.call(rbind, rows), "chained-linking")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
