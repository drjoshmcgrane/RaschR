# STUDY: misfit-both-channels
#
# Item misfit that is NOT frame-specific: items that over- or
# under-discriminate relative to the Rasch model, equally in both frames.
# This is ordinary misfit rather than DIF, and it is the common case in
# real instruments -- the wording case study's Q8 discriminates at 0.89
# where Q6 discriminates at 2.42 on the same scale.
#
# The distortion should partly cancel in a ratio of dispersions, since
# both frames carry the same aberrant items. It cannot cancel exactly:
# the frames apply different units to the same slope, so the effective
# discrimination rho_g * a_i differs and the Rasch calibration absorbs it
# differently in each.
#
# Both channels are measured on the same data:
#   item-side    ratio of the SDs of the item locations (Humphry 2005)
#   person-side  corrected true-score variance ratio over common persons
#                (what this package uses for item-set units)
#
# 12 items and 980 persons, matching the published design; planted unit
# ratio 1.30. Serial.
#   Rscript tools/simval/studies/misfit-both-channels.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "misfit-both-channels"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()

ratio <- 1.30
rho <- c(ratio^-0.5, ratio^0.5)
lt <- log(ratio)
K <- 12L
N <- 980L
delta <- seq(-2, 2, length.out = K)
R <- 100L

run_cell <- function(label, disc_fun) {
  m <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("item", "person")))
  n_ref <- 0L
  for (r in seq_len(R)) {
    set.seed(180000 + nchar(label) * 977 + r)
    a <- disc_fun(K)                       # SAME discriminations in both frames
    th <- rnorm(N, 0, 1.3)
    fits <- lapply(1:2, function(g) {
      X <- vapply(seq_len(K), function(i)
        rbinom(N, 1, plogis(rho[g] * a[i] * (th - delta[i]))), numeric(N))
      colnames(X) <- sprintf("I%02d", seq_len(K))
      tryCatch(rasch(X), error = function(e) NULL)
    })
    if (any(vapply(fits, is.null, TRUE))) { n_ref <- n_ref + 1L; next }
    m[r, "item"] <- log(sd(fits[[2]]$items$location) /
                        sd(fits[[1]]$items$location))
    ok <- !fits[[1]]$person$extreme & !fits[[2]]$person$extreme
    if (sum(ok) < 30) { n_ref <- n_ref + 1L; next }
    v <- vapply(1:2, function(g) {
      f <- fits[[g]]
      uh <- f$person$theta[ok]
      lm <- .person_link_moments(as.matrix(f$X), f$tau_list)
      (var(uh) - mean(lm$w[ok])) / mean(lm$g[ok])^2
    }, 0)
    if (any(!is.finite(v)) || any(v <= 0)) { n_ref <- n_ref + 1L; next }
    m[r, "person"] <- 0.5 * (log(v[2]) - log(v[1]))
  }
  for (ch in colnames(m)) {
    ok <- is.finite(m[, ch])
    add(label, sprintf("log unit ratio, %s channel", ch), sum(ok),
        bias = mean(m[ok, ch]) - lt, emp_sd = sd(m[ok, ch]),
        n_attempted = R, n_refused = n_ref,
        notes = sprintf("recovered %.3f vs planted %.3f",
                        exp(mean(m[ok, ch])), ratio))
  }
  oi <- is.finite(m[, "item"]); op <- is.finite(m[, "person"])
  cat(sprintf("[%s] %-38s item %+.4f (%.3f) | person %+.4f (%.3f)\n",
              format(Sys.time(), "%H:%M"), label,
              mean(m[oi, "item"]) - lt, exp(mean(m[oi, "item"])),
              mean(m[op, "person"]) - lt, exp(mean(m[op, "person"]))))
}

flat <- function(K) rep(1, K)
two_over <- function(K) { a <- rep(1, K); a[c(3, 8)] <- 2.0; a }
two_under <- function(K) { a <- rep(1, K); a[c(3, 8)] <- 0.5; a }
mixed <- function(K) { a <- rep(1, K); a[c(3, 8)] <- 2.0; a[c(5, 10)] <- 0.5; a }
four_over <- function(K) { a <- rep(1, K); a[c(2, 5, 8, 11)] <- 2.0; a }
four_under <- function(K) { a <- rep(1, K); a[c(2, 5, 8, 11)] <- 0.5; a }
scatter <- function(K) exp(rnorm(K, 0, 0.25))     # every item slightly off

run_cell("clean (all discriminations 1)", flat)
run_cell("2 of 12 over-discriminating (x2)", two_over)
run_cell("2 of 12 under-discriminating (x0.5)", two_under)
run_cell("2 over and 2 under", mixed)
run_cell("4 of 12 over-discriminating (x2)", four_over)
run_cell("4 of 12 under-discriminating (x0.5)", four_under)
run_cell("all items scattered (log-normal sd 0.25)", scatter)

sv_write(do.call(rbind, rows), "misfit-both-channels")
cat(sprintf("TOTAL elapsed: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))
