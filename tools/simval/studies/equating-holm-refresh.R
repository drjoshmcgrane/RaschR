# STUDY: equating-holm-refresh
#
# equate_tests() switched its drift adjustment from BH to Holm in the
# multiplicity standardisation, and the familywise calibration on record
# (equating-multiplicity.csv) predates the switch: its rows are stamped
# with the BH-adjusted construction and must not be cited for the current
# function. This study re-runs the null familywise cells under the code as
# it stands, mirroring the original run_null design exactly: two 12-item
# forms, k anchors named A1..Ak with identical true locations on
# seq(-1.2, 1.2), form-specific unique items with distinct names, persons
# N(0, 1), N = 500 per form, independent calibrations declared. Familywise
# error = any anchor flagged.
#   Rscript tools/simval/studies/equating-holm-refresh.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "equating-holm-refresh"
rows <- list()
t0 <- Sys.time()

N <- 500L
R <- 2000L
sim_form <- function(delta_named, N) {
  I <- length(delta_named)
  theta <- rnorm(N, 0, 1)
  X <- matrix(rbinom(N * I, 1, plogis(outer(theta, delta_named, "-"))), N, I)
  colnames(X) <- names(delta_named)
  X
}

for (k in c(3L, 5L, 10L)) {
  anchors_true <- setNames(seq(-1.2, 1.2, length.out = k), paste0("A", seq_len(k)))
  n_uniq <- 12L - k
  uniq1 <- setNames(seq(-2, 2, length.out = n_uniq), paste0("U1_", seq_len(n_uniq)))
  uniq2 <- setNames(seq(-2.3, 2.3, length.out = n_uniq), paste0("U2_", seq_len(n_uniq)))
  d1 <- c(anchors_true, uniq1); d2 <- c(anchors_true, uniq2)
  anames <- names(anchors_true)

  fam <- rep(NA, R)
  n_refused <- 0L
  for (r in seq_len(R)) {
    set.seed(910000 + k * 100000 + r)
    f1 <- tryCatch(suppressWarnings(rasch(sim_form(d1, N))), error = function(e) NULL)
    f2 <- tryCatch(suppressWarnings(rasch(sim_form(d2, N))), error = function(e) NULL)
    if (is.null(f1) || is.null(f2) ||
        !isTRUE(f1$est$converged) || !isTRUE(f2$est$converged)) {
      n_refused <- n_refused + 1L; next
    }
    eq <- tryCatch(equate_tests(f1, f2, independent = TRUE),
                   error = function(e) NULL)
    if (is.null(eq) || !isTRUE(eq$inferential)) { n_refused <- n_refused + 1L; next }
    fl <- eq$table$drift[match(anames, eq$table$item)]
    if (anyNA(fl)) { n_refused <- n_refused + 1L; next }
    fam[r] <- any(fl)
  }
  ok <- !is.na(fam)
  rows[[length(rows) + 1L]] <- sv_row(STUDY,
    sprintf("null k=%d (Holm)", k), "familywise flag", sum(ok),
    familywise = mean(fam[ok]), n_attempted = R, n_refused = n_refused,
    notes = "Holm-adjusted alpha=0.05, current equate_tests; anchors truth identical both forms")
  cat(sprintf("[%s] k=%-2d familywise %.4f (n=%d, refused %d)\n",
              format(Sys.time(), "%H:%M"), k, mean(fam[ok]), sum(ok), n_refused))
}
sv_write(do.call(rbind, rows), "equating-holm-refresh")
cat(sprintf("TOTAL %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
