# Classical versus hybrid DIF inference when the classical assumptions hold
# exactly: balanced crossed cells, independent normal errors, and one common
# residual variance. This isolates the finite-sample cost of using HC3 for the
# uniform factor term when no heteroskedasticity is present.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "5000"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "1")))
NITEM <- 8L

one_replicate <- function(n_cell, levels_group, effect = 0) {
  group <- factor(rep(seq_len(levels_group), each = 4L * n_cell))
  ci <- factor(rep(rep(seq_len(4L), each = n_cell), levels_group))
  ci_mean <- c(-0.45, -0.15, 0.15, 0.45)[ci]
  Y <- matrix(rnorm(length(group) * NITEM), nrow = length(group),
              ncol = NITEM)
  Y <- Y + ci_mean
  if (effect != 0)
    Y[group == levels_group, 3L] <-
      Y[group == levels_group, 3L] + effect

  pc <- ph <- matrix(NA_real_, NITEM, 2L,
    dimnames = list(NULL, c("group", "group:ci")))
  for (i in seq_len(NITEM)) {
    d <- data.frame(z = Y[, i], group = group, ci = ci)
    fc <- .dif_type2(d, c("group", "ci", "group:ci"),
                     variance = "classical")
    fh <- .dif_type2(d, c("group", "ci", "group:ci"),
                     variance = "hc3", robust_terms = "group")
    pc[i, ] <- fc$p[match(colnames(pc), fc$term)]
    ph[i, ] <- fh$p[match(colnames(ph), fh$term)]
  }
  ac <- stats::p.adjust(as.vector(pc), "holm")
  ah <- stats::p.adjust(as.vector(ph), "holm")
  c(classical_raw = mean(pc[, "group"] < 0.05),
    hybrid_raw = mean(ph[, "group"] < 0.05),
    classical_fwer = any(ac < 0.05),
    hybrid_fwer = any(ah < 0.05),
    classical_power = ac[3L] < 0.05,
    hybrid_power = ah[3L] < 0.05,
    p_disagree = (ac[3L] < 0.05) != (ah[3L] < 0.05))
}

run_condition <- function(n_cell, levels_group, effect = 0) {
  seeds <- sample.int(.Machine$integer.max, NREP)
  z <- do.call(rbind, parallel::mclapply(seq_len(NREP), function(r) {
    set.seed(seeds[r])
    one_replicate(n_cell, levels_group, effect)
  }, mc.cores = NCORE, mc.preschedule = TRUE))
  scenario <- sprintf("%d group levels; %d per group-by-interval cell; N=%d",
                      levels_group, n_cell,
                      levels_group * 4L * n_cell)
  mk <- function(quantity, col, field, note) {
    x <- z[, col]
    args <- list(study = "dif-hc3-homoskedastic", scenario = scenario,
      quantity = quantity, n_reps = NREP, n_attempted = NREP,
      n_refused = 0L, n_nonconv = 0L, effect = effect,
      notes = paste("balanced cells; iid N(0,1) errors;", note))
    args[[field]] <- mean(x)
    if (field == "type1") args$mc_override <- list(
      type1 = stats::sd(x) / sqrt(NREP))
    do.call(sv_row, args)
  }
  if (effect == 0) {
    rbind(
      mk("classical uniform Type I", "classical_raw", "type1",
         "mean item-wise rejection"),
      mk("hybrid uniform Type I", "hybrid_raw", "type1",
         "HC3 for uniform terms"),
      mk("classical Holm FWER", "classical_fwer", "familywise",
         "one family over uniform and non-uniform terms"),
      mk("hybrid Holm FWER", "hybrid_fwer", "familywise",
         "HC3 uniform; classical non-uniform"),
      mk("adjusted-decision disagreement", "p_disagree", "power",
         "I3 decision differs between methods under the null"))
  } else {
    rbind(
      mk("classical planted-item power", "classical_power", "power",
         "Holm-adjusted uniform term for I3"),
      mk("hybrid planted-item power", "hybrid_power", "power",
         "Holm-adjusted HC3 uniform term for I3"),
      mk("adjusted-decision disagreement", "p_disagree", "power",
         "I3 decision differs between methods under the alternative"))
  }
}

set.seed(8.27e7)
designs <- list(c(2L, 10L), c(2L, 30L), c(2L, 75L), c(2L, 150L),
                c(3L, 10L), c(3L, 50L))
rows <- do.call(rbind, lapply(designs, function(x) rbind(
  run_condition(n_cell = x[2L], levels_group = x[1L]),
  run_condition(n_cell = x[2L], levels_group = x[1L], effect = 0.5))))
sv_write(rows, "dif-hc3-homoskedastic")
