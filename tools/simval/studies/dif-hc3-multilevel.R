# Three-level check of the ordinary DIF HC3 term. The main study covers binary
# person factors; this study verifies the multi-degree-of-freedom Wald test.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "500"))

simulate_data <- function(effect = 0) {
  ns <- c(100L, 200L, 300L)
  grp <- factor(rep(c("A", "B", "C"), ns))
  theta <- rnorm(sum(ns)) + c(A = -0.6, B = 0, C = 0.6)[grp]
  delta <- seq(-1.4, 1.4, length.out = 8L)
  shift <- matrix(0, length(theta), 8L)
  shift[grp == "C", 3L] <- effect
  X <- matrix(rbinom(length(theta) * 8L, 1,
    plogis(outer(theta, delta, "-") - shift)), length(theta), 8L)
  colnames(X) <- paste0("I", 1:8)
  data.frame(X, group = grp)
}

term_p <- function(fit, variance) {
  grp <- factor(fit$factors$group)
  ci <- .dif_class_intervals(fit,
    .dif_n_groups(fit, grp))
  ans <- matrix(NA_real_, ncol(fit$residuals), 2L,
                dimnames = list(colnames(fit$residuals), c("f1", "f1:ci")))
  for (i in seq_len(nrow(ans))) {
    d <- data.frame(z = fit$residuals[, i], f1 = grp, ci = ci)
    d <- d[stats::complete.cases(d), , drop = FALSE]
    z <- .dif_type2(d, c("f1", "ci", "f1:ci"), variance = variance)
    j <- match(colnames(ans), z$term)
    ans[i, ] <- z$p[j]
  }
  ans
}

run_condition <- function(effect) {
  z <- matrix(NA_real_, NREP, 5L,
    dimnames = list(NULL, c("classical_raw", "hc3_raw", "classical_fwer",
                            "full_hc3_fwer", "hybrid_fwer")))
  power <- matrix(NA_real_, NREP, 3L,
    dimnames = list(NULL, c("classical", "full_hc3", "hybrid")))
  for (r in seq_len(NREP)) {
    fit <- tryCatch(rasch(simulate_data(effect), factors = "group"),
                    error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$est$converged)) next
    pc <- term_p(fit, "classical"); ph <- term_p(fit, "hc3")
    py <- pc; py[, "f1"] <- ph[, "f1"]
    ac <- stats::p.adjust(as.vector(pc), "holm")
    ah <- stats::p.adjust(as.vector(ph), "holm")
    ay <- stats::p.adjust(as.vector(py), "holm")
    z[r, ] <- c(mean(pc[, "f1"] < 0.05), mean(ph[, "f1"] < 0.05),
                 any(ac < 0.05), any(ah < 0.05), any(ay < 0.05))
    power[r, ] <- c(ac[3L] < 0.05, ah[3L] < 0.05, ay[3L] < 0.05)
  }
  if (effect == 0) {
    mk <- function(q, col, field) {
      ok <- is.finite(z[, col]); n <- sum(ok)
      args <- list(study = "dif-hc3-multilevel", scenario =
        "three levels; 1:2:3 groups; ability gaps", quantity = q,
        n_reps = n, n_attempted = NREP, n_refused = NREP - n,
        n_nonconv = 0L, effect = effect,
        notes = "three-level person factor; Holm family spans uniform and non-uniform item terms")
      args[[field]] <- mean(z[ok, col])
      if (field == "type1") args$mc_override <- list(type1 =
        stats::sd(z[ok, col]) / sqrt(n))
      do.call(sv_row, args)
    }
    rbind(mk("classical uniform Type I", "classical_raw", "type1"),
          mk("HC3 uniform Type I", "hc3_raw", "type1"),
          mk("classical Holm FWER", "classical_fwer", "familywise"),
          mk("full-HC3 Holm FWER", "full_hc3_fwer", "familywise"),
          mk("hybrid Holm FWER", "hybrid_fwer", "familywise"))
  } else {
    do.call(rbind, lapply(colnames(power), function(nm) {
      ok <- is.finite(power[, nm]); n <- sum(ok)
      sv_row("dif-hc3-multilevel",
        "three levels; 1:2:3 groups; planted C shift 0.6",
        paste(nm, "planted-item power"), n,
        power = mean(power[ok, nm]), effect = effect,
        n_attempted = NREP, n_refused = NREP - n, n_nonconv = 0L,
        notes = "Holm-adjusted uniform term for I3")
    }))
  }
}

set.seed(8.26e7)
sv_write(rbind(run_condition(0), run_condition(0.6)), "dif-hc3-multilevel")
