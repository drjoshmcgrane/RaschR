# Null calibration and power of ordinary between-person DIF with classical
# and HC3 covariance. The reconstruction below follows dif_anova() for a
# single between-person factor, including person aggregation, class intervals,
# and one Holm family over uniform and non-uniform item terms.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "1000"))

simulate_panel <- function(n_a, n_b, repeats = c(1L, 1L), ability_gap = 0,
                           effect = 0) {
  np <- n_a + n_b
  pid <- sprintf("P%05d", seq_len(np))
  grp0 <- factor(c(rep("A", n_a), rep("B", n_b)))
  theta0 <- rnorm(np) + ifelse(grp0 == "B", ability_gap, 0)
  nr <- ifelse(grp0 == "A", repeats[1L], repeats[2L])
  take <- rep(seq_len(np), nr)
  grp <- grp0[take]
  theta <- theta0[take]
  delta <- seq(-1.4, 1.4, length.out = 8L)
  shift <- matrix(0, length(take), length(delta))
  shift[grp == "B", 3L] <- effect
  X <- matrix(rbinom(length(take) * length(delta), 1,
    plogis(outer(theta, delta, "-") - shift)), length(take), length(delta))
  colnames(X) <- paste0("I", seq_along(delta))
  data.frame(X, group = grp, id = pid[take], check.names = FALSE)
}

dif_p <- function(fit, variance) {
  grp <- factor(fit$factors$group)
  id <- as.character(fit$person$id)
  repeated <- anyDuplicated(id) > 0L
  ng <- .dif_n_groups(fit, grp, id = if (repeated) id else NULL)
  ci <- if (repeated) .dif_person_ci(fit, id, ng) else
    .dif_class_intervals(fit, ng)
  out <- matrix(NA_real_, ncol(fit$residuals), 2L,
                dimnames = list(colnames(fit$residuals), c("f1", "f1:ci")))
  for (i in seq_len(ncol(fit$residuals))) {
    d <- data.frame(z = fit$residuals[, i], ci = ci, pid = id, f1 = grp)
    d <- d[stats::complete.cases(d), , drop = FALSE]
    key <- factor(d$pid)
    pz <- tapply(d$z, key, mean)
    first <- which(!duplicated(key))
    pdat <- d[first[match(levels(key), as.character(key[first]))],
              c("pid", "ci", "f1"), drop = FALSE]
    pdat$z <- as.numeric(pz)
    ft <- .dif_type2(pdat, c("f1", "ci", "f1:ci"), variance = variance)
    if (!is.null(ft)) {
      j <- match(colnames(out), ft$term)
      out[i, !is.na(j)] <- ft$p[j[!is.na(j)]]
    }
  }
  out
}

run_scenario <- function(label, n_a, n_b, repeats = c(1L, 1L),
                         ability_gap = 0, effect = 0) {
  per_rep <- matrix(NA_real_, NREP, 8L,
    dimnames = list(NULL, c("raw_classical", "raw_hc3", "fwer_classical",
                            "fwer_hc3", "fwer_hybrid", "power_classical",
                            "power_hc3", "power_hybrid")))
  refused <- nonconv <- 0L
  for (r in seq_len(NREP)) {
    dat <- simulate_panel(n_a, n_b, repeats, ability_gap, effect)
    fit <- tryCatch(rasch(dat, factors = "group", id = "id"),
                    error = function(e) NULL)
    if (is.null(fit)) { refused <- refused + 1L; next }
    if (!isTRUE(fit$est$converged)) { nonconv <- nonconv + 1L; next }
    pc <- dif_p(fit, "classical"); ph <- dif_p(fit, "hc3")
    py <- pc; py[, "f1"] <- ph[, "f1"]
    per_rep[r, "raw_classical"] <- mean(pc[, "f1"] < 0.05, na.rm = TRUE)
    per_rep[r, "raw_hc3"] <- mean(ph[, "f1"] < 0.05, na.rm = TRUE)
    ac <- stats::p.adjust(as.vector(pc), "holm")
    ah <- stats::p.adjust(as.vector(ph), "holm")
    ay <- stats::p.adjust(as.vector(py), "holm")
    per_rep[r, "fwer_classical"] <- any(ac < 0.05, na.rm = TRUE)
    per_rep[r, "fwer_hc3"] <- any(ah < 0.05, na.rm = TRUE)
    per_rep[r, "fwer_hybrid"] <- any(ay < 0.05, na.rm = TRUE)
    per_rep[r, "power_classical"] <- ac[3L] < 0.05
    per_rep[r, "power_hc3"] <- ah[3L] < 0.05
    per_rep[r, "power_hybrid"] <- ay[3L] < 0.05
  }
  ok <- stats::complete.cases(per_rep)
  n <- sum(ok)
  mk <- function(quantity, col, field, note) {
    args <- list(study = "dif-hc3", scenario = label, quantity = quantity,
      n_reps = n, n_attempted = NREP, n_refused = refused,
      n_nonconv = nonconv, effect = effect, notes = note)
    args[[field]] <- mean(per_rep[ok, col])
    if (field == "type1") args$mc_override <- list(
      type1 = stats::sd(per_rep[ok, col]) / sqrt(n))
    do.call(sv_row, args)
  }
  if (effect == 0) rbind(
      mk("classical uniform Type I", "raw_classical", "type1",
         "mean item-wise rejection; MCSE clustered by replicate"),
      mk("HC3 uniform Type I", "raw_hc3", "type1",
         "mean item-wise rejection; MCSE clustered by replicate"),
      mk("classical Holm FWER", "fwer_classical", "familywise",
         "one family over uniform and non-uniform item terms"),
      mk("HC3 Holm FWER", "fwer_hc3", "familywise",
         "diagnostic only: HC3 applied to uniform and class-interval interaction terms"),
      mk("hybrid Holm FWER", "fwer_hybrid", "familywise",
         "public method: HC3 uniform terms and residual-ANOVA non-uniform terms"))
  else rbind(
      mk("classical planted-item power", "power_classical", "power",
         "Holm-adjusted uniform term for I3"),
      mk("HC3 planted-item power", "power_hc3", "power",
         "diagnostic all-HC3 family; Holm-adjusted uniform term for I3"),
      mk("hybrid planted-item power", "power_hybrid", "power",
         "public method; Holm-adjusted uniform term for I3"))
}

set.seed(8.22e7)
rows <- rbind(
  run_scenario("balanced independent persons", 300L, 300L),
  run_scenario("1:4 groups; ability gap 0.8", 120L, 480L,
               ability_gap = 0.8),
  run_scenario("equal groups; 1 versus 4 observations per person", 250L, 250L,
               repeats = c(1L, 4L)),
  run_scenario("1:4 groups; planted DIF 0.6 logits", 120L, 480L,
               ability_gap = 0.8, effect = 0.6))
sv_write(rows, "dif-hc3")
