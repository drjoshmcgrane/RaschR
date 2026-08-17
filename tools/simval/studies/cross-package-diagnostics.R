# STUDY: cross-package-diagnostics
#
# Compares this package's diagnostics against independent implementations
# on the same datasets, at the level each comparison honestly supports:
#   A  item fit (infit/outfit) vs eRm::itemfit -- value-level after
#      aligning the mean-square convention (we divide by sum E[z^2] with
#      the information-share correction; eRm divides by n), plus
#      flag-direction agreement on a planted over-discriminating item
#   B  reliability: alpha vs psych::alpha (formula parity) and PSI vs
#      eRm::SepRel (convention-aligned rebuild: MLE persons, extremes
#      excluded, our calibration)
#   C  uniform DIF vs eRm::Waldtest and difR::difMH -- decision-level,
#      BH-aligned (their defaults are unadjusted)
#   D  residual correlations q3/q3_star vs TAM::tam.modelfit Q3/aQ3 --
#      value-level native, plus the WLE-vs-EAP-aligned recomputation
#   E  person fit residuals vs PerFit lzstar and U3 -- rank-level, with
#      planted careless-person detection
#   F  dimensionality_test vs eRm::NPtest (Tmd, T11) and sirt::conf.detect
#      -- decision-level only (different principles; no comparable values)
# Citation rows name the diagnostics with NO external parallel, which
# remain validated by simulation only.
# Serial. Rscript tools/simval/studies/cross-package-diagnostics.R

suppressWarnings(suppressMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(eRm)
}))
source("tools/simval/harness.R")
STUDY <- "cross-package-diagnostics"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M")),
                          sprintf(...), "\n")
R <- 25L

## no-parallel citation rows -------------------------------------------------
for (cc in list(
  c("dependence_magnitude (Andrich-Kreiner d)", "results/dependence-magnitude-fix.csv"),
  c("spread_test (Andrich 1985 binomial bound)", "results/round1-checks.csv"),
  c("tailored_analysis bootstrap", "results/tailored-bootstrap.csv"),
  c("comparative judgement diagnostics (transitivity, judge fit, guards)", "results/round1-checks.csv; btl-clustered.csv"),
  c("equate_tests multiplicity", "results/equate-familywise.csv")))
  add(cc[1], "no external parallel: simulation-validated only", 0L,
      notes = sprintf("see %s", cc[2]))

## ---- A: item fit vs eRm ---------------------------------------------------
plain_ms <- function(fit) {
  ne <- !fit$person$extreme
  Z <- fit$residuals[ne, , drop = FALSE]
  V <- fit$moments$V[ne, , drop = FALSE]; V[is.na(Z)] <- NA
  list(outfit = colMeans(Z^2, na.rm = TRUE),
       infit  = colSums(Z^2 * V, na.rm = TRUE) / colSums(V, na.rm = TRUE))
}
ai <- ao <- ri <- ro <- flag_agree <- rep(NA_real_, R)
for (r in seq_len(R)) {
  d <- simulate_rasch(500, 12, model = "dichotomous",
                      discrimination = c(3, rep(1, 11)),
                      difficulty = c(0, seq(-2.5, 2.5, length.out = 11)),
                      seed = 2500 + r)
  X <- as.matrix(d[grep("^I", names(d))])
  fit <- rasch(X); pl <- plain_ms(fit)
  it <- eRm::itemfit(eRm::person.parameter(eRm::RM(X)))
  keep <- fit$items$item != "I01"
  ai[r] <- cor(pl$infit[keep],  unname(it$i.infitMSQ)[keep])
  ao[r] <- cor(pl$outfit[keep], unname(it$i.outfitMSQ)[keep])
  ri[r] <- max(abs(pl$infit[keep]  - unname(it$i.infitMSQ)[keep]))
  ro[r] <- max(abs(pl$outfit[keep] - unname(it$i.outfitMSQ)[keep]))
  flag_agree[r] <- (fit$items$fit_resid[1] < -2 ||
                    fit$items$infit_ms[1] < 0.8) &&
                   unname(it$i.infitMSQ)[1] < 0.85
}
add("item fit (infit/outfit) vs eRm::itemfit, planted disc-3 item",
    "value agreement, aligned convention", R,
    notes = sprintf("mean cor infit %.4f / outfit %.4f; mean max|diff| %.4f / %.4f (non-planted items); both flag the planted item in the same direction in %d/%d replicates; reported MSQs sit 8-10%% above eRm's by the E[z^2] convention",
      mean(ai), mean(ao), mean(ri), mean(ro), sum(flag_agree), R))
tick("A: infit cor %.4f outfit cor %.4f flags %d/%d",
     mean(ai), mean(ao), sum(flag_agree), R)

## ---- B: reliability -------------------------------------------------------
psi_formula <- function(theta, se) (var(theta) - mean(se^2)) / var(theta)
da <- ds <- rep(NA_real_, R)
for (r in seq_len(R)) {
  d <- simulate_rasch(500, 12, model = "dichotomous", seed = 2700 + r)
  X <- as.matrix(d[grep("^I", names(d))])
  fit <- rasch(X)
  a_p <- suppressWarnings(psych::alpha(as.data.frame(X), check.keys = FALSE,
                                       warnings = FALSE)$total$raw_alpha)
  da[r] <- abs(fit$alpha$alpha - a_p)
  sr <- eRm::SepRel(eRm::person.parameter(eRm::RM(X)))$sep.rel
  st <- score_table(fit, method = "mle")
  raw <- rowSums(X)
  th <- st$theta[match(raw, st$score)]; se <- st$se[match(raw, st$score)]
  ok <- !is.na(th) & !is.na(se)
  ds[r] <- psi_formula(th[ok], se[ok]) - sr
}
add("alpha vs psych::alpha", "formula parity", R,
    notes = sprintf("max |diff| over %d replicates: %.2e (identical formula)", R, max(da)))
add("PSI vs eRm::SepRel (convention-aligned: MLE persons, extremes excluded)",
    "value agreement", R,
    notes = sprintf("mean diff %+.5f, max |diff| %.5f; native fit$psi$PSI differs by construction (WLE persons, extremes included) and floors at 0 where SepRel can go negative",
      mean(ds), max(abs(ds))))
tick("B: alpha %.1e, PSI aligned max %.5f", max(da), max(abs(ds)))

## ---- C: uniform DIF, decision-level --------------------------------------
det <- matrix(0, 2, 3, dimnames = list(c("dif","null"), c("rasch","eRm","MH")))
fp  <- matrix(0, 2, 3, dimnames = dimnames(det))
for (arm in 1:2) for (r in seq_len(R)) {
  d <- if (arm == 1)
    simulate_rasch(600, 10, dif = list(items = "I05", uniform = 0.8),
                   n_groups = 2, seed = 2900 + r)
  else simulate_rasch(600, 10, n_groups = 2, seed = 3100 + r)
  items <- grep("^I[0-9]+$", names(d), value = TRUE)
  X <- as.matrix(d[items]); grp <- d$group
  s <- dif_anova(rasch(d, id = "id", factors = "group"))$summary
  s <- s[s$term == "group", ]
  f_r <- s$uniform_DIF[match(items, s$item)]
  wt <- tryCatch(eRm::Waldtest(eRm::RM(X), splitcr = grp), error = function(e) NULL)
  f_e <- if (is.null(wt)) rep(NA, 10) else {
    ct <- wt$coef.table
    p <- setNames(ct[, "p-value"], sub("^beta ", "", rownames(ct)))
    p.adjust(p[items], "BH") < 0.05
  }
  mh <- difR::difMH(Data = X, group = as.character(grp),
                    focal.name = levels(grp)[2])
  f_m <- p.adjust(mh$p.value, "BH") < 0.05
  hit <- function(fl) if (arm == 1) isTRUE(fl[5]) else FALSE
  fpn <- function(fl) sum(fl[if (arm == 1) -5 else TRUE], na.rm = TRUE)
  det[arm, ] <- det[arm, ] + c(hit(f_r), hit(f_e), hit(f_m))
  fp[arm, ]  <- fp[arm, ]  + c(fpn(f_r), fpn(f_e), fpn(f_m))
}
add("uniform DIF vs eRm::Waldtest and difR::difMH (BH-aligned)",
    "decision agreement", 2L * R,
    notes = sprintf("planted-item detection over %d DIF replicates: rasch %d, eRm %d, MH %d; false positives per replicate (both arms pooled): rasch %.3f, eRm %.3f, MH %.3f",
      R, det[1,1], det[1,2], det[1,3],
      sum(fp[,1])/(2*R), sum(fp[,2])/(2*R), sum(fp[,3])/(2*R)))
tick("C: detection %d/%d/%d of %d", det[1,1], det[1,2], det[1,3], R)

## ---- D: Q3 vs TAM ---------------------------------------------------------
q3r <- q3ar <- alr <- top_r <- top_t <- rep(NA_real_, R)
for (r in seq_len(R)) {
  d <- simulate_rasch(500, 10, model = "dichotomous",
                      dependence = list(pairs = list(c("I04", "I05")),
                                        strength = 0.6),
                      seed = 3300 + r)
  X <- as.matrix(d[, grep("^I\\d", names(d))])
  fit <- rasch(X)
  rc <- residual_correlations(fit)
  ours <- rc$pairs
  ours$key <- paste(pmin(ours$item_a, ours$item_b),
                    pmax(ours$item_a, ours$item_b))
  mod <- TAM::tam.mml(resp = X, verbose = FALSE)
  mf <- TAM::tam.modelfit(mod, progress = FALSE)
  idx <- which(upper.tri(mf$Q3.matr), arr.ind = TRUE)
  tam <- data.frame(key = paste(pmin(rownames(mf$Q3.matr)[idx[,1]],
                                     colnames(mf$Q3.matr)[idx[,2]]),
                                pmax(rownames(mf$Q3.matr)[idx[,1]],
                                     colnames(mf$Q3.matr)[idx[,2]])),
                    Q3 = mf$Q3.matr[idx], aQ3 = mf$aQ3.matr[idx])
  m <- merge(ours[, c("key", "q3", "q3_star")], tam, by = "key")
  q3r[r] <- cor(m$q3, m$Q3); q3ar[r] <- cor(m$q3_star, m$aQ3)
  delta <- vapply(fit$tau_list, function(t) t[1], 0)
  Ee <- plogis(outer(mod$person$EAP, delta, "-")); Ve <- Ee * (1 - Ee)
  q3e <- cor((X - Ee) / sqrt(Ve))[idx]
  alr[r] <- cor(q3e, tam$Q3)
  pl <- m$key == "I04 I05"
  top_r[r] <- rank(-m$q3_star)[pl] == 1
  top_t[r] <- rank(-m$aQ3)[pl] == 1
}
add("residual correlations q3/q3_star vs TAM Q3/aQ3",
    "value agreement + planted-pair detection", R,
    notes = sprintf("mean cor q3 %.4f, q3_star %.4f native (WLE vs EAP residuals); %.4f when recomputed at TAM's EAP points; planted pair top-1: rasch (q3_star) %d/%d, TAM (aQ3) %d/%d",
      mean(q3r), mean(q3ar), mean(alr), sum(top_r), R, sum(top_t), R))
tick("D: q3 cor %.4f aligned %.4f top1 %d/%d", mean(q3r), mean(alr), sum(top_r), R)

## ---- E: person fit vs PerFit ----------------------------------------------
rlz <- ru3 <- ovl <- det_r <- det_lz <- rep(NA_real_, R)
for (r in seq_len(R)) {
  d <- simulate_rasch(500, 15, model = "dichotomous", careless = 0.10,
                      seed = 3500 + r)
  X <- as.matrix(d[, grep("^I\\d", names(d))])
  cid <- attr(d, "truth")$careless_idx
  fit <- rasch(X)
  pf <- fit$person$fit_resid
  lz <- unlist(PerFit::lzstar(X, IRT.PModel = "1PL")$PFscores)
  u3 <- unlist(PerFit::U3(X)$PFscores)
  ok <- is.finite(pf) & is.finite(lz) & is.finite(u3)
  rlz[r] <- cor(pf[ok], lz[ok], method = "spearman")
  ru3[r] <- cor(pf[ok], u3[ok], method = "spearman")
  nfl <- round(0.1 * sum(ok))
  fl_r  <- order(-pf[ok])[seq_len(nfl)]
  fl_lz <- order(lz[ok])[seq_len(nfl)]
  ovl[r] <- length(intersect(fl_r, fl_lz)) / nfl
  ids <- which(ok)
  det_r[r]  <- mean(cid %in% ids[fl_r])
  det_lz[r] <- mean(cid %in% ids[fl_lz])
}
add("person fit residuals vs PerFit lzstar/U3",
    "rank agreement + careless detection", R,
    notes = sprintf("mean |Spearman| vs lzstar %.3f, vs U3 %.3f; worst-decile overlap %.2f; planted careless caught in worst decile: rasch %.2f vs lzstar %.2f",
      mean(abs(rlz)), mean(abs(ru3)), mean(ovl), mean(det_r), mean(det_lz)))
tick("E: |rho| %.3f overlap %.2f", mean(abs(rlz)), mean(ovl))

## ---- F: dimensionality, decision-level ------------------------------------
dim_cell <- function(n_items, dim2, twod, seed0, Rn) {
  dec <- matrix(0, 4, 1, dimnames = list(c("rasch","Tmd","T11","DETECT"), NULL))
  for (r in seq_len(Rn)) {
    d <- if (twod)
      simulate_rasch(500, n_items, model = "dichotomous",
                     second_dim = list(items = dim2, rho = 0.3),
                     seed = seed0 + r)
    else simulate_rasch(500, n_items, model = "dichotomous", seed = seed0 + r)
    items <- grep("^I[0-9]", names(d), value = TRUE)
    X <- as.matrix(d[, items]); fit <- rasch(d, id = "id")
    dt <- dimensionality_test(fit, items_positive = dim2,
                              items_negative = setdiff(items, dim2))
    Xn <- X[rowSums(X) > 0 & rowSums(X) < ncol(X), ]
    rso <- eRm::rsampler(Xn, eRm::rsctrl(burn_in = 256, n_eff = 200,
                                         step = 32, seed = seed0 + r))
    p_tmd <- eRm::NPtest(rso, method = "Tmd",
                         idx1 = match(setdiff(items, dim2), items),
                         idx2 = match(dim2, items))$prop
    p_t11 <- eRm::NPtest(rso, method = "T11")$prop
    cl <- ifelse(items %in% dim2, 2L, 1L)
    cd <- NULL
    invisible(capture.output(
      cd <- sirt::conf.detect(data = X, score = fit$person$theta,
                              itemcluster = cl, progress = FALSE)))
    dec <- dec + c(isTRUE(dt$multidimensional), p_tmd < 0.05, p_t11 < 0.05,
                   cd$detect["DETECT", "unweighted"] >= 0.2)
  }
  dec / Rn
}
Rd <- 15L
n10 <- dim_cell(10, paste0("I0", 6:9), FALSE, 3700, Rd)
d10 <- dim_cell(10, paste0("I0", 6:9), TRUE,  3800, Rd)
n12 <- dim_cell(12, sprintf("I%02d", 7:12), FALSE, 3900, Rd)
d12 <- dim_cell(12, sprintf("I%02d", 7:12), TRUE,  4000, Rd)
add("dimensionality vs eRm::NPtest (Tmd/T11) and sirt::conf.detect",
    "decision rates (flag proportion per condition)", 4L * Rd,
    notes = sprintf("null 10 items: rasch %.2f Tmd %.2f T11 %.2f DETECT %.2f | 2D 4/6 split: %.2f %.2f %.2f %.2f (short-subtest caution fires; quiet verdicts inconclusive) | null 12: %.2f %.2f %.2f %.2f | 2D 6/6 split: %.2f %.2f %.2f %.2f",
      n10[1], n10[2], n10[3], n10[4], d10[1], d10[2], d10[3], d10[4],
      n12[1], n12[2], n12[3], n12[4], d12[1], d12[2], d12[3], d12[4]))
tick("F: 2D 6/6 rasch %.2f Tmd %.2f T11 %.2f DETECT %.2f",
     d12[1], d12[2], d12[3], d12[4])

sv_write(do.call(rbind, rows), "cross-package-diagnostics")
cat(sprintf("TOTAL elapsed: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
