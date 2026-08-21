# STUDY: cross-package-validation
#
# Compares this package's estimates against independent accepted
# implementations on the same datasets:
#   A  rasch() dichotomous locations  vs sirt::rasch.pairwise (pairwise
#      minimum chi-square) and eRm::RM (full conditional ML)
#   B  rasch(model="PCM") locations + thresholds vs eRm::PCM (full CML)
#   C  rasch_efrm item-set unit ratio vs TAM::tam.mml.2pl with
#      irtmodel="2PL.groups" (slopes constrained equal within set; MML)
#   D  btl() dichotomous vs BradleyTerry2::BTm (identical likelihood ->
#      machine-precision parity expected)
#   E  btl() polytomous (free and pc thresholds) vs VGAM::vglm acat with
#      btl's symmetric threshold constraint imposed (identical model ->
#      machine-precision parity expected)
#   F  btl_efrm panel-unit ratio vs per-panel intercept-free
#      vglm(binomialff) fits (through-origin LS slope of location vectors)
#   G  rasch_mfrm item locations, rater severities, and PCM thresholds vs
#      TAM::tam.mml.mfr (MML; formulaA = ~ item + rater + item:step)
#   H  rasch_efrm person-group unit (phi) vs per-group lme4::glmer Rasch
#      fits: through-origin slope of centred item coefficients (primary,
#      population-robust) and random-intercept SD ratio (secondary)
#   I  crossed 2 sets x 2 groups: rasch_efrm alpha and phi vs the glmer
#      per-cell anchors (per-set coefficient slopes for phi; 2x2
#      log-linear decomposition of random-intercept SDs for alpha)
#   J  polytomous rasch_efrm item-set unit ratio vs TAM GPCM with
#      set-constrained slopes (irtmodel="GPCM.design", set-indicator E;
#      "2PL.groups" does NOT impose proportional category loadings on
#      polytomous data and is unsound there)
#
# This is an AGREEMENT study: cross-package differences are computed per
# dataset, so simulator truths vary across replicates while the target
# unit ratios (cells C and F) are fixed by design. Serial, single process.
#   Rscript tools/simval/studies/cross-package-validation.R

suppressWarnings(suppressMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(VGAM)
}))
source("tools/simval/harness.R")
STUDY <- "cross-package-validation"
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- sv_row(STUDY, ...)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M")),
                          sprintf(...), "\n")

## ---- A: dichotomous locations vs sirt and eRm ----------------------------
R <- 25L
a_sirt <- a_erm <- matrix(NA_real_, R, 2)  # cor, maxdiff
rmse <- matrix(NA_real_, R, 3, dimnames = list(NULL, c("rasch","sirt","eRm")))
for (r in seq_len(R)) {
  d <- simulate_rasch(500, 10, model = "dichotomous", seed = 500 + r)
  X <- as.matrix(d[, grep("^I\\d", names(d))])
  tr <- attr(d, "truth")$difficulty; tr <- tr - mean(tr)
  b_r <- local({f <- rasch(X); v <- setNames(f$items$location, f$items$item); v - mean(v)})
  b_s <- local({f <- sirt::rasch.pairwise(X, zerosum = TRUE, progress = FALSE)
                v <- setNames(f$item$b, f$item$item)[names(b_r)]; v - mean(v)})
  b_e <- local({f <- eRm::RM(X); v <- -f$betapar
                names(v) <- sub("^beta ", "", names(v)); v <- v[names(b_r)]; v - mean(v)})
  a_sirt[r, ] <- c(cor(b_r, b_s), max(abs(b_r - b_s)))
  a_erm[r, ]  <- c(cor(b_r, b_e), max(abs(b_r - b_e)))
  rmse[r, ] <- c(sqrt(mean((b_r - tr)^2)), sqrt(mean((b_s - tr)^2)),
                 sqrt(mean((b_e - tr)^2)))
}
add("dichotomous locations vs sirt::rasch.pairwise", "agreement over items",
    R, notes = sprintf(paste0("mean cor %.6f; mean max|diff| %.4f, worst %.4f ",
      "logits; truth RMSE rasch %.3f vs sirt %.3f (sirt is the MINCHI ",
      "pairwise variant, not Zwinderman PCML)"),
      mean(a_sirt[,1]), mean(a_sirt[,2]), max(a_sirt[,2]),
      mean(rmse[,"rasch"]), mean(rmse[,"sirt"])))
add("dichotomous locations vs eRm::RM (full CML)", "agreement over items",
    R, notes = sprintf("mean cor %.6f; mean max|diff| %.4f, worst %.4f logits; truth RMSE rasch %.3f vs eRm %.3f",
      mean(a_erm[,1]), mean(a_erm[,2]), max(a_erm[,2]),
      mean(rmse[,"rasch"]), mean(rmse[,"eRm"])))
tick("A done: sirt cor %.6f, eRm cor %.6f", mean(a_sirt[,1]), mean(a_erm[,1]))

## ---- B: PCM vs eRm::PCM --------------------------------------------------
loc_m <- thr_m <- matrix(NA_real_, R, 2); n_ref <- 0L
for (r in seq_len(R)) {
  d <- simulate_rasch(400, 8, model = "PCM", n_categories = 4,
                      difficulty = c(-1.5, 1.5), seed = 700 + r)
  X <- as.matrix(d[, grep("^I\\d", names(d))])
  if (!all(apply(X, 2, function(x) length(unique(x))) == 4)) { n_ref <- n_ref + 1L; next }
  f_r <- rasch(X, model = "PCM")
  loc_r <- setNames(f_r$items$location, f_r$items$item)
  tau_r <- do.call(rbind, f_r$tau_list)
  f_e <- eRm::PCM(X)
  tt <- eRm::thresholds(f_e)$threshtable[[1]]
  sh <- mean(tt[, "Location"])
  loc_m[r, ] <- c(cor(loc_r, tt[, "Location"] - sh),
                  max(abs(loc_r - (tt[, "Location"] - sh))))
  thr_m[r, ] <- c(cor(as.vector(tau_r), as.vector(tt[, -1] - sh)),
                  max(abs(tau_r - (tt[, -1] - sh))))
}
ok <- stats::complete.cases(loc_m)
add("PCM locations vs eRm::PCM (full CML)", "agreement over items", sum(ok),
    n_attempted = R, n_refused = n_ref,
    notes = sprintf("mean cor %.6f; mean max|diff| %.4f, worst %.4f logits",
      mean(loc_m[ok,1]), mean(loc_m[ok,2]), max(loc_m[ok,2])))
add("PCM thresholds vs eRm::PCM", "agreement over 24 thresholds", sum(ok),
    n_attempted = R, n_refused = n_ref,
    notes = sprintf("mean cor %.6f; mean max|diff| %.4f, worst %.4f logits (pairwise vs full CML, largest at sparse extreme thresholds)",
      mean(thr_m[ok,1]), mean(thr_m[ok,2]), max(thr_m[ok,2])))
tick("B done: loc cor %.6f", mean(loc_m[ok,1]))

## ---- C: EFRM item-set unit ratio vs TAM ----------------------------------
lr <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("efrm","tam"))); n_refC <- 0L
for (r in seq_len(R)) {
  s <- simulate_efrm(n_per_group = 1000, items_per_set = 8, n_sets = 2,
                     n_groups = 1, set_unit_ratio = 1.4, seed = 900 + r)
  tr <- attr(s, "truth"); its <- unlist(tr$item_sets, use.names = FALSE)
  fe <- tryCatch(rasch_efrm(s, item_sets = tr$item_sets, groups = "group",
                            id = "id", boot_reps = 50), error = function(e) NULL)
  if (is.null(fe)) { n_refC <- n_refC + 1L; next }
  al <- fe$alpha_table$alpha[match(c("set1","set2"), fe$alpha_table$set)]
  mt <- TAM::tam.mml.2pl(resp = as.matrix(s[, its]), irtmodel = "2PL.groups",
                         est.slopegroups = rep(1:2, each = 8),
                         control = list(progress = FALSE))
  a <- mt$B[its, 2, 1]
  lr[r, ] <- c(log(al[2]/al[1]), log(mean(a[9:16])/mean(a[1:8])))
}
okC <- stats::complete.cases(lr); lt <- log(1.4)
add("EFRM item-set log unit ratio (truth log 1.4 = 0.3365)",
    "bias of rasch_efrm", sum(okC), bias = mean(lr[okC,"efrm"]) - lt,
    emp_sd = sd(lr[okC,"efrm"]), n_attempted = R, n_refused = n_refC,
    notes = paste("semiparametric finite-grid link with the within-set",
                  "calibrations held fixed; one latent distribution is",
                  "estimated for each observed person group"))
add("EFRM item-set log unit ratio (truth log 1.4 = 0.3365)",
    "bias of TAM 2PL.groups anchor", sum(okC), bias = mean(lr[okC,"tam"]) - lt,
    emp_sd = sd(lr[okC,"tam"]),
    notes = sprintf("external MML anchor (normal population satisfied by the simulator); cor(efrm, TAM) across replicates %.3f, mean |log-ratio diff| %.4f",
      cor(lr[okC,"efrm"], lr[okC,"tam"]), mean(abs(lr[okC,1] - lr[okC,2]))))
tick("C done: efrm bias %+.4f, TAM bias %+.4f", mean(lr[okC,1]) - lt, mean(lr[okC,2]) - lt)

## ---- D: btl dichotomous vs BradleyTerry2::BTm ----------------------------
mx <- se_lo <- se_hi <- rep(NA_real_, R)
for (r in seq_len(R)) {
  d <- simulate_btl(n_objects = 8, n_judges = 12, reps_per_pair = 12,
                    model = "dichotomous", seed = 1100 + r)
  f <- btl(d, "object_a", "object_b", winner = "winner")   # unclustered for SE parity
  b1 <- setNames(f$objects$location, f$objects$object); b1 <- b1 - mean(b1)
  objs <- sort(names(b1))
  key <- paste(pmin(d$object_a, d$object_b), pmax(d$object_a, d$object_b), sep = ":")
  first <- pmin(d$object_a, d$object_b); second <- pmax(d$object_a, d$object_b)
  agg <- data.frame(
    player1 = factor(tapply(first, key, `[`, 1), levels = objs),
    player2 = factor(tapply(second, key, `[`, 1), levels = objs),
    win1 = as.vector(tapply(d$winner == first, key, sum)),
    win2 = as.vector(tapply(d$winner == second, key, sum)))
  bt <- BradleyTerry2::BTm(outcome = cbind(win1, win2), player1 = player1,
                           player2 = player2, data = agg)
  ab <- BradleyTerry2::BTabilities(bt)
  b2 <- ab[, "ability"]; b2 <- b2 - mean(b2)
  mx[r] <- max(abs(b1 - b2[names(b1)]))
  # BTm SEs are reference-scale; project its covariance to sum-zero
  # contrasts before comparing with btl's sum-zero SEs
  K <- length(objs)
  V <- matrix(0, K, K, dimnames = list(objs, objs))
  vc <- vcov(bt); nm <- sub("^\\.\\.", "", rownames(vc))
  V[nm, nm] <- vc
  Cm <- diag(K) - 1/K
  se_btm <- setNames(sqrt(diag(Cm %*% V %*% t(Cm))), objs)
  se_btl <- setNames(f$objects$se, f$objects$object)
  sr <- se_btl[objs] / se_btm[objs]
  se_lo[r] <- min(sr); se_hi[r] <- max(sr)
}
add("btl dichotomous vs BradleyTerry2::BTm", "identical likelihood parity", R,
    notes = sprintf("worst max|diff| over %d replicates: %.2e logits (btl Newton tol 1e-8); unclustered-sandwich vs model-based sum-zero-projected SE ratios in [%.3f, %.3f]",
      R, max(mx), min(se_lo), max(se_hi)))
tick("D done: worst %.2e", max(mx))

## ---- E: btl polytomous vs constrained VGAM acat --------------------------
make_design <- function(a, b, objs) {
  K <- length(objs); ia <- match(a, objs); ib <- match(b, objs)
  Z <- matrix(0, length(a), K - 1L,
              dimnames = list(NULL, paste0("z", seq_len(K - 1L))))
  for (k in seq_len(K - 1L)) Z[, k] <- (ia == k) - (ib == k)
  Z[ia == K, ] <- Z[ia == K, ] - 1; Z[ib == K, ] <- Z[ib == K, ] + 1
  Z
}
beta_from_vglm <- function(fit, objs) {
  B <- coef(fit, matrix = TRUE)
  b <- B[setdiff(rownames(B), "(Intercept)"), 1]
  b <- c(b, -sum(b)); names(b) <- objs; b - mean(b)
}
mxF <- mxP <- llF <- rep(NA_real_, R); m5 <- 4L
v1 <- (seq_len(m5) - (m5 + 1)/2); v1 <- v1 / sqrt(sum(v1^2))
for (r in seq_len(R)) {
  d <- simulate_btl(n_objects = 8, n_judges = 10, reps_per_pair = 25,
                    model = "polytomous", n_categories = 5, seed = 1300 + r)
  objs <- sort(unique(c(d$object_a, d$object_b)))
  Z <- make_design(d$object_a, d$object_b, objs)
  dat <- data.frame(y = factor(d$response, levels = 0:m5, ordered = TRUE), Z)
  zc <- setNames(rep(list(matrix(1, m5, 1)), length(objs) - 1L), colnames(Z))
  bf <- btl(d, "object_a", "object_b", response = "response", thresholds = "free")
  bp <- btl(d, "object_a", "object_b", response = "response", thresholds = "pc")
  vf <- VGAM::vglm(y ~ ., family = VGAM::acat(reverse = FALSE), data = dat,
             constraints = c(list("(Intercept)" =
               cbind(c(-1,0,0,1), c(0,-1,1,0))), zc))
  vp <- VGAM::vglm(y ~ ., family = VGAM::acat(reverse = FALSE), data = dat,
             constraints = c(list("(Intercept)" = cbind(-v1)), zc))
  bfl <- setNames(bf$objects$location, bf$objects$object)[objs]
  bpl <- setNames(bp$objects$location, bp$objects$object)[objs]
  mxF[r] <- max(abs(bfl - mean(bfl) - beta_from_vglm(vf, objs)),
                abs(bf$thresholds$tau + coef(vf, matrix = TRUE)["(Intercept)", ]))
  mxP[r] <- max(abs(bpl - mean(bpl) - beta_from_vglm(vp, objs)),
                abs(bp$thresholds$tau + coef(vp, matrix = TRUE)["(Intercept)", ]))
  llF[r] <- abs(bf$loglik - VGAM::logLik(vf))
}
add("btl polytomous free thresholds vs VGAM acat (symmetry-constrained)",
    "identical model parity (locations + thresholds)", R,
    notes = sprintf("worst max|diff| %.2e logits; worst |logLik diff| %.2e (5 categories, 8 objects)",
      max(mxF), max(llF)))
add("btl polytomous pc thresholds vs VGAM acat (spread-constrained)",
    "identical model parity", R,
    notes = sprintf("worst max|diff| %.2e logits", max(mxP)))
tick("E done: free %.2e, pc %.2e", max(mxF), max(mxP))

## ---- F: btl_efrm panel-unit ratio vs per-panel vglm ----------------------
rat <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("efrm","vglm"))); n_refF <- 0L
for (r in seq_len(R)) {
  db <- simulate_btl_efrm(n_objects_per_set = 8, n_sets = 1, n_panels = 2,
                          n_judges_per_panel = 8, reps_within = 60,
                          panel_units = c(1, 1.3), seed = 1500 + r)
  trb <- attr(db, "truth")
  ef <- tryCatch(btl_efrm(db, "object_a", "object_b", winner = "winner",
                          judge = "judge", panels = "panel",
                          object_sets = trb$object_sets,
                          se_method = "conditional"), error = function(e) NULL)
  if (is.null(ef)) { n_refF <- n_refF + 1L; next }
  phi <- setNames(ef$phi_table$phi, ef$phi_table$panel)
  objsE <- sort(unique(c(db$object_a, db$object_b)))
  fit_bin <- function(dp) {
    Zp <- make_design(dp$object_a, dp$object_b, objsE)
    yb <- as.integer(dp$winner == dp$object_a)
    fb <- VGAM::vglm(yb ~ . - 1, family = VGAM::binomialff,
                     data = data.frame(yb = yb, Zp))
    b <- coef(fb); b <- c(b, -sum(b)); names(b) <- objsE; b - mean(b)
  }
  e1 <- fit_bin(db[db$panel == "panel1", ])
  e2 <- fit_bin(db[db$panel == "panel2", ])
  rat[r, ] <- c(phi["panel2"]/phi["panel1"], sum(e1*e2)/sum(e1^2))
}
okF <- stats::complete.cases(rat)
add("btl_efrm panel-unit ratio (truth 1.30)", "bias of btl_efrm phi ratio",
    sum(okF), bias = mean(rat[okF,"efrm"]) - 1.3, emp_sd = sd(rat[okF,"efrm"]),
    n_attempted = R, n_refused = n_refF,
    notes = sprintf("per-panel intercept-free vglm LS-slope anchor: bias %+.4f, sd %.4f; mean |ratio diff between methods| %.4f",
      mean(rat[okF,"vglm"]) - 1.3, sd(rat[okF,"vglm"]),
      mean(abs(rat[okF,1] - rat[okF,2]))))
tick("F done: efrm %+.4f, vglm %+.4f", mean(rat[okF,1]) - 1.3, mean(rat[okF,2]) - 1.3)

## ---- G: rasch_mfrm vs TAM::tam.mml.mfr -----------------------------------
# TAM's default constraint="cases" leaves the item facet uncentred and the
# rater facet centred; centre every facet in both packages before
# comparing. formulaA must use item:step (PCM thresholds); ~ step would
# impose a rating-scale structure. Facet levels matched by name.
it_m <- sv_m <- th_m <- matrix(NA_real_, R, 2)
rms <- matrix(NA_real_, R, 4,
              dimnames = list(NULL, c("it_r","it_t","sv_r","sv_t")))
ctr <- function(x) x - mean(x)
for (r in seq_len(R)) {
  d <- simulate_mfrm(n_persons = 200, n_items = 5, n_raters = 5,
                     n_categories = 4, theta_sd = 1.2, item_sd = 1,
                     rater_severity_sd = 0.8, seed = 1700 + r)
  truth <- attr(d, "truth")
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater")
  r_item <- ctr(setNames(mf$item_effects$location, mf$item_effects$item))
  r_sev  <- ctr(setNames(mf$facet_effects$rater$severity,
                         mf$facet_effects$rater$level))
  r_thr <- mf$item_thresholds[order(mf$item_thresholds$item,
                                    mf$item_thresholds$k), ]
  wide <- reshape(as.data.frame(d), idvar = c("person", "rater"),
                  timevar = "item", direction = "wide")
  iids <- sort(unique(d$item))
  resp <- wide[, paste0("score.", iids)]; colnames(resp) <- iids
  tam <- TAM::tam.mml.mfr(resp = resp,
                          facets = data.frame(rater = wide$rater),
                          formulaA = ~ item + rater + item:step,
                          pid = as.integer(factor(wide$person)),
                          verbose = FALSE)
  xf <- tam$xsi.facets
  t_im <- setNames(xf$xsi[xf$facet == "item"], xf$parameter[xf$facet == "item"])
  t_sv <- setNames(xf$xsi[xf$facet == "rater"],
                   sub("^rater", "", xf$parameter[xf$facet == "rater"]))
  st <- xf[xf$facet == "item:step", ]
  st_item <- sub(":.*$", "", st$parameter)
  t_thr <- t_im[st_item] + st$xsi
  o <- order(st_item, as.integer(sub("^.*step", "", st$parameter)))
  t_thr <- t_thr[o]
  t_im <- ctr(t_im[names(r_item)]); t_sv <- ctr(t_sv[names(r_sev)])
  it_m[r, ] <- c(cor(r_item, t_im), max(abs(r_item - t_im)))
  sv_m[r, ] <- c(cor(r_sev, t_sv), max(abs(r_sev - t_sv)))
  th_m[r, ] <- c(cor(ctr(r_thr$tau), ctr(t_thr)),
                 max(abs(ctr(r_thr$tau) - ctr(t_thr))))
  tr_i <- ctr(truth$difficulty[names(r_item)])
  tr_s <- ctr(truth$severity[names(r_sev)])
  rms[r, ] <- c(sqrt(mean((r_item - tr_i)^2)), sqrt(mean((t_im - tr_i)^2)),
                sqrt(mean((r_sev - tr_s)^2)),  sqrt(mean((t_sv - tr_s)^2)))
}
add("MFRM item locations vs TAM::tam.mml.mfr", "agreement over items", R,
    notes = sprintf("mean cor %.6f; mean max|diff| %.4f, worst %.4f logits; truth RMSE rasch %.3f vs TAM %.3f",
      mean(it_m[,1]), mean(it_m[,2]), max(it_m[,2]),
      mean(rms[,"it_r"]), mean(rms[,"it_t"])))
add("MFRM rater severities vs TAM::tam.mml.mfr", "agreement over raters", R,
    notes = sprintf("mean cor %.6f; mean max|diff| %.4f, worst %.4f logits; truth RMSE rasch %.3f vs TAM %.3f",
      mean(sv_m[,1]), mean(sv_m[,2]), max(sv_m[,2]),
      mean(rms[,"sv_r"]), mean(rms[,"sv_t"])))
add("MFRM PCM thresholds vs TAM::tam.mml.mfr", "agreement over 15 thresholds", R,
    notes = sprintf("mean cor %.6f; mean max|diff| %.4f, worst %.4f logits (pairwise CML vs MML)",
      mean(th_m[,1]), mean(th_m[,2]), max(th_m[,2])))
tick("G done: items cor %.6f, raters cor %.6f", mean(it_m[,1]), mean(sv_m[,1]))

## ---- H/I: EFRM person-side and crossed units vs lme4::glmer ---------------
cell_long <- function(d, item_cols, group_val) {
  di <- d[as.character(d$group) == group_val, , drop = FALSE]
  data.frame(person = factor(rep(di$id, times = length(item_cols))),
             item = factor(rep(item_cols, each = nrow(di)), levels = item_cols),
             resp = as.vector(as.matrix(di[, item_cols])))
}
fit_cell_glmer <- function(dat) {
  m <- suppressWarnings(lme4::glmer(resp ~ 0 + item + (1 | person),
        data = dat, family = binomial, nAGQ = 1,
        control = lme4::glmerControl(optimizer = "bobyqa")))
  list(coef = lme4::fixef(m),
       sd = sqrt(unname(lme4::VarCorr(m)$person[1, 1])))
}
ls_slope <- function(b1, b2) {
  c1 <- b1 - mean(b1); c2 <- b2 - mean(b2); sum(c1 * c2) / sum(c1 * c1)
}

## H: person side (1 set, 2 groups, phi ratio 1.3). Third anchor: TAM on
## item-by-person-group VIRTUAL COLUMNS (group g's persons respond to
## their group's block, NA elsewhere) with est.slopegroups indexed by the
## person group -- the blocks link only through the fixed N(0,1)
## population per group, so this anchor is normality-dependent but
## properly integrated (quadrature), unlike glmer's Laplace SD route.
ph <- matrix(NA_real_, R, 4,
             dimnames = list(NULL, c("efrm","slope","sdr","tamv")))
for (r in seq_len(R)) {
  d <- simulate_efrm(n_per_group = 400, items_per_set = 8, n_sets = 1,
                     n_groups = 2, set_unit_ratio = 1, group_unit_ratio = 1.3,
                     seed = 1900 + r)
  tr <- attr(d, "truth"); its <- unlist(tr$item_sets)
  f1 <- fit_cell_glmer(cell_long(d, its, "g1"))
  f2 <- fit_cell_glmer(cell_long(d, its, "g2"))
  X <- as.matrix(d[, its]); grp <- as.character(d$group)
  V <- matrix(NA_integer_, nrow(X), 16,
              dimnames = list(NULL, paste0(rep(its, 2), "_",
                                           rep(c("g1","g2"), each = 8))))
  V[grp == "g1", 1:8]  <- X[grp == "g1", ]
  V[grp == "g2", 9:16] <- X[grp == "g2", ]
  mt <- tryCatch(TAM::tam.mml.2pl(resp = V, irtmodel = "2PL.groups",
                                  est.slopegroups = rep(1:2, each = 8),
                                  control = list(progress = FALSE)),
                 error = function(e) NULL)
  ef <- tryCatch(rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                            id = "id", boot_reps = 30), error = function(e) NULL)
  if (is.null(ef) || is.null(mt)) next
  a <- mt$B[, 2, 1]
  ph[r, ] <- c(ef$phi_table$phi[2] / ef$phi_table$phi[1],
               ls_slope(f1$coef, f2$coef), f2$sd / f1$sd,
               mean(a[9:16]) / mean(a[1:8]))
}
okH <- stats::complete.cases(ph)
add("EFRM person-group unit phi (truth ratio 1.30) vs glmer and TAM anchors",
    "log phi ratio bias", sum(okH),
    bias = mean(log(ph[okH,"efrm"])) - log(1.3),
    emp_sd = sd(log(ph[okH,"efrm"])),
    notes = sprintf("glmer coefficient-slope anchor bias %+.4f (population-robust); glmer ranef-SD-ratio %+.4f (normality-dependent, upward at short tests); TAM virtual-column slope-group anchor %+.4f (normality-dependent, quadrature)",
      mean(log(ph[okH,"slope"])) - log(1.3), mean(log(ph[okH,"sdr"])) - log(1.3),
      mean(log(ph[okH,"tamv"])) - log(1.3)))
tick("H done: efrm %+.4f, slope anchor %+.4f, TAM virtual %+.4f",
     mean(log(ph[okH,"efrm"])) - log(1.3), mean(log(ph[okH,"slope"])) - log(1.3),
     mean(log(ph[okH,"tamv"])) - log(1.3))

## I: crossed 2 sets x 2 groups (alpha 1.4, phi 1.3)
cr <- matrix(NA_real_, R, 4,
             dimnames = list(NULL, c("e_phi","e_al","g_phi","g_al")))
for (r in seq_len(R)) {
  d <- simulate_efrm(n_per_group = 400, items_per_set = 8, n_sets = 2,
                     n_groups = 2, set_unit_ratio = 1.4, group_unit_ratio = 1.3,
                     seed = 2100 + r)
  tr <- attr(d, "truth")
  cells <- list()
  for (s in 1:2) for (g in c("g1", "g2"))
    cells[[paste0("s", s, g)]] <- fit_cell_glmer(
      cell_long(d, tr$item_sets[[s]], g))
  ef <- tryCatch(rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                            id = "id", boot_reps = 30), error = function(e) NULL)
  if (is.null(ef)) next
  L <- log(matrix(c(cells$s1g1$sd, cells$s1g2$sd,
                    cells$s2g1$sd, cells$s2g2$sd), 2, 2, byrow = TRUE))
  cr[r, ] <- c(ef$phi_table$phi[2] / ef$phi_table$phi[1],
               ef$alpha_table$alpha[2] / ef$alpha_table$alpha[1],
               sqrt(ls_slope(cells$s1g1$coef, cells$s1g2$coef) *
                    ls_slope(cells$s2g1$coef, cells$s2g2$coef)),
               exp(mean(L[2, ] - L[1, ])))
}
okI <- stats::complete.cases(cr)
add("crossed EFRM 2x2 (alpha 1.40, phi 1.30) vs glmer per-cell anchors",
    "log unit-ratio biases", sum(okI),
    bias = mean(log(cr[okI,"e_al"])) - log(1.4),
    emp_sd = sd(log(cr[okI,"e_al"])),
    notes = sprintf("efrm phi bias %+.4f (anchor %+.4f via per-set coefficient slopes); efrm alpha bias in `bias` (anchor %+.4f via ranef-SD 2x2 decomposition)",
      mean(log(cr[okI,"e_phi"])) - log(1.3),
      mean(log(cr[okI,"g_phi"])) - log(1.3),
      mean(log(cr[okI,"g_al"])) - log(1.4)))
tick("I done: efrm alpha %+.4f phi %+.4f",
     mean(log(cr[okI,"e_al"])) - log(1.4), mean(log(cr[okI,"e_phi"])) - log(1.3))

## J: polytomous EFRM vs TAM GPCM.design
pj <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("efrm","tam")))
Emat <- cbind(as.integer(rep(c(1, 0), each = 8)),
              as.integer(rep(c(0, 1), each = 8)))
for (r in seq_len(R)) {
  d <- simulate_efrm(n_per_group = 800, items_per_set = 8, n_sets = 2,
                     n_groups = 1, set_unit_ratio = 1.4, n_categories = 4,
                     seed = 2300 + r)
  tr <- attr(d, "truth"); its <- unlist(tr$item_sets)
  ef <- tryCatch(rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                            id = "id", boot_reps = 30), error = function(e) NULL)
  mt <- tryCatch(TAM::tam.mml.2pl(resp = as.matrix(d[, its]),
                                  irtmodel = "GPCM.design", E = Emat,
                                  control = list(progress = FALSE)),
                 error = function(e) NULL)
  if (is.null(ef) || is.null(mt)) next
  al <- ef$alpha_table$alpha[match(c("set1","set2"), ef$alpha_table$set)]
  a <- mt$B[its, 2, 1]
  pj[r, ] <- c(log(al[2]/al[1]), log(mean(a[9:16])/mean(a[1:8])))
}
okJ <- stats::complete.cases(pj)
add("polytomous EFRM unit ratio (truth log 1.4) vs TAM GPCM.design",
    "log unit ratio bias", sum(okJ),
    bias = mean(pj[okJ,"efrm"]) - log(1.4), emp_sd = sd(pj[okJ,"efrm"]),
    notes = sprintf("TAM GPCM set-constrained-slope anchor bias %+.4f (sd %.4f); proportional category loadings verified; irtmodel='2PL.groups' is unsound for polytomous data (free category loadings, ratio 1.19 vs truth 1.40 single-seed)",
      mean(pj[okJ,"tam"]) - log(1.4), sd(pj[okJ,"tam"])))
tick("J done: efrm %+.4f, TAM GPCM %+.4f",
     mean(pj[okJ,"efrm"]) - log(1.4), mean(pj[okJ,"tam"]) - log(1.4))

sv_write(do.call(rbind, rows), "cross-package-validation")
cat(sprintf("TOTAL elapsed: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
