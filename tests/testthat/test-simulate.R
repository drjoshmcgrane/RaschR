# Data simulation: each planted departure must trip its matching diagnostic
# (the plant -> detect loop), and the truth is attached for recovery checks.

test_that("simulate_rasch plants misfit the Rasch diagnostics detect", {
  # discrimination: a central over-discriminating item overfits (low outfit)
  d <- simulate_rasch(600, 11, discrimination = c(rep(1, 5), 3, rep(1, 5)),
                      seed = 1)
  f <- rasch(d)
  # the over-discriminating item overfits: its outfit is the lowest and
  # clearly below expectation (extreme persons are excluded from item fit,
  # so the mean-square is not deflated by their boundary residuals)
  expect_lt(f$items$outfit_ms[6], 0.75)
  expect_equal(which.min(f$items$outfit_ms), 6L)
  expect_s3_class(d, "rasch_sim")

  # DIF flags the planted item and (essentially) nothing else
  d <- simulate_rasch(800, 10, dif = list(items = "I05", uniform = 1.2),
                      n_groups = 2, seed = 2)
  s <- dif_anova(rasch(d, factors = "group"))$summary
  expect_true(s$uniform_DIF[s$item == "I05"])
  expect_equal(sum(s$uniform_DIF[s$item != "I05"]), 0L)

  # local dependence flags the planted pair on Q3*
  d <- simulate_rasch(1000, 10,
                      dependence = list(pairs = list(c("I03", "I04")),
                                        strength = 2.5), seed = 4)
  h <- residual_correlations(rasch(d), flag = 0.2)$flagged
  expect_true(any((h$item_a == "I03" & h$item_b == "I04") |
                  (h$item_a == "I04" & h$item_b == "I03")))

  # careless responders inflate person outfit
  d <- simulate_rasch(600, 15, careless = 0.12, seed = 6)
  po <- rasch(d)$person$outfit_ms; ci <- attr(d, "truth")$careless_idx
  expect_gt(mean(po[ci], na.rm = TRUE), mean(po[-ci], na.rm = TRUE) + 0.5)

  expect_output(print(d), "careless")
})

test_that("simulate_btl plants misfit the paired-comparison diagnostics detect", {
  # erratic judges carry large fit residuals and low consistency
  d <- simulate_btl(8, 12, erratic_judges = 0.17, seed = 1)
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  er <- attr(d, "truth")$erratic
  expect_gt(mean(bt$judges$fit_resid[bt$judges$judge %in% er]),
            mean(bt$judges$fit_resid[!bt$judges$judge %in% er]) + 1)

  # graded comparisons recover the object locations
  d <- simulate_btl(8, 12, model = "graded", n_categories = 4, seed = 4)
  bt <- btl(d, "object_a", "object_b", response = "response", judge = "judge")
  expect_gt(cor(bt$objects$location,
                attr(d, "truth")$location[bt$objects$object]), 0.95)

  # a planted carry-over dependence is recovered
  d <- simulate_btl(6, 10, reps_per_pair = 40,
                    dependence = list(exposure = 0, carry_over = 1.2), seed = 3)
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
            order = "order")
  co <- bt$dependence$estimate[bt$dependence$effect == "carry_over"]
  expect_gt(co, 0.5)
})

test_that("simulate_btl accepts explanatory object locations", {
  loc <- c(O1 = -1.2, O2 = -0.1, O3 = 0.2, O4 = 1.8)
  d <- simulate_btl(4, 8, 20, object_locations = loc, seed = 81)
  expect_equal(attr(d, "truth")$location, loc - mean(loc))
  expect_error(simulate_btl(4, object_locations = c(1, 2, 3)),
               "length 4")
  expect_error(simulate_btl(4, object_locations = c(A = 1, B = 2, C = 3, D = 4)),
               "every generated object")
  expect_error(simulate_btl(4, object_locations = rep(1, 4),
                            second_attribute = list(rho = 0.2)),
               "positive spread")
})

test_that("simulate_mfrm plants rater severity, misfit, and interaction", {
  d <- simulate_mfrm(120, 5, 6, rater_severity_sd = 0.7, erratic_raters = 0.17,
                     interaction = list(rater = "R3", item = "I2", bias = 1.8),
                     seed = 1)
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater", interaction = "rater")
  tr <- attr(d, "truth"); fe <- mf$facet_effects$rater
  rec <- fe$severity[match(names(tr$severity), fe$level)]
  # severities are recovered for the well-behaved raters (the erratic one's
  # true severity is meaningless once it rates at random)
  keep <- !(names(tr$severity) %in% tr$erratic)
  expect_gt(abs(cor(rec[keep], tr$severity[keep])), 0.9)
  er <- tr$erratic
  expect_gt(mean(fe$fit_resid[fe$level %in% er]),
            mean(fe$fit_resid[!fe$level %in% er]) + 1)      # erratic rater misfits
  ie <- mf$interaction_effects                             # interaction at R3xI2
  expect_equal(ie[which.max(abs(ie$gamma)), c("item", "level")],
               data.frame(item = "I2", level = "R3"), ignore_attr = TRUE)
})

test_that("simulate_efrm plants a frame-unit ratio rasch_efrm recovers", {
  d <- simulate_efrm(300, 8, set_unit_ratio = 1.35, seed = 2)
  tr <- attr(d, "truth")
  ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group")
  ratio <- max(ef$alpha_table$alpha) / min(ef$alpha_table$alpha)
  expect_gt(ratio, 1.2); expect_lt(ratio, 1.55)          # ~1.35 recovered
  expect_output(print(d), "set-unit ratio")
})

test_that("simulate_efrm generates partial credit items on request", {
  d <- simulate_efrm(250, 6, n_sets = 2, n_groups = 1, set_unit_ratio = 1.3,
                     n_categories = 4, seed = 5)
  tr <- attr(d, "truth")
  X <- as.matrix(d[, unlist(tr$item_sets)])
  expect_setequal(sort(unique(as.vector(X))), 0:3)
  expect_length(tr$thresholds, ncol(X))
  expect_true(all(vapply(tr$thresholds, length, 1L) == 3L))
  # thresholds centre on the item locations
  expect_equal(unname(vapply(tr$thresholds, mean, 0)),
               unname(tr$difficulty), tolerance = 1e-10)
  # the dichotomous draw stream is untouched by the generalisation
  d2 <- simulate_efrm(50, 4, n_sets = 2, n_groups = 2,
                      set_unit_ratio = 1.3, seed = 1)
  expect_identical(sum(as.matrix(d2[, 2:9])), 402L)
})

test_that("EFRM simulators plant frame-specific and judge departures", {
  clean <- simulate_efrm(600, 5, n_sets = 2, n_groups = 2,
                         set_unit_ratio = 1, seed = 710)
  drift <- simulate_efrm(600, 5, n_sets = 2, n_groups = 2,
    set_unit_ratio = 1,
    item_drift = list(items = "S1I03", group = "g2", shift = 1.2),
    careless = 0.05, missing = 0.03, seed = 710)
  tr <- attr(drift, "truth")
  expect_identical(tr$item_drift$group, "g2")
  expect_identical(tr$item_drift$items, "S1I03")
  expect_length(tr$careless_idx, 60L)
  expect_equal(length(tr$missing_cells), round(0.03 * 1200 * 10))
  # The same seed isolates the shifted item before the later contamination.
  g2 <- clean$group == "g2"
  expect_lt(mean(drift$S1I03[g2], na.rm = TRUE),
            mean(clean$S1I03[g2], na.rm = TRUE) - 0.08)
  expect_true(any(grepl("item drift", tr$planted)))

  b <- simulate_btl_efrm(5, 2, n_judges_per_panel = 5, n_panels = 2,
                         reps_within = 4, reps_cross = 4,
                         erratic_judges = 0.2, seed = 711)
  expect_length(attr(b, "truth")$erratic, 2L)
  expect_true(any(grepl("erratic judge", attr(b, "truth")$planted)))
  expect_setequal(unique(b$judge), sprintf("J%03d", 1:10))
  expect_setequal(unique(b$panel), c("panel1", "panel2"))
  expect_true(all(attr(b, "truth")$erratic %in% b$judge))
  expect_equal(length(unique(b$panel[b$judge %in%
                                    attr(b, "truth")$erratic])), 2L)
  expect_lte(diff(range(table(b$judge))), 1L)
  set_of <- attr(b, "truth")$set_of
  stratum <- paste(pmin(set_of[b$object_a], set_of[b$object_b]),
                   pmax(set_of[b$object_a], set_of[b$object_b]), sep = ":")
  for (ss in unique(stratum))
    expect_lte(diff(range(table(factor(b$panel[stratum == ss],
                                       levels = c("panel1", "panel2"))))), 1L)
})

test_that("paired-comparison simulators retain every declared judge", {
  d <- simulate_btl(n_objects = 3, n_judges = 20, reps_per_pair = 7,
                    erratic_judges = 0.2, seed = 713)
  expect_setequal(unique(d$judge), sprintf("J%d", 1:20))
  expect_lte(diff(range(table(d$judge))), 1L)
  expect_true(all(attr(d, "truth")$erratic %in% d$judge))
  expect_error(simulate_btl(n_objects = 3, n_judges = 20,
                            reps_per_pair = 1, seed = 714),
               "fewer than the 20 requested judges")
})

test_that("the extra misfit types plant detectable signals", {
  # extreme response style: style persons over-use the end categories
  d <- simulate_rasch(600, 12, model = "PCM", n_categories = 4,
                      response_style = list(type = "extreme", prop = 0.3), seed = 1)
  si <- attr(d, "truth")$style_idx; cats <- as.matrix(d[, grep("^I", names(d))])
  expect_gt(mean(cats[si, ] %in% c(0, 3)), mean(cats[-si, ] %in% c(0, 3)) + 0.1)

  # speededness: a missing tail growing toward the last item
  d <- simulate_rasch(800, 15, speeded = 0.5, seed = 2)
  miss <- colMeans(is.na(as.matrix(d[, grep("^I", names(d))])))
  expect_lt(miss[8], 0.02)
  expect_gt(miss[15], 0.3)
  expect_true(miss[15] > miss[13] && miss[13] > miss[11])   # monotone gradient

  # MFRM halo: halo raters barely differentiate items -> large interaction
  d <- simulate_mfrm(140, 6, 6, rater_severity_sd = 0.5, item_sd = 1.3,
                     halo = 0.17, seed = 5)
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater", interaction = "rater")
  hr <- attr(d, "truth")$halo; ie <- mf$interaction_effects
  expect_gt(mean(abs(ie$gamma[ie$level %in% hr])),
            2 * mean(abs(ie$gamma[!ie$level %in% hr])))
})

test_that("sim_replicate and sim_recovery support Monte Carlo and recovery", {
  b <- sim_replicate(simulate_rasch, 6, n_persons = 300, n_items = 8, seed = 1)
  expect_s3_class(b, "rasch_sim_batch")
  expect_length(b, 6)
  expect_false(identical(b[[1]]$I01, b[[2]]$I01))          # different datasets

  # recovery: a clean fit gets its planted parameters back
  d <- simulate_rasch(600, 12, seed = 1)
  rec <- sim_recovery(rasch(d), d)
  expect_s3_class(rec, "rasch_recovery")
  s <- rec$summary
  expect_gt(s$correlation[s$parameter == "item difficulty"], 0.95)
  # person ability is noisier (WLE precision from only 12 items limits it)
  expect_gt(s$correlation[s$parameter == "person ability"], 0.75)
  # bias is not identifiable for an origin-centred location parameter, so it
  # is reported NA rather than a structurally-zero value
  expect_true(is.na(s$bias[s$parameter == "item difficulty"]))
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot_recovery(rec))

  # recovery across the other layouts
  d <- simulate_btl(8, 12, seed = 2)
  rb <- sim_recovery(btl(d, "object_a", "object_b", winner = "winner",
                         judge = "judge"), d)
  expect_gt(rb$summary$correlation[1], 0.9)

  # paired-comparison frames use their fitted identifying conventions: the
  # common object scale, geometric-mean-one panel units, and the first set as
  # the unit and origin reference
  d <- simulate_btl_efrm(5, 2, 5, 2, 8, 8,
    panel_units = c(1, 1.25), set_units = c(1, 1.3),
    set_origins = c(0, 0.4), seed = 71)
  tr <- attr(d, "truth")
  bf <- btl_efrm(d, "object_a", "object_b", winner = "winner",
                 judge = "judge", panels = "panel",
                 object_sets = tr$object_sets, se_method = "conditional",
                 boot_reps = 0)
  rf <- sim_recovery(bf, d)
  expect_setequal(rf$summary$parameter,
    c("object location", "panel unit (log)", "set unit (log)", "set origin"))
  expect_gt(rf$summary$correlation[
    rf$summary$parameter == "object location"], 0.9)
})

test_that("person recovery matches generating truth by ID", {
  d <- simulate_rasch(500, 12, seed = 721)
  ord <- sample.int(nrow(d))
  shuffled <- d[ord, , drop = FALSE]
  fit <- rasch(shuffled, id = "id")
  rec <- sim_recovery(fit, shuffled)
  pp <- rec$pieces[["person ability"]]
  expect_identical(pp$label, as.character(fit$person$id)[
    is.finite(fit$person$theta)])
  expect_gt(rec$summary$correlation[
    rec$summary$parameter == "person ability"], 0.7)
})

test_that("audit fixes hold: PCM structure, truth honesty, recovery centring", {
  # PCM and RSM now genuinely differ: per-item threshold patterns for PCM,
  # one common pattern for RSM; PCM thresholds stay ordered
  d1 <- simulate_rasch(200, 8, model = "PCM", n_categories = 4, seed = 42)
  d2 <- simulate_rasch(200, 8, model = "RSM", n_categories = 4, seed = 42)
  expect_false(identical(as.matrix(d1[, 2:9]), as.matrix(d2[, 2:9])))
  rel1 <- lapply(attr(d1, "truth")$thresholds, function(t) round(t - mean(t), 3))
  expect_gt(length(unique(rel1)), 1)                       # PCM varies
  rel2 <- lapply(attr(d2, "truth")$thresholds, function(t) round(t - mean(t), 3))
  expect_length(unique(rel2), 1)                           # RSM common
  expect_true(all(vapply(attr(d1, "truth")$thresholds,
                         function(t) !is.unsorted(t), TRUE)))

  # dif without groups is an error, not a silent no-op with a false truth
  expect_error(simulate_rasch(50, 6, dif = list(items = "I03", uniform = 1)),
               "n_groups")
  # polytomous guessing warns and is not recorded as planted
  expect_warning(d <- simulate_rasch(50, 4, model = "PCM", n_categories = 4,
                                     guessing = 0.3, seed = 1), "dichotomous")
  expect_false(any(grepl("guessing", attr(d, "truth")$planted)))
  # disordered thresholds work at 3 categories and preserve the location
  t2 <- rasch:::.sim_thresholds(0.5, 2, 1.2, disordered = TRUE)
  expect_true(is.unsorted(t2)); expect_equal(mean(t2), 0.5)
  # careless overwrite is not double-counted in the truth
  d <- simulate_rasch(300, 10, model = "PCM", n_categories = 4, careless = 0.5,
                      response_style = list(type = "extreme", prop = 0.5),
                      seed = 3)
  tr <- attr(d, "truth")
  expect_length(intersect(tr$style_idx, tr$careless_idx), 0)
  # halo raters never overflow into NA when erratic raters shrink the pool
  d <- simulate_mfrm(30, 4, 5, erratic_raters = 0.4, halo = 0.8, seed = 1)
  expect_false(anyNA(attr(d, "truth")$halo))

  # person ability is centred in recovery: an asymmetric difficulty range
  # must not masquerade as person-ability bias. Bias is not identifiable
  # up to the origin, so it is reported NA; the alignment shows instead as
  # a high correlation with no residual scale error
  d <- simulate_rasch(400, 10, difficulty = c(0, 3), seed = 2)
  r <- sim_recovery(rasch(d), d)
  expect_true(is.na(r$summary$bias[r$summary$parameter == "person ability"]))
  expect_gt(r$summary$correlation[r$summary$parameter == "person ability"], 0.75)
  # MFRM recovery reports item difficulties from the item margins
  d <- simulate_mfrm(60, 5, 5, seed = 1)
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater")
  r <- sim_recovery(mf, d)
  expect_true("item difficulty" %in% r$summary$parameter)
  expect_gt(r$summary$correlation[r$summary$parameter == "item difficulty"], 0.9)
})

test_that("misfit layers compose: dependence and style respect DIF / 2nd dim", {
  # dependence regeneration keeps the target item's DIF
  d <- simulate_rasch(2000, 6, dif = list(items = "I02", uniform = 1.5),
                      n_groups = 2,
                      dependence = list(pairs = list(c("I01", "I02")),
                                        strength = 1), seed = 7)
  g <- attr(d, "truth")$groups
  expect_gt(mean(d$I02[g != "g2"]) - mean(d$I02[g == "g2"]), 0.12)
  # ...and the target item's second dimension
  d <- simulate_rasch(2000, 6, second_dim = list(items = "I02", rho = 0.2),
                      dependence = list(pairs = list(c("I01", "I02")),
                                        strength = 1), seed = 8)
  expect_lt(cor(d$I02, attr(d, "truth")$theta), 0.2)
  # response style keeps DIF for style-affected persons
  d <- simulate_rasch(3000, 6, model = "PCM", n_categories = 4,
                      dif = list(items = "I02", uniform = 2.5), n_groups = 2,
                      response_style = list(type = "extreme", prop = 0.5),
                      seed = 10)
  tr <- attr(d, "truth"); g <- tr$groups
  sty <- seq_len(nrow(d)) %in% tr$style_idx
  expect_gt(mean(d$I02[sty & g != "g2"]) - mean(d$I02[sty & g == "g2"]), 0.8)
})

test_that("btl_dimensionality reference honours fitted dependence effects", {
  skip_on_cran()   # heavy simulation; verified locally and on CI
  # one-dimensional data whose only structure is within-judge order effects,
  # fitted WITH order: the dependence-aware reference must not read the
  # order structure as a second attribute
  flags <- vapply(1:4, function(s) {
    d <- simulate_btl(8, 12, reps_per_pair = 30,
                      dependence = list(exposure = 1, carry_over = 0.8),
                      seed = 300 + s)
    bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
              order = "order")
    isTRUE(btl_dimensionality(bt, reps = 50)$leading_structured)
  }, TRUE)
  expect_lte(sum(flags), 1L)   # was ~36% false-positive before the fix
})

test_that("a seeded simulator call leaves the caller's RNG stream alone", {
  set.seed(99); before <- runif(3)
  set.seed(99); invisible(simulate_rasch(60, 5, seed = 7)); after <- runif(3)
  expect_equal(before, after)                     # stream not commandeered
  expect_identical(simulate_rasch(60, 5, seed = 7)$I01,
                   simulate_rasch(60, 5, seed = 7)$I01)   # still reproducible
  # and an absent stream is left absent rather than seeded behind the caller
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    rm(".Random.seed", envir = globalenv())
  invisible(simulate_btl(n_objects = 4, n_judges = 3, reps_per_pair = 2,
                         seed = 3))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("every positive planted proportion affects at least one observation", {
  styled <- simulate_rasch(
    10, 4, model = "PCM", n_categories = 3,
    response_style = list(type = "extreme", prop = 0.01, strength = 1),
    seed = 41)
  expect_length(attr(styled, "truth")$style_idx, 1L)

  careless <- simulate_rasch(10, 4, careless = 0.01, seed = 42)
  expect_length(attr(careless, "truth")$careless_idx, 1L)

  incomplete <- simulate_rasch(10, 4, speeded = 0.01,
                               missing = 0.001, seed = 43)
  truth <- attr(incomplete, "truth")
  expect_length(truth$speeded_idx, 1L)
  expect_length(truth$missing_cells, 1L)
  expect_true(all(is.na(as.matrix(incomplete[, paste0("I0", 1:4)]))[
    truth$missing_cells]))
  expect_true(any(grepl("1 response missing", truth$planted, fixed = TRUE)))

  erratic <- simulate_mfrm(10, 3, 2, erratic_raters = 0.01, seed = 44)
  expect_length(attr(erratic, "truth")$erratic, 1L)
  halo <- simulate_mfrm(10, 3, 2, halo = 0.01, seed = 45)
  expect_length(attr(halo, "truth")$halo, 1L)
  expect_error(simulate_mfrm(10, 3, 2, erratic_raters = 1,
                             halo = 0.01, seed = 45),
               "every rater is erratic")

  expect_warning(
    fully_careless <- simulate_rasch(
      10, 4, model = "PCM", n_categories = 3, careless = 1,
      response_style = list(type = "extreme", prop = 0.01, strength = 1),
      seed = 46),
    "response_style could not remain")
  expect_length(attr(fully_careless, "truth")$style_idx, 0L)

  framed <- simulate_efrm(n_per_group = 2, items_per_set = 2,
                          n_sets = 2, n_groups = 2,
                          careless = 0.01, missing = 0.001, seed = 47)
  truth <- attr(framed, "truth")
  expect_length(truth$careless_idx, 1L)
  expect_length(truth$missing_cells, 1L)
})
