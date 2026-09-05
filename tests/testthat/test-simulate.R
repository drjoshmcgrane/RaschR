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

test_that("response-style probabilities remain stable at large strengths", {
  extreme <- simulate_rasch(
    40, 5, model = "PCM", n_categories = 3,
    response_style = list(type = "extreme", prop = 0.5, strength = 1000),
    seed = 91)
  sx <- attr(extreme, "truth")$style_idx
  expect_true(all(as.matrix(extreme[sx, paste0("I0", 1:5)]) %in% c(0L, 2L)))

  middle <- simulate_rasch(
    40, 5, model = "PCM", n_categories = 3,
    response_style = list(type = "middle", prop = 0.5, strength = 1000),
    seed = 92)
  sm <- attr(middle, "truth")$style_idx
  expect_true(all(as.matrix(middle[sm, paste0("I0", 1:5)]) == 1L))
})

test_that("a requested threshold disorder is present in the generating truth", {
  # PCM threshold spans vary by item. The disorder must therefore be imposed
  # as an adjacent reversal, not as a subtraction that a wide random gap can
  # sometimes absorb.
  for (seed in 1:30) {
    d <- simulate_rasch(
      40, 4, model = "PCM", n_categories = 3,
      disordered = "I01", seed = seed)
    tr <- attr(d, "truth")
    expect_named(tr$thresholds, sprintf("I%02d", 1:4))
    expect_true(any(diff(tr$thresholds$I01) < 0))
    expect_equal(mean(tr$thresholds$I01), tr$difficulty[["I01"]],
                 tolerance = 1e-12)
  }
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

test_that("explanatory simulation departures cannot be absorbed by the design", {
  p <- data.frame(exposure = seq(-1, 1, length.out = 5),
                  type = rep(c("A", "B"), length.out = 5))
  B <- stats::model.matrix(~ exposure * type, p)
  d <- .sim_explanatory_departure(B, 1.25)
  expect_equal(unname(drop(crossprod(B, d$values))), rep(0, ncol(B)),
               tolerance = 1e-10)
  expect_equal(d$values[d$index], 1.25, tolerance = 1e-12)

  saturated <- stats::model.matrix(
    ~ exposure * type,
    data.frame(exposure = seq(-1, 1, length.out = 4),
               type = rep(c("A", "B"), length.out = 4)))
  expect_error(.sim_explanatory_departure(saturated, 1), "saturated")
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

test_that("BTL-EFRM simulator unit arguments are plain vectors", {
  expect_error(
    simulate_btl_efrm(n_panels = 2, panel_units = matrix(c(1, 1.2), 1)),
    "panel_units")
  expect_error(
    simulate_btl_efrm(n_sets = 2, set_units = matrix(c(1, 1.2), 1)),
    "set_units")
  expect_error(
    simulate_btl_efrm(n_sets = 2, set_origins = matrix(c(0, 0.2), 1)),
    "set_origins")
  expect_error(
    simulate_btl_efrm(n_sets = 2, set_units = c(1e-300, 1e300)),
    "numerically representable")
  expect_error(
    simulate_btl_efrm(n_sets = 2, set_origins = c(-1e308, 1e308)),
    "numerically representable")
})

test_that("EFRM simulator normalises very small ratios without underflow", {
  d <- simulate_efrm(n_per_group = 2, items_per_set = 2, n_sets = 2,
                     n_groups = 1, set_unit_ratio = 1e-300, seed = 713)
  expect_true(all(is.finite(attr(d, "truth")$alpha)))
  expect_true(all(attr(d, "truth")$alpha > 0))
})

test_that("planted-misfit selectors are plain vectors", {
  expect_error(
    simulate_rasch(40, 4, disordered = matrix("I01", 1),
                   model = "PCM"),
    "plain vectors")
  expect_error(
    simulate_mfrm(10, 3, 3,
      interaction = list(rater = matrix("R1", 1), item = "I1", bias = 1)),
    "each name one")
  expect_error(
    simulate_efrm(10, 3,
      item_drift = list(items = matrix("S1I01", 1), group = "g1", shift = 1)),
    "plain vector")
  expect_error(
    simulate_efrm(10, 3,
      item_drift = list(items = "S1I01", group = matrix("g1", 1), shift = 1)),
    "name one generated group")
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
  fit <- rasch(d)
  rec <- sim_recovery(fit, d)
  expect_s3_class(rec, "rasch_recovery")
  s <- rec$summary
  expect_gt(s$correlation[s$parameter == "item difficulty"], 0.95)
  # person ability is noisier (WLE precision from only 12 items limits it)
  expect_gt(s$correlation[s$parameter == "person ability"], 0.75)
  # bias is not identifiable for an origin-centred location parameter, so it
  # is reported NA rather than a structurally-zero value
  expect_true(is.na(s$bias[s$parameter == "item difficulty"]))
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot_recovery(rec))

  bad <- fit
  bad$est$converged <- FALSE
  expect_error(sim_recovery(bad, d), "did not converge")

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

  # Set names are arbitrary labels: recovery follows the object partition,
  # not the simulator's conventional set1/set2 spelling.
  renamed <- bf
  set_names <- c(set1 = "Form A", set2 = "Form B")
  renamed$objects$set <- unname(set_names[as.character(renamed$objects$set)])
  renamed$alpha_table$set <- unname(set_names[as.character(
    renamed$alpha_table$set)])
  renamed$kappa_table$set <- unname(set_names[as.character(
    renamed$kappa_table$set)])
  rr <- sim_recovery(renamed, d)
  expect_setequal(rr$summary$parameter,
    c("object location", "panel unit (log)", "set unit (log)", "set origin"))
})

test_that("EFRM recovery follows set membership rather than set labels", {
  d <- simulate_efrm(20, 4, n_sets = 2, n_groups = 2, seed = 72)
  tr <- attr(d, "truth")
  set_of <- stats::setNames(
    rep(c("Form A", "Form B"), lengths(tr$item_sets)),
    unlist(tr$item_sets, use.names = FALSE))
  fit <- list(
    est = list(converged = TRUE),
    linking = list(alpha_edges = data.frame(converged = TRUE)),
    set_of = set_of,
    alpha_table = data.frame(set = c("Form A", "Form B"), alpha = tr$alpha),
    phi_table = data.frame(group = names(tr$phi), phi = tr$phi))
  class(fit) <- c("rasch_efrm", "rasch")
  r <- sim_recovery(fit, d)
  expect_setequal(r$summary$parameter,
                  c("set unit (log)", "group unit (log)"))
  expect_identical(r$pieces[["set unit (log)"]]$label,
                   c("Form A", "Form B"))

  repartitioned <- fit
  repartitioned$set_of[tr$item_sets[[1]][1]] <- "Form B"
  expect_error(sim_recovery(repartitioned, d),
               "partition does not match")
})

test_that("sim_replicate preserves a seed at the integer boundary", {
  seen <- integer(0)
  generator <- function(seed) {
    seen <<- c(seen, seed)
    structure(data.frame(x = 1), truth = list(layout = "test"))
  }
  out <- sim_replicate(generator, 1, seed = .Machine$integer.max)
  expect_identical(seen, .Machine$integer.max)
  expect_length(out, 1L)
  expect_error(sim_replicate(generator, 2, seed = .Machine$integer.max),
               "exceeds the integer range")
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

test_that("recovery retains an externally anchored origin", {
  d <- simulate_rasch(600, 8, seed = 720)
  tr <- attr(d, "truth")
  anchors <- data.frame(item = names(tr$difficulty)[c(1, 8)], k = 1L,
                        tau = unname(tr$difficulty[c(1, 8)]))
  fit <- rasch(d, id = "id", anchors = anchors)
  rec <- sim_recovery(fit, d)
  ip <- rec$pieces[["item difficulty"]]
  pp <- rec$pieces[["person ability"]]
  expect_equal(ip$estimated,
               fit$items$location[match(ip$label, fit$items$item)])
  expect_equal(pp$true,
               unname(tr$theta[match(pp$label, tr$person_id)]))
  expect_equal(pp$estimated,
               fit$person$theta[match(pp$label, fit$person$id)])
  expect_true(all(is.finite(rec$summary$bias)))
  expect_equal(rec$summary$bias[rec$summary$parameter == "item difficulty"],
               mean(fit$items$location - unname(tr$difficulty[
                 match(fit$items$item, names(tr$difficulty))])))

  d <- simulate_btl(7, 20, 10, seed = 721)
  tr <- attr(d, "truth")
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
             anchors = tr$location[c(1, 7)])
  rec <- sim_recovery(fit, d)
  expect_true(is.finite(rec$summary$bias[
    rec$summary$parameter == "object location"]))
})

test_that("the recovery plot refuses an unrelated object directly", {
  expect_error(plot_recovery(list()), "recovery result")
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
  # requested proportions are either realised exactly or refused; one
  # mechanism cannot silently erase part of another
  expect_error(simulate_mfrm(30, 4, 5, erratic_raters = 0.4, halo = 0.8,
                             seed = 1), "cannot coexist")
  d <- simulate_mfrm(30, 4, 5, erratic_raters = 0.4, halo = 0.6, seed = 1)
  expect_length(attr(d, "truth")$erratic, 2L)
  expect_length(attr(d, "truth")$halo, 3L)

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

test_that("MFRM recovery identifies the simulated rater facet", {
  d <- simulate_mfrm(12, 3, 4, seed = 73)
  tr <- attr(d, "truth")
  # A fitted object may contain additional facets and need not put the rater
  # facet first.  Its levels, rather than its position, identify the planted
  # severity parameter when the rater column has also been renamed.
  fit <- structure(list(
    est = list(converged = TRUE),
    facet_effects = list(
      occasion = data.frame(level = c("O1", "O2"), severity = c(-0.2, 0.2)),
      judge = data.frame(level = names(tr$severity),
                         severity = unname(tr$severity))
    ),
    item_effects = data.frame(item = names(tr$difficulty),
                              location = unname(tr$difficulty)),
    person = data.frame(id = tr$person_id, theta = tr$theta)
  ), class = c("rasch_mfrm", "rasch"))
  r <- sim_recovery(fit, d)
  rs <- r$pieces[["rater severity"]]
  expect_identical(rs$label, names(tr$severity))
  expect_equal(rs$estimated, as.numeric(tr$severity - mean(tr$severity)))

  ambiguous <- fit
  ambiguous$facet_effects$occasion <- ambiguous$facet_effects$judge
  expect_error(sim_recovery(ambiguous, d), "cannot be matched to one")

  misleading <- fit
  misleading$facet_effects$rater <- data.frame(
    level = c(names(tr$severity)[1], "unrelated"),
    severity = c(99, -99))
  r2 <- sim_recovery(misleading, d)
  expect_identical(r2$pieces[["rater severity"]]$label,
                   names(tr$severity))

  ordinary <- rasch(simulate_rasch(80, 4, seed = 74))
  expect_error(sim_recovery(ordinary, d), "does not match")
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

  # Response style redraws the source response. The target must therefore be
  # regenerated from that styled source, not from the response it replaced.
  d <- simulate_rasch(
    4000, 6, model = "PCM", n_categories = 4,
    dependence = list(pairs = list(c("I01", "I02")), strength = 2),
    response_style = list(type = "extreme", prop = 0.5, strength = 2),
    seed = 771)
  tr <- attr(d, "truth")
  sty <- seq_len(nrow(d)) %in% tr$style_idx
  dep_coef <- function(rows) unname(stats::coef(stats::lm(
    d$I02[rows] ~ d$I01[rows] +
      stats::poly(tr$theta[rows], 5, raw = TRUE)))[2L])
  expect_gt(dep_coef(sty), 0.25)
  expect_gt(dep_coef(!sty), 0.20)
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

  speeded_only <- simulate_rasch(20, 5, speeded = 1, seed = 430)
  both <- simulate_rasch(20, 5, speeded = 1, missing = 0.2, seed = 430)
  speeded_cells <- which(is.na(as.matrix(speeded_only[, paste0("I0", 1:5)])))
  both_truth <- attr(both, "truth")
  expect_length(both_truth$missing_cells, 20L)
  expect_length(intersect(speeded_cells, both_truth$missing_cells), 0L)
  expect_error(simulate_rasch(10, 3, speeded = 1, missing = 1, seed = 431),
               "finite value")

  erratic <- simulate_mfrm(10, 3, 2, erratic_raters = 0.01, seed = 44)
  expect_length(attr(erratic, "truth")$erratic, 1L)
  halo <- simulate_mfrm(10, 3, 2, halo = 0.01, seed = 45)
  expect_length(attr(halo, "truth")$halo, 1L)
  expect_error(simulate_mfrm(10, 3, 2, erratic_raters = 1,
                             halo = 0.01, seed = 45),
               "cannot coexist")

  expect_error(simulate_rasch(
    10, 4, model = "PCM", n_categories = 3, careless = 1,
    response_style = list(type = "extreme", prop = 0.01, strength = 1),
    seed = 46), "cannot coexist")

  framed <- simulate_efrm(n_per_group = 2, items_per_set = 2,
                          n_sets = 2, n_groups = 2,
                          careless = 0.01, missing = 0.001, seed = 47)
  truth <- attr(framed, "truth")
  expect_length(truth$careless_idx, 1L)
  expect_length(truth$missing_cells, 1L)

  no_drift <- simulate_efrm(n_per_group = 4, items_per_set = 2,
                            n_sets = 2, n_groups = 2, seed = 48)
  zero_drift <- simulate_efrm(
    n_per_group = 4, items_per_set = 2, n_sets = 2, n_groups = 2,
    item_drift = list(items = "S1I01", group = "g2", shift = 0), seed = 48)
  expect_identical(as.data.frame(zero_drift), as.data.frame(no_drift))
  expect_null(attr(zero_drift, "truth")$item_drift)
})

test_that("response-data simulators cannot remove every observation", {
  expect_error(simulate_rasch(20, 5, missing = 1), "finite value")
  expect_error(simulate_efrm(10, 3, missing = 1), "finite value")
  expect_error(simulate_rasch(2, 2, missing = 0.99), "every remaining response")
  expect_error(simulate_efrm(2, 2, missing = 0.99), "every response")
})
