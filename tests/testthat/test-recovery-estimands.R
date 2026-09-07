test_that("BTL-EFRM recovery refuses a different cross-set reference restriction", {
  d <- simulate_btl_efrm(5, 2, 5, 2, 8, 8,
    panel_units = c(1, 1.25), set_units = c(1, 1.3),
    set_origins = c(0, .4), seed = 71)
  tr <- attr(d, "truth")
  sets <- tr$object_sets
  names(sets) <- c("Z form", "A form")
  fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
    judge = "judge", panels = "panel", object_sets = sets,
    se_method = "conditional", boot_reps = 0)
  expect_identical(fit$reference_set, "A form")
  expect_equal(fit$alpha_table$alpha[fit$alpha_table$set == "A form"], 1)
  expect_equal(unname(tr$alpha[2L]), 1.3)
  expect_error(sim_recovery(fit, d), "different cross-set restriction")

  unchanged <- btl_efrm(d, "object_a", "object_b", winner = "winner",
    judge = "judge", panels = "panel", object_sets = tr$object_sets,
    se_method = "conditional", boot_reps = 0)
  expect_no_error(sim_recovery(unchanged, d))
})

test_that("BTL-EFRM recovery aligns an equivalent reference origin", {
  d <- simulate_btl_efrm(5, 2, 5, 2, 8, 8,
    panel_units = c(1, 1.25), set_units = c(1, 1),
    set_origins = c(0, .4), seed = 71)
  tr <- attr(d, "truth")
  sets <- tr$object_sets
  names(sets) <- c("Z form", "A form")
  fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
    judge = "judge", panels = "panel", object_sets = sets,
    se_method = "conditional", boot_reps = 0)
  rec <- sim_recovery(fit, d)
  expect_null(rec$note)
  expect_equal(rec$pieces[["object location"]]$true,
    unname((tr$v - .4)[rec$pieces[["object location"]]$label]))
  expect_equal(rec$pieces[["set origin"]]$true, c(-.4, 0))

  # Exact generating values on the fitted origin must have zero recovery
  # error. This isolates alignment from the sampling error in a real fit.
  exact <- fit
  exact$objects$v <- unname((tr$v - .4)[exact$objects$object])
  exact$objects$location <- exact$objects$v
  exact$phi_table$phi <- unname(tr$phi[exact$phi_table$panel])
  exact$alpha_table$alpha <- 1
  shifted_origins <- setNames(c(-.4, 0), names(sets))
  exact$kappa_table$kappa <- unname(shifted_origins[exact$kappa_table$set])
  exact_rec <- sim_recovery(exact, d)
  expect_equal(exact_rec$summary$rmse, rep(0, 4))
  expect_equal(exact_rec$summary$bias, rep(0, 4))
})

test_that("recovery refuses category compression in an actual PCM fit", {
  d <- simulate_rasch(25, 5, model = "PCM", n_categories = 5,
    difficulty = c(-1, 1), seed = 1)
  fit <- rasch(d)
  expect_true(fit$est$converged)
  expect_identical(ncol(fit$X), 5L)
  expect_true(any(lengths(fit$tau_list) != lengths(attr(d, "truth")$thresholds)))
  expect_error(sim_recovery(fit, d), "after category removal or merging")

  complete <- simulate_rasch(300, 5, model = "PCM", n_categories = 5,
    difficulty = c(-1, 1), seed = 1)
  complete <- complete[rev(seq_len(nrow(complete))), ]
  ordinary <- rasch(complete, items = rev(names(attr(complete, "truth")$thresholds)))
  expect_no_error(sim_recovery(ordinary, complete))
})

test_that("frame recovery checks the generating scale by response cell", {
  d <- simulate_efrm(30, 3, n_sets = 1, n_groups = 2,
    n_categories = 3, missing = .01, seed = 17)
  tr <- attr(d, "truth")
  fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
    id = "id", boot_reps = 0)
  expect_no_error(sim_recovery(fit, d))
  expect_no_error(.recovery_check_categories(fit, d, tr))
  compressed <- fit
  compressed$tau_list[[1L]] <- compressed$tau_list[[1L]][1L]
  expect_error(.recovery_check_categories(compressed, d, tr),
    "generating scale")

  dm <- simulate_mfrm(30, 3, 3, n_categories = 3, seed = 18)
  tm <- attr(dm, "truth")
  expect_identical(tm$n_categories, 3L)
  fm <- rasch_mfrm(dm, person = "person", item = "item", score = "score",
    facets = "rater")
  expect_no_error(sim_recovery(fm, dm))
  expect_null(sim_recovery(fm, dm)$note)
  legacy <- dm
  attr(legacy, "truth")$n_categories <- NULL
  expect_match(sim_recovery(fm, legacy)$note,
    "original category count was not recorded.*cannot be verified")
  fm$tau_list[[1L]] <- fm$tau_list[[1L]][1L]
  expect_error(.recovery_check_categories(fm, dm, tm), "generating scale")
  tm$n_categories <- NULL
  expect_error(.recovery_check_categories(fm, dm, tm), "generating scale")
})
