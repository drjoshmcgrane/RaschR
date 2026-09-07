test_that("fixed-origin CJ equating cannot promote conditional frame errors", {
  d <- simulate_btl_efrm(n_objects_per_set = 4, n_sets = 2, n_panels = 2,
    n_judges_per_panel = 10, reps_within = 15, reps_cross = 15, seed = 73)
  truth <- attr(d, "truth")
  conditional <- btl_efrm(d, "object_a", "object_b", winner = "winner",
    judge = "judge", panels = "panel", object_sets = truth$object_sets,
    se_method = "conditional", boot_reps = 0)
  bank <- data.frame(object = names(truth$v), location = unname(truth$v), se = 0)
  expect_true(conditional$converged)
  expect_true(all(is.finite(conditional$objects$se)))
  expect_true(all(is.na(conditional$alpha_table$p)))
  expect_null(conditional$cov_beta)

  for (shift in c("mean", "none")) {
    eq <- btl_equate(conditional, bank, shift = shift)
    expect_false(eq$inferential)
    expect_true(all(is.na(eq$table$p_adj)))
    expect_true(all(is.na(eq$table$se_diff)))
    expect_true(all(is.finite(eq$equated$location)))
    expect_match(paste(eq$notes, collapse = " "), "conditional frame errors")
    if (shift == "none") {
      # The independently fixed bank remains fixed when no estimated shift
      # from the conditional calibration is applied to it.
      expect_identical(eq$equated$se, bank$se)
      expect_equal(eq$shift_se, 0)
    } else {
      expect_true(all(is.na(eq$equated$se)))
    }
  }

  other_data <- simulate_btl_efrm(n_objects_per_set = 4, n_sets = 2,
    n_panels = 2, n_judges_per_panel = 10, reps_within = 15,
    reps_cross = 15, seed = 74)
  ordinary <- btl(other_data, "object_a", "object_b", winner = "winner",
                   judge = "judge")
  expect_true(ordinary$converged)
  expect_true(ordinary$cl$inference_available)
  expect_false(identical(ordinary$thr_structure, conditional$thr_structure))
  valid <- btl_equate(ordinary, bank, shift = "none")
  expect_true(valid$inferential)
  expect_true(all(is.finite(valid$table$p_adj)))

  for (shift in c("mean", "none")) {
    eq <- btl_equate(ordinary, conditional, independent = TRUE, shift = shift)
    expect_false(eq$inferential)
    expect_true(all(is.na(eq$table$p_adj)))
    expect_true(all(is.na(eq$equated$se)))
    expect_null(attr(eq$equated, "cov_location"))
    expect_null(attr(eq$equated, "df_location"))
    expect_equal(eq$equated$location,
                 conditional$objects$location + eq$shift)
    expect_match(paste(eq$notes, collapse = " "),
                 "Equated location standard errors.*conditional frame errors")
    # Reusing the returned bank cannot silently recover inference after
    # removal of the fitted-model class and its conditional-method metadata.
    reused <- btl_equate(ordinary, eq$equated, shift = "none")
    expect_false(reused$inferential)
    expect_true(all(is.na(reused$table$p_adj)))
  }
  reverse <- btl_equate(conditional, ordinary, independent = TRUE, shift = "none")
  expect_false(reverse$inferential)
  expect_identical(reverse$equated$se, ordinary$objects$se)
  expect_equal(attr(reverse$equated, "cov_location"), ordinary$cov_beta)
})

test_that("parametric-bootstrap frame errors retain fixed-origin inference", {
  d <- simulate_btl_efrm(n_objects_per_set = 4, n_sets = 1, n_panels = 1,
    n_judges_per_panel = 12, reps_within = 25, seed = 75)
  truth <- attr(d, "truth")
  bootstrap <- btl_efrm(d, "object_a", "object_b", winner = "winner",
    judge = "judge", panels = "panel", object_sets = truth$object_sets,
    se_method = "bootstrap", boot_reps = 30, seed = 76)
  bank <- data.frame(object = names(truth$v), location = unname(truth$v), se = 0)
  expect_true(bootstrap$converged)
  expect_true(all(is.finite(bootstrap$objects$se)))
  eq <- btl_equate(bootstrap, bank, shift = "none")
  expect_true(eq$inferential)
  expect_true(all(is.finite(eq$table$p_adj)))
  expect_true(all(is.infinite(eq$table$df)))
  expect_equal(eq$table$p, 2 * pnorm(-abs(eq$table$t)), tolerance = 1e-12)
})

test_that("polytomous equating still requires compatible threshold structures", {
  d <- simulate_btl(5, 20, reps_per_pair = 20, model = "polytomous",
                     n_categories = 4, seed = 77)
  free <- btl(d, "object_a", "object_b", response = "response", judge = "judge",
                thresholds = "free")
  pc <- btl(d, "object_a", "object_b", response = "response", judge = "judge",
              thresholds = "pc")
  expect_identical(free$m, 3L)
  expect_true(free$converged && pc$converged)
  expect_false(identical(free$thr_structure, pc$thr_structure))
  expect_error(btl_equate(free, pc), "incompatible response scales or threshold")
})
