test_that("saved contrasts with superseded support rules are not restored", {
  d <- simulate_rasch(200, 6, n_groups = 2, seed = 814)
  fit <- rasch(d, factors = "group")
  current <- dif_contrasts(fit, items = "I01")
  expect_identical(current$algorithm, "complete-contrast-cells-1")
  old <- current
  old$algorithm <- NULL
  old$table$estimate <- old$table$estimate + 1
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    model_type = "rasch", data = as.data.frame(d), base_fit = fit,
    rasch_steps = list(), btl_steps = list(), kept_fits = list(base = fit),
    settings = list(pc_items = "I01"),
    results = list(contrasts = old, display = list(item = "I01"))))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)
  expect_error(.save_app_project(project, path), "superseded calculation")
  expect_warning(restored <- .read_app_project(path), "complete-cell.*omitted")
  expect_null(restored$results$contrasts)
  for (field in c("data", "base_fit", "rasch_steps", "btl_steps", "kept_fits",
                  "settings"))
    expect_identical(restored[[field]], project[[field]])
  expect_identical(restored$results$display, project$results$display)
  expect_match(attr(restored, "rasch_project_legacy_dropped"), "planned DIF")
  expect_no_error(.validate_app_project(restored))
  expect_no_error(.save_app_project(restored, path))
  expect_no_warning(.read_app_project(path))

  changed <- project
  changed$results$contrasts$table$p_adj <- 0
  saveRDS(changed, path)
  expect_error(.read_app_project(path), "changed since they were saved")

  project$results$contrasts <- current
  project <- .seal_app_project(project)
  expect_no_error(.save_app_project(project, path))
  expect_no_warning(restored <- .read_app_project(path))
  expect_identical(restored$results$contrasts, current)

  project$schema <- 1L
  project$binding <- NULL
  project$results$contrasts <- old
  saveRDS(project, path)
  expect_warning(restored <- .read_app_project(path), "schema-1")
  expect_null(restored$results$contrasts)
  expect_identical(restored$base_fit, fit)
  expect_no_error(.validate_app_project(restored))
})
