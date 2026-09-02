.app_test_path <- function() {
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  app
}

test_that("the Shiny application resolves in source and installed layouts", {
  expect_true(file.exists(.app_test_path()))
})

test_that("app formula and code symbols preserve exact source names", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))

  nms <- c("if", ".2x", "a b", "a:b", "x`y")
  quoted <- e$.app_quote_name(nms)
  expect_no_error(lapply(quoted, function(z)
    parse(text = paste0("dat$", z))))

  spec <- e$.app_explanatory_interactions(c("a:b", "if", "c"))
  expect_length(unique(unname(spec$choices)), 3L)
  expect_identical(spec$map[["interaction_0001"]], c("a:b", "if"))
  f <- e$.app_explanatory_formula(
    c("a:b", "if"), list(c("a:b", "if")))
  d <- data.frame(row = seq_len(6L)); d$row <- NULL
  d[["a:b"]] <- seq_len(6L)
  d[["if"]] <- rep(0:1, 3L)
  expect_no_error(mm <- model.matrix(f, d))
  expect_setequal(colnames(mm),
                  c("(Intercept)", "`a:b`", "`if`", "`a:b`:`if`"))
})

test_that("sourcing the in-tree app reuses the active source namespace", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")
  before <- get(".rasch_boot_apply", envir = asNamespace("rasch"),
                inherits = FALSE)
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  after <- get(".rasch_boot_apply", envir = asNamespace("rasch"),
               inherits = FALSE)
  expect_identical(after, before)
})

test_that("saved app analyses make a validated round trip", {
  set.seed(81)
  theta <- rnorm(120)
  X <- sapply(seq(-1, 1, length.out = 5), function(d)
    rbinom(length(theta), 1, plogis(theta - d)))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(X)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project",
    schema = 2L,
    package_version = "test",
    created = "2026-08-16T00:00:00+1000",
    data = as.data.frame(X),
    model_type = "rasch",
    base_fit = fit,
    rasch_steps = list(),
    btl_steps = list(),
    rcode = "fit <- rasch(X)",
    kept_fits = list(reference = fit),
    kept_fit_code = list(reference = list(code = "fit <- rasch(X)",
                                          value = "fit")),
    settings = list(model_type = "rasch", item_cols = colnames(X),
                    ng_auto = TRUE, maxit = 75),
    resources = list(key = NULL),
    simulation = list(),
    results = list()
  ))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))

  expect_identical(.save_app_project(project, path), path)
  restored <- .read_app_project(path)
  expect_identical(restored$format, "rasch-shiny-project")
  expect_equal(restored$data, project$data)
  expect_equal(restored$base_fit$items, fit$items)
  expect_identical(restored$rcode, project$rcode)
  expect_identical(restored$kept_fit_code, project$kept_fit_code)
  expect_identical(restored$settings, project$settings)
  expect_identical(restored$resources, project$resources)

  # Schema-2 projects written with the earlier text signature remain valid
  # across the LF/CRLF boundary.
  legacy <- project
  unsigned <- legacy
  unsigned$binding <- NULL
  legacy$binding <- .fit_boot_md5_legacy_candidates(unsigned)[2L]
  expect_no_error(.validate_app_project(legacy))
})

test_that("invalid app analysis files are refused", {
  bad <- tempfile(fileext = ".rasch")
  on.exit(unlink(bad))
  saveRDS(list(format = "other", schema = 1L), bad)
  expect_error(.read_app_project(bad), "not a rasch analysis file")

  set.seed(2)
  X <- matrix(rbinom(300, 1, .5), 60, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X)
  future <- list(format = "rasch-shiny-project", schema = 3L,
                 data = as.data.frame(X), base_fit = fit)
  saveRDS(future, bad)
  expect_error(.read_app_project(bad), "unsupported")

  malformed <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), base_fit = fit, settings = "not a list"))
  saveRDS(malformed, bad)
  expect_error(.read_app_project(bad), "invalid settings")

  legacy <- list(format = "rasch-shiny-project", schema = 1L,
                 data = as.data.frame(X), model_type = "rasch",
                 base_fit = fit, rasch_steps = list(), btl_steps = list(),
                 results = list())
  saveRDS(legacy, bad)
  expect_warning(upgraded <- .read_app_project(bad), "schema-1")
  expect_identical(upgraded$schema, 2L)
  expect_true(isTRUE(attr(upgraded, "rasch_project_legacy")))
  expect_no_error(.validate_app_project(upgraded))

  empty_fit <- structure(list(), class = "rasch")
  malformed_fit <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = empty_fit,
    rasch_steps = list(), btl_steps = list(), results = list()))
  expect_error(.validate_app_project(malformed_fit), "saved base fit")

  bound <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), results = list()))
  expect_no_error(.validate_app_project(bound))
  matrix_bound <- bound
  matrix_bound$data <- X
  matrix_bound <- .seal_app_project(matrix_bound)
  expect_no_error(.validate_app_project(matrix_bound))
  bound$data <- data.frame(unrelated = seq_len(nrow(X)))
  expect_error(.validate_app_project(bound), "changed since they were saved")
})

test_that("saved bootstrap results must belong to the active saved fit", {
  set.seed(83)
  X1 <- matrix(rbinom(600, 1, 0.5), 120, 5,
               dimnames = list(NULL, paste0("I", 1:5)))
  X2 <- X1
  X2[1:10, 1] <- 1L - X2[1:10, 1]
  fit1 <- rasch(X1)
  fit2 <- rasch(X2)
  bs <- suppressWarnings(fit_bootstrap(
    fit1, B = 1L, workers = 1L, seed = 83L))
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X1), model_type = "rasch", base_fit = fit1,
    rasch_steps = list(), btl_steps = list(), results = list(
      bootstrap = list(bs = bs, B = 1L, seed = 1L, kind = "rasch"))))
  expect_no_error(.validate_app_project(project))

  history_entry <- function(type, fit) list(
    type = type, label = type, fit = fit, details = list(), code = NULL,
    created = "2026-08-16T00:00:00+1000")
  project$rasch_steps <- list(history_entry("test", fit2))
  expect_error(.validate_app_project(project),
               "does not belong to the active fit")
  project$rasch_steps <- list(
    history_entry("corrupt", "not a fitted model"),
    history_entry("test", fit1))
  expect_error(.validate_app_project(project), "fitted-model history")
  wrong_family <- fit1
  class(wrong_family) <- c("rasch_efrm", class(wrong_family))
  project$rasch_steps <- list(history_entry("wrong", wrong_family))
  expect_error(.validate_app_project(project), "fitted-model history")
  project$rasch_steps <- list()
  project$results$bootstrap <- list(bs = list())
  expect_error(.validate_app_project(project), "saved bootstrap result")

  project$results <- list()
  project$base_fit <- wrong_family
  expect_error(.validate_app_project(project), "declares model type")
})

test_that("opening a project retains results tied to its active fit", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  set.seed(82)
  X <- matrix(rbinom(500, 1, .5), 100, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), rcode = "fit <- rasch(dat)",
    kept_fits = list(), kept_fit_code = list(),
    settings = list(model_type = "rasch", id_col = "(none)",
                    factor_cols = character(0), item_cols = colnames(X),
                    thr_structure = "pcm", ng_auto = FALSE, ng = "6",
                    maxit = 75, tol = 1e-7),
    resources = list(anchors = NULL, key = NULL), simulation = list(),
    results = list(lr = list(marker = 17L))))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  .save_app_project(project, path)

  shiny::testServer(e$server, {
    session$setInputs(project_file = list(
      datapath = path, name = "test.rasch", size = file.info(path)$size,
      type = "application/octet-stream"))
    session$flushReact()
    expect_identical(lr_res()$marker, 17L)
    expect_s3_class(fit(), "rasch")
    expect_equal(raw_data(), as.data.frame(X))
    expect_no_error(src <- data_source_code())
    expect_match(src, 'readRDS\\("test[.]rasch"\\)')
    expect_match(src, "dat <- project[$]data")
    expect_identical(restored_project_settings(), project$settings)
    expect_identical(restored_project_resources(), project$resources)
    html <- paste(as.character(output$data_main), collapse = " ")
    expect_false(grepl("Welcome to rasch", html, fixed = TRUE))
    expect_match(html, "Data preview", fixed = TRUE)

    # Saving again captures the current run-defining controls rather than
    # relying on column-name guesses at the next opening.
    session$setInputs(id_col = "(none)", item_cols = colnames(X),
                      factor_cols = character(0), thr_structure = "pcm",
                      ng_auto = FALSE, ng = "6", maxit = 75, tol = 1e-7)
    saved <- project_state()
    expect_identical(saved$settings$item_cols, colnames(X))
    expect_identical(saved$settings$ng, "6")
    expect_identical(saved$settings$maxit, 75)
  })
})

test_that("CJ results and background work remain tied to their launching fit", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  d1 <- simulate_btl(5, 12, 8, seed = 91)
  d2 <- simulate_btl(5, 12, 8, seed = 92)
  bt1 <- btl(d1, "object_a", "object_b", "winner", judge = "judge")
  bt2 <- btl(d2, "object_a", "object_b", "winner", judge = "judge")
  judges <- unique(bt1$comparisons$judge)
  group <- setNames(rep(c("A", "B"), length.out = length(judges)), judges)
  bdif1 <- btl_dif(bt1, group, objects = "O3", min_n = 2)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(d1), model_type = "btl", base_fit = bt1,
    rasch_steps = list(), btl_steps = list(), rcode = "bt <- btl(dat)",
    kept_fits = list(), kept_fit_code = list(),
    settings = list(model_type = "btl", bt_a = "object_a",
                    bt_b = "object_b", bt_win = "winner",
                    bt_judge = "judge"),
    resources = list(), simulation = list(),
    results = list(btl_dif = bdif1,
                   btl_dif_meta = list(judge_col = "judge"))))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  .save_app_project(project, path)

  shiny::testServer(e$server, {
    # A successful replacement owns none of the preceding fit's requested
    # results, even when the replacement is another CJ calibration.
    bdif_res(list(marker = 12L))
    btlef_res(list(marker = 11L))
    complete_fit(bt2, "bt <- btl(dat)", character(0), "dat <- replacement")
    session$flushReact()
    expect_null(bdif_res())
    expect_null(btlef_res())

    # The fit observer must not erase results serialised with a reopened fit.
    session$setInputs(project_file = list(
      datapath = path, name = "cj.rasch", size = file.info(path)$size,
      type = "application/octet-stream"))
    session$flushReact()
    expect_identical(bdif_res()$result_signature, bdif1$result_signature)
    expect_identical(bdif_meta()$judge_col, "judge")
    expect_no_error(rasch:::.validate_btl_dif_result(bdif_res(), bfit()))
    expect_no_error(.validate_app_project(project_state()))
    expect_identical(btl_fit()$objects$object, bt1$objects$object)

    # Opening a project cancels either kind of worker before restoring state.
    fake_process <- function() {
      p <- new.env(parent = emptyenv())
      p$alive <- TRUE
      p$is_alive <- function() p$alive
      p$kill_tree <- function() p$alive <- FALSE
      p$kill <- function() p$alive <- FALSE
      p
    }
    ep <- fake_process(); bp <- fake_process()
    fake_job <- function(process) list(
      process = process, progress = NULL,
      progress_file = tempfile(), log_file = tempfile())
    efrm_job(fake_job(ep)); btlef_job(fake_job(bp))
    session$setInputs(project_file = list(
      datapath = path, name = "cj-again.rasch", size = file.info(path)$size,
      type = "application/octet-stream"))
    session$flushReact()
    expect_false(ep$alive)
    expect_false(bp$alive)
    expect_null(efrm_job())
    expect_null(btlef_job())
    expect_identical(bdif_res()$result_signature, bdif1$result_signature)
    expect_identical(bdif_meta()$judge_col, "judge")

    # A worker result is current only while its context, data and CJ base fit
    # are the same as at launch.
    st <- list(context = analysis_context(), data = raw_data(),
               check_btl_fit = TRUE, base_fit = btl_fit())
    expect_true(background_job_is_current(st))
    advance_analysis_context()
    expect_false(background_job_is_current(st))

    st$context <- analysis_context()
    st$data <- raw_data()
    expect_true(background_job_is_current(st))
    btl_fit(bt2)
    expect_false(background_job_is_current(st))

    btl_fit(bt1)
    st$base_fit <- bt1
    sim_data(as.data.frame(d2))
    expect_false(background_job_is_current(st))

    # Both completion observers apply the guard before storing their result.
    completed_process <- function(value) {
      p <- new.env(parent = emptyenv())
      p$is_alive <- function() FALSE
      p$get_result <- function() list(value = value, warnings = character(0))
      p
    }
    common <- list(progress = NULL, progress_file = tempfile(),
                   log_file = tempfile(), context = analysis_context(),
                   data = raw_data())
    ej <- c(list(process = completed_process(bt2), code_call = NULL,
                 src_line = "dat <- old", se_method = "hybrid",
                 simulation_stamp = NULL), common)
    advance_analysis_context()
    efrm_job(ej)
    session$flushReact()
    expect_null(efrm_job())
    expect_identical(btl_fit()$objects$location, bt1$objects$location)

    btlef_res(NULL); clear_btl_analysis_steps()
    bj <- c(list(process = completed_process(list(marker = 99L)),
                 check_btl_fit = TRUE, base_fit = bt1), common)
    btlef_job(bj)
    session$flushReact()
    expect_null(btlef_job())
    expect_null(btlef_res())
    expect_length(btl_analysis_steps(), 0L)
  })
})

test_that("a failed replacement leaves the current app analysis intact", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  set.seed(83)
  X <- matrix(rbinom(500, 1, .5), 100, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  current <- rasch(X)

  shiny::testServer(e$server, {
    fit_val(current)
    push_analysis_step("test", "Existing change", current)
    complete_fit(simpleError("replacement failed"), NULL, character(0),
                 "dat <- existing")
    expect_identical(fit_val(), current)
    expect_length(analysis_steps(), 1L)
  })
})

test_that("the app calculates and restores external person weights", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  current <- rasch(simulate_rasch(80, 5, seed = 921))
  w <- data.frame(item = colnames(current$X),
                  weight = c(2, 2, 1, 1, 0.5))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  write.csv(w, path, row.names = FALSE)

  shiny::testServer(e$server, {
    fit_val(current)
    session$setInputs(
      person_weight_level = "item",
      person_weights_file = list(datapath = path, name = "weights.csv",
        size = file.info(path)$size, type = "text/csv"),
      person_weights_go = 1)
    session$flushReact()
    z <- person_weight_state()
    expect_equal(nrow(z$table), nrow(current$X))
    expect_identical(z$by, "item")
    expect_match(person_weight_code(), "weighted_person_estimates")
  })
})

test_that("the app simulator preserves explanatory metadata and frame calls", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  shiny::testServer(e$server, {
    session$setInputs(
      sim_layout = "rasch_exp", sim_seed = 17,
      sr_persons = 100, sr_items = 6, sr_model = "dichotomous", sr_cats = 4,
      sr_mean = 0, sr_sd = 1, sr_dist = "normal", sr_diff = c(-2, 2),
      sr_over = 0, sr_under = 0, sr_guess = FALSE, sr_2d = FALSE,
      sr_rho = 0.3, sr_dep = FALSE, sr_dif = FALSE, sr_difmag = 1,
      sr_style = FALSE, sr_styletype = "extreme", sr_speeded = 0,
      sr_careless = 0, sr_missing = 0, sx_cont = 1, sx_cat = 0.5,
      sx_interaction = TRUE, sx_int = 0.4, sx_depart = 0.7, sim_go = 1)
    session$flushReact()
    expect_equal(nrow(sim_predictors_val()), 6L)
    expect_identical(sim_interactions_val(), "exposure:type")
    regenerated <- eval(parse(text = sim_code_val()))
    expect_s3_class(regenerated, "rasch_sim")
    expect_equal(nrow(attr(regenerated, "predictors")), 6L)
    expect_identical(attr(regenerated, "truth")$explanatory_formula,
                     "~ exposure + type + exposure:type")
    expect_true(any(grepl("fixed explanatory departure",
                          attr(regenerated, "truth")$planted)))

    bundle <- tempfile(fileext = ".zip")
    unpacked <- tempfile("rasch-simulation-test-")
    dir.create(unpacked)
    on.exit(unlink(c(bundle, unpacked), recursive = TRUE, force = TRUE),
            add = TRUE)
    write_sim_bundle(bundle)
    bundle_files <- utils::unzip(bundle, list = TRUE)$Name
    expect_setequal(bundle_files,
      c("data.csv", "simulation.R", "truth.rds", "predictors.csv",
        "README.txt"))
    utils::unzip(bundle, exdir = unpacked)
    expect_equal(nrow(read.csv(file.path(unpacked, "data.csv"))), 100L)
    expect_identical(readRDS(file.path(unpacked, "truth.rds"))$layout,
                     "rasch")
    expect_equal(nrow(read.csv(file.path(unpacked, "predictors.csv"))), 6L)

    session$setInputs(
      sim_layout = "btl_efrm", sim_seed = 18,
      sbf_objects = 4, sbf_sets = 2, sbf_judges = 4, sbf_panels = 2,
      sbf_within = 5, sbf_cross = 5, sbf_objsd = 1,
      sbf_setratio = 1.3, sbf_panelratio = 1.2, sbf_origin = 0.5,
      sbf_erratic = 0.25, sim_go = 2)
    session$flushReact()
    regenerated <- eval(parse(text = sim_code_val()))
    tr <- attr(regenerated, "truth")
    expect_identical(tr$layout, "btl_efrm")
    expect_length(tr$erratic, 2L)

    # Recovery must use the active frame-adjusted comparison fit, not the
    # equal-unit base fit. A truth-valued stand-in isolates that app routing
    # from the estimator, which is covered by test-simulate.R.
    framed <- structure(list(
      objects = data.frame(object = names(tr$v), v = unname(tr$v)),
      phi_table = data.frame(panel = names(tr$phi), phi = unname(tr$phi)),
      alpha_table = data.frame(set = names(tr$alpha), alpha = unname(tr$alpha)),
      kappa_table = data.frame(set = names(tr$kappa), kappa = unname(tr$kappa))),
      class = c("rasch_btl_efrm", "rasch_btl"))
    btl_fit(framed)
    fitted_sim_gen(sim_gen())
    rr <- sim_recovery_val()
    expect_s3_class(rr, "rasch_recovery")
    expect_setequal(rr$summary$parameter,
      c("object location", "panel unit (log)", "set unit (log)", "set origin"))
  })
})

test_that("the app exposes reproducible WrightMap panel controls", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")
  skip_if_not_installed("WrightMap")
  skip_if_not("item.groups" %in% names(formals(WrightMap::wrightMap)),
              "WrightMap item panels require version 1.5")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))

  d <- simulate_efrm(n_per_group = 55, items_per_set = 4, n_sets = 2,
                     n_groups = 2, n_categories = 2, seed = 73)
  truth <- attr(d, "truth")
  current <- rasch_efrm(d, item_sets = truth$item_sets, groups = "group",
                        id = "id", boot_reps = 0)

  shiny::testServer(e$server, {
    fit_val(current)
    session$flushReact()
    session$setInputs(wright_renderer = "wrightmap",
                      wright_type = "thresholds",
                      wright_person_panels = "groups",
                      wright_item_panels = "sets_groups",
                      wright_person_style = "histogram",
                      tg_bins = 35,
                      tg_axis_mode = "standard")
    session$flushReact()
    expect_match(output$wright_code, "wright_map\\(fit")
    expect_match(output$wright_code, 'person_panels = "groups"', fixed = TRUE)
    expect_match(output$wright_code,
                 'item_panels = c("sets", "groups")', fixed = TRUE)
    expect_type(output$wright, "list")
  })
})

test_that("uploaded frame maps refuse conflicting or unknown assignments", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  duplicate <- tempfile(fileext = ".csv")
  unknown <- tempfile(fileext = ".csv")
  blank <- tempfile(fileext = ".csv")
  spaced <- tempfile(fileext = ".csv")
  on.exit(unlink(c(duplicate, unknown, blank, spaced)), add = TRUE)
  write.csv(data.frame(item = c("I1", "I1"), set = c("A", "B")),
            duplicate, row.names = FALSE)
  write.csv(data.frame(item = c("I1", "wrong"), set = c("A", "B")),
            unknown, row.names = FALSE)
  write.csv(data.frame(item = c("I1", "I2"), set = c("A", " ")),
            blank, row.names = FALSE)
  write.csv(data.frame(item = c(" I1 ", "I2"), set = c(" A ", "B")),
            spaced, row.names = FALSE)

  shiny::testServer(e$server, {
    expect_error(read_frame_map(list(datapath = duplicate), "item",
                                c("I1", "I2"), "item"),
                 "assigns item values more than once")
    expect_error(read_frame_map(list(datapath = unknown), "item",
                                c("I1", "I2"), "item"),
                 "unknown item values")
    expect_error(read_frame_map(list(datapath = blank), "item",
                                c("I1", "I2"), "item"),
                 "non-blank set")
    expect_identical(
      read_frame_map(list(datapath = spaced), "item", c("I1", "I2"),
                     "item"),
      c(I1 = "A", I2 = "B"))
  })
})

test_that("supplied app anchors fail closed", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  bad_rasch <- tempfile(fileext = ".csv")
  bad_btl <- tempfile(fileext = ".csv")
  on.exit(unlink(c(bad_rasch, bad_btl)), add = TRUE)
  write.csv(data.frame(item = "I1", value = 0), bad_rasch,
            row.names = FALSE)
  write.csv(data.frame(object = "A", value = 0), bad_btl,
            row.names = FALSE)

  shiny::testServer(e$server, {
    session$setInputs(anchor_file = list(
      datapath = bad_rasch, name = "bad.csv", size = file.info(bad_rasch)$size,
      type = "text/csv"))
    expect_error(anchors_in(), "needs columns item, k, tau")
    session$setInputs(bt_anchor_file = list(
      datapath = bad_btl, name = "bad.csv", size = file.info(bad_btl)$size,
      type = "text/csv"))
    expect_error(bt_anchors_in(), "needs columns object, location")
  })
})

test_that("app BTL DIF refuses a factor that varies within judge", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  pair <- as.data.frame(t(utils::combn(LETTERS[1:5], 2)),
                        stringsAsFactors = FALSE)
  names(pair) <- c("object_a", "object_b")
  beta <- setNames(seq(-1, 1, length.out = 5), LETTERS[1:5])
  set.seed(902)
  d <- do.call(rbind, lapply(sprintf("J%02d", 1:12), function(j) {
    z <- pair
    z$judge <- j
    p <- plogis(beta[z$object_a] - beta[z$object_b])
    z$winner <- ifelse(runif(nrow(z)) < p, z$object_a, z$object_b)
    z$group <- if (j <= "J06") "G1" else "G2"
    z
  }))
  d$group[d$judge == "J01" & seq_len(nrow(d)) == 2L] <- "changed"
  bt_current <- btl(d, "object_a", "object_b", winner = "winner",
                    judge = "judge")
  path <- tempfile(fileext = ".csv")
  stable_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(path, stable_path)), add = TRUE)
  write.csv(d, path, row.names = FALSE)
  stable <- d
  stable$group[stable$judge == "J01"] <- "G1"
  write.csv(stable, stable_path, row.names = FALSE)

  shiny::testServer(e$server, {
    btl_fit(bt_current)
    session$setInputs(
      file = list(datapath = path, name = "comparisons.csv",
                  size = file.info(path)$size, type = "text/csv"),
      bt_judge = "judge", bdif_factors = "group")
    session$flushReact()
    expect_error(bdif_factor_maps(), "varies within judge")
    expect_match(bdif_code_grp(), "judge_factor <- function", fixed = TRUE)
    expect_silent(parse(text = bdif_code_grp()))
    session$setInputs(file = list(
      datapath = stable_path, name = "stable.csv",
      size = file.info(stable_path)$size, type = "text/csv"))
    session$flushReact()
    maps <- bdif_factor_maps()
    expect_identical(unname(maps$group["J01"]), "G1")
  })
})
