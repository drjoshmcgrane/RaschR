test_that("saved app analyses make a validated round trip", {
  set.seed(81)
  theta <- rnorm(120)
  X <- sapply(seq(-1, 1, length.out = 5), function(d)
    rbinom(length(theta), 1, plogis(theta - d)))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(X)
  project <- list(
    format = "rasch-shiny-project",
    schema = 1L,
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
    simulation = list(),
    results = list()
  )
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))

  expect_identical(.save_app_project(project, path), path)
  restored <- .read_app_project(path)
  expect_identical(restored$format, "rasch-shiny-project")
  expect_equal(restored$data, project$data)
  expect_equal(restored$base_fit$items, fit$items)
  expect_identical(restored$rcode, project$rcode)
  expect_identical(restored$kept_fit_code, project$kept_fit_code)
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
  future <- list(format = "rasch-shiny-project", schema = 2L,
                 data = as.data.frame(X), base_fit = fit)
  saveRDS(future, bad)
  expect_error(.read_app_project(bad), "unsupported")
})

test_that("opening a project retains results tied to its active fit", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  set.seed(82)
  X <- matrix(rbinom(500, 1, .5), 100, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X)
  project <- list(
    format = "rasch-shiny-project", schema = 1L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), rcode = "fit <- rasch(dat)",
    kept_fits = list(), kept_fit_code = list(), simulation = list(),
    results = list(lr = list(marker = 17L)))
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
  })
})

test_that("a failed replacement leaves the current app analysis intact", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
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
