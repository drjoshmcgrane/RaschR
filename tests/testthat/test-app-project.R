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
