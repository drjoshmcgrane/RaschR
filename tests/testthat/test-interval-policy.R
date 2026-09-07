interval_fixture <- function(missing = FALSE, tied = FALSE, poly = FALSE,
                             n_groups = NULL) {
  set.seed(if (missing) 744 else if (poly) 740 else 738)
  if (tied) {
    scores <- if (poly) c(0L, 0L, 1L, 1L, 2L, 2L) else
      c(0L, 0L, 0L, 1L, 1L, 1L)
    X <- t(replicate(if (poly) 240L else 200L, sample(scores)))
  } else {
    X <- matrix(rbinom(400L * 8L, 1,
      plogis(outer(rnorm(400), seq(-1, 1, length.out = 8), "-"))), 400, 8)
    if (missing) X[121:400, 7:8] <- NA_integer_
  }
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  rasch(X, n_groups = n_groups)
}

test_that("bootstrap refits repeat automatic and explicit interval policies", {
  for (ng in list(NULL, 4L)) {
    fit <- interval_fixture(missing = TRUE, n_groups = ng)
    expect_identical(.refit_n_groups(fit), ng)
    refit <- .fit_refit(fit$X, fit$model, .refit_n_groups(fit),
                        fit$refit_spec$anchors, fit$m, 60, 1e-8)
    expect_equal(unname(refit$chisq), fit$items$chisq, tolerance = 1e-10)
    if (is.null(ng)) {
      counts <- vapply(fit$ci_item, function(x) length(unique(na.omit(x))), 0L)
      expect_equal(unname(counts), c(rep(6L, 6L), 2L, 2L))
      imposed <- .fit_refit(fit$X, fit$model, fit$n_groups,
                             fit$refit_spec$anchors, fit$m, 60, 1e-8)
      expect_gt(max(abs(imposed$chisq - fit$items$chisq)), 4)
    }
  }
})

test_that("public bootstrap uses the request and rejects superseded results", {
  original <- .fit_refit
  seen <- list()
  local_mocked_bindings(.fit_refit = function(X, model, n_groups, ...) {
    seen[length(seen) + 1L] <<- list(n_groups)
    original(X, model, n_groups, ...)
  })
  for (ng in list(NULL, 4L)) {
    fit <- interval_fixture(missing = TRUE, n_groups = ng)
    seen <- list()
    bs <- suppressWarnings(fit_bootstrap(fit, B = 20L, seed = 745,
                                         workers = 1L))
    expect_length(seen, 20L)
    expect_true(all(vapply(seen, identical, logical(1), ng)))
    expect_equal(bs$B_used, 20)
    expect_identical(bs$algorithm, "loo-maxt-2")
    expect_no_error(.validate_fit_bootstrap(bs, fit))
    old <- unclass(bs)
    old$algorithm <- "loo-maxt-1"
    old$result_signature <- NULL
    old$result_signature <- .fit_boot_md5(old)
    class(old) <- class(bs)
    expect_error(.validate_fit_bootstrap(old, fit), "recompute")
    if (is.null(ng)) {
      project <- .seal_app_project(list(
        format = "rasch-shiny-project", schema = 2L,
        data = as.data.frame(fit$X), model_type = "rasch", base_fit = fit,
        rasch_steps = list(), btl_steps = list(), settings = list(),
        results = list(bootstrap = list(bs = old, B = 20L, seed = 745L,
                                        kind = "rasch"))))
      path <- tempfile(fileext = ".rasch")
      on.exit(unlink(path), add = TRUE)
      saveRDS(project, path)
      expect_warning(restored <- .read_app_project(path), "class-interval")
      expect_null(restored$results$bootstrap)
      expect_equal(restored$base_fit, fit)
      expect_no_error(.validate_app_project(restored))
    }
  }
})

test_that("one realised interval is valid in plots and downstream refits", {
  fit <- interval_fixture(tied = TRUE)
  expect_equal(fit$n_groups, 1L)
  expect_null(.refit_n_groups(fit))
  expect_null(drop_items(fit, "I6")$refit_spec$n_groups)
  legacy <- fit
  legacy$refit_spec <- NULL
  expect_null(.refit_n_groups(legacy))
  expect_no_error(drop_items(legacy, "I6"))
  explicit <- interval_fixture(tied = TRUE, n_groups = 4L)
  expect_equal(explicit$n_groups, 1L)
  expect_identical(.refit_n_groups(explicit), 4L)
  expect_identical(drop_items(explicit, "I6")$refit_spec$n_groups, 4L)
  reference <- .scree_reference(fit, k = 2, reps = 20, seed = 739)
  expect_equal(attr(reference, "n_used"), 20L)
  expect_equal(attr(reference, "n_errors"), 0L)
  # With identical complete-response scores, tailoring removes whole items.
  # Reach that substantive refusal, not an invalid interval-count error.
  expect_error(tailored_analysis(fit, chance = 0.49),
                "tailoring removed an item entirely")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (observed in c(FALSE, TRUE)) {
    expect_no_error(plot_ccc(fit, "I1", observed = observed))
    expect_no_error(plot_threshold_prob(fit, "I1", observed = observed))
  }
  expect_error(plot_ccc(fit, "I1", n_groups = 1), "n_groups")
  expect_error(plot_threshold_prob(fit, "I1", n_groups = 1), "n_groups")
  path <- tempfile(fileext = ".pdf")
  on.exit(unlink(path), add = TRUE)
  expect_no_error(save_item_plots(fit, "ccc", path, items = "I1"))
  expect_gt(file.info(path)$size, 0)

  poly <- interval_fixture(tied = TRUE, poly = TRUE)
  expect_equal(poly$n_groups, 1L)
  comparison <- lr_test(poly)
  direct <- rasch(poly$X, model = "RSM")
  expect_equal(comparison$loglik_rsm, direct$est$loglik)
})

test_that("plot exports retain per-item automatic allocations", {
  fit <- interval_fixture(missing = TRUE)
  seen <- list()
  local_mocked_bindings(plot_icc = function(fit, item, n_groups, ...) {
    seen[length(seen) + 1L] <<- list(n_groups)
    invisible(NULL)
  })
  path <- tempfile(fileext = ".pdf")
  on.exit(unlink(path), add = TRUE)
  save_item_plots(fit, "icc", path)
  expect_length(seen, 8L)
  expect_true(all(vapply(seen, is.null, logical(1))))
})

test_that("explanatory refits and relaxation retain automatic intervals", {
  base <- interval_fixture(missing = TRUE)
  predictors <- data.frame(item = colnames(base$X), x = rep(0:3, 2L))
  fit <- rasch_explanatory(base$X, predictors, ~ x)
  refit <- .explanatory_refit_modified(fit, fit$X)
  expect_null(refit$refit_spec$n_groups)
  expect_equal(refit$items$chisq, fit$items$chisq, tolerance = 1e-8)
  relaxed <- relax_explanatory(fit, "I1")
  expect_null(relaxed$refit_spec$n_groups)
  counts <- vapply(relaxed$ci_item,
                    function(x) length(unique(na.omit(x))), 0L)
  expect_equal(unname(counts), c(rep(6L, 6L), 2L, 2L))
})

test_that("frame and distractor displays accept one realised interval", {
  X <- interval_fixture(tied = TRUE)$X
  frame <- rasch_efrm(X, item_sets = list(all = colnames(X)),
                       groups = rep("g", nrow(X)), boot_reps = 0)
  raw <- matrix(ifelse(X == 1, "A", "B"), nrow(X), dimnames = dimnames(X))
  keyed <- rasch(raw, key = setNames(rep("A", ncol(X)), colnames(X)))
  expect_equal(frame$n_groups, 1L)
  expect_equal(keyed$n_groups, 1L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_icc_frames(frame, "I1"))
  expect_no_error(plot_distractors(keyed, "I1"))
  expect_error(plot_icc_frames(frame, "I1", n_groups = 1), "n_groups")
  expect_error(plot_distractors(keyed, "I1", n_groups = 1), "n_groups")
})
