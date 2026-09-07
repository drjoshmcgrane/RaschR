.partial_information_fixture <- function(kind) {
  m <- c(1L, 2L, 3L, 1L)
  vm <- data.frame(vkey = paste0("v", 1:4), item = paste0("I", 1:4),
                   set = rep("core", 4), group = rep("G", 4),
                   rater = rep("A", 4))
  structure(list(est = list(converged = TRUE),
    alpha_table = data.frame(set = "core", alpha = 1), m = m,
    tau_list = lapply(m, function(k) rep(0, k)),
    disc = c(0.5, 1, 1.5, 2), virtual_map = vm, facet_spec = "rater",
    items = data.frame(item = vm$vkey),
    X = matrix(c(0, 1, NA, NA, NA, 1, 2, NA, 1, 0, 2, 1,
                 NA, NA, NA, NA, 1, 2, NA, NA), 5, byrow = TRUE,
               dimnames = list(NULL, vm$vkey))),
    class = c(paste0("rasch_", kind), "rasch"))
}

test_that("structural information uses exact observed item patterns", {
  for (kind in c("efrm", "mfrm")) {
    fit <- .partial_information_fixture(kind)
    blocks <- .design_blocks(fit)
    expect_setequal(unname(blocks), list(1:2, 2:3, 1:4))
    expect_length(unique(names(blocks)), 3L)
    expect_error(test_information(fit, grid = 0, items = 1e100), "lie in")
    expect_error(test_information(fit, grid = 0, items = 1 + 1i), "whole numbers")
    ti <- test_information(fit, grid = 0)
    # Equal thresholds at zero give a uniform score on 0:m, with variance
    # m(m+2)/12. This reference does not call the package's moment code.
    contribution <- fit$disc^2 * fit$m * (fit$m + 2) / 12
    expected <- vapply(blocks, function(ii) sum(contribution[ii]), 0)
    expect_equal(ti$info, unname(expected), tolerance = 1e-12)
    expect_equal(ti$sem, unname(1 / sqrt(expected)), tolerance = 1e-12)
    subset <- test_information(fit, grid = 0, items = c("v1", "v2"))
    expect_equal(subset$info, unname(vapply(blocks, function(ii)
      sum(contribution[intersect(ii, 1:2)]), 0)), tolerance = 1e-12)

    # The same correction must reach the plotted expected totals.
    curves <- list()
    grDevices::pdf(NULL)
    tryCatch(with_mocked_bindings(plot_tcc(fit, grid = c(0, 0.1)),
      lines = function(x, y, ...) curves[[length(curves) + 1L]] <<- y,
      .package = "rasch"), finally = grDevices::dev.off())
    expect_equal(vapply(curves, `[`, 0, 1L),
      unname(vapply(blocks, function(ii) sum(fit$m[ii]) / 2, 0)),
      tolerance = 1e-12)
  }
})

test_that("partial administration labels cannot merge distinct designs", {
  for (kind in c("efrm", "mfrm")) {
    fit <- .partial_information_fixture(kind)
    fit$virtual_map$item <- c("a+b", "c", "a", "b+c")
    fit$X <- rbind(c(0, 1, NA, NA), c(NA, NA, 2, 1))
    blocks <- .design_blocks(fit)
    expect_setequal(unname(blocks), list(1:2, 3:4))
    expect_length(unique(names(blocks)), 2L)
    expect_length(unique(test_information(fit, grid = 0)$design), 2L)
  }
})

test_that("EFRM information agrees with fitted partial-pattern score curves", {
  d <- simulate_efrm(n_per_group = 120, items_per_set = 6, n_sets = 1,
                     n_groups = 2, seed = 91731)
  sets <- attr(d, "truth")$item_sets
  for (g in unique(d$group))
    d[head(which(d$group == g), 30), sets[[1]][1:2]] <- NA
  fit <- rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                    boot_reps = 0, workers = 1)
  ti <- test_information(fit, grid = 0)
  sc <- fit$score_curves[abs(fit$score_curves$theta) < 1e-8, ]
  expect_equal(nrow(ti), 4L)
  expect_equal(sort(ti$sem), sort(sc$sem), tolerance = 1e-8)
})

test_that("person-item information excludes designs outside the selected subgroup", {
  for (kind in c("efrm", "mfrm")) {
    fit <- .partial_information_fixture(kind)
    fit$person <- data.frame(theta = seq(-1, 1, length.out = 5))
    fit$factors <- data.frame(subgroup = c("short", "other", "full", NA, "short"),
                              frame = rep("G", 5))
    fit$frame_group <- "frame"
    fit$thresholds <- data.frame(item = rep(1:4, fit$m),
                                 tau = unlist(fit$tau_list))
    original <- fit
    captured <- NULL
    information_fun <- test_information
    grDevices::pdf(NULL)
    tryCatch(with_mocked_bindings(
      plot_pimap(fit, group = "subgroup: short", information = TRUE),
      test_information = function(fit, grid, items = NULL) {
        captured <<- information_fun(fit, grid, items)
        captured
      }, .package = "rasch"), finally = grDevices::dev.off())
    expect_length(unique(captured$design), 1L)
    expected <- vapply(captured$theta, function(th)
      sum(vapply(1:2, function(i) fit$disc[i]^2 *
        item_moments(th, fit$tau_list[[i]], fit$disc[i])$V, 0)), 0)
    expect_equal(captured$info, expected, tolerance = 1e-12)
    expect_identical(fit, original)
  }
})
