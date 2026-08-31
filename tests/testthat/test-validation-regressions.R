test_that("numeric NaN is refused as corrupt response data", {
  set.seed(1201)
  x <- matrix(rbinom(120, 1, 0.5), 20, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  x[2, 3] <- NaN
  expect_error(rasch(x), "non-finite score.*I3")
})

test_that("structured simulator arguments fail at the public boundary", {
  expect_error(simulate_rasch(second_dim = 0),
               "second_dim must be NULL or a named list")
  expect_error(simulate_rasch(dependence = list(strength = 1)),
               "dependence must contain pairs")
  expect_error(simulate_rasch(dif = list(items = "I01", typo = 1)),
               "dif has unknown component.*typo")
  expect_error(simulate_rasch(response_style = list(1, 2)),
               "response_style must be a named list")
  expect_warning(
    x <- simulate_rasch(model = "dichotomous",
                        response_style = list(type = "extreme"), seed = 1),
    "polytomous items only")
  expect_false(any(grepl("response style", attr(x, "truth")$planted)))

  expect_error(simulate_btl(second_attribute = 0),
               "second_attribute must be NULL or a named list")
  expect_error(simulate_btl(dependence = list(unknown = 1)),
               "dependence has unknown component.*unknown")
  expect_error(simulate_mfrm(interaction = list(rater = "R1", item = "I1")),
               "interaction must contain bias")
})

test_that("the BTL second attribute has its requested finite correlation", {
  for (k in c(3L, 5L, 8L, 15L)) for (rho in c(-0.75, 0, 0.6, 1)) {
    reps <- max(2L, ceiling(8 / choose(k, 2)))
    d <- simulate_btl(n_objects = k, n_judges = 8, reps_per_pair = reps,
                      second_attribute = list(rho = rho), seed = 1300 + k)
    tr <- attr(d, "truth")
    expect_equal(unname(stats::cor(tr$location, tr$location2)), rho,
                 tolerance = 1e-12)
    expect_equal(tr$attribute_correlation, rho, tolerance = 1e-12)
  }
  expect_error(simulate_btl(object_sd = 0,
                            second_attribute = list(rho = 0.5)),
               "positive spread")
})

test_that("the Rasch second dimension has its requested finite correlation", {
  for (n in c(3L, 5L, 20L, 500L)) for (rho in c(-0.75, 0, 0.5, 1)) {
    d <- simulate_rasch(n_persons = n, n_items = 4,
                        second_dim = list(items = "I01", rho = rho),
                        seed = 1350 + n)
    tr <- attr(d, "truth")
    expect_equal(unname(stats::cor(tr$theta, tr$theta2)), rho,
                 tolerance = 1e-12)
    expect_equal(mean(tr$theta2), 0, tolerance = 1e-12)
    expect_equal(stats::sd(tr$theta2), 1, tolerance = 1e-12)
  }
  expect_error(simulate_rasch(n_persons = 2, n_items = 3,
                              second_dim = list(items = "I01", rho = 0.5)),
               "at least three persons")
  for (dist in c("uniform", "skew", "bimodal")) {
    d <- simulate_rasch(n_persons = 50, n_items = 4, theta_dist = dist,
                        second_dim = list(items = "I01", rho = 0.35),
                        seed = 1390)
    tr <- attr(d, "truth")
    expect_equal(unname(stats::cor(tr$theta, tr$theta2)), 0.35,
                 tolerance = 1e-12)
    expect_equal(stats::sd(tr$theta2), 1, tolerance = 1e-12)
  }
})

test_that("BTL dependence decisions use a single Holm family", {
  d <- simulate_btl(n_objects = 8, n_judges = 30, reps_per_pair = 12,
                    dependence = list(exposure = 0.25, carry_over = 0.35),
                    seed = 1401)
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
           order = "order")
  expect_true(all(c("p", "p_adj", "significant") %in%
                    names(f$dependence)))
  ok <- is.finite(f$dependence$p)
  expect_equal(f$dependence$p_adj[ok],
               stats::p.adjust(f$dependence$p[ok], method = "holm"))
  expect_identical(f$dependence$significant[ok],
                   f$dependence$p_adj[ok] < 0.05)
  expect_match(paste(capture.output(print(f)), collapse = "\n"), "Holm p")
})

test_that("a one-judge named panel map is not mistaken for a data column", {
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 1, n_panels = 1,
                         n_judges_per_panel = 1, reps_within = 80, seed = 1501)
  panel_map <- stats::setNames("mapped-panel", unique(d$judge))
  f <- btl_efrm(d, "object_a", "object_b", winner = "winner",
                judge = "judge", panels = panel_map,
                object_sets = attr(d, "truth")$object_sets,
                se_method = "conditional")
  expect_identical(f$phi_table$panel, "mapped-panel")
})

test_that("split_items omits levels without an estimable item", {
  set.seed(1601)
  n <- 300
  group <- factor(rep(c("A", "B", "C"), each = n / 3))
  theta <- rnorm(n)
  x <- sapply(seq(-1.2, 1.2, length.out = 6), function(b)
    rbinom(n, 1, plogis(theta - b)))
  colnames(x) <- paste0("I", 1:6)
  x[group == "C", "I1"] <- NA
  f <- rasch(data.frame(x, group), factors = "group")
  s <- split_items(f, "I1", by = "group")
  expect_setequal(grep("^I1 \\(", s$items$item, value = TRUE),
                  c("I1 (A)", "I1 (B)"))
  expect_false("I1 (C)" %in% s$items$item)
  expect_true(any(grepl("unavailable level.*C", s$notes)))

  x[group != "A", "I2"] <- NA
  f2 <- rasch(data.frame(x, group), factors = "group")
  expect_error(split_items(f2, "I2", by = "group"),
               "fewer than two factor levels")
})

test_that("EFRM downstream refits avoid item-name collisions", {
  d <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                     n_groups = 2, seed = 13)
  tr <- attr(d, "truth")
  old <- tr$item_sets[[1]][1]
  names(d)[names(d) == old] <- ".rasch_id"
  tr$item_sets[[1]][tr$item_sets[[1]] == old] <- ".rasch_id"
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 0)
  out <- drop_items(f, tr$item_sets[[1]][2], boot_reps = 0)
  expect_s3_class(out, "rasch_efrm")
  expect_true(".rasch_id" %in% names(out$set_of))
})

test_that("wide MFRM conversion avoids generated-name collisions", {
  d <- simulate_mfrm(50, 3, 4, seed = 1801)
  wide <- reshape(as.data.frame(d), idvar = c("person", "rater"),
                  timevar = "item", direction = "wide")
  names(wide)[names(wide) == "score.I1"] <- "..person"
  names(wide)[names(wide) == "score.I2"] <- "Q2"
  names(wide)[names(wide) == "score.I3"] <- "Q3"
  factors <- data.frame("..item" = factor(rep(c("a", "b"),
                                               length.out = nrow(wide))),
                        check.names = FALSE)
  f <- rasch_mfrm(wide, person = "person",
                  items = c("..person", "Q2", "Q3"), facets = "rater",
                  factors = factors)
  expect_s3_class(f, "rasch_mfrm")
  expect_identical(names(f$factors), "..item")
})

test_that("EFRM reports nuisance-mass non-convergence and bootstrap accounting", {
  # The documented contract is at least 30 usable draws and a majority of
  # those requested.  At the permitted minimum, all 30 therefore have to
  # succeed.
  expect_identical(rasch:::.rasch_min_boot_success(30L), 30L)
  expect_identical(rasch:::.rasch_min_boot_success(50L), 30L)
  expect_identical(rasch:::.rasch_min_boot_success(60L), 31L)
  expect_identical(rasch:::.rasch_min_boot_success(100L), 51L)
  expect_identical(rasch:::.rasch_min_boot_success(200L), 101L)
  expect_identical(rasch:::.rasch_min_boot_success(300L), 151L)
  expect_identical(rasch:::.efrm_min_boot_success(30L), 30L)
  expect_identical(rasch:::.efrm_min_boot_success(50L), 30L)
  expect_identical(rasch:::.efrm_min_boot_success(60L), 31L)
  expect_identical(rasch:::.efrm_min_boot_success(100L), 51L)
  expect_identical(rasch:::.efrm_min_boot_success(200L), 101L)
  expect_identical(rasch:::.efrm_min_boot_success(300L), 151L)

  d <- simulate_efrm(n_per_group = 120, items_per_set = 6, n_sets = 2,
                     n_groups = 1, set_unit_ratio = 1.25, seed = 1901)
  tr <- attr(d, "truth")
  old_mass <- rasch:::efrm_fit_weights_cpp
  testthat::local_mocked_bindings(
    efrm_fit_weights_cpp = function(...) {
      z <- old_mass(...)
      z$converged <- FALSE
      z
    },
    .package = "rasch")
  expect_warning(
    f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                    boot_reps = 0, workers = 1),
    "nuisance masses")
  expect_true(f$est$stage1_converged)
  expect_false(f$est$converged)
  expect_true(any(f$linking$alpha_edges$converged %in% FALSE))
  expect_identical(f$boot_reps_requested, 0L)
  expect_identical(f$boot_reps_used, 0L)
  expect_identical(f$boot_reps_failed, 0L)
})

test_that("a failed EFRM full bootstrap retains its accounting", {
  d <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 1,
                     n_groups = 1, seed = 1904)
  tr <- attr(d, "truth")
  testthat::local_mocked_bindings(
    .efrm_boot_apply = function(n, fun, ...) rep(list(NULL), n),
    .package = "rasch")
  expect_warning(
    f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                    se_method = "bootstrap", boot_reps = 30, workers = 1),
    "0 of 30 replicates were usable; at least 30 are required")
  expect_identical(f$se_method, "hybrid")
  expect_identical(f$full_boot_reps_requested, 30L)
  expect_identical(f$full_boot_reps_attempted, 30L)
  expect_identical(f$full_boot_reps_used, 0L)
  expect_identical(f$full_boot_reps_failed, 30L)
  expect_true(any(grepl("hybrid standard errors were returned", f$notes,
                        fixed = TRUE)))
})

test_that("fixed-iteration NPML skips convergence checks without changing EM", {
  set.seed(1905)
  L <- matrix(rnorm(80), 8, 10)
  logw <- matrix(-log(10), 2, 10)
  mix_idx <- rep(1:2, each = 4)
  count <- sample(1:4, 8, replace = TRUE)
  rr <- rasch:::.efrm_npml_fit_weights_r(
    L, logw, mix_idx, count, maxit = 4, tol = 0)
  cc <- rasch:::efrm_fit_weights_cpp(
    L, logw, mix_idx, count, maxit = 4, tol = 0)
  expect_false(rr$converged)
  expect_false(cc$converged)
  expect_identical(rr$iterations, 4L)
  expect_identical(cc$iterations, 4L)
  expect_true(is.finite(rr$loglik))
  expect_true(is.finite(cc$loglik))
  expect_equal(rr$logw, cc$logw, tolerance = 1e-12)
  expect_equal(rr$loglik, cc$loglik, tolerance = 1e-12)
})

test_that("resolve_dif uses the adjusted omnibus rather than pairwise rejection", {
  set.seed(1902)
  n <- 800
  g <- factor(rep(c("a", "b"), each = n / 2))
  difficulty <- seq(-2, 2, length.out = 10)
  shift <- matrix(0, n, 10)
  shift[g == "b", 3] <- 1.4
  X <- matrix(rbinom(n * 10, 1,
                     plogis(outer(rnorm(n), difficulty, "-") - shift)),
              n, 10)
  colnames(X) <- sprintf("I%02d", seq_len(10))
  f <- rasch(data.frame(X, group = g), factors = "group")

  # Deliberately make every follow-up non-significant. The already adjusted
  # omnibus remains the decision; the follow-up supplies only the magnitude.
  testthat::local_mocked_bindings(
    dif_posthoc = function(...) list(table = data.frame(
      estimate = 1.2, se = 0.3, statistic = 4, p = 0.5, p_adj = 1)),
    .package = "rasch")
  rr <- resolve_dif(f, max_splits = 1)
  expect_identical(rr$n_splits, 1L)
  expect_identical(rr$splits$item, "I03")
})

test_that("BTL-DIF contrast df do not pool multi-cell effective counts", {
  # The validated ordinary two-level rule is retained.
  expect_equal(rasch:::.btl_dif_contrast_df(c(9.31, 9.31)), 16.62)
  # A four-cell difference-in-differences is referenced to its weakest cell,
  # not to the invalid pooled-count extension of the two-cell formula.
  expect_equal(rasch:::.btl_dif_contrast_df(c(12, 10, 9, 8.5)), 7.5)
})

test_that("an unidentified BTL-EFRM set unit is not counted as estimated", {
  set.seed(1903)
  set1 <- paste0("A", 1:4)
  set2 <- paste0("B", 1:4)
  pr1 <- t(utils::combn(set1, 2))
  rank1 <- stats::setNames(4:1, set1)
  within1 <- do.call(rbind, lapply(seq_len(nrow(pr1)), function(i) {
    a <- pr1[i, 1]; b <- pr1[i, 2]
    p <- plogis(1.5 * (rank1[[a]] - rank1[[b]]))
    data.frame(a = a, b = b,
               winner = ifelse(runif(40) < p, a, b))
  }))
  pr2 <- t(utils::combn(set2, 2))
  within2 <- do.call(rbind, lapply(seq_len(nrow(pr2)), function(i) {
    a <- pr2[i, 1]; b <- pr2[i, 2]
    data.frame(a = a, b = b,
               winner = ifelse(seq_len(40) %% 2L == 0L, a, b))
  }))
  cross <- expand.grid(a = set1, b = set2, stringsAsFactors = FALSE)
  cross <- cross[rep(seq_len(nrow(cross)), each = 20), ]
  cross$winner <- ifelse(runif(nrow(cross)) < 0.7, cross$a, cross$b)
  d <- rbind(within1, within2, cross)
  d$judge <- sample(sprintf("J%02d", 1:24), nrow(d), replace = TRUE)
  panel <- stats::setNames(rep("panel1", 24), sprintf("J%02d", 1:24))

  f <- btl_efrm(d, "a", "b", "winner", "judge", panel,
                list(set1 = set1, set2 = set2),
                se_method = "conditional")
  expect_true(is.na(f$alpha_table$alpha[f$alpha_table$set == "set2"]))
  expect_identical(f$n_parameters, 7L)
  expect_identical(f$equal_unit$parameters_frames, 7L)
  expect_identical(f$total_df, 21L)
})

test_that("frame bootstrap keeps a singleton stratum in its own group", {
  set.seed(20260828)
  ii <- rasch:::.frame_stratified_resample(
    list(main = 1:10, sparse = 17L, second = 21:25))

  expect_length(ii, 16L)
  expect_identical(ii[11L], 17L)
  expect_true(all(ii[1:10] %in% 1:10))
  expect_true(all(ii[12:16] %in% 21:25))
})
