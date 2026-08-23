test_that("WrightMap data are prepared for dichotomous and polytomous fits", {
  set.seed(71)
  n <- 180
  X <- matrix(rbinom(n * 6, 1, 0.5), n, 6)
  colnames(X) <- paste0("D", seq_len(ncol(X)))
  fd <- rasch(X)

  dd <- rasch:::.wright_map_data(fd, "thresholds")
  expect_equal(dim(dd$persons), c(n, 1L))
  expect_equal(dim(dd$items), c(6L, 1L))
  expect_equal(rownames(dd$items), colnames(X))
  expect_equal(unname(drop(dd$items)), fd$thresholds$tau)

  sim_p <- function(th, tau) {
    score <- 0:length(tau)
    z <- exp(score * th - c(0, cumsum(tau)))
    z / sum(z)
  }
  theta <- rnorm(n)
  Xp <- sapply(seq_len(5), function(j)
    vapply(theta, function(th) sample(0:3, 1, prob = sim_p(th, c(-1, 0, 1))),
           integer(1)))
  colnames(Xp) <- paste0("P", seq_len(ncol(Xp)))
  fp <- rasch(Xp)
  dp <- rasch:::.wright_map_data(fp, "thresholds")
  expect_equal(dim(dp$items), c(5L, 3L))
  for (i in seq_len(5)) {
    expect_equal(unname(dp$items[i, ]),
                 fp$thresholds$tau[fp$thresholds$item == i])
  }
  expect_equal(unname(drop(rasch:::.wright_map_data(fp, "locations")$items)),
               fp$items$location)
})

test_that("person and item panels preserve their labels and ordering", {
  set.seed(72)
  n <- 160
  group <- rep(c("Reference", "Focal"), each = n / 2)
  occasion <- rep(c("First", "Second"), length.out = n)
  X <- matrix(rbinom(n * 6, 1, 0.5), n, 6)
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(data.frame(X, group, occasion),
               factors = c("group", "occasion"))

  expect_equal(ncol(rasch:::.wright_map_data(fit)$persons), 1L)
  expect_null(rasch:::.wright_map_data(fit)$item_panels)

  d1 <- rasch:::.wright_map_data(
    fit, person_panels = "group",
    item_panels = setNames(c("A", "A", "B", "B", "C", "C"),
                           rev(colnames(X)))
  )
  expect_equal(colnames(d1$persons), c("Reference", "Focal"))
  expect_equal(unname(colSums(is.finite(d1$persons))), c(80, 80))
  expect_equal(as.character(d1$item_panels),
               c("C", "C", "B", "B", "A", "A"))

  dl <- rasch:::.wright_map_data(
    fit, item_panels = list(First = c("I1", "I3", "I5"),
                            Second = c("I2", "I4", "I6")))
  expect_equal(as.character(dl$item_panels),
               rep(c("First", "Second"), 3))

  d2 <- rasch:::.wright_map_data(fit, person_panels = c("group", "occasion"))
  expect_equal(ncol(d2$persons), 4L)
  expect_setequal(colnames(d2$persons),
                  c("Reference x First", "Focal x First",
                    "Reference x Second", "Focal x Second"))

  supplied <- cbind(Cohort1 = fit$person$theta,
                    Cohort2 = fit$person$theta + 0.2)
  d3 <- rasch:::.wright_map_data(fit, person_panels = supplied)
  expect_equal(d3$persons, supplied)
})

test_that("EFRM maps recover person groups and item sets from the fitted design", {
  set.seed(73)
  d <- simulate_efrm(n_per_group = 55, items_per_set = 4, n_sets = 2,
                     n_groups = 2, n_categories = 2, seed = 73)
  truth <- attr(d, "truth")
  fit <- rasch_efrm(d, item_sets = truth$item_sets, groups = "group", id = "id",
                    boot_reps = 0)

  ds <- rasch:::.wright_map_data(
    fit, person_panels = "groups", item_panels = "sets"
  )
  expect_equal(colnames(ds$persons), fit$phi_table$group)
  expect_equal(table(ds$item_panels), table(fit$virtual_map$set))

  df <- rasch:::.wright_map_data(fit, item_panels = c("sets", "groups"))
  expect_equal(nlevels(df$item_panels), nrow(fit$frames))
  expect_equal(length(df$item_panels), nrow(fit$virtual_map))
})

test_that("wright_map calls the installed WrightMap interface", {
  skip_if_not_installed("WrightMap")
  set.seed(74)
  X <- matrix(rbinom(160 * 6, 1, 0.5), 160, 6)
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(X)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  out <- wright_map(fit, main.title = NULL)
  expect_equal(out, rasch:::.wright_map_data(fit)$items)

  if ("item.groups" %in% names(formals(WrightMap::wrightMap))) {
    expect_no_error(wright_map(
      fit,
      person_panels = rep(c("Sample A", "Sample B"), each = 80),
      item_panels = rep(c("Domain 1", "Domain 2"), each = 3),
      main.title = NULL
    ))
  }
})

test_that("wright_map rejects inappropriate objects and malformed panels", {
  set.seed(75)
  X <- matrix(rbinom(100 * 5, 1, 0.5), 100, 5)
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(X)

  expect_error(rasch:::.wright_map_data(structure(list(), class = "rasch_btl")),
               "not judge locations")
  expect_error(rasch:::.wright_map_data(fit, person_panels = 1:4),
               "one value per person")
  expect_error(rasch:::.wright_map_data(fit, item_panels = 1:4),
               "one value per item")
  expect_error(rasch:::.wright_map_data(
    fit, item_panels = setNames(1:5, paste0("bad", 1:5))),
    "each item name")
})
