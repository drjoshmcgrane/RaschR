# Planned DIF contrasts: the auto-derived family, logit-scale estimation on
# resolved locations, and person-level residual tests for stacked
# repeated-measures designs (Maxwell & Delaney 2004; Hagquist & Andrich 2017).

test_that("the auto family follows the factor structure and finds planted DIF", {
  set.seed(42); n <- 800
  gender <- rep(c("m", "f"), n / 2)
  age <- sample(c("20", "30", "40", "50"), n, replace = TRUE)
  d <- seq(-1.5, 1.5, length.out = 8)
  th <- rnorm(n)
  shift <- matrix(0, n, 8)
  shift[gender == "f", 3] <- 0.8
  shift[, 5] <- (as.numeric(age) - 35) / 10 * 0.4
  X <- matrix(rbinom(n * 8, 1, plogis(outer(th, d, "-") - shift)), n, 8)
  colnames(X) <- paste0("I", 1:8)
  fit <- rasch(data.frame(X, gender = gender, age = age),
               factors = c("gender", "age"))
  dc <- dif_contrasts(fit, items = c("I3", "I5", "I7"))

  # the derived family: difference, linear + quadratic trends, interaction
  expect_setequal(dc$family$contrast,
                  c("gender: m - f", "age: linear", "age: quadratic",
                    "gender(m - f) x age(linear)"))
  t <- dc$table
  # planted uniform gender DIF on I3, in the planted direction
  g3 <- t[t$item == "I3" & t$contrast == "gender: m - f", ]
  expect_true(g3$significant && g3$practical && g3$estimate < -0.5)
  # planted linear age drift on I5
  a5 <- t[t$item == "I5" & t$contrast == "age: linear", ]
  expect_true(a5$significant && a5$estimate > 0.5)
  # the clean item stays clean
  expect_false(any(t$significant[t$item == "I7"]))
  # every contrast reads as a difference of weighted averages (cells strings
  # are rounded to 2 dp for display, hence the loose tolerance)
  expect_true(all(abs(vapply(strsplit(dc$family$cells, ", "), function(cc)
    sum(abs(as.numeric(sub("^.* ", "", cc)))), 0) - 2) < 0.05))
  # z-statistics equal estimate/se in the independent-rows case
  expect_equal(t$statistic, t$estimate / t$se, tolerance = 1e-10)
})

test_that("stacked designs use person-level scores and detect drift over time", {
  set.seed(9); n <- 400
  d <- seq(-1.5, 1.5, length.out = 8)
  th <- rnorm(n)
  gender <- rep(c("m", "f"), n / 2)
  gen <- function(shift4) {
    s <- matrix(0, n, 8); s[, 4] <- shift4
    matrix(rbinom(n * 8, 1, plogis(outer(th, d, "-") - s)), n, 8)
  }
  X <- rbind(gen(0), gen(0.6)); colnames(X) <- paste0("I", 1:8)
  dat <- data.frame(X, time = rep(c("1", "2"), each = n),
                    gender = rep(gender, 2))
  id <- rep(sprintf("P%03d", 1:n), 2)
  fit <- rasch(dat, factors = c("time", "gender"), id = id)
  dc <- dif_contrasts(fit, items = c("I4", "I6"), id = id)

  # time varies within id, so it is detected as within-subject
  expect_identical(dc$within, "time")
  expect_true(dc$paired)
  t <- dc$table
  w4 <- t[t$item == "I4" & t$contrast == "time: 2 - 1", ]
  # paired t (df = persons - 1), significant, positive drift, sign-aligned
  expect_equal(w4$df, n - 1)
  expect_true(w4$significant && w4$estimate > 0.4 && w4$statistic > 0)
  expect_true(all(is.na(t$se)))
  expect_true(all(is.na(t$lower)) && all(is.na(t$upper)))
  # clean item, null gender effect, null interaction all stay quiet
  expect_false(any(t$significant[t$item == "I6"]))
  expect_false(t$significant[t$item == "I4" & t$contrast == "gender: m - f"])
  expect_false(t$significant[t$item == "I4" &
                             t$contrast == "time(2 - 1) x gender(m - f)"])
  # The resolved point size remains available, but its row-independent
  # calibration covariance is not a repeated-person sampling covariance.
  ds <- dif_size(fit, "I4", by = "time")
  expect_true(any(abs(ds$pairs$difference) > 0.3))
  expect_true(all(is.na(ds$levels$se)))
  expect_true(all(is.na(ds$pairs$se)) && all(is.na(ds$pairs$significant)))
  expect_match(paste(ds$notes, collapse = " "), "sampling SEs")
  # within requires id
  expect_error(dif_contrasts(fit, items = "I4", within = "time"), "id")
})

test_that("custom cell-weight contrasts are accepted and normalised", {
  set.seed(5); n <- 500
  g <- sample(c("a", "b", "c"), n, replace = TRUE)
  d <- seq(-1, 1, length.out = 6)
  sh <- matrix(0, n, 6); sh[g == "c", 2] <- 0.9   # DIF on I2 only
  X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, g = g), factors = "g")
  dc <- dif_contrasts(fit, items = "I2",
                      contrasts = list("c vs rest" = c(a = -1, b = -1, c = 2)))
  r <- dc$table
  expect_equal(nrow(r), 1L)
  expect_true(r$significant && r$estimate > 0.5)
  expect_error(dif_contrasts(fit, items = "I2",
                             contrasts = list(bad = c(x = 1, y = -1))),
               "design cells")
})

test_that("DIF post-hocs give marginal pairs and pure interaction magnitudes", {
  set.seed(31); n <- 1800
  g1 <- factor(rep(c("a", "b", "c"), each = n / 3))
  g2 <- factor(rep(rep(c("x", "y", "z"), each = n / 9), 3))
  d <- seq(-1.3, 1.3, length.out = 6)
  shift <- ifelse(g1 == "c" & g2 == "z", 1.0, 0)
  X <- matrix(rbinom(n * 6, 1,
    plogis(outer(rnorm(n), d, "-") - outer(shift, c(0, 1, 0, 0, 0, 0)))),
    n, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, g1 = g1, g2 = g2),
               factors = c("g1", "g2"))

  main <- dif_posthoc(fit, "I2", "g1")
  expect_s3_class(main, "rasch_dif_posthoc")
  expect_equal(nrow(main$table), choose(3, 2))
  expect_true(all(main$table$p_adj >= main$table$p - 1e-12))

  intr <- dif_posthoc(fit, "I2", "g1:g2")
  expect_equal(nrow(intr$table), choose(3, 2)^2)
  expect_identical(intr$type, "interaction magnitude")

  # The c-a by z-x row is the difference-in-differences of the jointly
  # resolved cell locations, not merely the largest pair of cell means.
  cells <- dif_size(fit, "I2", by = c("g1", "g2"))$levels
  loc <- setNames(cells$location, cells$level)
  manual <- (loc[["c:z"]] - loc[["a:z"]]) -
    (loc[["c:x"]] - loc[["a:x"]])
  row <- intr$table[intr$table$contrast == "c - a x z - x", ]
  expect_equal(row$estimate, manual, tolerance = 1e-8)
  expect_gt(abs(row$estimate), 0.6)
})
