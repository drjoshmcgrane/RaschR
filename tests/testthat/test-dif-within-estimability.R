test_that("within effects require an estimable adjusted intercept", {
  set.seed(910)
  pd <- expand.grid(g1 = factor(c("A", "B")),
                    g2 = factor(c("C", "D")), person = seq_len(20))
  Y <- cbind(rnorm(nrow(pd)), rnorm(nrow(pd)) + 1)
  check_design <- function(keep, estimable) {
    dd <- droplevels(pd[keep, , drop = FALSE])
    yy <- Y[keep, , drop = FALSE]
    X <- model.matrix(~ g1 * g2, dd,
                      contrasts.arg = list(g1 = "contr.sum", g2 = "contr.sum"))
    full <- lm.fit(X, drop(yy %*% c(-1, 1) / sqrt(2)))
    reduced <- lm.fit(X[, -1, drop = FALSE],
                      drop(yy %*% c(-1, 1) / sqrt(2)))
    expect_equal(full$rank - reduced$rank, as.integer(estimable))
    z <- rasch:::.dif_within_tests(yy, dd, "occasion", list(occasion = 2L),
                                   "occasion", c("g1", "g2", "g1:g2"))
    main <- z[z$term == "occasion", ]
    expect_equal(nrow(main), 1L)
    if (estimable) {
      f <- (sum(reduced$residuals^2) - sum(full$residuals^2)) /
        (sum(full$residuals^2) / (nrow(dd) - full$rank))
      expect_equal(main$df, 1)
      expect_equal(main$F_value, f)
      expect_equal(main$p, pf(f, 1, nrow(dd) - full$rank, lower.tail = FALSE))
    } else {
      expect_true(all(is.na(main[c("df", "F_value", "p")])))
    }
  }
  check_design(rep(TRUE, nrow(pd)), TRUE)
  check_design(!(pd$g1 == "B" & pd$g2 == "D"), FALSE)
  check_design((pd$g1 == "A") == (pd$g2 == "C"), FALSE)
})

test_that("public factorial DIF retains unavailable within tests in its family", {
  d <- simulate_rasch(160, 6, seed = 9102)
  d$id <- rep(seq_len(80), each = 2)
  d$occasion <- rep(c("pre", "post"), 80)
  d$g1 <- rep(c("A", "B"), each = 80)
  d$g2 <- rep(c("C", "D"), each = 80)
  fit <- rasch(d, id = "id", factors = c("occasion", "g1", "g2"))
  z <- dif_anova(fit, within = "occasion", n_groups = 2, effects = "factorial")
  tab <- as.data.frame(z$terms)
  main <- tab[tab$term == "occasion", ]
  expect_equal(nrow(main), 6L)
  expect_true(all(is.na(main$df)))
  expect_true(all(is.na(main$F_value)))
  expect_true(all(is.na(main$p)))
  expect_true(all(is.na(main$p_adj)))
  family <- !tab$term %in% c("ci", "Residuals")
  usable <- family & is.finite(tab$p)
  expect_true(any(usable))
  expect_equal(tab$p_adj[usable],
               p.adjust(tab$p[usable], "holm", n = sum(family)))
})

test_that("absent crossed within cells do not abort DIF for other items", {
  d <- simulate_rasch(320, 6, seed = 9104)
  d$id <- rep(seq_len(80), each = 4)
  d$time <- rep(c("pre", "pre", "post", "post"), 80)
  d$condition <- rep(c("A", "B", "A", "B"), 80)
  analyse <- function(dat) {
    fit <- rasch(dat, id = "id", factors = c("time", "condition"))
    dif_anova(fit, within = c("time", "condition"),
              n_groups = 2, effects = "factorial")
  }
  control <- analyse(d)
  expect_true(all(is.finite(control$terms$p[
    control$terms$term == "time:condition"])))
  absent <- d$time == "pre" & d$condition == "B"
  d$I01[absent] <- NA
  z <- analyse(d)
  tab <- as.data.frame(z$terms)
  family <- !tab$term %in% c("ci", "Residuals")
  expect_equal(sum(family), 6L * 6L)
  unsupported <- tab$item == "I01" & family
  expect_equal(sum(unsupported), 6L)
  expect_true(all(is.na(tab$p[unsupported])))
  expect_true(all(is.na(tab$p_adj[unsupported])))
  expect_true(all(is.finite(tab$p[
    tab$item != "I01" & tab$term == "time:condition"])))
  usable <- family & is.finite(tab$p)
  expect_equal(tab$p_adj[usable],
               p.adjust(tab$p[usable], "holm", n = sum(family)))
  expect_true(any(grepl("missing a within-subject cell", z$notes)))

  # The same absent combination for every item still leaves the between
  # stratum available; no within test has a complete factorial panel.
  global <- analyse(d[!absent, ])
  w <- !global$terms$term %in% c("ci", "Residuals")
  expect_equal(sum(w), 36L)
  expect_true(all(is.na(global$terms$p[w])))
  expect_true(any(is.finite(global$terms$p[global$terms$term == "ci"])))
})

test_that("within tests do not silently omit a level absent for one item", {
  d <- simulate_rasch(240, 6, seed = 9105)
  d$id <- rep(seq_len(80), each = 3)
  d$occasion <- rep(c("first", "second", "third"), 80)
  d$I01[d$occasion == "third"] <- NA
  fit <- rasch(d, id = "id", factors = "occasion")
  z <- dif_anova(fit, within = "occasion", n_groups = 2)
  main <- z$terms[z$terms$term == "occasion", ]
  expect_true(is.na(main$p[main$item == "I01"]))
  expect_true(all(is.finite(main$p[main$item != "I01"])))
  expect_equal(main$df[main$item != "I01"], rep(2, 5))
})
