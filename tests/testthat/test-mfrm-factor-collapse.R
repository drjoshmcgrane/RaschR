test_that("row-level MFRM factor frames retain observed person values", {
  d <- simulate_mfrm(100, 4, 3, n_categories = 3, seed = 9471)
  group <- ifelse(as.integer(factor(d$person)) %% 2L, "A", "B")
  group[!duplicated(d$person)] <- NA_character_
  d$group <- group
  fit <- function(data, factors) rasch_mfrm(
    data, person = "person", item = "item", score = "score",
    facets = "rater", factors = factors)
  named <- fit(d, "group")
  external <- fit(d, data.frame(group = group))

  expect_identical(sum(is.na(named$factors$group)), 0L)
  expect_identical(external$factors, named$factors)
  expect_identical(external$X, named$X)
  expect_equal(external$thresholds, named$thresholds)
  expect_equal(external$person$theta, named$person$theta)
  expect_equal(dif_anova(external, "group")$summary,
               dif_anova(named, "group")$summary)

  # Reversing the response rows makes an observed factor value appear first.
  # Match person and item identities before comparing the reordered fit.
  ii <- rev(seq_len(nrow(d)))
  reversed <- fit(d[ii, ], data.frame(group = group[ii]))
  rows <- match(named$person$id, reversed$person$id)
  cols <- match(colnames(named$X), colnames(reversed$X))
  expect_identical(reversed$factors$group[rows], named$factors$group)
  expect_equal(unname(reversed$X[rows, cols]), unname(named$X))
  expect_equal(reversed$person$theta[rows], named$person$theta,
               tolerance = 1e-9)
  expect_equal(dif_anova(reversed, "group")$summary,
               dif_anova(named, "group")$summary, tolerance = 1e-8)
})

test_that("wide MFRM factor frames collapse independently by column", {
  d <- simulate_mfrm(80, 4, 3, n_categories = 3, seed = 9472)
  wide <- reshape(d, idvar = c("person", "rater"), timevar = "item",
                   direction = "wide")
  items <- names(wide)[startsWith(names(wide), "score.")]
  wide$group <- ifelse(as.integer(factor(wide$person)) %% 2L, "A", "B")
  wide$band <- ifelse(as.integer(factor(wide$person)) %% 3L, "X", "Y")
  wide$group[!duplicated(wide$person)] <- NA_character_
  wide$band[!duplicated(wide$person, fromLast = TRUE)] <- NA_character_
  fit <- function(factors) rasch_mfrm(
    wide, person = "person", items = items, facets = "rater", factors = factors)
  named <- fit(c("group", "band"))
  external <- fit(wide[c("group", "band")])
  expect_identical(external$factors, named$factors)
  expect_identical(external$X, named$X)
  expect_equal(external$thresholds, named$thresholds)
  expect_equal(dif_anova(external, "group")$summary,
               dif_anova(named, "group")$summary)
})

test_that("MFRM factor collapse keeps missing persons and rejects conflicts", {
  d <- simulate_mfrm(60, 3, 2, n_categories = 3, seed = 9473)
  people <- unique(d$person)
  group <- ifelse(as.integer(factor(d$person)) %% 2L, "A", "B")
  group[!duplicated(d$person)] <- NA_character_
  group[d$person == people[1L]] <- NA_character_
  d$group <- group
  fit <- function(factors) rasch_mfrm(
    d, person = "person", item = "item", score = "score",
    facets = "rater", factors = factors)
  named <- fit("group")
  external <- fit(data.frame(group = factor(group, levels = c("A", "B"))))
  expect_identical(as.character(external$factors$group), named$factors$group)
  expect_true(is.factor(external$factors$group))
  expect_identical(which(is.na(external$factors$group)),
                   match(people[1L], external$person$id))

  # With no observation in a factor column, all values remain missing.
  empty <- fit(data.frame(group = rep(NA_character_, nrow(d))))
  expect_true(all(is.na(empty$factors$group)))

  conflict <- group
  ii <- which(d$person == people[2L])
  conflict[ii[1L]] <- if (group[ii[2L]] == "A") "B" else "A"
  expect_error(fit(data.frame(group = conflict)), "varies within person")
  d$group <- conflict
  expect_error(fit("group"), "varies within person")
})
