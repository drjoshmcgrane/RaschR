
test_that("drop_items removes items and refits, for rasch and efrm", {
  d <- simulate_rasch(300, 8, seed = 1)
  f <- rasch(d, id = "id")
  f2 <- drop_items(f, "I03")
  expect_equal(nrow(f2$items), nrow(f$items) - 1L)
  expect_false("I03" %in% f2$items$item)
  expect_match(tail(f2$notes, 1), "dropped item")
  expect_s3_class(f2, "rasch")

  e <- simulate_efrm(n_per_group = 250, items_per_set = 8, n_sets = 2,
                     n_groups = 2, set_unit_ratio = 1.3, seed = 7)
  tr <- attr(e, "truth")
  fe <- rasch_efrm(e, item_sets = tr$item_sets, groups = "group", id = "id",
                   boot_reps = 0)
  fe2 <- drop_items(fe, "S1I03")
  expect_s3_class(fe2, "rasch_efrm")
  expect_false("S1I03" %in% names(fe2$set_of))
  expect_equal(length(names(fe2$set_of)), length(names(fe$set_of)) - 1L)
  expect_setequal(unique(fe2$set_of), unique(fe$set_of))   # both sets survive
  expect_equal(nrow(fe2$person), nrow(fe$person))          # persons preserved

  # guards
  expect_error(drop_items(fe, "NOPE"), "not in the fit")
  expect_error(drop_items(fe, names(fe$set_of)[1:8]), "empty set")
  expect_error(drop_items(f, character()), "at least one item")
})
