test_that("reshape roles use literal column names, never numeric positions", {
  d <- data.frame(row_tag = paste0("row", 1:4),
                  person = c("P1", "P1", "P2", "P2"),
                  visit = c("pre", "post", "pre", "post"),
                  Q = c(0, 1, 1, 0))
  names(d)[2:3] <- c("1", "2")
  r <- rack_data(d, person = "1", time = "2", items = "Q")
  expect_identical(r$id, c("P1", "P2"))
  expect_equal(r[["Q@post"]], c(1, 0))
  expect_equal(r[["Q@pre"]], c(0, 1))
  s <- stack_data(d, person = "1", time = "2", items = "Q")
  expect_identical(s$id, d[["1"]])
  expect_identical(as.character(s$time), d[["2"]])
  for (fun in list(rack_data, stack_data)) {
    expect_error(fun(d, person = 1, time = "2", items = "Q"),
                 "`person` must name exactly one column")
    expect_error(fun(d, person = "1", time = 2, items = "Q"),
                 "`time` must name exactly one column")
    expect_error(fun(d, person = "1", time = "1", items = "Q"),
                 "distinct")
  }
})

test_that("the shared role validator refuses numeric-looking selectors", {
  d <- data.frame(a = 1, b = 2, c = 3, d = 4, check.names = FALSE)
  names(d) <- c("1", "1.5", "Inf", "-1")
  for (bad in c(1, 1.5, Inf, -1))
    expect_error(rasch:::.check_reshape_column(d, bad, "role"),
                 "must name exactly one column")
  for (name in names(d))
    expect_identical(rasch:::.check_reshape_column(d, name, "role"), name)
})
