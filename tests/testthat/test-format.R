
test_that("results tables print without scientific notation", {
  # base R prints a p of 4e-83 as "4.00e-83"; a large sample produces
  # probabilities that small routinely, and they are unreadable that way
  d <- data.frame(item = c("I1", "I2"), n = c(10L, 20L),
                  location = c(1.23456, -0.5), p_adj = c(4e-83, 0.42),
                  flagged = c(TRUE, FALSE))
  out <- capture.output(print(.tag_tables(d)))
  expect_false(any(grepl("e[-+][0-9]", out)))
  expect_true(any(grepl("< 0.001", out, fixed = TRUE)))
  expect_true(any(grepl("1.235", out, fixed = TRUE)))   # fixed decimals
  expect_true(any(grepl("\\b10\\b", out)))              # integers stay integers

  # tagging keeps the frame a data frame and survives subsetting
  td <- .tag_tables(d)
  expect_s3_class(td, "data.frame")
  expect_true(is.data.frame(td))
  expect_s3_class(td[, c("item", "p_adj")], "rasch_table")
  expect_identical(as.data.frame(td)$p_adj, d$p_adj)   # values untouched

  # a table that already carries a class of its own keeps its own printing
  keep <- d
  class(keep) <- c("something_else", "data.frame")
  expect_false(inherits(.tag_tables(keep), "rasch_table"))
})

test_that("observed and expected proportions are not read as probabilities", {
  # obs_p and est_p are category proportions: a category nobody chose has an
  # observed proportion of exactly zero, which "< 0.001" would misreport
  expect_false(.is_pcol("obs_p"))
  expect_false(.is_pcol("est_p"))
  expect_true(.is_pcol("p"))
  expect_true(.is_pcol("p_adj"))
  expect_true(.is_pcol("total_p"))
  # an empty category beside a populated one: the zero must read as zero
  expect_equal(.fmt_df(data.frame(obs_p = c(0, 0.2273)))$obs_p,
               c("0.000", "0.227"))
  expect_equal(.fmt_df(data.frame(p_adj = c(0, 0.2273)))$p_adj,
               c("< 0.001", "0.227"))
})

test_that("a fitted object's tables print through the shared formatter", {
  set.seed(1)
  N <- 1200; K <- 6
  th <- stats::rnorm(N); dl <- seq(-1.2, 1.2, length.out = K)
  X <- vapply(seq_len(K), function(i)
    stats::rbinom(N, 1, stats::plogis(th - dl[i])), numeric(N))
  colnames(X) <- sprintf("I%02d", seq_len(K))
  f <- rasch(data.frame(id = seq_len(N), X), id = "id")
  expect_s3_class(f$items, "rasch_table")
  expect_s3_class(f$item_trait, "rasch_table")
  expect_false(any(grepl("e[-+][0-9]", capture.output(print(f$items)))))
  expect_false(any(grepl("e[-+][0-9]", capture.output(print(f)))))
})

test_that("a saved table carries plain decimals rather than exponents", {
  p <- tempfile(fileext = ".csv")
  on.exit(unlink(p))
  .write_csv_plain(data.frame(item = "I1", p_adj = 2.3e-17, loc = 1.2345678), p)
  txt <- readLines(p)
  expect_false(any(grepl("e[-+][0-9]", txt)))
  # and the precision is kept, not rounded away
  expect_equal(utils::read.csv(p)$p_adj, 2.3e-17)
  expect_equal(utils::read.csv(p)$loc, 1.2345678)
})
