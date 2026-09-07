test_that("large crossed designs preserve exact factor-cell identity", {
  x <- as.data.frame(matrix("b", 4, 54))
  x[1, ] <- "a"
  x[3, 1] <- "a"
  x[4, 2] <- "a"
  expect_identical(nrow(unique(x)), 4L)
  cells <- .factor_cells(x)
  expect_identical(nlevels(cells), 4L)
  expect_identical(as.integer(cells), c(1L, 4L, 3L, 2L))
  names(x)[1:5] <- c("sep", "collapse", "method", "decreasing", "na.last")
  expect_identical(.factor_cells(x), cells)

  repeated <- rbind(x, x[3, ], x[1, ])
  repeated[6, 2] <- NA_character_
  result <- .factor_cells(repeated)
  expect_identical(as.integer(result[5]), as.integer(result[3]))
  expect_true(is.na(result[6]))
  expect_identical(nlevels(result), 4L)

  # Labels are also kept distinct when pasted values coincide in this path.
  x[[1]] <- c("a:b", "b", "a", "b")
  x[[2]] <- c("c", "b", "b:c", "a")
  x[3, 3:54] <- x[1, 3:54]
  expect_identical(nlevels(.factor_cells(x)), 4L)
  expect_false(anyDuplicated(levels(.factor_cells(x))) > 0L)
})

test_that("factor keys do not interpret column names as paste arguments", {
  x <- data.frame(A = c("a", "b", "a", NA),
                  B = c("g", "h", "g", "i"))
  expected <- .factor_keys(x)
  for (nm in c("collapse", "recycle0", "sep")) {
    names(x)[1] <- nm
    expect_identical(.factor_keys(x), expected)
  }

  fit <- structure(list(
    tau_list = rep(list(0), 4), facet_spec = "collapse",
    virtual_map = data.frame(collapse = c("a", "a", "b", "b")),
    X = matrix(c(1, 0, NA, NA, NA, NA, 1, 0), 2, byrow = TRUE)
  ), class = c("rasch_mfrm", "rasch"))
  blocks <- .design_blocks(fit)
  expect_length(blocks, 2L)
  expect_setequal(unname(blocks), list(1:2, 3:4))
})
