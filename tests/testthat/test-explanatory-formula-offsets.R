test_that("explanatory designs refuse formula offsets without confusing names", {
  p <- data.frame(item = paste0("I", 1:8), x = rep(c(-1, 1), 4),
                  z = c(-1.4, -.8, .2, 1.1, 1.4, -.2, .8, -1.1))
  X <- matrix(rep(0:2, 8), 3, 8, dimnames = list(NULL, p$item))
  q <- p; names(q)[1L] <- "object"
  formulas <- list(~ x + offset(z), ~ x + offset(z) + offset(x),
                   ~ x + offset(0 * z), ~ offset(z))
  for (formula in formulas) {
    expect_error(.explanatory_metadata(p, formula, X),
                  "formula offsets are not supported")
    expect_error(.btl_explanatory_design(q, formula, q$object),
                  "formula offsets are not supported")
  }
  # A column called offset is an ordinary estimated predictor, not an offset
  # term. Unused metadata and arithmetic inside I() must also remain valid.
  p$offset <- p$z; q$offset <- q$z
  r <- .explanatory_metadata(p, ~ x + offset, X)
  b <- .btl_explanatory_design(q, ~ x + offset, q$object)
  expect_identical(colnames(r$B), c("x", "offset"))
  expect_identical(colnames(b$B), c("x", "offset"))
  expect_equal(.explanatory_metadata(p, ~ x + I(z + 3), X)$B,
               .explanatory_metadata(p, ~ x + z, X)$B,
               ignore_attr = TRUE)
  expect_equal(.btl_explanatory_design(q, ~ x + I(z + 3), q$object)$B,
               .btl_explanatory_design(q, ~ x + z, q$object)$B,
               ignore_attr = TRUE)
})

test_that("public Rasch and CJ fits cannot silently drop an offset", {
  set.seed(712)
  p <- data.frame(item = paste0("I", 1:8), x = rep(c(-1, 1), 4),
                  z = c(-1.4, -.8, .2, 1.1, 1.4, -.2, .8, -1.1))
  theta <- rnorm(200)
  for (m in c(1L, 2L)) {
    X <- vapply(.5 * p$x + p$z, function(b) vapply(theta, function(th)
      sample.int(m + 1L, 1L,
                 prob = item_moments(th, b + seq(-.5, .5, length.out = m))$P) - 1L,
      integer(1)), integer(length(theta)))
    colnames(X) <- p$item
    expect_error(rasch_explanatory(X, p, ~ x + offset(z)),
                  "formula offsets are not supported")
    pt <- p[rep(seq_len(nrow(p)), each = m), ]
    pt$threshold <- rep(seq_len(m), nrow(p))
    expect_error(rasch_explanatory(X, pt, ~ x + offset(z), level = "threshold"),
                  "formula offsets are not supported")
  }
  q <- p; names(q)[1L] <- "object"
  pairs <- t(combn(q$object, 2))
  d <- data.frame(a = rep(pairs[, 1L], each = 40),
                  b = rep(pairs[, 2L], each = 40))
  beta <- setNames(.5 * q$x + q$z, q$object)
  eta <- beta[d$a] - beta[d$b]
  d$winner <- ifelse(runif(nrow(d)) < plogis(eta), d$a, d$b)
  d$grade <- vapply(eta, function(v)
    sample.int(3L, 1L, prob = item_moments(v, c(-.5, .5))$P) - 1L, integer(1))
  expect_error(btl_explanatory(d, q, ~ x + offset(z), "a", "b", winner = "winner"),
                "formula offsets are not supported")
  expect_error(btl_explanatory(d, q, ~ x + offset(z), "a", "b", response = "grade"),
                "formula offsets are not supported")
})
