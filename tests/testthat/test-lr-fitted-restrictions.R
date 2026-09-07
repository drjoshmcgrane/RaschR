test_that("LR comparisons retain fitted restrictions without replay metadata", {
  set.seed(9671)
  theta <- rnorm(350L)
  tau <- lapply(seq_len(5L), function(i) c(-1, 0, 1) + (i - 3) * 0.2)
  X <- sapply(tau, function(tt) vapply(theta, function(th)
    sample(0:3, 1L, prob = item_moments(th, tt)$P), 0L))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  unrestricted <- rasch(X)
  anchored <- rasch(X, anchors = data.frame(
    item = c("I1", "I2"), k = c(1L, 2L), tau = c(-1, 0.3)))
  components <- rasch(X, pc_components = 2L)
  expect_true(anchored$est$converged)
  expect_true(components$est$converged)
  expect_lt(anchored$est$n_parameters, unrestricted$est$n_parameters)
  expect_lt(components$est$n_parameters, unrestricted$est$n_parameters)

  for (metadata in list(NULL, list())) {
    a <- anchored
    a$refit_spec <- metadata
    expect_error(lr_test(a), "fixed threshold anchors")
    pc <- components
    pc$refit_spec <- metadata
    expect_error(lr_test(pc), "principal-component threshold constraints")
  }
  # The actual threshold flags are enough to identify anchoring even when
  # the separate anchor table is not present in an older saved fit.
  anchored$refit_spec <- NULL
  anchored$est$anchors <- NULL
  expect_true(any(anchored$est$thr$anchored))
  expect_error(lr_test(anchored), "fixed threshold anchors")

  expected <- lr_test(unrestricted)
  unrestricted$refit_spec <- NULL
  actual <- lr_test(unrestricted)
  for (name in c("chisq", "df", "p", "chisq_adj", "p_adj", "lambda"))
    expect_equal(actual[[name]], expected[[name]])
  expect_identical(actual$fit_rsm$X, expected$fit_rsm$X)
})
