sim_lltm <- function(n = 600, seed = 4101) {
  set.seed(seed)
  predictors <- data.frame(
    item = paste0("I", 1:10),
    operation = rep(0:1, each = 5),
    format = rep(c("A", "B"), 5),
    stringsAsFactors = FALSE)
  delta <- 0.75 * predictors$operation +
    0.35 * (predictors$format == "B")
  theta <- rnorm(n)
  X <- sapply(delta, function(d) rbinom(n, 1, plogis(theta - d)))
  colnames(X) <- predictors$item
  list(X = X, predictors = predictors, theta = theta, delta = delta)
}

sim_lpcm_explanatory <- function(n = 700, seed = 4102) {
  set.seed(seed)
  predictors <- data.frame(
    item = paste0("P", 1:10),
    format = rep(c("A", "B"), 5),
    stringsAsFactors = FALSE)
  theta <- rnorm(n)
  X <- matrix(0L, n, nrow(predictors))
  for (i in seq_len(ncol(X))) {
    fb <- as.integer(predictors$format[i] == "B")
    tau <- c(-0.7 + 0.25 * fb, 0.65 + 0.55 * fb)
    X[, i] <- vapply(theta, function(b)
      sample.int(3L, 1L, prob = item_moments(b, tau)$P) - 1L, integer(1))
  }
  colnames(X) <- predictors$item
  list(X = X, predictors = predictors)
}

test_that("LLTM recovers item-feature effects and retains Rasch scoring", {
  d <- sim_lltm()
  f <- rasch_explanatory(d$X, d$predictors,
                         ~ operation + format, level = "item")

  expect_s3_class(f, "rasch_explanatory")
  expect_identical(f$explanatory_model, "LLTM")
  expect_true(f$est$converged)
  expect_equal(f$est$n_parameters, 2L)
  expect_lt(abs(f$est$coefficients["operation", "estimate"] - 0.75), .25)
  expect_lt(abs(f$est$coefficients["formatB", "estimate"] - 0.35), .25)
  expect_equal(f$person$theta, f$score_table$theta[f$person$raw + 1L],
               tolerance = 1e-10)

  z <- explanatory_test(f)
  expect_equal(z$df, f$reference_fit$est$n_parameters - f$est$n_parameters)
  expect_true(is.finite(z$p_kent))
})

test_that("LPCM accepts threshold effects and selected interactions", {
  d <- sim_lpcm_explanatory()
  f <- rasch_explanatory(d$X, d$predictors,
    ~ format + threshold + format:threshold, level = "item")

  expect_identical(f$explanatory_model, "LPCM")
  expect_true(f$est$converged)
  expect_true(all(c("formatB", "threshold2", "formatB:threshold2") %in%
                    f$est$coefficients$term))
  expect_true(all(is.finite(f$est$coefficients$se)))
  expect_equal(nrow(f$explanatory$metadata), sum(f$m))
})

test_that("Rasch explanatory predictors may be continuous, categorical or ordinal", {
  set.seed(4104)
  p <- expand.grid(category = c("A", "B"),
                   ordinal = c("low", "moderate", "high"),
                   replicate = 1:3, stringsAsFactors = FALSE)
  p$item <- paste0("T", seq_len(nrow(p)))
  p$continuous <- scale(rnorm(nrow(p)))[, 1]
  p$ordinal <- ordered(p$ordinal,
                       levels = c("low", "moderate", "high"))
  delta <- .35 * p$continuous + .4 * (p$category == "B") +
    .3 * (as.integer(p$ordinal) - 2)
  theta <- rnorm(700)
  X <- sapply(delta, function(d) rbinom(length(theta), 1,
                                         plogis(theta - d)))
  colnames(X) <- p$item
  fit <- rasch_explanatory(
    X, p, ~ continuous + category + ordinal + category:ordinal)

  terms <- fit$est$coefficients$term
  expect_true("continuous" %in% terms)
  expect_true("categoryB" %in% terms)
  ordinal_terms <- terms[grepl("^ordinal", terms) & !grepl(":", terms)]
  expect_equal(length(ordinal_terms), 2L)
  expect_false(any(grepl("ordinal\\.[LQCTest]", terms)))
  expect_equal(sum(grepl("ordinal", terms) & grepl("categoryB", terms)), 2L)
  expect_true(all(abs(fit$est$coefficients[ordinal_terms, "estimate"] - .3) < .2))
  expect_true(is.ordered(fit$explanatory$metadata$ordinal))
})

test_that("threshold metadata and unidentified formulae are checked", {
  d <- sim_lpcm_explanatory(n = 350)
  p <- merge(d$predictors,
             data.frame(threshold = 1:2), by = NULL)
  p <- p[order(match(p$item, d$predictors$item), p$threshold), ]
  p$scoring <- p$threshold * ifelse(p$format == "B", 2, 1)
  f <- rasch_explanatory(d$X, p, ~ format + scoring,
                         level = "threshold")
  expect_s3_class(f, "rasch_explanatory")

  expect_error(rasch_explanatory(d$X, p[-1, ], ~ scoring,
                                 level = "threshold"), "missing")
  q <- d$predictors; q$constant <- 1
  expect_error(rasch_explanatory(d$X, q, ~ constant),
               "no estimable variation")
  q$bad <- seq_len(nrow(q)); q$bad[1] <- Inf
  expect_error(rasch_explanatory(d$X, q, ~ bad), "non-finite")
})

test_that("diagnostics use one Holm family and fixed departures refit all outputs", {
  d <- sim_lltm(n = 500)
  f <- rasch_explanatory(d$X, d$predictors, ~ operation + format)
  dg <- explanatory_diagnostics(f)
  expect_equal(dg$p_adj, p.adjust(dg$p, "holm"))
  expect_equal(sort(unique(dg$component)), "Item location")

  before <- f$person$theta
  g <- relax_explanatory(f, dg$item[1], "location")
  expect_equal(nrow(g$explanatory$relaxations), 1L)
  expect_gt(g$est$n_parameters, f$est$n_parameters)
  expect_false(isTRUE(all.equal(before, g$person$theta)))
  expect_equal(g$residuals,
               (g$X - g$moments$E) / sqrt(g$moments$V),
               tolerance = 1e-10)
})

test_that("item changes preserve and refit explanatory calibrations", {
  d <- sim_lltm(n = 450)
  groups <- data.frame(group = rep(c("A", "B"), length.out = nrow(d$X)))
  f <- rasch_explanatory(d$X, d$predictors, ~ operation + format,
                         factors = groups)

  dropped <- drop_items(f, "I10")
  expect_s3_class(dropped, "rasch_explanatory")
  expect_equal(ncol(dropped$X), 9L)

  split <- split_items(f, "I1", by = "group")
  expect_s3_class(split, "rasch_explanatory")
  expect_equal(ncol(split$X), 11L)
  expect_true(any(split$explanatory$relaxations$component == "Item location"))
  expect_false(isTRUE(all.equal(f$person$theta, split$person$theta)))

  combined <- combine_items(f, list(c("I1", "I2")))
  expect_s3_class(combined, "rasch_explanatory")
  expect_equal(ncol(combined$X), 9L)
  expect_true(nrow(combined$explanatory$relaxations) >= 1L)
  expect_false(isTRUE(all.equal(f$person$theta, combined$person$theta)))

  dependence <- dependence_magnitude(f, dependent = "I2", independent = "I1")
  expect_s3_class(dependence$refit, "rasch_explanatory")
  expect_true(all(grepl("I2\\|I1=", tail(dependence$refit$items$item, 2L))))
  expect_false(isTRUE(all.equal(f$person$theta, dependence$refit$person$theta)))
})

test_that("multiple-choice data remain scored and inspectable after item changes", {
  set.seed(4106)
  n <- 450
  pred <- data.frame(item = paste0("M", 1:10),
                     feature = rep(0:1, each = 5))
  theta <- rnorm(n)
  raw <- sapply(.65 * pred$feature, function(d) {
    correct <- rbinom(n, 1, plogis(theta - d)) == 1L
    ifelse(correct, "A", sample(c("B", "C", "D"), n, replace = TRUE))
  })
  colnames(raw) <- pred$item
  group <- factor(rep(c("A", "B"), length.out = n))
  fit <- rasch_explanatory(
    raw, pred, ~ feature, factors = data.frame(group = group),
    key = stats::setNames(rep("A", ncol(raw)), colnames(raw)))

  split <- split_items(fit, "M1", by = "group")
  expect_type(split$X, "integer")
  expect_true(all(c("M1 (A)", "M1 (B)") %in% colnames(split$mc$raw)))
  expect_true(all(is.na(split$mc$raw[group == "B", "M1 (A)"])))
  expect_true(all(is.na(split$mc$raw[group == "A", "M1 (B)"])))
  expect_s3_class(distractor_analysis(split, "M1 (A)"), "data.frame")

  relaxed <- relax_explanatory(split, "M2", "location")
  expect_identical(relaxed$mc$raw, split$mc$raw)

  dropped <- drop_items(fit, "M10")
  expect_equal(ncol(dropped$mc$raw), 9L)
  combined <- combine_items(fit, list(c("M1", "M2")))
  expect_equal(ncol(combined$mc$raw), 8L)
  expect_true(all(colnames(combined$mc$raw) %in% paste0("M", 3:10)))
})

sim_explanatory_cj <- function(reps = 30, polytomous = FALSE, seed = 4103) {
  set.seed(seed)
  predictors <- data.frame(
    object = LETTERS[1:8],
    domain = rep(0:1, each = 4),
    format = rep(c("A", "B"), 4), stringsAsFactors = FALSE)
  beta <- 0.8 * predictors$domain + 0.35 * (predictors$format == "B")
  names(beta) <- predictors$object
  pairs <- t(combn(predictors$object, 2L))
  d <- data.frame(a = rep(pairs[, 1], each = reps),
                  b = rep(pairs[, 2], each = reps),
                  stringsAsFactors = FALSE)
  if (!polytomous) {
    p <- plogis(beta[d$a] - beta[d$b])
    d$winner <- ifelse(runif(nrow(d)) < p, d$a, d$b)
  } else {
    tau <- c(-0.8, 0.8)
    d$response <- vapply(seq_len(nrow(d)), function(r)
      sample.int(3L, 1L,
        prob = item_moments(beta[d$a[r]] - beta[d$b[r]], tau)$P) - 1L,
      integer(1))
  }
  list(data = d, predictors = predictors)
}

test_that("explanatory CJ recovers object effects and compares with a free fit", {
  d <- sim_explanatory_cj(reps = 35)
  f <- btl_explanatory(d$data, d$predictors, ~ domain + format,
                       object_a = "a", object_b = "b", winner = "winner")
  expect_s3_class(f, "rasch_btl_explanatory")
  expect_true(f$converged)
  expect_lt(abs(f$object_coefficients["domain", "estimate"] - .8), .3)
  expect_lt(abs(f$object_coefficients["formatB", "estimate"] - .35), .3)
  z <- explanatory_test(f)
  expect_equal(z$df, 5L)
  expect_true(is.finite(z$p_kent))

  dg <- explanatory_diagnostics(f)
  expect_equal(dg$p_adj, p.adjust(dg$p, "holm"))
  g <- relax_btl_explanatory(f, dg$object[1])
  expect_s3_class(g, "rasch_btl_explanatory")
  expect_equal(nrow(g$explanatory$relaxations), 1L)
  expect_false(isTRUE(all.equal(f$objects$location, g$objects$location)))
})

test_that("ordered explanatory CJ retains its response-threshold model", {
  d <- sim_explanatory_cj(reps = 45, polytomous = TRUE)
  f <- btl_explanatory(d$data, d$predictors,
                       ~ domain + format + domain:format,
                       object_a = "a", object_b = "b",
                       response = "response", thresholds = "pc")
  expect_s3_class(f, "rasch_btl_explanatory")
  expect_equal(f$m, 2L)
  expect_identical(f$thr_structure, "pc")
  expect_true("domain:formatB" %in% f$object_coefficients$term)
  expect_equal(f$thresholds$tau, f$reference_fit$thresholds$tau,
               tolerance = .5)
})

test_that("comparative judgement predictors retain all three variable types", {
  set.seed(4105)
  p <- expand.grid(category = c("A", "B"),
                   ordinal = c("low", "moderate", "high"),
                   replicate = 1:2, stringsAsFactors = FALSE)
  p$object <- paste0("O", seq_len(nrow(p)))
  p$continuous <- scale(rnorm(nrow(p)))[, 1]
  p$ordinal <- ordered(p$ordinal,
                       levels = c("low", "moderate", "high"))
  beta <- .3 * p$continuous + .45 * (p$category == "B") +
    .25 * (as.integer(p$ordinal) - 2)
  names(beta) <- p$object
  pairs <- t(combn(p$object, 2L))
  d <- data.frame(a = rep(pairs[, 1], each = 18),
                  b = rep(pairs[, 2], each = 18))
  prob <- plogis(beta[d$a] - beta[d$b])
  d$winner <- ifelse(runif(nrow(d)) < prob, d$a, d$b)
  fit <- btl_explanatory(
    d, p, ~ continuous + category + ordinal + category:ordinal,
    "a", "b", winner = "winner")

  terms <- fit$object_coefficients$term
  expect_true(all(c("continuous", "categoryB") %in% terms))
  ordinal_terms <- terms[grepl("^ordinal", terms) & !grepl(":", terms)]
  expect_equal(length(ordinal_terms), 2L)
  expect_false(any(grepl("ordinal\\.[LQCTest]", terms)))
  expect_equal(sum(grepl("ordinal", terms) & grepl("categoryB", terms)), 2L)
  expect_true(all(abs(fit$object_coefficients[ordinal_terms, "estimate"] - .25) < .2))
  expect_true(is.ordered(fit$explanatory$metadata$ordinal))
})

test_that("saturated explanatory designs reproduce free calibrations", {
  d <- sim_lltm(n = 450)
  fr <- rasch(d$X)
  fe <- rasch_explanatory(d$X, d$predictors, ~ item)
  expect_equal(fe$thresholds$tau, fr$thresholds$tau, tolerance = 1e-7)
  expect_equal(explanatory_test(fe)$df, 0L)

  cj <- sim_explanatory_cj(reps = 30)
  br <- btl(cj$data, "a", "b", winner = "winner")
  be <- btl_explanatory(cj$data, cj$predictors, ~ object,
                        "a", "b", winner = "winner")
  expect_equal(be$objects$location, br$objects$location, tolerance = 1e-7)
  expect_equal(explanatory_test(be)$df, 0L)
})
