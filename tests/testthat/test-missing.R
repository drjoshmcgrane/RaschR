# Robustness under missing data: available-case traditional statistics that
# reduce exactly to the textbook complete-case forms, and paired-comparison
# semantics where an unmatched winner is missing, never a tie.

test_that("CTT statistics support complete- and available-case summaries", {
  set.seed(1)
  simP <- function(t, tau) {
    x <- 0:length(tau); p <- exp(x * t - c(0, cumsum(tau))); p / sum(p)
  }
  tau <- lapply(seq(-1, 1, length.out = 8), function(d) d + c(-0.7, 0.7))
  th <- rnorm(400)
  X <- sapply(tau, function(tt) vapply(th, function(t)
    sample(0:2, 1, prob = simP(t, tt)), 0L))
  colnames(X) <- sprintf("I%02d", 1:8)
  f0 <- rasch(X); c0 <- ctt_table(f0)
  tot <- rowSums(X)
  thirds <- .class_intervals(tot, rep(FALSE, length(tot)), 3)
  i <- 3
  expect_equal(c0$table$item_total[i], cor(X[, i], tot), tolerance = 1e-10)
  expect_equal(c0$table$item_rest[i], cor(X[, i], tot - X[, i]),
               tolerance = 1e-10)

  # Equal total scores remain together. The discrimination index therefore
  # cannot change when the same response rows are presented in another order.
  ord <- sample.int(nrow(X))
  c1 <- ctt_table(rasch(X[ord, , drop = FALSE]))
  expect_equal(c1$table$di, c0$table$di, tolerance = 1e-12)
  expect_equal(c0$table$di[i],
               (mean(X[thirds == 3, i]) - mean(X[thirds == 1, i])) / 2,
               tolerance = 1e-10)
  Xr <- X[, -i]
  expect_equal(c0$table$alpha_drop[i],
               (8 - 1) / (8 - 2) * (1 - sum(apply(Xr, 2, var)) /
                                      var(rowSums(Xr))), tolerance = 1e-10)

  # 30% missing plus a linked design with ZERO complete cases
  Xm <- X
  Xm[matrix(runif(length(X)) < 0.3, nrow(X))] <- NA
  Xm[1:200, 1:2] <- NA; Xm[201:400, 7:8] <- NA
  cc <- ctt_table(rasch(Xm))
  expect_equal(cc$n, 0L)
  expect_true(all(is.na(cc$table$facility)))

  cm <- ctt_table(rasch(Xm), missing = "available")
  expect_equal(cm$n, 0L)
  expect_true(all(is.finite(cm$table$facility)))
  expect_true(all(is.finite(cm$table$item_total)))
  expect_true(all(cm$table$item_rest > 0.15))
  # the linked booklets leave item pairs with NO respondents in common:
  # their pairwise covariance does not exist, so alpha and alpha-if-deleted
  # are withheld with the reason -- treating the unobservable covariances
  # as zero would fabricate an alpha from invented independence
  expect_true(all(is.na(cm$table$alpha_drop)))
  expect_true(is.na(cm$alpha))
  expect_match(cm$note, "alpha withheld")
  expect_true(is.na(cm$mean))          # honest: no complete responders
  expect_no_error(print(cm))
})

test_that("complete-case CTT alpha is recomputed on the reported sample", {
  set.seed(13)
  X <- matrix(rbinom(300 * 6, 1, 0.5), 300, 6,
              dimnames = list(NULL, paste0("A", 1:6)))
  Xm <- X
  miss <- matrix(FALSE, 300, 6)
  miss[101:300, ] <- matrix(runif(200 * 6) < 0.35, 200, 6)
  Xm[miss] <- NA
  fit <- rasch(Xm)
  ct <- ctt_table(fit)
  Xc <- Xm[complete.cases(Xm), , drop = FALSE]
  C <- cov(Xc)
  manual <- ncol(Xc) / (ncol(Xc) - 1) *
    (1 - sum(diag(C)) / sum(C))
  expect_equal(ct$n, nrow(Xc))
  expect_equal(ct$alpha, manual, tolerance = 1e-12)
})

test_that("an unmatched winner is missing, not a tie", {
  set.seed(4)
  beta <- c(A = -1, B = 0, C = 1)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 80), b = rep(pr[, 2], each = 80))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  # corrupt some winners: blanks and stray strings are MISSING; "tie" is a tie
  d$win[1:10] <- ""
  d$win[11:20] <- "unsure"
  d$win[21:35] <- "tie"
  f <- btl(d, "a", "b", winner = "win", ties = "drop")
  expect_true(any(grepl("treated as missing", f$notes)))
  expect_true(any(grepl("15 tie\\(s\\) dropped", f$notes)))
  expect_equal(f$n_comparisons, nrow(d) - 35)
  # margin path: same semantics
  d$margin <- factor(ifelse(d$win %in% c(d$a, d$b), "much", NA),
                     levels = "much", ordered = TRUE)
  d$margin[21:35] <- NA
  fm <- btl(d, "a", "b", winner = "win", margin = "margin")
  expect_true(any(grepl("treated as missing", fm$notes)))
  expect_true(any(grepl("tie\\(s\\) placed in the middle", fm$notes)))
  expect_equal(fm$m, 2L)   # much-worse / tie / much-better
})
