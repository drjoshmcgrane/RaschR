# Multiple-choice scoring: double keying and polytomous option scoring
# (Andrich & Styles 2011), with the distractor_rescore proposal.

sim_mc_partial <- function(N = 800, L = 6, seed = 4) {
  # trait-driven 3-category items relabelled as options: A = full credit,
  # B = informative distractor, D = uninformative wrong, (C = rare noise)
  set.seed(seed)
  th <- rnorm(N)
  d0 <- seq(-0.8, 0.8, length.out = L)
  raw <- sapply(d0, function(d) {
    x <- vapply(th, function(b) {
      p <- item_moments(b, c(d - 0.8, d + 0.8))$P
      sample(0:2, 1, prob = p)
    }, 0L)
    o <- c("D", "B", "A")[x + 1]
    swap <- runif(N) < 0.05 & o == "D"
    o[swap] <- "C"
    o
  })
  colnames(raw) <- paste0("M", seq_len(L))
  list(raw = raw, th = th)
}

test_that("double keying credits every listed option", {
  expect_error(
    .resolve_key(c(M1 = "A/a")),
    "each credited option must be named once"
  )
  expect_error(.resolve_key(c(M1 = "A/")), "empty credited option")
  expect_error(.resolve_key(c(M1 = "A//C")), "empty credited option")
  set.seed(9); N <- 400
  th <- rnorm(N)
  raw <- sapply(seq(-1, 1, length.out = 5), function(d) {
    ok <- rbinom(N, 1, plogis(th - d))
    ifelse(ok == 1, sample(c("A", "C"), N, replace = TRUE),
           sample(c("B", "D"), N, replace = TRUE))
  })
  colnames(raw) <- paste0("M", 1:5)
  fit <- rasch(raw, key = setNames(rep("A/C", 5), colnames(raw)))
  expect_true(all(fit$m == 1))
  # both A and C score 1
  expect_equal(unname(fit$X[, 1]),
               as.integer(raw[, 1] %in% c("A", "C")))
  da <- distractor_analysis(fit)
  expect_true(all(da$keyed[da$option %in% c("A", "C")]))
  expect_true(all(!da$keyed[da$option %in% c("B", "D")]))
  expect_equal(unname(fit$mc$key[1]), "A/C")
  rp <- .person_estimates(fit$X[, -1L, drop = FALSE], fit$tau_list[-1L])
  expect_true(any(rp$extreme))
  expected_n <- sum(!is.na(fit$mc$raw[, 1L]) & is.finite(rp$theta) &
                      !rp$extreme)
  expect_equal(sum(da$n[da$item == "M1"]), expected_n)

  # Full-credit options are one scored category. The miskey reference is
  # therefore their pooled takers, not the highest of their separate means.
  ok <- !is.na(fit$mc$raw[, "M1"]) & is.finite(rp$theta) & !rp$extreme
  pooled <- mean(rp$theta[ok & fit$mc$raw[, "M1"] %in% c("A", "C")])
  expected_flag <- with(da[da$item == "M1", ],
    !keyed & n >= 10L & mean_location > pooled)
  expect_identical(da$flag[da$item == "M1"], expected_flag)

  failed <- fit
  failed$est$converged <- FALSE
  expect_error(distractor_analysis(failed), "did not converge")
  expect_error(distractor_rescore(failed), "did not converge")
  expect_error(plot_distractors(failed, "M1"), "did not converge")
})

test_that("multiple-choice scoring honours declared and negative missing codes", {
  set.seed(96)
  raw <- matrix(sample(c("A", "B"), 900, replace = TRUE), 300, 3,
                dimnames = list(NULL, paste0("M", 1:3)))
  raw[c(1, 4), 1] <- c("9", "09")
  raw[2, 2] <- "-2"
  fit <- rasch(raw, key = setNames(rep("A", 3), colnames(raw)),
               na_codes = 9)

  expect_true(all(is.na(fit$mc$raw[cbind(c(1, 4, 2), c(1, 1, 2))])))
  expect_true(all(is.na(fit$X[cbind(c(1, 4, 2), c(1, 1, 2))])))
  expect_false(any(fit$mc$raw == "9", na.rm = TRUE))

  expect_error(rasch(raw, key = setNames(rep("A", 3), colnames(raw)),
                     na_codes = matrix(9)), "plain numeric or character")
  expect_error(rasch(raw, key = setNames(rep("A", 3), colnames(raw)),
                     na_codes = 9.5), "integer score values")
  expect_error(rasch(raw, key = setNames(rep("A", 3), colnames(raw)),
                     na_codes = NA_character_), "without missing values")
})

test_that("polytomous option scoring fits credited distractors as categories", {
  s <- sim_mc_partial()
  os <- do.call(rbind, lapply(colnames(s$raw), function(it)
    data.frame(item = it, option = c("A", "B"), score = c(2, 1))))
  fit <- rasch(s$raw, key = os)
  expect_true(all(fit$m == 2))
  expect_true(any(grepl("polytomous option-score maps", fit$notes,
                        fixed = TRUE)))
  expect_false(any(grepl("scored 0/1 against the key", fit$notes,
                         fixed = TRUE)))
  expect_equal(unname(fit$X[, 3]),
               unname(c(A = 2L, B = 1L, C = 0L, D = 0L)[s$raw[, 3]]))
  # the polytomous scoring recovers the trait better than binary scoring
  bin <- rasch(s$raw, key = setNames(rep("A", 6), colnames(s$raw)))
  expect_gt(cor(fit$person$theta, s$th, use = "complete.obs"),
            cor(bin$person$theta, s$th, use = "complete.obs"))
  # display forms
  expect_equal(unname(fit$mc$key[1]), "A=2, B=1")
  da <- distractor_analysis(fit)
  expect_equal(da$score[da$item == "M1" & da$option == "B"], 1L)
  expect_true(da$keyed[da$item == "M1" & da$option == "A"])
})

test_that("distractor_rescore proposes credit for the informative distractor", {
  s <- sim_mc_partial()
  bin <- rasch(s$raw, key = setNames(rep("A", 6), colnames(s$raw)))
  pr <- distractor_rescore(bin)
  os <- pr$option_scores
  # B credited on (nearly) all items; C and D never; A keeps top score
  b <- os[os$option == "B", "score"]
  expect_gte(mean(b > 0), 5 / 6)
  expect_true(all(os[os$option %in% c("C", "D"), "score"] == 0))
  for (it in unique(os$item)) {
    d <- os[os$item == it, ]
    expect_equal(d$score[d$option == "A"], max(d$score))
  }
  # the proposal covers every observed option and feeds rasch directly
  refit <- rasch(s$raw, key = os)
  expect_s3_class(refit, "rasch")
  expect_true(all(refit$m >= 1))
  expect_true(is.finite(refit$psi$PSI))
  # evidence table carries the separation statistic
  expect_true(all(c("z_sep", "proposed", "se_location") %in% names(pr$evidence)))

  # Reviewing one item still returns a complete key for the original raw
  # dataset; items outside the review retain their fitted scoring.
  one <- distractor_rescore(bin, items = "M1")
  expect_setequal(unique(one$option_scores$item), colnames(s$raw))
  expect_equal(one$option_scores$score[
    one$option_scores$item == "M2" & one$option_scores$option == "A"], 1L)
  expect_no_error(rasch(s$raw, key = one$option_scores))
})

test_that("rescoring uses all other non-keyed options as the baseline", {
  set.seed(190)
  n <- 300L
  ability <- seq(-2.5, 2.5, length.out = n)
  raw <- matrix("B", n, 6L,
                dimnames = list(NULL, paste0("M", seq_len(6L))))
  ord <- order(ability, decreasing = TRUE)
  raw[ord[seq_len(90L)], "M1"] <- "A"
  raw[ord[91:170], "M1"] <- "B"
  raw[ord[171:235], "M1"] <- "C"
  raw[ord[236:300], "M1"] <- "D"
  difficulty <- seq(-1, 1, length.out = 5L)
  for (j in 2:6)
    raw[, j] <- ifelse(runif(n) < plogis(ability - difficulty[j - 1L]),
                         "A", "B")
  fit <- rasch(raw, key = setNames(rep("A", 6L), colnames(raw)))
  proposal <- distractor_rescore(fit, items = "M1", min_n = 75L, z = 1)
  expect_gt(proposal$option_scores$score[
    proposal$option_scores$item == "M1" &
      proposal$option_scores$option == "B"], 0L)
})

test_that("distractor summaries refuse dependent repeated rows", {
  s <- sim_mc_partial(N = 250)
  raw <- rbind(s$raw, s$raw)
  fit <- rasch(raw, id = rep(seq_len(nrow(s$raw)), 2L),
               key = setNames(rep("A", ncol(raw)), colnames(raw)))
  expect_error(distractor_analysis(fit), "one response row per person")
  expect_error(plot_distractors(fit, "M1"), "one response row per person")
  expect_error(distractor_rescore(fit), "one response row per person")
})

test_that("distractor rescoring refuses an unobserved full-credit option", {
  set.seed(94)
  raw <- matrix(sample(c("B", "C"), 800, replace = TRUE), 200, 4,
                dimnames = list(NULL, paste0("M", 1:4)))
  key <- do.call(rbind, lapply(colnames(raw), function(item)
    data.frame(item = item, option = c("A", "B", "C"),
               score = c(2L, 1L, 0L))))
  fit <- rasch(raw, key = key)
  expect_error(distractor_rescore(fit, min_n = 5),
               "no observed full-credit option")
})

test_that("a one-person rest-measure cell withholds its point-biserial", {
  s <- sim_mc_partial(N = 300)
  fit <- rasch(s$raw, key = setNames(rep("A", 6), colnames(s$raw)))
  rp <- .person_estimates(fit$X[, -1L, drop = FALSE], fit$tau_list[-1L])
  ie <- which(rp$extreme)[1L]
  ii <- which(!rp$extreme & is.finite(rp$theta))[1L]
  expect_true(all(is.finite(c(ie, ii))))
  # A routed design can leave one option represented only by an extreme rest
  # score and another by a single usable rest score. Neither point-biserial
  # is defined, and the keyed option cannot anchor a rescoring proposal.
  fit$mc$raw[, "M1"] <- NA_character_
  fit$mc$raw[ie, "M1"] <- "A"
  fit$mc$raw[ii, "M1"] <- "B"
  out <- distractor_analysis(fit, "M1", min_n = 1)
  expect_setequal(out$option, c("A", "B"))
  expect_equal(out$n[out$option == "A"], 0L)
  expect_equal(out$n[out$option == "B"], 1L)
  expect_true(all(is.na(out$point_biserial)))
  expect_error(distractor_rescore(fit, "M1", min_n = 1),
               "no observed full-credit option with")
})

test_that("key validation guards remain informative", {
  raw <- matrix(sample(c("A", "B"), 60, TRUE), 20, 3,
                dimnames = list(NULL, paste0("M", 1:3)))
  expect_error(rasch(raw, key = data.frame(item = "M1", option = "A",
                                           score = -1)),
               "non-negative")
  expect_error(rasch(raw, key = data.frame(item = "M1", option = c("A", "A"),
                                           score = c(1, 2))),
               "duplicate")
  expect_error(rasch(raw, key = data.frame(item = "M1", option = "A",
                                           score = 0)),
               "credits no option")
  expect_error(rasch(raw, key = data.frame(item = "M9", key = "A")),
               "no key item matches")
})

test_that("data-frame keys preserve literal item-column names", {
  set.seed(95)
  raw <- matrix(sample(c("A", "B"), 900, replace = TRUE), 300, 3,
                dimnames = list(NULL, c(" M1 ", "M2", "M3")))
  key <- data.frame(item = colnames(raw), key = rep("A", 3),
                    check.names = FALSE)
  fit <- rasch(raw, key = key)
  expect_true(" M1 " %in% fit$items$item)

  option_key <- data.frame(item = rep(colnames(raw), each = 2L),
                           option = rep(c("A", "B"), 3L),
                           score = rep(c(1L, 0L), 3L),
                           check.names = FALSE)
  fit2 <- rasch(raw, key = option_key)
  expect_true(" M1 " %in% fit2$items$item)

  ordinary <- raw
  colnames(ordinary)[1L] <- "M1"
  padded_key <- key
  padded_key$item[1L] <- " M1 "
  expect_true("M1" %in% rasch(ordinary, key = padded_key)$items$item)
})

test_that("a factor item selector is read by its label in distractor plots", {
  s <- sim_mc_partial(N = 300)
  fit <- rasch(s$raw, key = setNames(rep("A", 6), colnames(s$raw)))
  fit$mc$raw[, "M1"] <- rep(c("A", "B"), length.out = nrow(fit$mc$raw))
  fit$mc$raw[, "M3"] <- rep(c("A", "C"), length.out = nrow(fit$mc$raw))
  labels <- NULL
  testthat::local_mocked_bindings(
    .rr_legend = function(pos, ...) labels <<- list(...)[[1L]],
    .package = "rasch")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  plot_distractors(fit, factor("M3"))
  expect_true(any(grepl("^C", labels)))
  expect_false(any(grepl("^B", labels)))
})

test_that("distractor legend colours cycle with more options than the palette", {
  s <- sim_mc_partial(N = 400)
  fit <- rasch(s$raw, key = setNames(rep("A", 6), colnames(s$raw)))
  fit$mc$raw[, "M3"] <- rep(LETTERS[1:10], length.out = nrow(fit$mc$raw))
  legend_cols <- NULL
  testthat::local_mocked_bindings(
    .rr_legend = function(pos, ...) legend_cols <<- list(...)$col,
    .package = "rasch")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  plot_distractors(fit, "M3")
  expect_length(legend_cols, 10L)
  expect_false(anyNA(legend_cols))
  expect_identical(legend_cols[9L], legend_cols[1L])
})
