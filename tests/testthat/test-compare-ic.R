# Composite-likelihood information criteria in compare_fits(): the penalty
# tr(H^-1 J) absorbs the pairwise over-counting, so CL-AIC/CL-BIC (Varin &
# Vidoni 2005; Gao & Song 2010) must select the planted structure where the
# nominal-count criteria could not.

gen_pcm <- function(N, btrue, tau_mat, seed) {
  set.seed(seed)
  L <- length(btrue); m <- ncol(tau_mat); th <- rnorm(N)
  X <- matrix(0L, N, L, dimnames = list(NULL, sprintf("I%02d", 1:L)))
  for (j in 1:L) for (i in 1:N) {
    d <- th[i] - btrue[j] - tau_mat[j, ]
    p <- c(1, exp(cumsum(d)))
    X[i, j] <- sample(0:m, 1, prob = p / sum(p))
  }
  X
}

test_that("CL-ICs select RSM when the rating structure is true", {
  L <- 8; btrue <- seq(-1.2, 1.2, length.out = L)
  tau_rsm <- matrix(rep(c(-0.9, 0, 0.9), each = L), L, 3)
  X <- gen_pcm(500, btrue, tau_rsm, 21)
  cmp <- compare_fits(PCM = rasch(X), RSM = rasch(X, model = "RSM"))
  expect_equal(cmp$label[which.min(cmp$cl_bic)], "RSM")
  expect_equal(cmp$label[which.min(cmp$cl_aic)], "RSM")
  # the effective count sits well above the nominal one: each response
  # enters every pair its item forms
  expect_true(all(cmp$eff_params > cmp$parameters))
  expect_true(all(cmp$eff_params < cmp$parameters * (L - 1)))
})

test_that("CL-ICs select PCM when the threshold spreads truly vary", {
  L <- 8; btrue <- seq(-1.2, 1.2, length.out = L)
  set.seed(5)
  tau_pcm <- t(sapply(1:L, function(j)
    sort(rnorm(3, 0, 1)) * runif(1, 0.4, 1.8)))
  X <- gen_pcm(500, btrue, tau_pcm, 22)
  cmp <- compare_fits(PCM = rasch(X), RSM = rasch(X, model = "RSM"))
  expect_equal(cmp$label[which.min(cmp$cl_bic)], "PCM")
  expect_equal(cmp$label[which.min(cmp$cl_aic)], "PCM")
})

test_that("CL-BIC counts independent persons when response rows repeat", {
  L <- 6L; N <- 260L
  X <- gen_pcm(N, seq(-1, 1, length.out = L),
               matrix(rep(c(-0.7, 0, 0.7), each = L), L, 3), 220)
  ids <- sprintf("P%03d", seq_len(N))
  fit <- rasch(rbind(X, X), id = rep(ids, 2L))
  eligible <- !fit$person$extreme & rowSums(!is.na(fit$X)) >= 2L
  n_person <- length(unique(fit$person$id[eligible]))
  eff <- abs(sum(diag(fit$est$cov_beta %*% fit$est$H_beta)))
  ic <- .cl_ic(fit)
  expect_equal(unname(ic["bic"]),
               -2 * fit$est$loglik + log(n_person) * eff)
  expect_false(isTRUE(all.equal(unname(ic["bic"]),
                                -2 * fit$est$loglik + log(sum(eligible)) * eff)))
})

sim_btl_pos <- function(pos_effect, seed) {
  set.seed(seed)
  K <- 8; beta <- seq(-1.5, 1.5, length.out = K)
  names(beta) <- paste0("O", K:1)
  n <- 1200
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  jd <- sample(sprintf("J%02d", 1:15), n, TRUE)
  win <- rbinom(n, 1, plogis(beta[ia] - beta[ib] + pos_effect))
  data.frame(object_a = names(beta)[ia], object_b = names(beta)[ib],
             winner = ifelse(win == 1, names(beta)[ia], names(beta)[ib]),
             judge = jd)
}

test_that("CL-ICs detect a planted BTL position effect and reject an absent one", {
  d1 <- sim_btl_pos(0.5, 35)
  b0 <- btl(d1, "object_a", "object_b", "winner", judge = "judge")
  b1 <- btl(d1, "object_a", "object_b", "winner", judge = "judge",
            position = TRUE)
  cmp <- compare_fits(plain = b0, position = b1)
  expect_equal(cmp$label[which.min(cmp$cl_bic)], "position")
  expect_true(all(c("judges", "objects", "OSI") %in% names(cmp)))
  expect_true(cmp$same_data[2])
  d0 <- sim_btl_pos(0, 30)
  c0 <- btl(d0, "object_a", "object_b", "winner", judge = "judge")
  c1 <- btl(d0, "object_a", "object_b", "winner", judge = "judge",
            position = TRUE)
  cmp0 <- compare_fits(plain = c0, position = c1)
  expect_equal(cmp0$label[which.min(cmp0$cl_bic)], "plain")
})

test_that("BTL covariance inference and CL-ICs are withheld with few judges", {
  set.seed(74)
  K <- 6; n <- 600; beta <- seq(-1, 1, length.out = K)
  ia <- sample(K, n, TRUE)
  ib <- (ia + sample(K - 1L, n, TRUE) - 1L) %% K + 1L
  win <- rbinom(n, 1, plogis(beta[ia] - beta[ib]))
  d <- data.frame(a = paste0("O", ia), b = paste0("O", ib),
                  winner = paste0("O", ifelse(win == 1L, ia, ib)),
                  judge = sample(paste0("J", 1:6), n, TRUE))
  f <- btl(d, "a", "b", "winner", judge = "judge")
  expect_false(f$cl$inference_available)
  expect_true(all(is.na(f$objects$se)))
  expect_true(is.na(f$osi$PSI))
  cmp <- compare_fits(a = f, b = f)
  expect_true(all(is.na(cmp$cl_aic)) && all(is.na(cmp$cl_bic)))
  expect_match(attr(cmp, "note"), "too few independent clusters")
})

test_that("mixtures are refused; MFRM fits get NA ICs with the reason noted", {
  X <- gen_pcm(150, seq(-1, 1, length.out = 5),
               matrix(rep(c(-0.5, 0.5), each = 5), 5, 2), 9)
  f <- rasch(X)
  d <- sim_btl_pos(0, 11)
  b <- btl(d, "object_a", "object_b", "winner")
  expect_error(compare_fits(f, b), "not a mixture")
  # MFRM: no Godambe matrices on the assembled est -> NA ICs, note says why
  long <- data.frame(person = rep(sprintf("P%03d", 1:120), each = 4),
                     rater = rep(c("R1", "R2"), 240),
                     item = rep(rep(c("A", "B"), each = 2), 120))
  set.seed(13)
  long$score <- rbinom(nrow(long), 2, 0.5)
  mf <- rasch_mfrm(long, person = "person", item = "item", score = "score",
                   facets = "rater")
  cmp <- compare_fits(a = mf, b = mf)
  expect_true(all(is.na(cmp$cl_aic)))
  expect_match(attr(cmp, "note"), "MFRM/EFRM")
})

test_that("compare_fits withholds ICs for an unconverged fit", {
  set.seed(50); N <- 300; L <- 5
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N), d, "-"))), N, L)
  colnames(X) <- paste0("I", 1:L)
  good <- rasch(as.data.frame(X))
  bad <- suppressWarnings(rasch(as.data.frame(X), maxit = 1L))  # non-convergence
  expect_false(isTRUE(bad$est$converged))
  expect_true(all(is.na(bad$thresholds$se)))
  expect_true(all(is.na(bad$items$se)))
  expect_true(all(is.na(bad$items$p_adj)))
  expect_true(all(is.na(bad$person$se)))
  expect_true(is.na(bad$total_chisq_p))
  expect_true(is.na(bad$psi$PSI))
  cmp <- compare_fits(good = good, bad = bad)
  expect_true("converged" %in% names(cmp))
  expect_true(is.na(cmp$cl_aic[cmp$label == "bad"]))
  expect_true(is.na(cmp$loglik[cmp$label == "bad"]))
})

test_that("compare_fits withholds ICs for invalid Godambe ingredients", {
  set.seed(501)
  X <- matrix(rbinom(300 * 5, 1, 0.5), 300, 5)
  colnames(X) <- paste0("I", 1:5)
  fit <- rasch(X)
  expect_true(is.finite(compare_fits(a = fit, b = fit)$cl_aic[1]))

  bad_cov <- fit
  bad_cov$est$cov_beta[1, 2] <- bad_cov$est$cov_beta[1, 2] + 1
  cmp_cov <- compare_fits(a = bad_cov, b = bad_cov)
  expect_true(all(is.na(cmp_cov$eff_params)))
  expect_true(all(is.na(cmp_cov$cl_aic)))
  expect_true(all(is.na(cmp_cov$cl_bic)))

  singular_h <- fit
  singular_h$est$H_beta[,] <- 0
  cmp_h <- compare_fits(a = singular_h, b = singular_h)
  expect_true(all(is.na(cmp_h$eff_params)))
  expect_true(all(is.na(cmp_h$cl_aic)))
  expect_true(all(is.na(cmp_h$cl_bic)))
})

test_that("compare_fits same_data compares actual responses", {
  set.seed(51); N <- 300; L <- 5
  d <- seq(-1.5, 1.5, length.out = L)
  mk <- function(s) { set.seed(s)
    X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N), d, "-"))), N, L)
    colnames(X) <- paste0("I", 1:L); rasch(as.data.frame(X)) }
  f1 <- mk(1); f2 <- mk(2)      # same items/max/N, different responses
  cmp <- compare_fits(a = f1, b = f2)
  expect_false(cmp$same_data[cmp$label == "b"])
  expect_true(is.na(cmp$two_delta_ll[cmp$label == "b"]))
  expect_true(is.na(cmp$cl_aic[cmp$label == "b"]))
  expect_true(is.na(cmp$cl_bic[cmp$label == "b"]))
})

test_that("same_data includes the independent-person allocation", {
  set.seed(510)
  X <- matrix(rbinom(480 * 5, 1, 0.5), 480, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  ids1 <- rep(sprintf("P%03d", 1:240), 2L)
  ids2 <- rep(sprintf("Q%03d", 1:240), 2L)
  same <- compare_fits(a = rasch(X, id = ids1),
                       b = rasch(X, id = ids2))
  expect_true(same$same_data[2L])

  ids2[c(1L, 241L)] <- ids2[c(2L, 242L)]
  changed <- compare_fits(a = rasch(X, id = ids1),
                          b = rasch(X, id = ids2))
  expect_false(changed$same_data[2L])
  expect_true(is.na(changed$cl_aic[2L]))
  expect_true(is.na(changed$cl_bic[2L]))
})

test_that("Rasch same_data is invariant to row, item and ID presentation", {
  set.seed(511)
  X <- matrix(rbinom(360 * 6, 1, 0.5), 360, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  id <- rep(sprintf("P%03d", 1:180), 2L)
  f1 <- rasch(X, id = id)
  ord_r <- sample.int(nrow(X)); ord_c <- sample.int(ncol(X))
  # The response rows, item columns and person labels all change presentation;
  # the multiset of rows within each independent person does not.
  relabel <- setNames(sprintf("Q%03d", sample.int(180)), unique(id))
  f2 <- rasch(X[ord_r, ord_c, drop = FALSE], id = relabel[id[ord_r]])
  same <- compare_fits(original = f1, reordered = f2)
  expect_true(same$same_data[2L])
  expect_true(is.finite(same$cl_aic[2L]))
  expect_equal(same$two_delta_ll[2L], 0, tolerance = 1e-8)
})

test_that("BTL same_data is exact and invariant to row order", {
  set.seed(3)
  pr <- t(combn(LETTERS[1:4], 2))
  d <- data.frame(a = rep(pr[, 1], each = 40),
                  b = rep(pr[, 2], each = 40))
  x <- rbinom(nrow(d), 1, 0.5); x[1:3] <- c(0, 0, 1)
  d$win <- ifelse(x == 1, d$a, d$b)
  f1 <- btl(d, "a", "b", "win")
  frev <- btl(d[rev(seq_len(nrow(d))), ], "a", "b", "win")
  same <- compare_fits(original = f1, reordered = frev)
  expect_true(same$same_data[2])
  expect_equal(same$two_delta_ll[2], 0, tolerance = 1e-8)

  # This outcome change has zero change under the old positional checksum:
  # +1 on rows 1 and 2, -1 on row 3. It must nevertheless be different data.
  d2 <- d; x2 <- x; x2[1:3] <- c(1, 1, 0)
  d2$win <- ifelse(x2 == 1, d2$a, d2$b)
  f2 <- btl(d2, "a", "b", "win")
  different <- compare_fits(original = f1, changed = f2)
  expect_false(different$same_data[2])
  expect_true(is.na(different$two_delta_ll[2]))

  # Equal total comparison count is insufficient when row weights differ.
  dw1 <- d; dw1$count <- 2
  # Move one replicate between rows with opposing outcomes; redistributing
  # counts between duplicate rows is merely another compression.
  dw2 <- dw1; dw2$count[c(1, 3)] <- c(3, 1)
  fw1 <- btl(dw1, "a", "b", "win", count = "count")
  fw2 <- btl(dw2, "a", "b", "win", count = "count")
  weighted <- compare_fits(original = fw1, reweighted = fw2)
  expect_false(weighted$same_data[2])

  # Judge allocation defines the independent clusters even when outcomes
  # and the number of judges are unchanged.
  dj1 <- d; dj1$judge <- rep(sprintf("J%02d", 1:12), length.out = nrow(d))
  dj2 <- dj1; dj2$judge[c(1, 3)] <- rev(dj2$judge[c(1, 3)])
  fj1 <- btl(dj1, "a", "b", "win", judge = "judge")
  fj2 <- btl(dj2, "a", "b", "win", judge = "judge")
  clustered <- compare_fits(original = fj1, reassigned = fj2)
  expect_false(clustered$same_data[2])

  # Judge names themselves do not alter the independent-cluster allocation.
  relabel <- setNames(paste0("K", sample.int(12)), unique(dj1$judge))
  dj3 <- dj1
  dj3$judge <- unname(relabel[dj1$judge])
  fj3 <- btl(dj3, "a", "b", "win", judge = "judge")
  renamed <- compare_fits(original = fj1, renamed = fj3)
  expect_true(renamed$same_data[2])
  expect_true(is.finite(renamed$cl_aic[2]))
})

test_that("BTL same_data follows independent units, not count representation", {
  set.seed(512)
  pr <- t(combn(LETTERS[1:5], 2L))
  d <- data.frame(a = rep(pr[, 1L], each = 30L),
                  b = rep(pr[, 2L], each = 30L))
  y <- rbinom(nrow(d), 1L, 0.5)
  d$win <- ifelse(y == 1L, d$a, d$b)
  expanded <- btl(d, "a", "b", "win")
  counted <- stats::aggregate(rep(1L, nrow(d)),
                              d[c("a", "b", "win")], sum)
  names(counted)[4L] <- "count"
  compressed <- btl(counted, "a", "b", "win", count = "count")
  same <- compare_fits(expanded = expanded, compressed = compressed)
  expect_true(same$same_data[2L])
  expect_equal(same$loglik[2L], same$loglik[1L], tolerance = 1e-10)
  expect_equal(same$cl_aic[2L], same$cl_aic[1L], tolerance = 1e-8)

  judged <- d
  judged$judge <- rep(sprintf("J%02d", 1:15), length.out = nrow(judged))
  judged_fit <- btl(judged, "a", "b", "win", judge = "judge")
  judged_count <- stats::aggregate(
    rep(1L, nrow(judged)), judged[c("a", "b", "win", "judge")], sum)
  names(judged_count)[5L] <- "count"
  judged_compressed <- btl(judged_count, "a", "b", "win",
                           judge = "judge", count = "count")
  judged_same <- compare_fits(expanded = judged_fit,
                              compressed = judged_compressed)
  expect_true(judged_same$same_data[2L])
  expect_equal(judged_same$cl_bic[2L], judged_same$cl_bic[1L],
               tolerance = 1e-8)

  # One half-scored tie is one independent unit with two score
  # contributions.  It is not two independent opposing judgements.
  base <- d[seq_len(120L), ]
  tie <- data.frame(a = base$a, b = base$b, win = base$win,
                    stringsAsFactors = FALSE)
  tie$win[1L] <- "tie"
  tie$count <- 1L
  tie$count[1L] <- 2L
  tie_fit <- btl(tie, "a", "b", "win", count = "count", ties = "half")
  opposite <- tie[-1L, c("a", "b", "win"), drop = FALSE]
  opposite <- rbind(data.frame(a = tie$a[1L], b = tie$b[1L],
                               win = c(tie$a[1L], tie$b[1L])), opposite)
  opposite_fit <- btl(opposite, "a", "b", "win")
  different <- compare_fits(tie = tie_fit, independent = opposite_fit)
  expect_false(different$same_data[2L])
  expect_true(is.na(different$cl_aic[2L]))
})
