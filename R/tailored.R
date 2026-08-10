# rasch :: tailored analysis for guessing
# ===========================================================================
# The tailored test-of-fit procedure for detecting and correcting
# guessing on multiple-choice proficiency items (Waller 1989 ARRG; Andrich,
# Marais and Humphry 2012; Andrich & Marais 2019 ch. 17). Where a person's
# probability of success is below the chance level (for example 0.25 with
# four options), a correct response carries more guessing than information;
# tailoring converts those responses to missing and re-estimates. The
# four-step procedure compares the initial analysis with the tailored one
# on a common origin: (3) the initial data re-analysed with the mean of the
# anchor items fixed at their tailored values (average anchoring aligns the
# origins), and (4) the initial data with every item anchored at its
# tailored difficulty, re-estimating only the persons. Guessing shows as
# difficult items becoming harder in the tailored analysis.
# ===========================================================================

#' Tailored analysis for guessing
#'
#' Runs the four-step tailored procedure of Andrich, Marais and Humphry
#' (2012) on a dichotomous analysis.
#' Step 1 is the supplied fit. Step 2 (tailored) sets to missing every
#' observed response whose modelled probability of success, at the step-1
#' person and item estimates, is below \code{chance}, and re-estimates
#' items and persons. Step 3 (origin-equated) re-analyses the
#' \emph{original} data with the mean location of the anchor items fixed at
#' their tailored values by average anchoring, so the two calibrations
#' share an origin. Step 4 (all-anchored) fixes every item at its tailored
#' difficulty and re-estimates persons on the original data. Guessing is
#' indicated when difficult items are estimated harder in the tailored
#' analysis than in the origin-equated one; the comparison table and
#' \code{\link{plot_equate}} on the two calibrations show it directly.
#'
#' @param fit A dichotomous fit from \code{\link{rasch}}.
#' @param chance The guessing floor: the probability of success by chance
#'   (1/number of options; default 0.25).
#' @param anchor_items Items whose mean location fixes the common origin in
#'   step 3. The default takes the third of the test (at least two items)
#'   least affected by tailoring -- fewest responses removed, ties broken
#'   towards the easier tailored location -- which are the easy items the
#'   procedure trusts.
#' @param se_method \code{"none"} (default) reports the item shifts
#'   descriptively. \code{"bootstrap"} resamples persons and repeats the
#'   complete four-step procedure, including automatic anchor selection, to
#'   obtain standard errors, percentile intervals, and Holm-adjusted tests.
#' @param boot_reps Person-bootstrap replicates when
#'   \code{se_method = "bootstrap"}; at least 50, default 999. The
#'   sign-count bootstrap p-value has resolution floor \code{2/(boot_reps
#'   + 1)}, so after the Holm adjustment across m items the smallest
#'   achievable adjusted p is \code{2m/(boot_reps + 1)}; a warning fires
#'   when that floor exceeds 0.05 (detection would be impossible).
#' @return A list of class \code{"rasch_tailored"}: \code{tailored},
#'   \code{origin_equated}, and \code{anchored} fits, the comparison
#'   \code{table} (initial, tailored, origin-equated locations, the
#'   tailored-minus-equated \code{shift}; bootstrap uncertainty columns when
#'   requested), the number of
#'   responses removed, the anchor items used, \code{se_method}, and the
#'   number of usable bootstrap replicates \code{boot_reps_used}.
#' @references Waller, M. I. (1989). Modeling guessing behavior: A
#'   comparison of two IRT models. Applied Psychological Measurement, 13,
#'   233-243. Andrich, D., Marais, I. and Humphry, S. (2012). Using a
#'   theorem by Andersen and the dichotomous Rasch model to assess the
#'   presence of random guessing in multiple choice items. Journal of
#'   Educational and Behavioral Statistics, 37, 417-442.
#' @examples
#' set.seed(1); N <- 800
#' d <- seq(-2, 2.5, length.out = 10); th <- rnorm(N)
#' P <- plogis(outer(th, d, "-"))
#' P <- 0.25 + 0.75 * P            # uniform guessing floor
#' X <- matrix(rbinom(N * 10, 1, P), N, 10)
#' colnames(X) <- paste0("I", 1:10)
#' ta <- tailored_analysis(rasch(X), chance = 0.25)
#' ta$table
#' @export
tailored_analysis <- function(fit, chance = 0.25, anchor_items = NULL,
                              se_method = c("none", "bootstrap"),
                              boot_reps = 999L) {
  se_method <- match.arg(se_method)
  if (!inherits(fit, "rasch")) stop("tailored_analysis needs a rasch fit")
  if (inherits(fit, "rasch_efrm") || inherits(fit, "rasch_mfrm"))
    stop("tailored_analysis currently requires an ordinary dichotomous Rasch fit")
  if (max(fit$m) > 1L)
    stop("tailored analysis applies to dichotomous (multiple-choice) items")
  if (length(chance) != 1L || !is.finite(chance) || chance <= 0 || chance >= 1)
    stop("chance must be one finite value strictly between 0 and 1")

  anchor_requested <- anchor_items
  # step 2: remove responses where the model gives success below chance
  P <- fit$moments$E                     # dichotomous: E = P(x = 1)
  Xt <- fit$X
  cut_cells <- !is.na(Xt) & !is.na(P) & P < chance
  if (!any(cut_cells))
    stop("no responses fall below the chance level; nothing to tailor")
  Xt[cut_cells] <- NA
  tailored <- rasch(Xt, model = fit$model, id = fit$person$id,
                    factors = fit$factors, n_groups = fit$n_groups)
  if (!identical(tailored$items$item, fit$items$item))
    stop("tailoring removed an item entirely; lower 'chance' or drop the item first")

  # step 3: common origin via average anchoring on the easy items -- those
  # least affected by tailoring (fewest responses removed, ties broken by
  # the lower tailored location), at least two and about a third of the test
  removed_per_item <- colSums(cut_cells)
  if (is.null(anchor_items)) {
    n_anchor <- max(2L, ceiling(ncol(fit$X) / 3))
    ord <- order(removed_per_item,
                 tailored$items$location[match(colnames(fit$X),
                                               tailored$items$item)])
    anchor_items <- colnames(fit$X)[ord[seq_len(n_anchor)]]
  }
  anchor_items <- intersect(anchor_items, fit$items$item)
  if (length(anchor_items) < 2)
    stop("need at least two anchor items for the common origin; ",
         "nominate easy items via anchor_items")
  ta_loc <- tailored$items$location[match(anchor_items, tailored$items$item)]
  a3 <- data.frame(item = anchor_items, k = NA, tau = ta_loc)
  origin_equated <- rasch(fit$X, model = fit$model, id = fit$person$id,
                          factors = fit$factors, n_groups = fit$n_groups,
                          anchors = a3)

  # step 4: original data, every item fixed at its tailored value, persons
  # free. With no free item parameter there is nothing for pcml() to do, so
  # the fit is assembled directly on the anchored thresholds.
  thr_t <- tailored$thresholds
  thr4 <- thr_t; thr4$se <- 0; thr4$anchored <- TRUE
  est4 <- list(model = fit$model, thr = thr4,
               cov_tau = matrix(0, nrow(thr4), nrow(thr4)),
               loglik = NA_real_, iterations = 0L, converged = TRUE,
               m = fit$m, anchors = thr4, n_parameters = 0L)
  anchored <- .assemble_fit(fit$model, fit$X, est4, fit$person$id,
                            fit$factors, fit$n_groups, NA,
                            c(fit$notes,
                              "all item parameters anchored at their tailored values; persons re-estimated"))

  idx_t <- match(fit$items$item, tailored$items$item)
  idx_o <- match(fit$items$item, origin_equated$items$item)
  shift <- tailored$items$location[idx_t] - origin_equated$items$location[idx_o]
  tab <- data.frame(item = fit$items$item,
                    initial = fit$items$location,
                    tailored = tailored$items$location[idx_t],
                    origin_equated = origin_equated$items$location[idx_o],
                    removed = removed_per_item,
                    shift = shift, se = NA_real_, ci_low = NA_real_,
                    ci_high = NA_real_, p = NA_real_, p_adj = NA_real_,
                    significant = NA)
  boot_used <- NA_integer_
  if (se_method == "bootstrap") {
    if (length(boot_reps) != 1L || !is.finite(boot_reps) || boot_reps < 50L ||
        boot_reps != floor(boot_reps))
      stop("boot_reps must be one whole number of at least 50 for tailored bootstrap inference")
    boot_reps <- as.integer(boot_reps)
    N <- nrow(fit$X); draws <- list()
    for (bb in seq_len(boot_reps)) {
      take <- sample.int(N, N, replace = TRUE)
      Xb <- fit$X[take, , drop = FALSE]
      fb <- tryCatch(suppressWarnings(rasch(
        Xb, model = fit$model, id = seq_len(N),
        factors = if (is.null(fit$factors)) NULL else
          fit$factors[take, , drop = FALSE], n_groups = fit$n_groups)),
        error = function(e) NULL)
      tb <- if (is.null(fb)) NULL else tryCatch(suppressWarnings(
        tailored_analysis(fb, chance = chance,
                          anchor_items = anchor_requested,
                          se_method = "none")), error = function(e) NULL)
      if (!is.null(tb) && setequal(tb$table$item, tab$item))
        draws[[length(draws) + 1L]] <-
          tb$table$shift[match(tab$item, tb$table$item)]
    }
    if (length(draws) < max(30L, ceiling(boot_reps / 2)))
      stop("fewer than half the tailored bootstrap replicates were estimable; ",
           "the tailoring or anchor design is too unstable for inference")
    B <- do.call(rbind, draws); boot_used <- nrow(B)
    tab$se <- apply(B, 2, stats::sd)
    tab$ci_low <- apply(B, 2, stats::quantile, probs = 0.025, names = FALSE)
    tab$ci_high <- apply(B, 2, stats::quantile, probs = 0.975, names = FALSE)
    tab$p <- vapply(seq_len(ncol(B)), function(j) {
      2 * min((sum(B[, j] <= 0) + 1) / (nrow(B) + 1),
              (sum(B[, j] >= 0) + 1) / (nrow(B) + 1), 0.5)
    }, 0)
    tab$p_adj <- stats::p.adjust(tab$p, method = "holm")
    tab$significant <- tab$p_adj < 0.05
    # the sign-count bootstrap p has resolution floor 2/(B+1); after the
    # Holm step the smallest achievable adjusted p is 2m/(B+1). If that
    # floor sits above alpha, NO item could ever be flagged at these
    # settings -- an inert test must say so, not sit quietly
    floor_adj <- 2 * ncol(B) / (nrow(B) + 1)
    if (floor_adj > 0.05)
      warning(sprintf(paste0(
        "with %d usable replicates and %d items the smallest achievable ",
        "Holm-adjusted p is %.3f (> 0.05): no shift can reach significance ",
        "at these settings -- use boot_reps >= %d"),
        nrow(B), ncol(B), floor_adj, ceiling(2 * ncol(B) / 0.05) - 1L),
        call. = FALSE)
  }
  out <- list(tailored = tailored, origin_equated = origin_equated,
              anchored = anchored, table = tab,
              n_removed = sum(cut_cells), chance = chance,
              anchor_items = anchor_items, se_method = se_method,
              boot_reps_used = boot_used)
  class(out) <- "rasch_tailored"
  out
}

#' @export
print.rasch_tailored <- function(x, ...) {
  cat(sprintf("Tailored analysis: %d response(s) below chance %.2f set to missing\n",
              x$n_removed, x$chance))
  cat(sprintf("Origin from average-anchored items: %s\n",
              paste(x$anchor_items, collapse = ", ")))
  tab <- x$table
  num <- vapply(tab, is.numeric, TRUE)
  tab[num] <- lapply(tab[num], round, 3)
  print(tab, row.names = FALSE)
  if (x$se_method == "bootstrap") {
    up <- sum(x$table$significant %in% TRUE & x$table$shift > 0)
    cat(sprintf(paste0("%d item(s) significantly harder after the ",
                       "Holm-adjusted person bootstrap (%d usable draws).\n"),
                up, x$boot_reps_used))
  } else cat("Item shifts are descriptive; use se_method = 'bootstrap' for inference.\n")
  invisible(x)
}
