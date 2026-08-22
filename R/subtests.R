# rasch :: structure amendments (subtests and item splitting)
# ===========================================================================
# Two standard remedies applied by re-analysis. Local dependence: combine the
# dependent items into a single polytomous super-item (the subtest) whose
# score is the sum of its members. Differential item functioning: split the
# offending item into one item per group level, each observed only by that
# group, so every group receives its own location and the invariance
# violation is resolved within the model.
# ===========================================================================

#' Combine items into subtests and re-analyse
#'
#' Replaces each nominated item group by a polytomous super-item whose score is
#' the sum of its members, then refits the model. The function is commonly used
#' to examine item groups identified by \code{\link{residual_correlations}}.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param groups A list of character vectors, each naming two or more items to
#'   combine; a single vector is also accepted.
#' @param model Model for the re-analysis; defaults to \code{"PCM"}, which is
#'   almost always required because subtests change the maximum scores.
#' @return A new \code{\link{rasch}} fit on the combined structure, with the
#'   combinations recorded in its notes. Person and item estimates are
#'   recalculated. Fit grouping, external anchors on unchanged items, keyed
#'   scoring, PCM constraints and optimisation controls are retained.
#'   Anchored items cannot be combined because the resulting superitem has no
#'   corresponding external anchor.
#' @examples
#' set.seed(1); Np <- 500; L <- 8
#' d <- seq(-2, 2, length.out = L)
#' X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
#' colnames(X) <- paste0("I", 1:L)
#' X[, 5] <- ifelse(runif(Np) < 0.9, X[, 4], X[, 5])   # dependent pair
#' fit <- rasch(X)
#' fit2 <- combine_items(fit, list(c("I4", "I5")))
#' fit2$items$item
#' @seealso \code{\link{drop_items}} to remove an item rather than combine
#'   it, and \code{\link{residual_correlations}} for the dependence that
#'   motivates combining.
#' @export
combine_items <- function(fit, groups, model = "PCM") {
  if (!inherits(fit, "rasch")) stop("combine_items needs a rasch fit")
  if (inherits(fit, "rasch_mfrm"))
    stop("combine items in the long-format data and refit rasch_mfrm instead")
  if (inherits(fit, "rasch_efrm"))
    stop("combine items in the source data and refit rasch_efrm instead")
  if (is.character(groups)) groups <- list(groups)
  X <- fit$X
  used <- unlist(groups)
  bad <- setdiff(used, colnames(X))
  if (length(bad)) stop("item(s) not in the fit: ", paste(bad, collapse = ", "))
  if (anyDuplicated(used)) stop("an item appears in more than one subtest")
  short <- vapply(groups, length, 1L) < 2
  if (any(short)) stop("each subtest needs at least two items")

  keep <- setdiff(colnames(X), used)
  anchors <- fit$refit_spec$anchors
  if (!is.null(anchors) && any(as.character(anchors$item) %in% used))
    stop("an anchored item cannot be combined into a subtest; refit with ",
         "anchors on items that remain unchanged")
  # a super-item requires all its members: summing whatever happens to be
  # present would let the maximum vary by person, so partial responses stay
  # missing -- under linked designs this can empty a subtest entirely
  newcols <- lapply(groups, function(g) rowSums(X[, g, drop = FALSE]))
  empty <- vapply(newcols, function(v) sum(!is.na(v)) < 10, TRUE)
  if (any(empty))
    stop("subtest(s) with fewer than 10 persons answering every member ",
         "item (missing data): ",
         paste(vapply(groups[empty], paste, "", collapse = "+"),
               collapse = ", "),
         "; choose items answered by the same persons")
  Xn <- cbind(X[, keep, drop = FALSE], do.call(cbind, newcols))
  super_names <- vapply(groups, paste, "", collapse = "+")
  if (anyDuplicated(super_names) || any(super_names %in% keep))
    stop("the generated subtest names duplicate an existing item name")
  colnames(Xn) <- c(keep, super_names)

  refit <- if (inherits(fit, "rasch_explanatory")) {
    inherit <- c(stats::setNames(keep, keep),
                 stats::setNames(vapply(groups, `[`, "", 1L), super_names))
    .explanatory_refit_modified(fit, Xn, inherit = inherit,
                                fully_relaxed = super_names)
  } else .rasch_refit(fit, Xn, model = model)
  if (!isTRUE(refit$est$converged))
    stop("the subtest calibration did not converge; the combined analysis is unavailable")
  old_map <- fit$subtest_map %||% list()
  old_binary <- fit$subtest_binary %||% logical(0)
  members <- lapply(groups, function(g) unlist(lapply(g, function(it)
    old_map[[it]] %||% it), use.names = FALSE))
  binary <- vapply(groups, function(g) all(vapply(g, function(it) {
    if (it %in% names(old_binary)) isTRUE(old_binary[[it]])
    else fit$m[match(it, fit$items$item)] == 1L
  }, logical(1))), logical(1))
  keep_map <- old_map[intersect(names(old_map), keep)]
  keep_binary <- old_binary[intersect(names(old_binary), keep)]
  refit$subtest_map <- c(keep_map, stats::setNames(members, super_names))
  refit$subtest_binary <- c(keep_binary,
                            stats::setNames(binary, super_names))
  refit$notes <- c(refit$notes,
                   vapply(groups, function(g)
                     sprintf("subtest formed from: %s", paste(g, collapse = ", ")), ""))
  refit
}

#' Split items by a person factor to resolve DIF
#'
#' Replaces each nominated item with one item per level of a person factor,
#' each carrying that level's responses only (other levels missing). Every
#' group then receives its own item location, which resolves the invariance
#' violation flagged by \code{\link{dif_anova}}; the distance between the
#' split locations estimates the DIF size. The model is refitted with the
#' same settings.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param items Character vector naming the item(s) to split.
#' @param by The name of a person factor nominated in the fit, or a grouping
#'   vector with one entry per person.
#' @return A new \code{\link{rasch}} fit in which each split item appears as
#'   \code{"item (level)"}, with the splits recorded in its notes. Person and
#'   item estimates are recalculated with the fitted grouping, keyed scoring,
#'   anchors on unchanged items, PCM constraints and optimisation controls.
#'   An anchored item cannot be split.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' fit2 <- split_items(fit, "I3", by = "grp")
#' fit2$items$item
#' @seealso \code{\link{resolve_dif}}, which applies this iteratively;
#'   \code{\link{drop_items}}, which removes an item rather than resolving
#'   it; and \code{\link{dif_anova}}, which identifies the items to split.
#' @export
split_items <- function(fit, items, by) {
  if (!inherits(fit, "rasch")) stop("split_items needs a rasch fit")
  if (inherits(fit, "rasch_mfrm"))
    stop("split items in the long-format data and refit rasch_mfrm instead")
  if (inherits(fit, "rasch_efrm"))
    stop("amend the source data and refit rasch_efrm instead")
  X <- fit$X
  bad <- setdiff(items, colnames(X))
  if (length(bad)) stop("item(s) not in the fit: ", paste(bad, collapse = ", "))
  if (is.character(by) && length(by) == 1L) {
    if (is.null(fit$factors) || !by %in% names(fit$factors))
      stop("'", by, "' is not a person factor nominated in the fit")
    grp <- fit$factors[[by]]
  } else {
    if (length(by) != nrow(X)) stop("'by' must have one entry per person")
    grp <- by
  }
  grp <- factor(grp)
  if (nlevels(grp) < 2) stop("the splitting factor needs at least two levels")

  keep <- setdiff(colnames(X), items)
  anchors <- fit$refit_spec$anchors
  if (!is.null(anchors) && any(as.character(anchors$item) %in% items))
    stop("an anchored item cannot be split; nominate invariant anchors and ",
         "refit before resolving DIF")
  Xn <- as.data.frame(X[, keep, drop = FALSE], check.names = FALSE)
  key_extra <- NULL
  made <- list()
  for (it in items) for (lv in levels(grp)) {
    col <- X[, it]
    col[is.na(grp) | grp != lv] <- NA
    nm <- paste0(it, " (", lv, ")")
    if (nm %in% names(Xn)) stop("generated split-item name already exists: ", nm)
    if (!inherits(fit, "rasch_explanatory") &&
        !is.null(fit$mc) && it %in% colnames(fit$mc$raw)) {
      col <- fit$mc$raw[, it]
      col[is.na(grp) | grp != lv] <- NA
      kr <- fit$refit_spec$key
      if (is.data.frame(kr)) {
        kr <- kr[as.character(kr$item) == it, , drop = FALSE]
        kr$item <- nm
        key_extra <- rbind(key_extra, kr)
      }
    }
    Xn[[nm]] <- col
    made[[it]] <- c(made[[it]], nm)
  }
  refit <- if (inherits(fit, "rasch_explanatory")) {
    inherit <- c(stats::setNames(keep, keep),
                 unlist(lapply(items, function(it)
                   stats::setNames(rep(it, length(made[[it]])), made[[it]]))))
    .explanatory_refit_modified(fit, Xn, inherit = inherit,
                                location_relaxed = unlist(made,
                                                          use.names = FALSE))
  } else .rasch_refit(fit, Xn, key_extra = key_extra)
  if (!isTRUE(refit$est$converged))
    stop("the split-item calibration did not converge; the resolved analysis is unavailable")
  if (length(fit$subtest_map)) {
    sm <- fit$subtest_map[intersect(names(fit$subtest_map), keep)]
    sb <- fit$subtest_binary[intersect(names(fit$subtest_binary), keep)]
    for (it in intersect(items, names(fit$subtest_map))) {
      sm[made[[it]]] <- rep(list(fit$subtest_map[[it]]), length(made[[it]]))
      sb[made[[it]]] <- rep(isTRUE(fit$subtest_binary[[it]]),
                            length(made[[it]]))
    }
    refit$subtest_map <- sm
    refit$subtest_binary <- sb
  }
  old_map <- fit$split_map %||%
    stats::setNames(colnames(fit$X), colnames(fit$X))
  refit$split_map <- c(old_map[keep], unlist(lapply(items, function(it)
    stats::setNames(rep(unname(old_map[it]), length(made[[it]])), made[[it]]))))
  refit$notes <- c(refit$notes,
                   vapply(items, function(it)
                     sprintf("item %s split by group into: %s", it,
                             paste(made[[it]], collapse = ", ")), ""))
  refit
}

#' Resolve differential item functioning by iterative item splitting
#'
#' Splits items with uniform DIF one at a time, beginning with the largest
#' estimated effect,
#' and refits after each split. This order addresses the artificial DIF that a
#' large departure can induce in otherwise invariant items (Andrich and
#' Hagquist 2012, 2015). Each split gives the item a separate location and
#' threshold structure in every factor cell. A location split does not model a
#' group-specific discrimination, so items with non-uniform DIF are left for
#' review rather than being made untestable by a split. The procedure stops
#' when no resolvable uniform DIF remains or the remaining unsplit reference
#' set reaches \code{min_anchors}. Items fixed by external anchors are not
#' split.
#'
#' @param fit A fitted object from \code{\link{rasch}} carrying person
#'   factors.
#' @param factors Person factors to test, as in \code{\link{dif_anova}};
#'   defaults to every nominated factor.
#' @param alpha Significance level for the adjusted probabilities.
#' @param p_adjust Multiplicity adjustment across items each round.
#' @param min_anchors Minimum number of original items to leave unsplit as the
#'   internal reference set. The procedure stops before this set becomes
#'   smaller; pervasive DIF is not artificial DIF. Default
#'   \code{max(3, items / 4)}.
#' @param max_splits Hard cap on the number of splits. Default: the number
#'   of items.
#' @return A list of class \code{"rasch_resolve_dif"}: the final resolved
#'   \code{fit}, the \code{splits} performed (order, item, factor, partial
#'   eta-squared, source item, DIF magnitude in logits), the \code{stopped}
#'   reason, and the residual \code{dif} table for the final fit.
#' @references Andrich, D., & Hagquist, C. (2012). Real and artificial
#'   differential item functioning. \emph{Journal of Educational and
#'   Behavioral Statistics}, 37(3), 387-416.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1.2      # one strong DIF item
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' resolve_dif(fit)$splits
#' @seealso \code{\link{split_items}} for a single split,
#'   \code{\link{drop_items}} to remove an item instead, and
#'   \code{\link{dif_anova}} for the test it resolves.
#' @export
resolve_dif <- function(fit, factors = NULL, alpha = 0.05, p_adjust = "holm",
                        min_anchors = NULL, max_splits = NULL) {
  if (!inherits(fit, "rasch") || inherits(fit, c("rasch_mfrm", "rasch_efrm")))
    stop("resolve_dif needs an ordinary rasch fit with person factors")
  fac0 <- .dif_factors(fit, factors)
  fnames <- names(fac0)
  L0 <- ncol(fit$X)
  if (L0 < 3L) stop("DIF resolution needs at least three fitted items")
  if (is.null(min_anchors))
    min_anchors <- min(L0 - 1L, max(3L, L0 %/% 4L))
  if (is.null(max_splits)) max_splits <- L0
  if (length(min_anchors) != 1L || !is.finite(min_anchors) ||
      min_anchors < 2L || min_anchors != floor(min_anchors) ||
      min_anchors >= L0)
    stop("min_anchors must be a whole number from 2 to one fewer than the ",
         "number of fitted items")
  if (length(max_splits) != 1L || !is.finite(max_splits) ||
      max_splits < 0L || max_splits != floor(max_splits))
    stop("max_splits must be a non-negative whole number")

  # significant, non-superseded group terms of the current fit, with the
  # factors to split by and the partial eta-squared to rank on
  flagged <- function(cur, resolvable_only = TRUE) {
    keep <- intersect(fnames, names(cur$factors))
    if (!length(keep)) return(NULL)
    da <- dif_anova(cur, factors = keep, p_adjust = p_adjust, alpha = alpha)
    s <- da$summary
    # Splitting supplies cell-specific locations and thresholds, not
    # cell-specific slopes. Restrict automatic resolution to uniform-only DIF;
    # otherwise the split removes the between-group overlap and merely makes a
    # non-uniform violation disappear from the DIF table.
    take <- if (resolvable_only)
      s$uniform_DIF & !s$nonuniform_DIF & !s$superseded
    else
      (s$uniform_DIF | s$nonuniform_DIF) & !s$superseded
    vars <- da$summary_factors[take]
    s <- s[take, , drop = FALSE]
    if (!nrow(s)) return(NULL)
    out <- data.frame(
      item = s$item, factor = vapply(vars, .dif_term_label, ""),
      eta2 = pmax(ifelse(s$uniform_DIF, s$eta2_uniform, 0),
                  ifelse(s$nonuniform_DIF, s$eta2_nonuniform, 0)),
      uniform = s$uniform_DIF, nonuniform = s$nonuniform_DIF,
      stringsAsFactors = FALSE)
    attr(out, "term_factors") <- vars
    out
  }
  cur <- fit
  splits <- list(); done <- character(0); skipped <- character(0)
  skipped_anchor <- character(0)
  stopped <- "no significant DIF remains"
  repeat {
    fl <- flagged(cur)
    if (is.null(fl) || !nrow(fl)) {
      remaining <- flagged(cur, resolvable_only = FALSE)
      if (!is.null(remaining) && any(remaining$nonuniform))
        stopped <- paste("no resolvable uniform DIF remains;",
                         "non-uniform DIF requires item review")
      break
    }
    vars <- attr(fl, "term_factors")
    ord <- order(-fl$eta2)
    fl <- fl[ord, , drop = FALSE]
    vars <- vars[ord]
    # first flagged item-factor not already split, ranked by effect size
    pick <- NULL
    for (r in seq_len(nrow(fl))) {
      vk <- paste0(nchar(vars[[r]], type = "bytes"), ":", vars[[r]], ";",
                   collapse = "")
      key <- paste0(nchar(fl$item[r], type = "bytes"), ":", fl$item[r], ";", vk)
      if (!key %in% done) { pick <- fl[r, ]; pick_vars <- vars[[r]]; break }
    }
    if (is.null(pick)) { stopped <- "remaining DIF cannot be resolved further"; break }
    if (length(splits) >= max_splits) { stopped <- "reached the split cap"; break }
    n_anchor <- L0 - length(unique(vapply(splits, `[[`, "", "base_item")))
    if (n_anchor <= min_anchors) {
      stopped <- sprintf("stopped to keep %d anchor items (pervasive DIF is not artificial DIF)",
                         min_anchors)
      break
    }
    by_vars <- pick_vars
    grp <- if (length(by_vars) == 1L) cur$factors[[by_vars]] else
      .factor_cells(cur$factors[by_vars], sep = ":")
    external <- cur$refit_spec$anchors
    if (!is.null(external) &&
        pick$item %in% as.character(external$item)) {
      skipped_anchor <- c(skipped_anchor, paste(pick$item, pick$factor))
      done <- c(done, key)
      next
    }
    # Confirm the uniform flag with dif_size before splitting: dif_size drops
    # thin cells (min_n) and withholds weak categories, so if it can no
    # longer see the DIF the ANOVA flag rested on a near-empty cell, not on
    # real differential functioning -- splitting on it chases noise.
    ds <- tryCatch(dif_size(cur, pick$item, by = grp), error = function(e) NULL)
    if (isTRUE(pick$uniform)) {
      # Failure to obtain a magnitude is itself a failed confirmation. This
      # most often means min_n removed a thin level and fewer than two usable
      # levels remain; proceeding would recreate exactly the thin-cell split
      # this gate is intended to prevent.
      if (is.null(ds)) {
        skipped <- c(skipped, paste(pick$item, pick$factor))
        done <- c(done, key)
        next
      }
      n_grp_lev <- nlevels(droplevels(as.factor(grp)))
      reduced <- nrow(ds$levels) < n_grp_lev || isTRUE(any(ds$levels$weak))
      if (reduced && !isTRUE(any(ds$pairs$significant, na.rm = TRUE))) {
        skipped <- c(skipped, paste(pick$item, pick$factor))
        done <- c(done, key)
        next
      }
    }
    # DIF magnitude in logits, over the trustworthy (non-weak) pairs only
    mag <- if (!is.null(ds)) {
      d <- abs(ds$pairs$difference[!is.na(ds$pairs$se)])
      if (length(d)) max(d) else NA_real_
    } else NA_real_
    refit <- tryCatch(split_items(cur, pick$item, by = grp),
                      error = function(e) NULL)
    if (is.null(refit)) { done <- c(done, key); next }
    base_map <- cur$split_map %||%
      stats::setNames(colnames(cur$X), colnames(cur$X))
    splits[[length(splits) + 1L]] <- list(
      order = length(splits) + 1L, item = pick$item, factor = pick$factor,
      base_item = unname(base_map[pick$item]), eta2 = pick$eta2,
      magnitude = mag)
    done <- c(done, key)
    cur <- refit
  }
  split_df <- if (length(splits))
    do.call(rbind, lapply(splits, function(s) as.data.frame(s,
                                                            stringsAsFactors = FALSE))) else
    data.frame(order = integer(), item = character(), factor = character(),
               base_item = character(), eta2 = numeric(), magnitude = numeric())
  rownames(split_df) <- NULL
  final_dif <- tryCatch(flagged(cur, resolvable_only = FALSE),
                        error = function(e) NULL)
  notes <- character(0)
  if (length(skipped)) notes <- c(notes,
    sprintf("%d flagged item-factor(s) not split because the resolved magnitude could not be confirmed with a near-empty category (%s)",
            length(skipped), paste(skipped, collapse = "; ")))
  if (length(skipped_anchor)) notes <- c(notes,
    sprintf("%d externally anchored item-factor(s) not split (%s)",
            length(skipped_anchor), paste(skipped_anchor, collapse = "; ")))
  out <- list(fit = cur, splits = split_df, n_splits = nrow(split_df),
              stopped = stopped, dif = final_dif, notes = notes,
              n_remaining_dif = if (is.null(final_dif)) 0L else nrow(final_dif),
              n_nonuniform = if (is.null(final_dif)) 0L else
                sum(final_dif$nonuniform %in% TRUE))
  out <- .tag_tables(out)
  class(out) <- "rasch_resolve_dif"
  out
}

#' @export
print.rasch_resolve_dif <- function(x, ...) {
  cat(sprintf("Iterative DIF resolution: %d split(s); %s\n",
              x$n_splits, x$stopped))
  if (x$n_splits) {
    d <- x$splits; d$eta2 <- round(d$eta2, 3); d$magnitude <- round(d$magnitude, 3)
    print(d, row.names = FALSE)
  }
  cat(sprintf("Remaining items with significant DIF: %d\n", x$n_remaining_dif))
  if (!is.null(x$n_nonuniform) && x$n_nonuniform)
    cat(sprintf("Non-uniform item-factor findings requiring review: %d\n",
                x$n_nonuniform))
  invisible(x)
}
