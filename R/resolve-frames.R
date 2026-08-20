# rasch :: resolving an item that does not hold across frames
# ===========================================================================
# The frame model gives each item one location, shared across the frames it
# appears in and scaled by the frame unit. That tie is what makes the frames
# comparable, and it is an assumption. When frame_invariance() finds an item
# the tie does not hold for, there are two remedies. Removing the item takes
# it out of every frame, so it stops measuring anyone, including in the
# frames it behaved perfectly well in. Resolving it gives it a separate
# location in each frame: it stops linking the frames, which is what the
# diagnosis says is wrong with it, and goes on measuring the person within
# their own frame.
#
# Resolving costs a parameter per extra frame and removes the item from the
# link, so the group units then rest on the items that remain common. It is
# the better remedy when the item measures well within frames and the worse
# one when it does not measure at all.
# ===========================================================================

#' Resolve items that do not hold across frames
#'
#' Gives named items a separate location in each frame and refits, so they
#' continue to measure persons within their own frame while no longer
#' constraining the comparison between frames. The result is an ordinary
#' frame fit, so every diagnostic applies to it unchanged.
#'
#' @details
#' The fitted model represents an item's threshold in a frame as the frame
#' unit times one location shared across frames, so an item that behaves
#' differently in one frame is misrepresented in all of them, and the group
#' units absorb part of the discrepancy.
#' \code{\link{frame_invariance}} tests for such items;
#' this function and \code{\link{drop_items}} are the two remedies.
#'
#' Resolving is the milder one. The item is replaced by one version per
#' frame, named \code{"item (frame)"}, each answered by that frame's persons
#' alone. Each version keeps its own location, so the item still contributes
#' to the person estimates of everyone who answered it, and it no longer
#' contributes to the link between frames. Dropping the item removes that
#' contribution as well, from every frame at once.
#'
#' The link is what pays for it. Person-group units are identified by the
#' items two frames have in common, so a resolved item leaves the units
#' resting on the items that remain shared, and resolving too many leaves
#' them unidentified. This function refuses when a set would be left with
#' fewer than two common items; the model's own connectivity check catches
#' the remaining cases.
#'
#' Prefer resolving when the item measures well inside each frame and only
#' its comparability is in doubt, and dropping when the item is a poor
#' measure wherever it appears. The distinction is empirical: compare the
#' unit estimates and the person standard errors the two remedies produce.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param items Item names to resolve.
#' @param boot_reps Bootstrap replicates for the refit, resolved as in
#'   \code{\link{drop_items}}: the refit keeps the character of the fit it
#'   came from unless a number is given.
#' @return A refitted object of class \code{"rasch_efrm"}, carrying a note
#'   for each item resolved. The resolved versions appear in the item table
#'   as \code{"item (frame)"}.
#' @seealso \code{\link{frame_invariance}}, which identifies the items to
#'   resolve; \code{\link{drop_items}}, which removes an item instead; and
#'   \code{\link{split_items}}, the equivalent for an ordinary fit.
#' @examples
#' d <- simulate_efrm(n_per_group = 200, items_per_set = 6, n_sets = 2,
#'                    n_groups = 2, set_unit_ratio = 1.3, seed = 5)
#' tr <- attr(d, "truth")
#' fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
#'                   id = "id", boot_reps = 0)
#' fit2 <- resolve_frames(fit, "S1I02", boot_reps = 0)
#' grep("S1I02", fit2$items$item, value = TRUE)
#' @export
resolve_frames <- function(fit, items, boot_reps = NULL) {
  if (!inherits(fit, "rasch"))
    stop("resolve_frames needs a fit from rasch_efrm()")
  if (!inherits(fit, "rasch_efrm"))
    stop("resolve_frames needs a frame model; for an ordinary fit with ",
         "person factors use split_items() instead")
  if (!length(items)) stop("name at least one item to resolve")
  items <- as.character(items)

  all_items <- names(fit$set_of)
  bad <- setdiff(items, all_items)
  if (length(bad))
    stop("item(s) not in the fit: ", paste(bad, collapse = ", "),
         "; the fit holds: ", paste(utils::head(all_items, 8), collapse = ", "),
         if (length(all_items) > 8) ", ..." else "")

  # a resolved item is common to no two frames, so it stops identifying the
  # group units; the items left shared have to carry that on their own
  common <- setdiff(all_items, items)
  by_set <- split(common, fit$set_of[common])
  short <- setdiff(unique(fit$set_of), names(by_set))
  thin <- c(short, names(by_set)[vapply(by_set, length, 0L) < 2L])
  if (length(thin))
    stop("resolving those items would leave set(s) ",
         paste(sort(thin), collapse = ", "),
         " with fewer than two items common to the frames, so the group ",
         "units would be unidentified; resolve fewer items, or use ",
         "drop_items() if the item does not measure well in any frame")

  grp <- fit$factors[[fit$frame_group]]
  glev <- levels(factor(grp))
  gvec <- as.character(grp)

  # the virtual columns of one item across frames are complementary, so the
  # source response is the non-missing value among them
  src <- vapply(all_items, function(it) {
    cols <- intersect(paste0(it, ":", glev), colnames(fit$X))
    v <- rep(NA_real_, nrow(fit$X))
    for (cc in cols) {
      w <- !is.na(fit$X[, cc])
      v[w] <- fit$X[w, cc]
    }
    v
  }, numeric(nrow(fit$X)))
  colnames(src) <- all_items

  cols <- list()
  new_set <- character()
  made <- list()
  for (it in all_items) {
    if (!it %in% items) {
      cols[[it]] <- src[, it]
      new_set[it] <- fit$set_of[[it]]
      next
    }
    for (lv in glev) {
      v <- src[, it]
      v[is.na(gvec) | gvec != lv] <- NA
      if (all(is.na(v))) next          # the item was not taken in this frame
      nm <- paste0(it, " (", lv, ")")
      cols[[nm]] <- v
      new_set[nm] <- fit$set_of[[it]]
      made[[it]] <- c(made[[it]], lv)
    }
  }

  d <- data.frame(id = fit$person$id, as.data.frame(cols, check.names = FALSE),
                  fit$factors, check.names = FALSE, stringsAsFactors = FALSE)
  extra <- setdiff(names(fit$factors), fit$frame_group)
  reps <- boot_reps
  if (is.null(reps)) {
    reps <- fit$boot_reps_used
    if (is.null(reps) || !is.finite(reps))
      reps <- if (any(is.finite(fit$alpha_table$se_log_alpha))) NULL else 0L
  }
  refit <- rasch_efrm(d, item_sets = split(names(new_set), new_set),
                      groups = fit$frame_group, id = "id",
                      factors = if (length(extra)) extra else NULL,
                      se_method = fit$se_method, boot_reps = reps)
  refit$notes <- c(refit$notes,
                   vapply(names(made), function(it)
                     sprintf("item %s resolved by frame into: %s", it,
                             paste0(it, " (",
                                    paste(made[[it]], collapse = "/"), ")")),
                     character(1)))
  refit
}
