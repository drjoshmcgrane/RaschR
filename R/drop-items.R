# rasch :: removing items from a fitted analysis
# ===========================================================================
# Item screening is a step in an analysis, not a preliminary to it: an item
# is judged by its behaviour in the fit, and judging it means refitting
# without it. For frame models the refit matters more than usual, because a
# set unit is estimated from the dispersion its own items produce, so an
# item that shares no unit with its set moves the very quantity used to
# decide whether the set differs from another.
# ===========================================================================

#' Drop items and refit
#'
#' Removes named items from a fitted analysis and refits it, keeping the
#' original model, person identifiers, factors, and (for frame models) the
#' set structure and standard-error method. The result is an ordinary fit of
#' the same class, so every diagnostic applies to it unchanged.
#'
#' @details
#' Screening items is part of an analysis of frames rather than a step before
#' it. A person-group unit is estimated from the same items in every frame,
#' so an item with differential item functioning distorts it; an item-set
#' unit is estimated from the dispersion of person estimates within each set,
#' so an item that fits its set badly distorts that set alone, with nothing
#' in the other set to offset it. In simulation the second effect is large:
#' at eight items per set, two under-discriminating items in one set moved a
#' planted unit ratio of 1.40 to 1.73, and four over-discriminating items
#' moved it to 1.02. Dropping such an item and comparing the units before
#' and after is therefore a substantive sensitivity analysis, not
#' housekeeping.
#'
#' An item that fits no set well is usually better removed than reassigned,
#' and it cannot be given a set of its own: a single item carries no
#' dispersion from which to estimate a unit.
#'
#' @param fit A fitted object from \code{\link{rasch}} or
#'   \code{\link{rasch_efrm}}.
#' @param items Item names to remove.
#' @param boot_reps For frame models, the number of linking bootstrap
#'   replicates for the refit. Defaults to the number the original fit used.
#' @return A refitted object of the same class as \code{fit}, carrying a note
#'   recording which items were dropped.
#' @seealso \code{\link{frame_invariance}}, which identifies the items a
#'   frame model's assumption does not hold for; \code{\link{split_items}}
#'   and \code{\link{resolve_dif}}, which resolve an item rather than
#'   remove it; and \code{\link{combine_items}}.
#' @examples
#' d <- simulate_rasch(300, 8, seed = 1)
#' fit <- rasch(d, id = "id")
#' fit2 <- drop_items(fit, "I03")
#' nrow(fit2$items)
#' @export
drop_items <- function(fit, items, boot_reps = NULL) {
  if (!inherits(fit, "rasch"))
    stop("drop_items needs a fit from rasch() or rasch_efrm()")
  if (inherits(fit, "rasch_mfrm"))
    stop("remove the item's rows from the long-format data and refit ",
         "rasch_mfrm() instead")
  if (!length(items)) stop("name at least one item to drop")
  items <- as.character(items)

  if (inherits(fit, "rasch_efrm")) {
    all_items <- names(fit$set_of)
    bad <- setdiff(items, all_items)
    if (length(bad))
      stop("item(s) not in the fit: ", paste(bad, collapse = ", "),
           "; the fit holds: ", paste(utils::head(all_items, 8), collapse = ", "),
           if (length(all_items) > 8) ", ..." else "")
    keep <- setdiff(all_items, items)
    if (!length(keep)) stop("dropping those items would leave no items")
    sets_left <- fit$set_of[keep]
    gone <- setdiff(unique(fit$set_of), unique(sets_left))
    if (length(gone))
      stop("dropping those items would empty set(s): ",
           paste(gone, collapse = ", "),
           "; a set with no items cannot carry a unit")
    if (length(unique(sets_left)) < 2L && length(unique(fit$set_of)) > 1L)
      stop("dropping those items would leave a single item set, so no set ",
           "unit is identified; drop fewer items or refit with rasch()")

    # rebuild the source responses: each person answered an item in exactly
    # one frame, so the item's virtual columns are complementary
    grp <- fit$factors[[fit$frame_group]]
    glev <- levels(factor(grp))
    src <- vapply(keep, function(it) {
      cols <- intersect(paste0(it, ":", glev), colnames(fit$X))
      v <- rep(NA_real_, nrow(fit$X))
      for (cc in cols) {
        w <- !is.na(fit$X[, cc])
        v[w] <- fit$X[w, cc]
      }
      v
    }, numeric(nrow(fit$X)))
    colnames(src) <- keep

    d <- data.frame(id = fit$person$id, src, fit$factors,
                    check.names = FALSE, stringsAsFactors = FALSE)
    extra <- setdiff(names(fit$factors), fit$frame_group)
    # boot_reps_used is NA for two different fits: one whose standard errors
    # came from the analytic route rather than a bootstrap, which is ordinary
    # for a single person group, and one asked for no standard errors at all.
    # Reading both as zero refits the first without the standard errors it
    # had, leaving no test on the set units. Tell them apart by whether the
    # source fit carries any, so the refit keeps the character of the fit it
    # came from. An explicit boot_reps is still honoured.
    reps <- boot_reps
    if (is.null(reps)) {
      reps <- fit$boot_reps_used
      if (is.null(reps) || !is.finite(reps))
        reps <- if (any(is.finite(fit$alpha_table$se_log_alpha))) NULL else 0L
    }
    refit <- rasch_efrm(d, item_sets = split(keep, sets_left),
                        groups = fit$frame_group, id = "id",
                        factors = if (length(extra)) extra else NULL,
                        se_method = fit$se_method, boot_reps = reps)
  } else {
    bad <- setdiff(items, colnames(fit$X))
    if (length(bad))
      stop("item(s) not in the fit: ", paste(bad, collapse = ", "))
    keep <- setdiff(colnames(fit$X), items)
    if (length(keep) < 2L)
      stop("dropping those items would leave fewer than two items")
    refit <- rasch(fit$X[, keep, drop = FALSE], model = fit$model,
                   id = fit$person$id, factors = fit$factors,
                   n_groups = fit$n_groups)
  }
  refit$notes <- c(refit$notes,
                   sprintf("dropped item(s): %s", paste(items, collapse = ", ")))
  refit
}
