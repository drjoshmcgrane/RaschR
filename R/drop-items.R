# rasch :: removing items from a fitted analysis
# ===========================================================================

.frame_group_vars <- function(fit) {
  fg <- fit$frame_group
  if (is.null(fg) || !length(fg)) stop("the frame variables are not recorded")
  if (length(fg) > 1L) fg[-1L] else fg[1L]
}

.frame_group_values <- function(fit) {
  fg <- fit$frame_group[1L]
  if (is.null(fit$factors) || !fg %in% names(fit$factors))
    stop("the fitted frame membership is not available")
  fit$factors[[fg]]
}

# Recover the unexpanded item responses through the explicit virtual map.
# Virtual names are deliberately not parsed: item and level labels may contain
# colons, regular-expression metacharacters, or other valid punctuation.
.efrm_source_matrix <- function(fit, items = names(fit$set_of)) {
  vm <- fit$virtual_map
  if (is.null(vm) || !all(c("item", "vkey") %in% names(vm)))
    stop("the fitted frame-to-item map is not available")
  bad <- setdiff(items, unique(vm$item))
  if (length(bad)) stop("item(s) absent from the frame-to-item map: ",
                        paste(bad, collapse = ", "))
  src <- vapply(items, function(it) {
    cols <- vm$vkey[vm$item == it & vm$vkey %in% colnames(fit$X)]
    v <- rep(NA_real_, nrow(fit$X))
    for (cc in cols) {
      take <- !is.na(fit$X[, cc])
      v[take] <- fit$X[take, cc]
    }
    v
  }, numeric(nrow(fit$X)))
  if (!is.matrix(src)) src <- matrix(src, ncol = length(items))
  colnames(src) <- items
  src
}

.efrm_refit <- function(fit, source, set_of, boot_reps = NULL,
                        ids = fit$person$id, factors = fit$factors,
                        se_method = NULL) {
  spec <- fit$refit_spec
  if (is.null(spec)) spec <- list()
  group_vars <- spec$groups
  if (is.null(group_vars) || !all(group_vars %in% names(factors)))
    group_vars <- .frame_group_vars(fit)
  extra <- spec$factors
  if (is.null(extra)) extra <- setdiff(names(factors), fit$frame_group)
  extra <- intersect(extra, names(factors))
  need <- unique(c(group_vars, extra))
  d <- data.frame(.rasch_id = ids, source, factors[, need, drop = FALSE],
                  check.names = FALSE, stringsAsFactors = FALSE)
  reps <- if (is.null(boot_reps)) spec$boot_reps else boot_reps
  if (is.null(reps)) {
    reps <- fit$boot_reps_used
    if (is.null(reps) || !is.finite(reps))
      reps <- if (any(is.finite(fit$alpha_table$se_log_alpha))) NULL else 0L
  }
  do.call(rasch_efrm, list(
    data = d, item_sets = split(names(set_of), unname(set_of)),
    groups = group_vars, id = ".rasch_id",
    factors = if (length(extra)) extra else NULL, items = colnames(source),
    n_groups = spec$n_groups %||% fit$n_groups,
    adjust_N = spec$adjust_N %||% NA_real_, na_codes = spec$na_codes %||% -1,
    maxit = spec$maxit %||% 50, tol = spec$tol %||% 1e-7,
    min_link_persons = spec$min_link_persons %||% 30,
    se_method = se_method %||% spec$se_method %||% fit$se_method,
    boot_reps = reps, workers = spec$workers %||% 1L,
    seed = spec$seed %||% NULL))
}

.rasch_refit <- function(fit, source, model = NULL, key_extra = NULL,
                         require_anchor = TRUE, na_codes = NULL) {
  spec <- fit$refit_spec
  if (is.null(spec)) spec <- list()
  if (!is.null(na_codes)) spec$na_codes <- na_codes
  source <- as.data.frame(source, check.names = FALSE,
                          stringsAsFactors = FALSE)
  keep <- names(source)
  if (inherits(fit, "rasch_explanatory") && all(keep %in% colnames(fit$X)))
    return(.explanatory_refit_modified(
      fit, source,
      inherit = stats::setNames(keep, keep)))
  key <- spec$key
  if (!is.null(fit$mc) && !is.null(key)) {
    raw_items <- intersect(colnames(fit$mc$raw), keep)
    for (it in raw_items) source[[it]] <- fit$mc$raw[, it]
    if (is.data.frame(key)) {
      key <- key[as.character(key$item) %in% keep, , drop = FALSE]
    } else {
      key <- key[intersect(names(key), keep)]
    }
  }
  if (!is.null(key_extra)) {
    if (!is.data.frame(key) && !is.null(key))
      stop("internal refit cannot combine incompatible key formats")
    key <- rbind(key, key_extra)
  }
  if (is.data.frame(key) && !nrow(key)) key <- NULL
  if (!is.data.frame(key) && !is.null(key) && !length(key)) key <- NULL
  anchors <- spec$anchors
  if (!is.null(anchors)) {
    anchor_items <- as.character(anchors$item)
    anchors <- anchors[anchor_items %in% keep, , drop = FALSE]
    if (!nrow(anchors) && require_anchor)
      stop("dropping those items would remove every anchor and lose the ",
           "fitted scale origin; retain an anchor or refit explicitly")
    if (!nrow(anchors)) anchors <- NULL
  }
  do.call(rasch, list(
    data = source, model = model %||% spec$model %||% fit$model,
    id = fit$person$id,
    factors = fit$factors, n_groups = spec$n_groups %||% fit$n_groups,
    adjust_N = spec$adjust_N %||% NA_real_, anchors = anchors,
    na_codes = spec$na_codes %||% -1, key = key,
    pc_components = spec$pc_components,
    maxit = spec$maxit %||% 60, tol = spec$tol %||% 1e-8))
}

.rasch_refit_after_drop <- function(fit, keep) {
  .rasch_refit(fit, fit$X[, keep, drop = FALSE])
}
# Item screening is a step in an analysis, not a preliminary to it: an item
# is judged by its behaviour in the fit, and judging it means refitting
# without it. For frame models the refit matters more than usual: an item that
# does not follow its set's unit changes both the within-frame calibration and
# the person-side link used to compare that set with another.
# ===========================================================================

#' Drop items and refit
#'
#' Removes named items and refits the analysis with the same model
#' specification.
#'
#' @details
#' The refit retains person identifiers and factors, class-interval settings,
#' optimisation controls, anchors, multiple-choice scoring and PCM component
#' constraints. An EFRM refit also retains the item-set and crossed-frame
#' design, linking controls and uncertainty method. The operation is refused
#' if it would remove every anchor, empty an item set or leave the model
#' unidentified.
#'
#' Item removal changes both the item calibration and the person estimates.
#' For an EFRM it can also change the estimated frame units. Compare the
#' original and refitted results as a sensitivity analysis.
#'
#' @param fit A fitted object from \code{\link{rasch}} or
#'   \code{\link{rasch_efrm}}. Many-facet fits are refused: remove the
#'   item's rows from the long-format data and refit
#'   \code{\link{rasch_mfrm}} instead.
#' @param items Item names to remove.
#' @param boot_reps Bootstrap replicates for an EFRM refit. The default retains
#'   the fitted specification; a number overrides it.
#' @return A refitted object of the same class as \code{fit}, carrying a note
#'   recording which items were dropped.
#' @seealso \code{\link{frame_invariance}} and
#'   \code{\link{resolve_frames}} for frame models;
#'   \code{\link{split_items}} and \code{\link{resolve_dif}} for DIF; and
#'   \code{\link{combine_items}} for response dependence.
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

    src <- .efrm_source_matrix(fit, keep)
    refit <- .efrm_refit(fit, src, sets_left, boot_reps = boot_reps)
  } else {
    bad <- setdiff(items, colnames(fit$X))
    if (length(bad))
      stop("item(s) not in the fit: ", paste(bad, collapse = ", "))
    keep <- setdiff(colnames(fit$X), items)
    if (length(keep) < 2L)
      stop("dropping those items would leave fewer than two items")
    refit <- .rasch_refit_after_drop(fit, keep)
  }
  if (!isTRUE(refit$est$converged))
    stop("the reduced calibration did not converge; the dropped-item analysis is unavailable")
  if (inherits(refit, "rasch_efrm") &&
      any(refit$linking$alpha_edges$converged %in% FALSE))
    stop("the reduced calibration's set-unit link did not converge; the dropped-item analysis is unavailable")
  if (length(fit$subtest_map)) {
    sk <- intersect(names(fit$subtest_map), refit$items$item)
    refit$subtest_map <- fit$subtest_map[sk]
    refit$subtest_binary <- fit$subtest_binary[sk]
  }
  refit$notes <- c(refit$notes,
                   sprintf("dropped item(s): %s", paste(items, collapse = ", ")))
  refit
}
