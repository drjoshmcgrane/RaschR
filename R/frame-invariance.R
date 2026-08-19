# rasch :: testing the invariance a frame model assumes
# ===========================================================================
# rasch_efrm() constrains each item to a single location shared across the
# frames it appears in, scaled by the frame unit. That constraint is the
# model's substantive claim -- the frames differ in unit, not in what the
# items measure -- and it cannot be tested from the fit itself, because the
# fit imposes it. Testing it means calibrating each frame on its own and
# comparing, which is what this function does.
# ===========================================================================

#' Test the item invariance a frame model assumes
#'
#' Calibrates each frame separately and compares each item across the frames
#' it appears in, on two counts: whether it keeps its location once the
#' frame units are accounted for, and whether it discriminates alike. A
#' frame model assumes both, differing only by the frame's unit; this
#' function tests the assumption rather than imposing it.
#'
#' @details
#' The comparison is possible only where an item set is taken by more than
#' one person group, since an item must appear in at least two frames to be
#' compared across them. Item sets partition the items, so there is no
#' equivalent test across sets: screen those with the ordinary item fit
#' statistics within each set instead.
#'
#' Each frame is refitted with \code{\link{rasch}} on its own persons and
#' items, giving locations in that frame's natural unit and centred on that
#' frame's origin. Dividing by the frame unit from the original fit puts
#' them on the common scale, where the model says they should agree. The
#' reported difference is between the two calibrations, its standard error
#' combines theirs (the person samples are disjoint, so the calibrations are
#' independent), and probabilities are Holm-adjusted across items.
#'
#' The summary compares the root mean squared difference with the root mean
#' squared standard error, following Humphry (2005). A root mean squared
#' difference materially larger than the root mean squared standard error
#' indicates item behaviour that the frame units do not account for, whether
#' or not individual items reach significance.
#'
#' A location comparison cannot detect a difference in discrimination: a
#' steeper item still crosses one half in the same place, so its location
#' survives intact. The second comparison uses the within-frame fit
#' statistics, which carry no unit because each is computed against its own
#' frame's model, and treats the difference of two independent standardised
#' infit statistics as having variance 2. It is conservative and needs a
#' reasonable sample: in simulation, with 8 items and two items
#' discriminating half again as steeply in one frame, it detected them in
#' 12\% of replicates at 500 persons per frame and 85\% at 2,000, with
#' false-positive rates near 1\% throughout. Treat a null result at a few
#' hundred persons per frame as uninformative rather than reassuring.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param alpha Significance level for flagging items.
#' @return A list of class \code{"rasch_frame_invariance"} with
#'   \code{locations} (one row per item and frame pair: locations on the
#'   common scale, their difference, its standard error, statistic, adjusted
#'   probability, and flag), \code{discrimination} (the same items compared
#'   on their within-frame infit statistics), and \code{summary} (per frame
#'   pair: the number of items, the root mean squared difference, the root
#'   mean squared standard error, their ratio, and the number of items
#'   flagged on each count).
#' @references
#' Humphry, S. M. (2005). \emph{Maintaining a Common Arbitrary Unit in Social
#' Measurement}. PhD thesis, Murdoch University.
#' @seealso \code{\link{drop_items}} to remove an item the test flags, and
#'   \code{\link{rasch_efrm}} for the model whose assumption is tested.
#' @examples
#' d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 1,
#'                    n_groups = 2, group_unit_ratio = 1.4, seed = 2)
#' tr <- attr(d, "truth")
#' fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
#'                   id = "id", boot_reps = 0)
#' frame_invariance(fit)
#' @export
frame_invariance <- function(fit, alpha = 0.05) {
  if (!inherits(fit, "rasch_efrm"))
    stop("frame_invariance needs a fit from rasch_efrm()")
  grp <- fit$factors[[fit$frame_group]]
  glev <- levels(factor(grp))
  if (length(glev) < 2L)
    stop("item invariance across frames needs at least two person groups: ",
         "with one group each item appears in a single frame, and item sets ",
         "partition the items, so use the item fit statistics within each ",
         "set instead")

  sets <- unique(fit$set_of)
  rho <- fit$frames
  out <- list(); dsc <- list()
  for (s in sets) {
    its <- names(fit$set_of)[fit$set_of == s]
    # calibrate the set separately within each group
    cal <- list()
    for (g in glev) {
      cols <- intersect(paste0(its, ":", g), colnames(fit$X))
      if (!length(cols)) next
      Xg <- fit$X[!is.na(fit$X[, cols[1]]), cols, drop = FALSE]
      colnames(Xg) <- sub(paste0(":", g, "$"), "", colnames(Xg))
      f <- tryCatch(rasch(Xg), error = function(e) NULL)
      if (is.null(f)) next
      r <- rho$rho[rho$set == s & rho$group == g]
      if (!length(r) || !is.finite(r) || r <= 0) next
      cal[[g]] <- data.frame(item = f$items$item,
                             loc = f$items$location / r,
                             se = f$items$se / r,
                             infit = f$items$infit_ms,
                             infit_z = f$items$infit_z,
                             stringsAsFactors = FALSE)
    }
    if (length(cal) < 2L) next
    gg <- names(cal)
    for (a in seq_len(length(gg) - 1L)) for (b in (a + 1L):length(gg)) {
      m <- merge(cal[[gg[a]]], cal[[gg[b]]], by = "item",
                 suffixes = c("_1", "_2"))
      if (!nrow(m)) next
      d <- m$loc_2 - m$loc_1
      d <- d - mean(d)                       # origins are separately centred
      se <- sqrt(m$se_1^2 + m$se_2^2)
      z <- d / se
      out[[length(out) + 1L]] <- data.frame(
        set = s, frame_1 = gg[a], frame_2 = gg[b], item = m$item,
        location_1 = m$loc_1, location_2 = m$loc_2,
        difference = d, se = se, statistic = z,
        p = 2 * stats::pnorm(-abs(z)), stringsAsFactors = FALSE)
      # discrimination: a location comparison cannot see a slope
      # difference, since a steeper item still crosses one half in the
      # same place. Fit statistics can: each is computed against its own
      # frame's model, so it carries no unit and the frames are
      # comparable directly. Two independent standardised statistics
      # differ with variance 2.
      zd <- (m$infit_z_1 - m$infit_z_2) / sqrt(2)
      dsc[[length(dsc) + 1L]] <- data.frame(
        set = s, frame_1 = gg[a], frame_2 = gg[b], item = m$item,
        infit_1 = m$infit_1, infit_2 = m$infit_2,
        statistic = zd, p = 2 * stats::pnorm(-abs(zd)),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out))
    stop("no item set is taken by two person groups, so no item appears in ",
         "two frames to be compared")
  cmp <- do.call(rbind, out)
  cmp$p_adj <- stats::p.adjust(cmp$p, method = "holm")
  cmp$flagged <- cmp$p_adj < alpha
  rownames(cmp) <- NULL
  dsc <- do.call(rbind, dsc)
  dsc$p_adj <- stats::p.adjust(dsc$p, method = "holm")
  dsc$flagged <- dsc$p_adj < alpha
  rownames(dsc) <- NULL

  key <- paste(cmp$set, cmp$frame_1, cmp$frame_2, sep = "|")
  keyd <- paste(dsc$set, dsc$frame_1, dsc$frame_2, sep = "|")
  smry <- do.call(rbind, lapply(unique(key), function(k) {
    z <- cmp[key == k, ]; y <- dsc[keyd == k, ]
    data.frame(set = z$set[1], frame_1 = z$frame_1[1], frame_2 = z$frame_2[1],
      n_items = nrow(z),
      rmsd = sqrt(mean(z$difference^2)),
      rmse = sqrt(mean(z$se^2)),
      ratio = sqrt(mean(z$difference^2)) / sqrt(mean(z$se^2)),
      n_location = sum(z$flagged), n_discrimination = sum(y$flagged),
      stringsAsFactors = FALSE)
  }))
  rownames(smry) <- NULL
  structure(list(locations = cmp, discrimination = dsc, summary = smry,
                 alpha = alpha), class = "rasch_frame_invariance")
}

#' @export
print.rasch_frame_invariance <- function(x, ...) {
  cat("Item invariance across frames (each frame calibrated separately)\n\n")
  print(.fmt_df(x$summary), row.names = FALSE)
  cat("\nrmsd/rmse above 1 indicates item behaviour the frame units do not account for\n")
  fl <- x$locations[x$locations$flagged, ]
  if (nrow(fl)) {
    cat(sprintf("\nLocation differs across frames for %d item(s) at alpha = %.2f (Holm-adjusted):\n",
                nrow(fl), x$alpha))
    print(.fmt_df(fl[, c("set", "frame_1", "frame_2", "item", "difference",
                         "se", "statistic", "p_adj")]), row.names = FALSE)
  } else {
    cat(sprintf("\nNo item's location differs across frames at alpha = %.2f.\n",
                x$alpha))
  }
  fd <- x$discrimination[x$discrimination$flagged, ]
  if (nrow(fd)) {
    cat(sprintf("\nDiscrimination differs across frames for %d item(s):\n",
                nrow(fd)))
    print(.fmt_df(fd[, c("set", "frame_1", "frame_2", "item", "infit_1",
                         "infit_2", "statistic", "p_adj")]), row.names = FALSE)
  } else {
    cat("\nNo item's discrimination differs across frames.\n")
  }
  invisible(x)
}
