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
#' false-positive rates near 1\% throughout.
#'
#' The two comparisons are therefore not equally sensitive, and the gap
#' matters because the two departures do comparable damage. At departures
#' that each move a unit ratio by six or seven per cent -- two items
#' shifted a logit, against two items discriminating half again as steeply
#' -- 500 persons per frame detected the shifted items about 97\% of the
#' time and the steeper items about 16\%. Only by around 2,000 persons per
#' frame do the two reach comparable power. A clean result at a few hundred
#' persons per frame is thus much stronger evidence against differential
#' item functioning than against differential discrimination: it has ruled
#' out one departure and barely tested the other.
#'
#' The discrimination table also reports a Winsteps-style index for each
#' frame (\code{disc_1}, \code{disc_2}) and their ratio, because a
#' standardised difference says only that something differs while the index
#' says how much and in which direction. It is fitted by maximum likelihood
#' on each item's own responses with the person measures and item location
#' held at their Rasch values, so it is relative to the frame's own model
#' and the unit cancels. Read it as description rather than estimate: the
#' measures are estimated including the item being scored, which biased the
#' index to about 1.19 for a true discrimination of 1.0 in simulation and
#' attenuated a true frame ratio of 1.5 to about 1.22. Tested on its own it
#' is also the weaker instrument, detecting a 1.5-fold difference in 63\%
#' of replicates at 2,000 persons per frame against the infit comparison's
#' 90\%, which is why the test column comes from the latter.
#'
#' Which item to flag is a screening decision, not a confirmatory one, and
#' the two call for different thresholds. \code{adjust} chooses: Holm across
#' every item and frame pair, or none. Both probabilities are reported
#' either way, so the choice changes only \code{flagged}. Screening 8 to 10
#' items with Holm costs between 20 and 60 points of sensitivity in
#' simulation, and that shows up in the repair: dropping the flagged items
#' from a planted unit ratio of 1.40 left the ratio at 1.479 under Holm and
#' 1.433 unadjusted, against 1.406 for dropping the items actually planted.
#' The loose screen is not free -- where misfit is strong it flags sound
#' items too, and \code{\link{drop_items}} then refuses drops that would
#' empty a set -- so use \code{"none"} to decide which items to examine and
#' \code{"holm"} to report which ones differ.
#'
#' Dropping is a complete cure where the item is found: removing the
#' planted items restored the ratio in every simulated departure, so what
#' limits the repair is detection rather than removal. Past roughly a fifth
#' of the items breaking invariance, no threshold rescues the ratio and the
#' item set itself is the problem.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param alpha Significance level for flagging items.
#' @param adjust Multiplicity adjustment used to flag items: \code{"holm"}
#'   across all comparisons, or \code{"none"} for screening. Both
#'   probabilities are reported regardless.
#' @return A list of class \code{"rasch_frame_invariance"} with
#'   \code{locations} (one row per item and frame pair: locations on the
#'   common scale, their difference, its standard error, statistic, both
#'   probabilities, and flag), \code{discrimination} (the same items compared
#'   on their within-frame infit statistics, with the Winsteps-style
#'   discrimination index for each frame and its ratio alongside), and
#'   \code{summary} (per frame pair: the number of items, the root mean
#'   squared difference, the root mean squared standard error, their ratio,
#'   and the number of items flagged on each count).
#' @references
#' Humphry, S. M. (2005). \emph{Maintaining a Common Arbitrary Unit in Social
#' Measurement}. PhD thesis, Murdoch University.
#' @seealso \code{\link{drop_items}} to remove an item the test flags, and
#'   \code{\link{rasch_efrm}} for the model whose assumption is tested.
#' @export
#' @examples
#' d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 1,
#'                    n_groups = 2, group_unit_ratio = 1.4, seed = 2)
#' tr <- attr(d, "truth")
#' fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
#'                   id = "id", boot_reps = 0)
#' frame_invariance(fit)

frame_invariance <- function(fit, alpha = 0.05, adjust = c("holm", "none")) {
  if (!inherits(fit, "rasch_efrm"))
    stop("frame_invariance needs a fit from rasch_efrm()")
  adjust <- match.arg(adjust)
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
                             disc = .winsteps_disc(f),
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
        infit_z = zd, p = 2 * stats::pnorm(-abs(zd)),
        disc_1 = m$disc_1, disc_2 = m$disc_2,
        disc_ratio = m$disc_2 / m$disc_1,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out))
    stop("no item set is taken by two person groups, so no item appears in ",
         "two frames to be compared")
  # both probabilities are always reported; adjust chooses which one flags
  cmp <- do.call(rbind, out)
  cmp$p_adj <- stats::p.adjust(cmp$p, method = "holm")
  cmp$flagged <- (if (adjust == "holm") cmp$p_adj else cmp$p) < alpha
  rownames(cmp) <- NULL
  dsc <- do.call(rbind, dsc)
  dsc$p_adj <- stats::p.adjust(dsc$p, method = "holm")
  dsc$flagged <- (if (adjust == "holm") dsc$p_adj else dsc$p) < alpha
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
                 alpha = alpha, adjust = adjust),
            class = "rasch_frame_invariance")
}

#' @export
print.rasch_frame_invariance <- function(x, ...) {
  adj <- if (is.null(x$adjust)) "holm" else x$adjust
  pcol <- if (adj == "holm") "p_adj" else "p"
  rule <- if (adj == "holm") "Holm-adjusted" else "unadjusted, screening"
  cat("Item invariance across frames (each frame calibrated separately)\n\n")
  print(.fmt_df(x$summary), row.names = FALSE)
  cat("\nrmsd/rmse above 1 indicates item behaviour the frame units do not account for\n")
  fl <- x$locations[x$locations$flagged, ]
  if (nrow(fl)) {
    cat(sprintf("\nLocation differs across frames for %d item(s) at alpha = %.2f (%s):\n",
                nrow(fl), x$alpha, rule))
    print(.fmt_df(fl[, c("set", "frame_1", "frame_2", "item", "difference",
                         "se", "statistic", pcol)]), row.names = FALSE)
  } else {
    cat(sprintf("\nNo item's location differs across frames at alpha = %.2f (%s).\n",
                x$alpha, rule))
  }
  fd <- x$discrimination[x$discrimination$flagged, ]
  if (nrow(fd)) {
    cat(sprintf("\nDiscrimination differs across frames for %d item(s):\n",
                nrow(fd)))
    print(.fmt_df(fd[, c("set", "frame_1", "frame_2", "item", "infit_1",
                         "infit_2", "infit_z", pcol, "disc_1", "disc_2",
                         "disc_ratio")]), row.names = FALSE)
    cat("The test is the infit comparison (infit_z, ", pcol, "). The disc\n",
        "columns describe size and direction only: they run high, and their\n",
        "ratio understates the difference.\n", sep = "")
  } else {
    cat("\nNo item's discrimination differs across frames.\n")
  }
  invisible(x)
}


# Winsteps-style item discrimination: the slope fitted by maximum
# likelihood on the item's own responses with the person measures and the
# item location held at their Rasch values. The frame unit is absorbed
# into those fixed measures, so the index is relative to the frame's own
# model. It is descriptive, not a two-parameter estimate: the measures are
# estimated including the item being scored, which biases the index upward
# (about 1.19 for a true 1.0 in simulation) and attenuates a ratio between
# frames (1.5 recovered as about 1.22).
.winsteps_disc <- function(f) {
  ok <- !f$person$extreme
  th <- f$person$theta[ok]
  X <- as.matrix(f$X)[ok, , drop = FALSE]
  d <- f$items$location
  vapply(seq_len(ncol(X)), function(i) {
    y <- X[, i]; g <- is.finite(y) & is.finite(th)
    if (length(unique(y[g])) < 2L) return(NA_real_)
    z <- th[g] - d[i]
    m <- tryCatch(suppressWarnings(
      stats::glm(y[g] ~ 0 + z, family = stats::binomial)),
      error = function(e) NULL)
    if (is.null(m)) NA_real_ else unname(stats::coef(m))
  }, 0)
}
