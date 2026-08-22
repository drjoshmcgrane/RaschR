# rasch :: testing the invariance a frame model assumes
# ===========================================================================
# rasch_efrm() constrains each item to a single location shared across the
# frames it appears in, scaled by the frame unit. That constraint is the
# model's substantive claim -- the frames differ in unit, not in what the
# items measure -- and it cannot be tested from the fit itself, because the
# fit imposes it. Testing it means calibrating each frame on its own and
# comparing, which is what this function does.
# ===========================================================================

.item_location_covariance <- function(fit) {
  thr <- fit$thresholds
  V <- fit$est$cov_tau
  if (is.null(V) || !nrow(thr) || nrow(V) != nrow(thr)) return(NULL)
  item_names <- fit$items$item
  A <- matrix(0, length(item_names), nrow(thr),
              dimnames = list(item_names, NULL))
  for (i in seq_along(item_names)) {
    rows <- which(thr$item == match(item_names[i], colnames(fit$X)))
    if (length(rows)) A[i, rows] <- 1 / length(rows)
  }
  A %*% V %*% t(A)
}

#' Test item invariance across frames
#'
#' Calibrates each frame separately and compares the locations and
#' discriminations of items administered in more than one frame.
#'
#' @details
#' Let \eqn{\hat\delta_{if}} be the location of item \eqn{i} from a separate
#' calibration of frame \eqn{f}, and let \eqn{\hat\rho_f} be that frame's
#' unit from the fitted EFRM. The common-scale location is
#' \eqn{\hat\delta_{if}^{*}=\hat\delta_{if}/\hat\rho_f}. Because each
#' separate calibration has its own origin, pairwise differences are centred
#' over the common thresholds before testing.
#'
#' The conditional method treats the fitted frame units as fixed.
#' Let \eqn{w_i=m_i/\sum_jm_j}, where \eqn{m_i} is the number of
#' thresholds for item \eqn{i}. If
#' \eqn{C=I-\mathbf{1}\mathbf{w}^{\mathsf T}} centres the common items on the
#' threshold-weighted origin, the covariance of the location differences is
#' \deqn{C\{V_1/\hat\rho_1^2+V_2/\hat\rho_2^2\}C^{\mathsf T}.}
#' This is fast and conditions on the estimated units. The discrimination
#' table gives the difference between the two standardised infit statistics,
#' divided by \eqn{\sqrt{2}}, together with fitted slopes and their ratio.
#' These quantities are descriptive under the conditional method; it does not
#' report discrimination probabilities.
#'
#' With \code{se_method = "bootstrap"}, persons are resampled within frame
#' and the EFRM and separate frame calibrations are refitted. Location tests
#' then use the empirical covariance of the centred differences. The
#' discrimination test uses the bootstrap standard error of the log slope
#' ratio. This includes uncertainty in the fitted frame units but is more
#' computationally demanding.
#'
#' Raw and Holm-adjusted probabilities are reported. With conditional
#' uncertainty, Holm adjustment covers the location comparisons. With
#' bootstrap uncertainty, it covers the combined family of location and
#' discrimination comparisons.
#' The summary gives the root mean squared location difference and root mean
#' squared standard error for each set and frame pair. Items from different
#' sets cannot be compared because the sets partition the items.
#' Location differences are relative to the mean difference of the common
#' items. Concentrated DIF can therefore produce non-zero centred contrasts
#' for items that were not themselves shifted. The table identifies the
#' pattern of relative departures; item content or external anchors are needed
#' to determine which items provide the defensible reference.
#' A compared set-by-frame cell must contain at least 50 persons with two or
#' more responses. Items with weakly determined standard errors in either
#' separate calibration are listed in \code{excluded} rather than tested.
#'
#' A flagged item may be resolved with \code{\link{resolve_frames}} when it
#' remains useful within frames, or removed with \code{\link{drop_items}}
#' when it fits poorly more generally. Either change requires a refit. The
#' invariance tests require a converged frame calibration.
#'
#' @name frame_invariance
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param alpha Significance level used for flags.
#' @param adjust Either \code{"holm"} or \code{"none"}. Both raw and
#'   adjusted probabilities are returned.
#' @param se_method \code{"conditional"} treats the estimated frame units as
#'   fixed; \code{"bootstrap"} refits the complete analysis to person
#'   resamples within frame.
#' @param boot_reps Number of bootstrap replicates. At least 30 are required.
#' @param seed Optional bootstrap seed.
#' @return An object of class \code{"rasch_frame_invariance"}. The
#'   \code{locations} and \code{discrimination} tables contain the pairwise
#'   item comparisons; \code{summary} contains set-level RMSD and RMSE
#'   summaries. Under the conditional method, discrimination \code{p},
#'   \code{p_adj}, and \code{flagged} are \code{NA}. \code{excluded} lists
#'   items whose observed category structures differed between calibrations
#'   or whose separate-frame estimate was weakly determined.
#'   The remaining components record the multiplicity and uncertainty settings.
#' @references
#' Humphry, S. M. (2005). \emph{Maintaining a Common Arbitrary Unit in Social
#' Measurement}. PhD thesis, Murdoch University.
#' @seealso \code{\link{resolve_frames}} to give a flagged item a location
#'   per frame, \code{\link{drop_items}} to remove it altogether, and
#'   \code{\link{rasch_efrm}} for the model whose assumption is tested.
#' @examples
#' d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 1,
#'                    n_groups = 2, group_unit_ratio = 1.4, seed = 2)
#' tr <- attr(d, "truth")
#' fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
#'                   id = "id", boot_reps = 0)
#' frame_invariance(fit)
NULL

.frame_invariance_conditional <- function(fit) {
  if (!isTRUE(fit$est$converged)) return(NULL)
  grp <- .frame_group_values(fit)
  glev <- levels(factor(grp))
  sets <- unique(fit$set_of)
  rho <- fit$frames
  vm <- fit$virtual_map
  spec <- fit$refit_spec
  if (is.null(spec)) spec <- list()
  out <- list(); dsc <- list(); excluded <- list()
  for (s in sets) {
    cal <- list()
    for (g in glev) {
      vr <- vm$set == s & vm$group == g & vm$vkey %in% colnames(fit$X)
      if (!any(vr)) next
      rows <- !is.na(grp) & as.character(grp) == g
      Xg <- fit$X[rows, vm$vkey[vr], drop = FALSE]
      colnames(Xg) <- vm$item[vr]
      Xg <- Xg[, colSums(!is.na(Xg)) > 0L, drop = FALSE]
      if (ncol(Xg) < 2L) next
      category_signature <- vapply(seq_len(ncol(Xg)), function(j)
        paste(sort(unique(Xg[!is.na(Xg[, j]), j])), collapse = ","), "")
      names(category_signature) <- colnames(Xg)
      f <- tryCatch(do.call(rasch, list(
        data = Xg, model = "PCM", n_groups = spec$n_groups,
        adjust_N = spec$adjust_N %||% NA_real_,
        maxit = spec$maxit %||% 50, tol = spec$tol %||% 1e-7)),
        error = function(e) NULL)
      if (is.null(f) || !isTRUE(f$est$converged)) next
      r <- rho$rho[rho$set == s & rho$group == g]
      if (length(r) != 1L || !is.finite(r) || r <= 0) next
      V <- .item_location_covariance(f)
      if (is.null(V) || any(!is.finite(V))) next
      good <- is.finite(f$items$location) & is.finite(f$items$se)
      weak_items <- f$items$item[!good]
      if (sum(good) < 2L) next
      V <- V[good, good, drop = FALSE]
      V <- V / r^2
      cal[[g]] <- list(
        table = data.frame(item = f$items$item[good],
          loc = f$items$location[good] / r, se = f$items$se[good] / r,
          n_thresholds = f$m[good],
          category_signature = unname(category_signature[f$items$item[good]]),
          infit = f$items$infit_ms[good], infit_z = f$items$infit_z[good],
          disc = f$items$disc[good], stringsAsFactors = FALSE),
        covariance = V, weak_items = weak_items)
    }
    if (length(cal) < 2L) next
    gg <- names(cal)
    for (a in seq_len(length(gg) - 1L)) for (b in (a + 1L):length(gg)) {
      weak_pair <- union(cal[[gg[a]]]$weak_items, cal[[gg[b]]]$weak_items)
      if (length(weak_pair)) excluded[[length(excluded) + 1L]] <- data.frame(
        set = s, frame_1 = gg[a], frame_2 = gg[b], item = weak_pair,
        reason = "weakly determined in a separate frame calibration",
        stringsAsFactors = FALSE)
      m <- merge(cal[[gg[a]]]$table, cal[[gg[b]]]$table, by = "item",
                 suffixes = c("_1", "_2"))
      if (!nrow(m)) next
      comparable <- m$n_thresholds_1 == m$n_thresholds_2 &
        m$category_signature_1 == m$category_signature_2
      if (any(!comparable)) {
        excluded[[length(excluded) + 1L]] <- data.frame(
          set = s, frame_1 = gg[a], frame_2 = gg[b],
          item = m$item[!comparable],
          reason = "different observed category structure",
          stringsAsFactors = FALSE)
        m <- m[comparable, , drop = FALSE]
      }
      if (nrow(m) < 2L) next
      raw_d <- m$loc_2 - m$loc_1
      w <- m$n_thresholds_1 / sum(m$n_thresholds_1)
      C <- diag(nrow(m)) - outer(rep(1, nrow(m)), w)
      d <- drop(C %*% raw_d)
      V1 <- cal[[gg[a]]]$covariance
      V2 <- cal[[gg[b]]]$covariance
      V1 <- V1[m$item, m$item, drop = FALSE]
      V2 <- V2[m$item, m$item, drop = FALSE]
      Vd <- C %*% (V1 + V2) %*% C
      se <- sqrt(pmax(diag(Vd), 0))
      z <- d / se
      out[[length(out) + 1L]] <- data.frame(
        set = s, frame_1 = gg[a], frame_2 = gg[b], item = m$item,
        location_1 = m$loc_1, location_2 = m$loc_2,
        difference = d, se = se, statistic = z,
        p = 2 * stats::pnorm(-abs(z)), stringsAsFactors = FALSE)
      zd <- (m$infit_z_1 - m$infit_z_2) / sqrt(2)
      dsc[[length(dsc) + 1L]] <- data.frame(
        set = s, frame_1 = gg[a], frame_2 = gg[b], item = m$item,
        infit_1 = m$infit_1, infit_2 = m$infit_2,
        infit_z = zd, p = 2 * stats::pnorm(-abs(zd)),
        disc_1 = m$disc_1, disc_2 = m$disc_2,
        disc_ratio = m$disc_2 / m$disc_1,
        disc_boundary = m$disc_1 <= 0.050001 | m$disc_1 >= 4.9999 |
          m$disc_2 <= 0.050001 | m$disc_2 >= 4.9999,
        stringsAsFactors = FALSE)
    }
  }
  excluded <- if (length(excluded)) do.call(rbind, excluded) else
    data.frame(set = character(), frame_1 = character(),
               frame_2 = character(), item = character(),
               reason = character())
  if (!length(out)) {
    if (!nrow(excluded)) return(NULL)
    return(list(
      locations = data.frame(
        set = character(), frame_1 = character(), frame_2 = character(),
        item = character(), location_1 = numeric(), location_2 = numeric(),
        difference = numeric(), se = numeric(), statistic = numeric(),
        p = numeric()),
      discrimination = data.frame(
        set = character(), frame_1 = character(), frame_2 = character(),
        item = character(), infit_1 = numeric(), infit_2 = numeric(),
        infit_z = numeric(), p = numeric(), disc_1 = numeric(),
        disc_2 = numeric(), disc_ratio = numeric(),
        disc_boundary = logical()),
      excluded = excluded))
  }
  list(locations = do.call(rbind, out), discrimination = do.call(rbind, dsc),
       excluded = excluded)
}

#' @rdname frame_invariance
#' @export
frame_invariance <- function(fit, alpha = 0.05, adjust = c("holm", "none"),
                             se_method = c("conditional", "bootstrap"),
                             boot_reps = 200, seed = NULL) {
  if (!inherits(fit, "rasch_efrm"))
    stop("frame_invariance needs a fit from rasch_efrm()")
  if (!isTRUE(fit$est$converged))
    stop("the frame calibration did not converge; invariance tests are unavailable")
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must be one number between 0 and 1")
  adjust <- match.arg(adjust)
  se_method <- match.arg(se_method)
  grp <- .frame_group_values(fit)
  glev <- levels(factor(grp))
  if (length(glev) < 2L)
    stop("item invariance across frames needs at least two person groups: ",
         "with one group each item appears in a single frame, and item sets ",
         "partition the items, so use the item fit statistics within each ",
         "set instead")

  # A separate-frame calibration supplies the covariance used by every item
  # comparison. Small frames produced valid-looking but unstable normal tests
  # in simulation, so require 50 informative persons in every observed
  # set-by-frame cell before reporting invariance probabilities.
  vm <- fit$virtual_map
  fr <- unique(vm[, c("set", "group")])
  n_frames_by_set <- table(fr$set)
  sparse <- vapply(seq_len(nrow(fr)), function(i) {
    if (n_frames_by_set[fr$set[i]] < 2L) return(FALSE)
    cc <- which(vm$set == fr$set[i] & vm$group == fr$group[i] &
                  vm$vkey %in% colnames(fit$X))
    if (length(cc) < 2L) return(FALSE)
    sum(rowSums(!is.na(fit$X[, vm$vkey[cc], drop = FALSE])) >= 2L) < 50L
  }, logical(1))
  if (any(sparse)) stop(
    "frame-invariance inference needs at least 50 persons with two or more ",
    "responses in every compared set-by-frame cell; sparse cell(s): ",
    paste(paste(fr$set[sparse], fr$group[sparse], sep = "/"), collapse = ", "))

  ans <- .frame_invariance_conditional(fit)
  if (is.null(ans))
    stop("no item set is taken by two person groups, so no item appears in ",
         "two frames to be compared")
  cmp <- ans$locations
  dsc <- ans$discrimination
  reps_used <- 0L
  if (se_method == "bootstrap") {
    if (length(boot_reps) != 1L || !is.finite(boot_reps) ||
        boot_reps < 30L || boot_reps != floor(boot_reps))
      stop("boot_reps must be a whole number of at least 30")
    if (!is.null(seed)) {
      if (length(seed) != 1L || !is.finite(seed)) stop("seed must be one number")
      old_seed <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        get(".Random.seed", .GlobalEnv) else NULL
      on.exit(if (is.null(old_seed)) {
        if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
      } else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
      set.seed(seed)
    }
    source <- .efrm_source_matrix(fit)
    strata <- split(seq_len(nrow(source)), as.character(grp), drop = TRUE)
    key4 <- function(x) .factor_keys(
      x[, c("set", "frame_1", "frame_2", "item"), drop = FALSE])
    lk <- key4(cmp)
    dk <- key4(dsc)
    bd <- matrix(NA_real_, boot_reps, nrow(cmp))
    ba <- matrix(NA_real_, boot_reps, nrow(dsc))
    for (b in seq_len(boot_reps)) {
      ii <- unlist(lapply(strata, sample, replace = TRUE), use.names = FALSE)
      fb <- tryCatch(.efrm_refit(
        fit, source[ii, , drop = FALSE], fit$set_of, boot_reps = 0,
        ids = sprintf("B%04dP%06d", b, seq_along(ii)),
        factors = fit$factors[ii, , drop = FALSE], se_method = "hybrid"),
        error = function(e) NULL)
      if (is.null(fb)) next
      if (!isTRUE(fb$est$converged) ||
          any(fb$linking$alpha_edges$converged %in% FALSE)) next
      ib <- .frame_invariance_conditional(fb)
      if (is.null(ib)) next
      xk <- key4(ib$locations)
      yk <- key4(ib$discrimination)
      bd[b, ] <- ib$locations$difference[match(lk, xk)]
      ba[b, ] <- log(ib$discrimination$disc_ratio[match(dk, yk)])
    }
    good <- rowSums(is.finite(bd)) == ncol(bd) &
      rowSums(is.finite(ba)) == ncol(ba)
    reps_used <- sum(good)
    if (reps_used < max(30L, ceiling(0.8 * boot_reps)))
      stop("only ", reps_used, " of ", boot_reps,
           " frame-invariance bootstrap refits succeeded")
    cmp$se <- apply(bd[good, , drop = FALSE], 2, stats::sd)
    cmp$statistic <- cmp$difference / cmp$se
    cmp$p <- 2 * stats::pnorm(-abs(cmp$statistic))
    dsc$log_disc_ratio <- log(dsc$disc_ratio)
    dsc$se_log_disc_ratio <- apply(ba[good, , drop = FALSE], 2, stats::sd)
    dsc$statistic <- dsc$log_disc_ratio / dsc$se_log_disc_ratio
    dsc$p <- 2 * stats::pnorm(-abs(dsc$statistic))
  } else {
    # The standardised-infit comparison was anti-conservative in the
    # validation study (7.1% combined Holm FWER over 2,000 null replicates,
    # with strong item-position dependence). Retain its descriptive value,
    # but do not attach an inferential probability to it.
    dsc$p <- NA_real_
  }
  # Conditional inference covers locations only. The validated bootstrap
  # adds discrimination and controls the two tables as one family.
  allp <- c(cmp$p, dsc$p)
  padj <- rep(NA_real_, length(allp))
  usable <- is.finite(allp)
  padj[usable] <- stats::p.adjust(allp[usable], method = "holm")
  cmp$p_adj <- padj[seq_len(nrow(cmp))]
  dsc$p_adj <- padj[nrow(cmp) + seq_len(nrow(dsc))]
  p_cmp <- if (adjust == "holm") cmp$p_adj else cmp$p
  cmp$flagged <- ifelse(is.finite(p_cmp), p_cmp < alpha, NA)
  rownames(cmp) <- NULL
  p_dsc <- if (adjust == "holm") dsc$p_adj else dsc$p
  dsc$flagged <- ifelse(is.finite(p_dsc), p_dsc < alpha, NA)
  rownames(dsc) <- NULL

  frame_pairs <- unique(cmp[c("set", "frame_1", "frame_2")])
  smry <- lapply(seq_len(nrow(frame_pairs)), function(i) {
    k <- frame_pairs[i, ]
    same <- function(d) d$set == k$set & d$frame_1 == k$frame_1 &
      d$frame_2 == k$frame_2
    z <- cmp[same(cmp), ]; y <- dsc[same(dsc), ]
    nx <- sum(same(ans$excluded))
    data.frame(set = z$set[1], frame_1 = z$frame_1[1], frame_2 = z$frame_2[1],
      n_items = nrow(z), n_excluded = nx,
      rmsd = sqrt(mean(z$difference^2)),
      rmse = sqrt(mean(z$se^2)),
      ratio = sqrt(mean(z$difference^2)) / sqrt(mean(z$se^2)),
      n_location = sum(z$flagged %in% TRUE),
      n_discrimination = if (identical(se_method, "bootstrap"))
        sum(y$flagged %in% TRUE) else NA_integer_,
      stringsAsFactors = FALSE)
  })
  smry <- if (length(smry)) do.call(rbind, smry) else data.frame(
    set = character(), frame_1 = character(), frame_2 = character(),
    n_items = integer(), n_excluded = integer(), rmsd = numeric(),
    rmse = numeric(), ratio = numeric(), n_location = integer(),
    n_discrimination = integer())
  rownames(smry) <- NULL
  structure(.tag_tables(list(locations = cmp, discrimination = dsc,
                             summary = smry, excluded = ans$excluded,
                             alpha = alpha, adjust = adjust,
                             se_method = se_method,
                             boot_reps_used = reps_used)),
            class = "rasch_frame_invariance")
}

#' @export
print.rasch_frame_invariance <- function(x, ...) {
  adj <- if (is.null(x$adjust)) "holm" else x$adjust
  pcol <- if (adj == "holm") "p_adj" else "p"
  rule <- if (adj == "holm") "Holm-adjusted" else "unadjusted, screening"
  cat("Item invariance across frames (each frame calibrated separately)\n\n")
  cat("Uncertainty:", if (identical(x$se_method, "bootstrap"))
    sprintf("person-within-frame bootstrap (%d successful refits)",
            x$boot_reps_used) else "conditional on the fitted frame units",
    "\n\n")
  print(.fmt_df(x$summary), row.names = FALSE)
  if (!is.null(x$excluded) && nrow(x$excluded))
    cat(sprintf("\n%d item comparison(s) excluded because the observed category structure differed between frames.\n",
                nrow(x$excluded)))
  cat("\nrmsd/rmse above 1 indicates item behaviour the frame units do not account for\n")
  fl <- x$locations[x$locations$flagged %in% TRUE, , drop = FALSE]
  if (nrow(fl)) {
    cat(sprintf("\nLocation differs across frames for %d item(s) at alpha = %.2f (%s):\n",
                nrow(fl), x$alpha, rule))
    print(.fmt_df(fl[, c("set", "frame_1", "frame_2", "item", "difference",
                         "se", "statistic", pcol)]), row.names = FALSE)
  } else {
    cat(sprintf("\nNo item's location differs across frames at alpha = %.2f (%s).\n",
                x$alpha, rule))
  }
  if (!identical(x$se_method, "bootstrap")) {
    cat("\nThe discrimination comparisons are descriptive:\n")
    cols <- c("set", "frame_1", "frame_2", "item", "infit_1", "infit_2",
              "infit_z", "disc_1", "disc_2", "disc_ratio", "disc_boundary")
    print(.fmt_df(x$discrimination[, intersect(cols,
                                                names(x$discrimination)),
                                           drop = FALSE]), row.names = FALSE)
    cat("Use se_method = \"bootstrap\" for discrimination probabilities.\n")
    return(invisible(x))
  }
  fd <- x$discrimination[x$discrimination$flagged %in% TRUE, , drop = FALSE]
  if (nrow(fd)) {
    cat(sprintf("\nDiscrimination differs across frames for %d item(s):\n",
                nrow(fd)))
    cols <- c("set", "frame_1", "frame_2", "item", "log_disc_ratio",
              "se_log_disc_ratio", "statistic", pcol, "disc_1", "disc_2",
              "disc_ratio", "disc_boundary")
    print(.fmt_df(fd[, intersect(cols, names(fd)), drop = FALSE]),
          row.names = FALSE)
  } else {
    cat("\nNo item's discrimination differs across frames.\n")
  }
  invisible(x)
}
