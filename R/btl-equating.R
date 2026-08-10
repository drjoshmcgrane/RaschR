# rasch :: common-object equating for paired comparisons
# ===========================================================================
# The paired-comparison analogue of equate_tests(). Two Bradley-Terry-Luce
# calibrations that share a set of common objects -- the same scripts,
# performances, or products judged by different panels, or by the same panel
# in different years -- each fix their own origin by the sum-zero constraint.
# Because the two constraints are imposed over DIFFERENT object sets, the two
# origins do not coincide even when the objects are unchanged: each scale is
# centred on the mean of a different collection. Equating therefore estimates
# the scale shift between the origins (the precision-weighted mean difference
# over the common objects) and then tests each common object against the
# shifted identity line. Its variance includes both the covariance induced by
# each sum-zero calibration and estimation of the shift from those same
# objects; it is not the naive sum of two marginal variances.
#
# Objects that survive define the equating link and carry the second panel's
# whole scale onto the first; objects that fail show drift -- a script the two
# panels valued differently, or a standard that moved between years -- and
# should be reviewed before the link is trusted. This is the standards-
# maintenance use of comparative judgement (Bramley 2007): a common set of
# anchor scripts lets panels judged apart be placed on one scale.
# ===========================================================================

# Coerce the second calibration to (object, location, se). It may be another
# btl fit or a "bank" -- a data frame of previously banked object locations,
# the paired-comparison counterpart of an item bank.
.btl_equate_ref <- function(reference) {
  if (inherits(reference, "rasch_btl"))
    return(data.frame(object = as.character(reference$objects$object),
                      location = reference$objects$location,
                      se = reference$objects$se,
                      stringsAsFactors = FALSE))
  if (is.data.frame(reference) || (is.list(reference) && !is.null(names(reference)))) {
    reference <- as.data.frame(reference, stringsAsFactors = FALSE)
    if (!all(c("object", "location") %in% names(reference)))
      stop("a bank needs columns 'object' and 'location' (and ideally 'se')")
    se_supplied <- "se" %in% names(reference)
    if (!"se" %in% names(reference)) reference$se <- NA_real_
    out <- data.frame(object = trimws(as.character(reference$object)),
                      location = suppressWarnings(as.numeric(reference$location)),
                      se = suppressWarnings(as.numeric(reference$se)),
                      stringsAsFactors = FALSE)
    if (anyNA(out$object) || any(!nzchar(out$object)))
      stop("bank object names must be non-missing and non-empty")
    if (anyDuplicated(out$object))
      stop("bank object names must be unique: ",
           paste(unique(out$object[duplicated(out$object)]), collapse = ", "))
    if (any(!is.finite(out$location)))
      stop("bank object locations must be finite")
    if (any(!is.na(out$se) & (!is.finite(out$se) | out$se < 0)))
      stop("bank standard errors must be non-negative finite values or NA")
    attr(out, "se_supplied") <- se_supplied
    return(out)
  }
  stop("`fit2` must be a btl fit or a bank data frame (object, location, se)")
}

.btl_equate_bank_cov <- function(reference, ids) {
  C <- attr(reference, "cov_location", exact = TRUE)
  if (is.null(C)) return(NULL)
  if (!is.matrix(C) || !is.numeric(C) || any(!is.finite(C)) ||
      !identical(dim(C), c(length(ids), length(ids))))
    stop("attr(fit2, 'cov_location') must be a finite numeric square matrix ",
         "with one row and column per bank object")
  if (!isTRUE(all.equal(C, t(C), tolerance = 1e-8)))
    stop("attr(fit2, 'cov_location') must be symmetric")
  if (!is.null(rownames(C)) || !is.null(colnames(C))) {
    if (is.null(rownames(C)) || is.null(colnames(C)) ||
        anyNA(match(ids, rownames(C))) || anyNA(match(ids, colnames(C))))
      stop("named bank covariance rows and columns must match every bank object")
    C <- C[ids, ids, drop = FALSE]
  }
  ev <- eigen((C + t(C)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -1e-8 * max(1, max(abs(ev))))
    stop("attr(fit2, 'cov_location') must be positive semidefinite")
  C
}

#' Equate two paired-comparison calibrations through their common objects
#'
#' Compares the object locations of a Bradley-Terry-Luce fit with those of a
#' second fit (or a banked table of object locations), matched by object name.
#' The use case is standards maintenance and comparative judgement across
#' panels or years: the same scripts, performances, or products are judged by
#' different panels -- or by one panel in successive years -- and a common set
#' of anchor objects is carried through so that the two rounds land on a single
#' scale.
#'
#' Each calibration is identified by the sum-zero constraint, but the two
#' constraints are imposed over \emph{different} object sets, so the origins do
#' not coincide even when the shared objects are unchanged: each scale is
#' centred on the mean of a different collection. A scale shift between the two
#' origins is therefore estimated by the precision-weighted mean difference
#' over the common objects, and each common object is then tested against the
#' shifted identity line. A flagged object shows drift -- a script the two
#' panels valued differently, or a standard that moved between years -- and
#' weakens the equating link; the surviving objects carry the second
#' calibration's whole scale onto the first (\code{loc2 + shift}).
#'
#' The single shift presumes the drifting objects are a \emph{minority}. When
#' most of the common objects have genuinely moved, the precision-weighted
#' shift is pulled toward the movers and the drift tests can invert -- the
#' stable anchors flag as the apparent drifters. Read wholesale flagging
#' (several objects, one direction) as a contaminated link, not as evidence
#' about the individual objects; equate through a vetted anchor subset
#' instead. The \code{shift_se} accounts for the covariance of the location
#' estimates within each calibration (each is sum-zero constrained, so its
#' locations are not independent).
#'
#' @param fit1 A fitted object from \code{\link{btl}}: the calibration whose
#'   scale (origin) the equating targets.
#' @param fit2 A second \code{\link{btl}} fit, or a bank: a data frame with
#'   columns \code{object}, \code{location}, and optionally \code{se}; object
#'   names must be unique and locations finite. Bank-based drift inference
#'   requires the joint location covariance as a square matrix in
#'   \code{attr(fit2, "cov_location")}, ordered like the bank rows (or named by
#'   object), unless the bank is treated as fixed with zero SEs. For a graded
#'   fit the bank must carry
#'   \code{attr(bank, "m")} matching the number of fitted score steps.
#' @param alpha Significance level for the (multiplicity-adjusted) drift tests.
#' @param p_adjust Adjustment for the common-object tests, passed to
#'   \code{stats::p.adjust}. The default is \code{"holm"}.
#' @param independent Whether the calibrations have independent judges and
#'   comparisons. For two fitted objects the default \code{NULL} withholds
#'   drift tests until independence is stated explicitly. Bank tables are
#'   treated as independent unless \code{FALSE} is supplied. Dependent
#'   calibrations require a joint or paired bootstrap for inference.
#' @return A list of class \code{"rasch_btl_equate"}: the comparison
#'   \code{table} (per common object: object, both locations and standard
#'   errors, their \code{difference}, the \code{shifted_difference} against the
#'   estimated origin, the pooled \code{se_diff}, \code{t}, raw and adjusted
#'   \code{p}, and the \code{drifting} flag); the estimated \code{shift} and
#'   its \code{shift_se}; \code{equated}, the second calibration's full object
#'   table re-expressed on \code{fit1}'s scale; the number of common objects
#'   \code{n_common}; the number usable for inference \code{n_inference};
#'   whether inference was available \code{inferential}; \code{alpha};
#'   \code{p_adjust}; and \code{notes}.
#' @references Bramley, T. (2007). Paired comparison methods. In P. Newton,
#'   J. Baird, H. Goldstein, H. Patrick, & P. Tymms (Eds.), \emph{Techniques
#'   for monitoring the comparability of examination standards} (pp. 246-294).
#'   London: Qualifications and Curriculum Authority.
#' @examples
#' set.seed(1)
#' beta <- setNames(seq(-2, 2, length.out = 8), paste0("O", 1:8))
#' sim <- function(objs) {
#'   pr <- t(utils::combn(objs, 2))
#'   d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#'   d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#'   btl(d, "a", "b", "win")
#' }
#' eq <- btl_equate(sim(paste0("O", 1:7)), sim(paste0("O", 2:8)),
#'                   independent = TRUE)
#' eq$table
#' @export
btl_equate <- function(fit1, fit2, alpha = 0.05, p_adjust = "holm",
                       independent = NULL) {
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must be one probability strictly between 0 and 1")
  if (length(p_adjust) != 1L || !p_adjust %in% stats::p.adjust.methods)
    stop("p_adjust must name a method in stats::p.adjust.methods")
  if (!is.null(independent) && (length(independent) != 1L ||
      is.na(independent) || !is.logical(independent)))
    stop("independent must be NULL, TRUE, or FALSE")
  if (!inherits(fit1, "rasch_btl"))
    stop("`fit1` must be a paired-comparison (btl) fit")
  # equating a non-converged calibration carries its unidentified locations
  # and understated covariance straight into the drift table: refuse it
  for (nm in c("fit1", "fit2")) {
    f <- get(nm)
    if (inherits(f, "rasch_btl") && isFALSE(f$converged))
      stop("`", nm, "` did not converge (its comparison design does not ",
           "identify some object locations); resolve that before equating -- ",
           "the drift statistics would inherit boundary estimates and ",
           "understated standard errors", call. = FALSE)
  }
  if (inherits(fit2, "rasch_btl")) {
    if (!identical(fit1$m, fit2$m) ||
        !identical(as.character(fit1$categories),
                   as.character(fit2$categories)) ||
        !identical(fit1$thr_structure, fit2$thr_structure))
      stop("the paired-comparison calibrations use incompatible response ",
           "scales or threshold structures; shift-only equating requires ",
           "the same model and category scale")
  } else if (fit1$m > 1L) {
    fm <- attr(fit2, "m", exact = TRUE)
    if (is.null(fm) || !identical(as.integer(fm), as.integer(fit1$m)))
      stop("a graded-comparison bank must carry attr(bank, 'm') matching ",
           "the fitted number of score steps; otherwise scale compatibility ",
           "cannot be established")
  }
  cur <- data.frame(object = as.character(fit1$objects$object),
                    location = fit1$objects$location,
                    se = fit1$objects$se, stringsAsFactors = FALSE)
  ref <- .btl_equate_ref(fit2)
  common <- intersect(cur$object, ref$object)
  if (length(common) < 3)
    stop("need at least three common objects to equate paired-comparison scales")
  a <- cur[match(common, cur$object), ]
  b <- ref[match(common, ref$object), ]
  bank_cov <- if (inherits(fit2, "rasch_btl")) NULL else
    .btl_equate_bank_cov(fit2, ref$object)
  if (!is.null(bank_cov)) {
    cov_se <- sqrt(pmax(diag(bank_cov), 0))
    stated <- is.finite(ref$se)
    if (any(stated & abs(ref$se - cov_se) >
            1e-6 * pmax(1, ref$se, cov_se)))
      stop("the bank standard errors must agree with the diagonal of ",
           "attr(fit2, 'cov_location')")
    if (!isTRUE(attr(ref, "se_supplied", exact = TRUE)))
      ref$se <- cov_se
    b <- ref[match(common, ref$object), ]
  }
  d <- a$location - b$location
  v <- a$se^2 + b$se^2
  usable <- is.finite(d) & is.finite(v)
  independent_ok <- if (is.null(independent)) !inherits(fit2, "rasch_btl")
                    else isTRUE(independent)
  joint_cov_ok <- inherits(fit2, "rasch_btl") || !is.null(bank_cov) ||
    all(b$se[usable] == 0)
  inferential <- independent_ok && sum(usable) >= 3L && joint_cov_ok
  w <- if (sum(usable) >= 3L) 1 / pmax(v[usable], 1e-10) else numeric(0)
  # precision-weighted mean difference: the shift between the two sum-zero
  # origins, best estimated where both calibrations are most certain
  c0 <- if (length(w)) sum(w * d[usable]) / sum(w)
        else mean(d[is.finite(d)])
  # the common objects' location estimates are CORRELATED within each
  # sum-zero calibration, so Var(c0) = u' (Sigma1 + Sigma2) u with
  # u = w / sum(w), taken from the stored sandwich covariances. A bank
  # contributes an explicitly attached joint covariance; without it,
  # inference is withheld unless the bank is fixed (all SEs zero).
  covsub <- function(fit, objs_c) {
    if (inherits(fit, "rasch_btl") && !is.null(fit$cov_beta)) {
      i <- match(objs_c, as.character(fit$objects$object))
      fit$cov_beta[i, i, drop = FALSE]
    } else {
      C <- .btl_equate_bank_cov(fit, ref$object)
      if (is.null(C)) NULL else {
        ii <- match(objs_c, ref$object)
        C[ii, ii, drop = FALSE]
      }
    }
  }
  shift_se <- NA_real_; se_diff <- t <- p <- p_adj <-
    rep(NA_real_, length(common)); drifting <- rep(NA, length(common))
  if (inferential) {
    u <- w / sum(w)
    S1all <- covsub(fit1, common)
    S1 <- if (is.null(S1all)) diag(a$se[usable]^2, sum(usable))
          else S1all[usable, usable, drop = FALSE]
    S2all <- covsub(fit2, common)
    S2 <- if (is.null(S2all)) diag(b$se[usable]^2, sum(usable))
          else S2all[usable, usable, drop = FALSE]
    shift_se <- sqrt(max(drop(t(u) %*% (S1 + S2) %*% u), 0))
  # each drift test compares d_i - c0, and c0 is estimated from the SAME
  # common objects: Var(d_i - c0) = [(I - 1u') Sigma (I - u 1')]_ii
  #   = Sigma_ii - 2 (Sigma u)_i + u' Sigma u,
  # not the naive Sigma_ii -- ignoring the estimated shift (and the
  # within-calibration covariance) mis-states every drift p-value
    Sg <- S1 + S2
    Su <- drop(Sg %*% u)
    var_d <- pmax(diag(Sg) - 2 * Su + drop(t(u) %*% Su), 1e-10)
    se_diff[usable] <- sqrt(var_d)
    t[usable] <- (d[usable] - c0) / se_diff[usable]
    p[usable] <- 2 * pnorm(-abs(t[usable]))
    p_adj[usable] <- p.adjust(p[usable], method = p_adjust)
    drifting[usable] <- p_adj[usable] < alpha
  }
  tab <- data.frame(object = common,
                    location_1 = a$location, se_1 = a$se,
                    location_2 = b$location, se_2 = b$se,
                    difference = d, shifted_difference = d - c0,
                    se_diff = se_diff, t = t, p = p, p_adj = p_adj,
                    drifting = drifting, stringsAsFactors = FALSE)
  rownames(tab) <- NULL
  # the second calibration, whole, carried onto fit1's scale
  equated <- data.frame(object = ref$object,
                        location = ref$location + c0,
                        se = ref$se, stringsAsFactors = FALSE)
  rownames(equated) <- NULL
  notes <- sprintf(paste0("Origins differ because each calibration is sum-zero ",
                          "over its own object set; a shift of %.3f logits ",
                          "aligns fit2 to fit1."), c0)
  if (any(!usable))
    notes <- c(notes, paste0(
      "Drift tests unavailable for objects without standard errors: ",
      paste(common[!usable], collapse = ", "), "."))
  if (is.null(independent) && inherits(fit2, "rasch_btl"))
    notes <- c(notes, paste(
      "Drift tests withheld because independence between fitted",
      "calibrations was not stated; set independent = TRUE only for",
      "independent judges and comparisons."))
  if (identical(independent, FALSE))
    notes <- c(notes, paste(
      "Drift tests withheld for dependent calibrations because cross-fit",
      "covariance is unavailable; use a joint or paired bootstrap."))
  if (independent_ok && sum(usable) >= 3L && !joint_cov_ok)
    notes <- c(notes, paste(
      "Drift tests withheld because a bank with non-zero marginal SEs",
      "needs its joint object-location covariance in",
      "attr(fit2, 'cov_location'); marginal SEs do not carry the",
      "calibration-origin covariance."))
  if (independent_ok && sum(usable) < 3L)
    notes <- c(notes, paste(
      "Drift tests withheld because at least three common objects with",
      "standard errors are required."))
  if (any(drifting %in% TRUE))
    notes <- c(notes, sprintf(
      "%d common object(s) drift beyond the shifted link: %s",
      sum(drifting %in% TRUE), paste(common[drifting %in% TRUE], collapse = ", ")))
  structure(class = "rasch_btl_equate",
            list(table = tab, shift = c0, shift_se = shift_se,
                 equated = equated, n_common = length(common),
                 n_inference = sum(usable), inferential = inferential,
                 alpha = alpha, p_adjust = p_adjust, notes = notes))
}

#' Plot a paired-comparison equating comparison
#'
#' Scatter of the two calibrations' common-object locations with the shifted
#' identity line, per-object 95 per cent error bars, and a dotted guide band
#' at the average pooled precision; objects that drift (after
#' the multiplicity adjustment) are highlighted and labelled. The counterpart
#' of \code{\link{plot_equate}} for Bradley-Terry-Luce scales.
#'
#' @param fit1 A fitted object from \code{\link{btl}}.
#' @param fit2 A second \code{\link{btl}} fit, or a bank data frame with columns
#'   \code{object}, \code{location}, and optionally \code{se}.
#' @param ... Passed to \code{\link{btl_equate}} (e.g. \code{alpha},
#'   \code{p_adjust}).
#' @return Called for its plotting side effect; invisibly the
#'   \code{\link{btl_equate}} result.
#' @examples
#' set.seed(1)
#' beta <- setNames(seq(-2, 2, length.out = 8), paste0("O", 1:8))
#' sim <- function(objs) {
#'   pr <- t(utils::combn(objs, 2))
#'   d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#'   d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#'   btl(d, "a", "b", "win")
#' }
#' plot_btl_equate(sim(paste0("O", 1:7)), sim(paste0("O", 2:8)),
#'                  independent = TRUE)
#' @export
plot_btl_equate <- function(fit1, fit2, ...) {
  eq <- btl_equate(fit1, fit2, ...)
  tab <- eq$table
  rng <- range(c(tab$location_1, tab$location_2)) + c(-0.4, 0.4)
  op <- .rr_canvas(rng, rng, "Calibration 2 location (logits)",
                   "Calibration 1 location (logits)",
                   sprintf("%d common objects, shift %.3f, r = %.3f",
                           eq$n_common, eq$shift,
                           stats::cor(tab$location_1, tab$location_2)),
                   grid_x = TRUE)
  on.exit(par(op))
  abline(eq$shift, 1, col = .rr$ink, lwd = 2)
  band <- 1.96 * sqrt(mean(tab$se_1^2 + tab$se_2^2, na.rm = TRUE))
  if (is.finite(band)) {
    abline(eq$shift + band, 1, lty = 3, col = .rr$soft)
    abline(eq$shift - band, 1, lty = 3, col = .rr$soft)
  }
  hs <- is.finite(tab$se_1)
  segments(tab$location_2[hs], tab$location_1[hs] - 1.96 * tab$se_1[hs],
           tab$location_2[hs], tab$location_1[hs] + 1.96 * tab$se_1[hs],
           col = paste0(.rr$soft, "88"))
  points(tab$location_2, tab$location_1, pch = 21, cex = 1.6,
         bg = ifelse(tab$drifting %in% TRUE, .rr$red, .rr$blue),
         col = "white", lwd = 1.2)
  dr <- tab$drifting %in% TRUE
  if (any(dr))
    text(tab$location_2[dr], tab$location_1[dr],
         tab$object[dr], pos = 3, offset = 0.5, cex = 0.75,
         col = .rr$red)
  invisible(eq)
}

#' @export
print.rasch_btl_equate <- function(x, ...) {
  tab <- x$table
  cat(sprintf(paste0("Common-object equating over %d object(s): shift %.3f ",
                     "(SE %s), correlation %.3f, RMSD %.3f\n"),
              x$n_common, x$shift,
              if (is.finite(x$shift_se)) sprintf("%.3f", x$shift_se) else "withheld",
              stats::cor(tab$location_1, tab$location_2),
              sqrt(mean(tab$shifted_difference^2))))
  core <- c("object", "location_1", "location_2", "shifted_difference", "t",
            "p_adj", "drifting")
  print(.fmt_df(tab[, intersect(core, names(tab))]), row.names = FALSE)
  if (isTRUE(x$inferential))
    cat(sprintf("%d object(s) drift beyond the %s-adjusted %.0f%% level.\n",
                sum(tab$drifting %in% TRUE), x$p_adjust, 100 * (1 - x$alpha)))
  else cat("Drift inference withheld; see $notes.\n")
  cat("(standard errors and unadjusted columns on $table; fit2 on fit1's scale in $equated)\n")
  invisible(x)
}
