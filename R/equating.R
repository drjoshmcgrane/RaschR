# rasch :: test equating
# ===========================================================================
# Comparison of common-item locations across two separately analysed
# datasets. Because each analysis fixes its own origin, the comparison
# allows for a scale shift (estimated by precision-weighted mean difference)
# and, when the calibrations are independent, tests each common item against
# the shifted identity line. The variance includes estimation of the shift
# from the same correlated item locations; it is not the naive sum of two
# marginal variances.
#
# Items that survive define the equating link; items that fail show item
# drift and should be dropped from the link (or anchored individually).
# ===========================================================================

# Item-location covariance over a set of items: block means of the fit's
# threshold covariance (recentred parameterisation), which carries the
# negative correlations induced by the identification constraint. A bank
# table contributes its attached joint covariance. Marginal standard errors
# alone cannot recover the correlations induced by its calibration origin;
# they are retained for weighting, but drift inference is then withheld.
.equate_loc_cov <- function(obj, items) {
  if (inherits(obj, "rasch") && !is.null(obj$est$cov_tau)) {
    thr <- obj$est$thr
    # thr$item holds integer item positions in column order, which is the
    # order of obj$items: the match below is by position, not name
    idx <- match(items, obj$items$item)
    rows <- lapply(idx, function(i) thr$id[thr$item == i])
    S <- matrix(0, length(items), length(items))
    for (i in seq_along(rows)) for (j in seq_along(rows))
      S[i, j] <- mean(obj$est$cov_tau[rows[[i]], rows[[j]]])
    return(S)
  }
  ref <- .equate_ref(obj)
  supplied <- .equate_bank_cov(obj, ref$item)
  if (!is.null(supplied)) {
    ii <- match(items, ref$item)
    return(supplied[ii, ii, drop = FALSE])
  }
  se <- ref$se[match(items, ref$item)]
  diag(se^2, length(items))
}

.equate_bank_cov <- function(reference, ids) {
  C <- attr(reference, "cov_location", exact = TRUE)
  if (is.null(C)) return(NULL)
  if (!is.matrix(C) || !is.numeric(C) || any(!is.finite(C)) ||
      !identical(dim(C), c(length(ids), length(ids))))
    stop("attr(reference, 'cov_location') must be a finite numeric square ",
         "matrix with one row and column per bank item")
  if (!isTRUE(all.equal(C, t(C), tolerance = 1e-8)))
    stop("attr(reference, 'cov_location') must be symmetric")
  if (!is.null(rownames(C)) || !is.null(colnames(C))) {
    if (is.null(rownames(C)) || is.null(colnames(C)) ||
        anyNA(match(ids, rownames(C))) || anyNA(match(ids, colnames(C))))
      stop("named bank covariance rows and columns must match every bank item")
    C <- C[ids, ids, drop = FALSE]
  }
  ev <- eigen((C + t(C)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -1e-8 * max(1, max(abs(ev))))
    stop("attr(reference, 'cov_location') must be positive semidefinite")
  C
}

.equate_ref <- function(reference) {
  if (inherits(reference, "rasch"))
    return(data.frame(item = reference$items$item,
                      location = reference$items$location,
                      se = reference$items$se,
                      max = reference$items$max))
  reference <- as.data.frame(reference)
  if (!all(c("item", "location") %in% names(reference)))
    stop("reference needs columns item, location (and ideally se)")
  se_supplied <- "se" %in% names(reference)
  if (!"se" %in% names(reference)) reference$se <- NA_real_
  if (!"max" %in% names(reference)) reference$max <- NA_integer_
  out <- reference[, c("item", "location", "se", "max")]
  out$item <- trimws(as.character(out$item))
  out$location <- suppressWarnings(as.numeric(out$location))
  out$se <- suppressWarnings(as.numeric(out$se))
  out$max <- suppressWarnings(as.numeric(out$max))
  if (anyNA(out$item) || any(!nzchar(out$item)))
    stop("reference item names must be non-missing and non-empty")
  if (anyDuplicated(out$item))
    stop("reference item names must be unique: ",
         paste(unique(out$item[duplicated(out$item)]), collapse = ", "))
  if (any(!is.finite(out$location)))
    stop("reference locations must be finite")
  if (any(!is.na(out$se) & (!is.finite(out$se) | out$se < 0)))
    stop("reference standard errors must be non-negative finite values or NA")
  if (any(!is.na(out$max) & (!is.finite(out$max) | out$max < 1 |
                             out$max != floor(out$max))))
    stop("reference max values must be positive whole numbers or NA")
  attr(out, "se_supplied") <- se_supplied
  out
}

#' Equate two test calibrations through their common items
#'
#' Places two calibrations on a common origin using their shared items, then
#' tests the shared items for drift. The reference may be a fitted model or an
#' item bank.
#'
#' @details
#' Let \eqn{d_j} be the location difference for common item \eqn{j} and
#' \eqn{v_j} its marginal variance. With \code{shift = "mean"} the scale
#' shift is the precision-weighted mean
#' \deqn{\hat s=\frac{\sum_j d_j/v_j}{\sum_j 1/v_j},}
#' and each item is tested using \eqn{d_j-\hat s} with a variance that
#' accounts for the estimated shift through the items' joint covariance.
#' Drift inference requires
#' independent calibrations and at least three common items with usable joint
#' covariance information. Otherwise the function returns a descriptive link.
#' Fitted calibrations must have converged.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param reference A second \code{\link{rasch}} fit, or a data frame with
#'   columns \code{item}, \code{location}, and optionally \code{se}. Item
#'   names must be unique and locations finite. For bank-based drift inference,
#'   attach the bank's joint item-location covariance as a square matrix in
#'   \code{attr(reference, "cov_location")}, ordered like the bank rows (or
#'   named by item); marginal SEs alone do not carry the centring covariance.
#'   A bank treated as fixed may instead have zero SEs. A polytomous bank must
#'   also include \code{max}, the maximum item score.
#' @param shift \code{"mean"} (default) allows a scale shift between the two
#'   analyses; \code{"none"} compares raw locations, appropriate when both
#'   analyses are already on a shared (anchored) scale.
#' @param independent Whether the two calibrations use independent sampling
#'   units. For two fitted objects this must be stated explicitly: the default
#'   \code{NULL} withholds drift tests because cross-fit covariance is
#'   otherwise unknown. A bank table is treated as independent unless
#'   \code{FALSE} is supplied. When \code{FALSE}, descriptive equating is
#'   returned but inferential drift columns are withheld.
#' @return A list with the comparison \code{table} (locations, standard
#'   errors, difference, t, raw and Holm-adjusted p, drift flag), the
#'   estimated \code{shift},
#'   the location \code{correlation}, the root mean square difference after
#'   shifting (\code{rmsd}), the number of common items \code{n_common}, the
#'   number with usable standard errors \code{n}, and whether drift inference
#'   was available (\code{inferential}). The \code{notes} component records
#'   exclusions and the reason inference was withheld, where applicable.
#' @examples
#' set.seed(1); d <- seq(-1.5, 1.5, length.out = 8)
#' mk <- function() {
#'   X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
#'   colnames(X) <- paste0("I", 1:8); rasch(X)
#' }
#' eq <- equate_tests(mk(), mk(), independent = TRUE)
#' eq$table
#' @export
equate_tests <- function(fit, reference, shift = c("mean", "none"),
                         independent = NULL) {
  shift <- match.arg(shift)
  if (!is.null(independent) && (length(independent) != 1L ||
      is.na(independent) || !is.logical(independent)))
    stop("independent must be NULL, TRUE, or FALSE")
  if (!inherits(fit, "rasch") ||
      inherits(fit, "rasch_efrm") || inherits(fit, "rasch_mfrm"))
    stop("fit must be an ordinary person-by-item Rasch calibration")
  if (inherits(reference, "rasch_efrm") || inherits(reference, "rasch_mfrm"))
    stop("reference must be an ordinary Rasch calibration or item bank")
  # under an explanatory design an item location is not an item parameter:
  # items sharing a design cell are forced to one location and carry the
  # design coefficients' uncertainty, so drift can be neither localised nor
  # tested at its stated size
  if (inherits(fit, "rasch_explanatory") ||
      inherits(reference, "rasch_explanatory"))
    stop("item drift is not defined for an explanatory calibration: the ",
         "item locations are functions of their predictors, so a drifted ",
         "item is smeared over every item sharing its design cell and the ",
         "standard errors are the design coefficients'. Equate the ",
         "unrestricted calibrations, or test the design coefficients")
  if (!isTRUE(fit$est$converged))
    stop("the current calibration did not converge; equating is unavailable")
  if (inherits(reference, "rasch") && !isTRUE(reference$est$converged))
    stop("the reference calibration did not converge; equating is unavailable")
  ref <- .equate_ref(reference)
  cur <- data.frame(item = fit$items$item, location = fit$items$location,
                    se = fit$items$se, max = fit$items$max)
  common <- intersect(cur$item, ref$item)
  if (length(common) < 2) stop("need at least two common items to equate")
  a <- cur[match(common, cur$item), ]
  b <- ref[match(common, ref$item), ]
  bank_cov <- if (inherits(reference, "rasch")) NULL else
    .equate_bank_cov(reference, ref$item)
  if (!is.null(bank_cov)) {
    cov_se <- sqrt(pmax(diag(bank_cov), 0))
    stated <- is.finite(ref$se)
    if (any(stated & abs(ref$se - cov_se) >
            1e-6 * pmax(1, ref$se, cov_se)))
      stop("the bank standard errors must agree with the diagonal of ",
           "attr(reference, 'cov_location')")
    if (!isTRUE(attr(ref, "se_supplied", exact = TRUE)))
      ref$se <- cov_se
    b <- ref[match(common, ref$item), ]
  }
  known_max <- is.finite(b$max)
  if (any(known_max & a$max != b$max))
    stop("common items use different maximum scores: ",
         paste(common[known_max & a$max != b$max], collapse = ", "))
  if (any(a$max > 1L) && !all(known_max))
    stop("a polytomous item bank must include a 'max' column so scoring-scale ",
         "compatibility can be checked")
  d <- a$location - b$location
  v <- a$se^2 + b$se^2
  # an item whose location or SE is unavailable (for example a weakly
  # determined item whose SE is honestly NA) cannot contribute to the
  # precision-weighted shift or the drift tests, but it must not poison
  # the remaining items either: exclude it, say so, and carry its row in
  # the table with NA test columns
  usable <- is.finite(d) & is.finite(v)
  notes <- character(0)
  if (any(!usable)) {
    notes <- c(notes, sprintf(
      "common item(s) excluded from the shift and drift tests (location or SE unavailable): %s",
      paste(common[!usable], collapse = ", ")))
  }
  independent_ok <- if (is.null(independent)) !inherits(reference, "rasch")
                    else isTRUE(independent)
  joint_cov_ok <- inherits(reference, "rasch") || !is.null(bank_cov) ||
    all(b$se[usable] == 0)
  if (is.null(independent) && inherits(reference, "rasch"))
    notes <- c(notes, paste(
      "drift tests withheld because independence between the two fitted",
      "calibrations was not stated; set independent = TRUE only for",
      "independent sampling units"))
  if (identical(independent, FALSE))
    notes <- c(notes, paste(
      "drift tests withheld for dependent calibrations because cross-fit",
      "covariance is not available; use a joint or paired bootstrap"))
  if (independent_ok && sum(usable) >= 3L && !joint_cov_ok)
    notes <- c(notes, paste(
      "drift tests withheld because a bank with non-zero marginal SEs needs",
      "its joint item-location covariance in attr(reference, 'cov_location');",
      "marginal SEs do not carry the calibration-origin covariance"))
  if (independent_ok && sum(usable) < 3L)
    notes <- c(notes, paste(
      "drift tests withheld because at least three common items with",
      "standard errors are needed to distinguish item drift from the link"))
  inferential <- independent_ok && sum(usable) >= 3L && joint_cov_ok
  # the shift c0 is estimated from the same common items the drift tests
  # then examine, and each calibration's locations are correlated through
  # its identification constraint: the drift denominators use
  # Var(d_i - c0) = [(I - 1u') Sigma (I - u 1')]_ii over the usable items,
  # with Sigma the sum of the two calibrations' item-location covariances
  var_d <- rep(NA_real_, length(common))
  if (shift == "mean") {
    if (sum(usable) >= 2L) {
      w <- 1 / pmax(v[usable], 1e-10)
      c0 <- sum(w * d[usable]) / sum(w)
    } else c0 <- mean(d[is.finite(d)])
  } else {
    c0 <- 0
  }
  if (inferential) {
    Sg <- .equate_loc_cov(fit, common) + .equate_loc_cov(reference, common)
    if (shift == "mean") {
      u <- w / sum(w)
      Suu <- Sg[usable, usable, drop = FALSE]
      Su <- drop(Suu %*% u)
      var_d[usable] <- pmax(diag(Suu) - 2 * Su +
                              drop(t(u) %*% Su), 1e-10)
    } else var_d[usable] <- pmax(diag(Sg)[usable], 1e-10)
  }
  t <- ifelse(usable, (d - c0) / sqrt(var_d), NA_real_)
  p <- 2 * pnorm(-abs(t))
  n <- sum(usable)
  p_adj <- rep(NA_real_, length(p))
  if (inferential) p_adj[usable] <- p.adjust(p[usable], method = "holm")
  tab <- data.frame(item = common,
                    location_1 = a$location, se_1 = a$se,
                    location_2 = b$location, se_2 = b$se,
                    difference = d, adj_difference = d - c0,
                    t = t, p = p, p_adj = p_adj,
                    drift = ifelse(is.na(p_adj), NA, p_adj < 0.05))
  rownames(tab) <- NULL
  finite_loc <- is.finite(d)
  cor_link <- if (sum(finite_loc) >= 2L)
    stats::cor(a$location[finite_loc], b$location[finite_loc]) else NA_real_
  structure(class = "rasch_equate", list(table = tab, shift = c0,
       correlation = cor_link,
       rmsd = sqrt(mean((d[finite_loc] - c0)^2)), n = n,
       n_common = length(common), inferential = inferential,
       note = if (length(notes)) paste(notes, collapse = "; ") else NULL))
}

#' Plot a test-equating comparison
#'
#' Scatter of the two calibrations' common-item locations with the shifted
#' identity line and per-item 95 per cent bands; drifting items
#' (Holm-adjusted) are highlighted and labelled.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param reference A second \code{\link{rasch}} fit, or a data frame with
#'   columns \code{item}, \code{location}, and optionally \code{se}; a
#'   polytomous bank also needs \code{max}.
#' @param shift Passed to \code{\link{equate_tests}}.
#' @param independent Passed to \code{\link{equate_tests}}.
#' @return Called for its plotting side effect; invisibly the
#'   \code{\link{equate_tests}} result.
#' @examples
#' set.seed(1); d <- seq(-1.5, 1.5, length.out = 8)
#' mk <- function() {
#'   X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
#'   colnames(X) <- paste0("I", 1:8); rasch(X)
#' }
#' plot_equate(mk(), mk(), independent = TRUE)
#' @export
plot_equate <- function(fit, reference, shift = c("mean", "none"),
                        independent = NULL) {
  eq <- equate_tests(fit, reference, shift, independent = independent)
  tab <- eq$table
  paired <- is.finite(tab$location_1) & is.finite(tab$location_2)
  if (!any(paired))
    .refuse("no common item has finite locations in both calibrations; ",
            "there is nothing to display")
  rng <- range(c(tab$location_1[paired], tab$location_2[paired])) +
    c(-0.4, 0.4)
  op <- .rr_canvas(rng, rng, "Reference location (logits)",
                   "Current location (logits)",
                   sprintf("%d common items, shift %.3f, r = %.3f",
                           eq$n_common, eq$shift, eq$correlation),
                   grid_x = TRUE)
  on.exit(par(op))
  abline(eq$shift, 1, col = .rr$ink, lwd = 2)
  # excluded items (NA SEs) still appear as points, but the average band
  # and the drift highlighting are computed over the usable rows only
  band <- 1.96 * sqrt(mean(tab$se_1^2 + tab$se_2^2, na.rm = TRUE))
  if (is.finite(band)) {
    abline(eq$shift + band, 1, lty = 3, col = .rr$soft)
    abline(eq$shift - band, 1, lty = 3, col = .rr$soft)
  }
  hs <- is.finite(tab$se_1)
  segments(tab$location_2[hs], tab$location_1[hs] - 1.96 * tab$se_1[hs],
           tab$location_2[hs], tab$location_1[hs] + 1.96 * tab$se_1[hs],
           col = paste0(.rr$soft, "88"))
  dr <- tab$drift %in% TRUE
  points(tab$location_2, tab$location_1, pch = 21, cex = 1.6,
         bg = ifelse(dr, .rr$red, .rr$blue), col = "white", lwd = 1.2)
  if (any(dr))
    text(tab$location_2[dr], tab$location_1[dr],
         tab$item[dr], pos = 3, offset = 0.5, cex = 0.75, col = .rr$red)
  invisible(eq)
}


#' @export
print.rasch_equate <- function(x, ...) {
  cat(sprintf("Common-item equating over %d item(s): shift %.3f, correlation %.3f, RMSD %.3f\n",
              x$n_common, x$shift, x$correlation, x$rmsd))
  core <- c("item", "location_1", "location_2", "adj_difference", "t",
            "p_adj", "drift")
  print(.fmt_df(x$table[, intersect(core, names(x$table))]), row.names = FALSE)
  if (!isTRUE(x$inferential)) cat("Drift inference withheld; see note below.\n")
  cat("(standard errors and unadjusted columns on $table)\n")
  if (!is.null(x$note)) cat("Note:", x$note, "\n")
  invisible(x)
}
