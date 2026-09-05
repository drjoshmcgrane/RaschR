# rasch :: traditional (classical test theory) statistics
# ===========================================================================
# Classical test theory companion statistics, computed on complete cases
# only (alpha and its relatives have no missing-data form): item facility,
# item-total and item-rest correlations,
# the upper-lower discrimination index over total-score thirds (Andrich &
# Marais 2019 ch. 5), coefficient alpha with alpha-if-item-deleted, and the
# classical summary (mean, SD, SEM) for comparison with the Rasch results.
# ===========================================================================

#' Traditional (classical test theory) statistics
#'
#' The classical companion table conventionally reported alongside a Rasch
#' analysis (Andrich and Marais 2019, chs. 3-5), on complete cases only: per item the facility (mean score over
#' maximum), the item-total and corrected item-rest correlations, the
#' discrimination index DI = PRU - PRL (mean proportion-of-maximum in the
#' upper third of total scores minus the lower third). Equal total scores
#' remain in the same third, so the group sizes can differ and DI is withheld
#' when three distinct score groups cannot be formed. The table also gives
#' alpha if the item is deleted; the summary gives coefficient alpha, the
#' raw-score mean, SD, and the
#' classical standard error of measurement \eqn{s\sqrt{1 - \alpha}}, which
#' unlike the Rasch SE is one value for all persons.
#'
#' @param fit A fitted object from \code{\link{rasch}} whose columns form one
#'   administered item set. Expanded EFRM and MFRM response-cell matrices are
#'   not accepted; one-cell-per-item reductions are.
#' @param missing \code{"complete"} (default) computes the classical table on
#'   respondents who answered every item, matching the textbook total-score
#'   definitions. \code{"available"} retains itemwise and pairwise available
#'   cases as an explicitly exploratory summary; respondents answering
#'   different item sets need not be comparable.
#' @return A list of class \code{"rasch_ctt"}: the per-item \code{table}
#'   (\code{item}, \code{n}, \code{min}, \code{max}, \code{facility},
#'   \code{item_total}, \code{item_rest}, \code{di}, \code{alpha_drop}), and
#'   the scalars
#'   \code{alpha}, \code{n} (complete cases), \code{mean}, \code{sd}, and
#'   \code{sem}.
#' @references
#' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
#' tests. Psychometrika, 16, 297--334.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
#' colnames(X) <- paste0("I", 1:8)
#' ctt_table(rasch(X))
#' @export
ctt_table <- function(fit, missing = c("complete", "available")) {
  if (!inherits(fit, "rasch")) stop("ctt_table needs a rasch fit")
  structural <- inherits(fit, c("rasch_efrm", "rasch_mfrm"))
  if (!.classical_design_applicable(fit))
    stop("traditional whole-test statistics are not defined when an item is ",
         "represented by several frame or facet response cells; ",
         "summarise an observable item design instead")
  missing <- match.arg(missing)
  Xall <- fit$X
  if (structural) {
    item_names <- fit$virtual_map$item[
      match(colnames(Xall), fit$virtual_map$vkey)]
    if (anyNA(item_names) || anyDuplicated(item_names))
      stop("the one-cell-per-item reduction could not be matched to its ",
           "item names")
    colnames(Xall) <- item_names
  }
  X <- if (missing == "complete")
    Xall[stats::complete.cases(Xall), , drop = FALSE] else Xall
  L <- ncol(X)
  n_i <- colSums(!is.na(X))
  if (all(n_i < 3)) {
    tab <- data.frame(item = colnames(X), n = n_i, min = NA_real_,
                      max = fit$m, facility = NA_real_, item_total = NA_real_,
                      item_rest = NA_real_, di = NA_real_,
                      alpha_drop = NA_real_)
    out <- list(table = tab, alpha = NA_real_,
                n = sum(stats::complete.cases(X)),
                n_range = range(n_i), mean = NA_real_, sd = NA_real_,
                sem = NA_real_, missing = missing,
                note = if (missing == "complete")
                  "fewer than 3 complete responders: classical statistics withheld"
                else "fewer than 3 usable responses per item: classical statistics withheld")
    out <- .tag_tables(out)
  class(out) <- "rasch_ctt"
    return(out)
  }
  # Under available-case mode, total proportions and pairwise covariances are
  # exploratory because different answered sets can differ in difficulty.
  Mmat <- matrix(rep(fit$m, each = nrow(X)), nrow(X), L)
  Mmat[is.na(X)] <- NA
  tot_p <- rowSums(X, na.rm = TRUE) / rowSums(Mmat, na.rm = TRUE)
  thirds <- rep(NA_integer_, nrow(X))
  has_total <- is.finite(tot_p)
  if (sum(has_total) >= 3L)
    thirds[has_total] <- .class_intervals(
      tot_p[has_total], rep(FALSE, sum(has_total)), 3)
  # pairwise covariance carries the alpha-if-deleted computation under
  # missing data (equal to the variance form when data are complete). A
  # pairwise-complete covariance need not be positive semidefinite, and a
  # near-zero total covariance turns the alpha formula into an absurd
  # number (an observed case printed alpha = -257): alpha is a ratio of
  # variances, so it is only reported when C is a valid covariance matrix
  # with a non-degenerate total
  C <- suppressWarnings(stats::cov(X, use = "pairwise.complete.obs"))
  csum <- sum(C, na.rm = TRUE)
  dsum <- sum(diag(C), na.rm = TRUE)
  # Pairwise deletion can produce a genuinely indefinite matrix; this is not
  # numerical rounding and cannot define coefficient alpha. Use the same
  # scale-relative PSD tolerance as the package's other covariance checks.
  C_ok <- !anyNA(C) && .covariance_is_psd(C) &&
    is.finite(csum) && csum > 0.05 * dsum
  alpha <- if (L > 1L && C_ok)
    L / (L - 1L) * (1 - dsum / csum) else NA_real_
  min_i <- suppressWarnings(vapply(seq_len(L), function(i)
    if (n_i[i] == 0) NA_real_ else min(X[, i], na.rm = TRUE), 0))
  facility <- colMeans(X, na.rm = TRUE) / fit$m
  facility[n_i == 0L] <- NA_real_
  tab <- data.frame(item = colnames(X), n = n_i, min = min_i, max = fit$m,
                    facility = facility,
                    item_total = NA_real_, item_rest = NA_real_,
                    di = NA_real_, alpha_drop = NA_real_)
  for (i in seq_len(L)) {
    x <- X[, i]; ok <- !is.na(x)
    if (sum(ok) >= 3 && stats::sd(x[ok]) > 0) {
      rest_max <- rowSums(Mmat, na.rm = TRUE) - ifelse(ok, fit$m[i], 0)
      rest_p <- (rowSums(X, na.rm = TRUE) - ifelse(ok, x, 0)) / rest_max
      # A respondent who answered only the focal item has no rest score.
      # Treating the empty rest set as zero biases the corrected correlation
      # in sparse or routed designs.
      rest_p[!is.finite(rest_p) | rest_max <= 0] <- NA_real_
      tab$item_total[i] <- .safe_cor(x[ok], tot_p[ok])
      tab$item_rest[i] <- .safe_cor(x[ok], rest_p[ok])
      hi <- ok & thirds == 3; lo <- ok & thirds == 1
      if (sum(hi) >= 2 && sum(lo) >= 2)
        tab$di[i] <- mean(x[hi]) / fit$m[i] - mean(x[lo]) / fit$m[i]
    }
    if (L > 2 && C_ok) {
      Cr <- C[-i, -i, drop = FALSE]
      sr <- sum(Cr, na.rm = TRUE)
      if (is.finite(sr) && sr > 0.05 * sum(diag(Cr)))
        tab$alpha_drop[i] <- (L - 1) / (L - 2) *
          (1 - sum(diag(Cr), na.rm = TRUE) / sr)
    }
  }
  rownames(tab) <- NULL
  cc <- stats::complete.cases(X)
  tot_cc <- rowSums(X[cc, , drop = FALSE])
  out <- list(table = tab, alpha = alpha, n = sum(cc),
              n_range = range(n_i),
              mean = if (sum(cc) >= 3) mean(tot_cc) else NA_real_,
              sd = if (sum(cc) >= 3) stats::sd(tot_cc) else NA_real_,
              sem = if (is.finite(alpha) && alpha <= 1 && sum(cc) >= 3)
                stats::sd(tot_cc) * sqrt(1 - alpha) else NA_real_,
              missing = missing,
              note = {
                nt <- c(
                  if (missing == "available") paste(
                    "available-case item statistics are exploratory; persons",
                    "answering different item sets are not necessarily comparable"),
                  if (missing == "available" && !C_ok) paste(
                    "alpha withheld: the pairwise covariance under this",
                    "missingness is not a valid (positive semidefinite,",
                    "non-degenerate) covariance matrix"))
                if (length(nt)) paste(nt, collapse = "; ") else NULL
              })
  out <- .tag_tables(out)
  class(out) <- "rasch_ctt"
  out
}

#' @export
print.rasch_ctt <- function(x, ...) {
  cat(sprintf("Traditional statistics (%s cases; item n %d-%d; %d complete)\n",
              x$missing, x$n_range[1], x$n_range[2], x$n))
  if (is.finite(x$mean))
    cat(sprintf("Raw score mean %.2f, SD %.2f (complete responders); alpha %.3f; SEM %.2f\n",
                x$mean, x$sd, x$alpha, x$sem))
  else
    cat(sprintf("Too few complete responders for total-score summaries; alpha %.3f\n",
                x$alpha))
  y <- x$table
  num <- vapply(y, is.numeric, TRUE)
  y[num] <- lapply(y[num], round, 3)
  print(y, row.names = FALSE)
  if (!is.null(x$note)) cat("Note:", x$note, "\n")
  invisible(x)
}

#' Reshape repeated measurements for racked or stacked analysis
#'
#' Repeated measurements (the same persons and items at two or more time
#' points) enter a Rasch analysis in one of two designs (Andrich & Marais
#' 2019, ch. 26). \emph{Racking} keeps one row per person and duplicates
#' the items per time point (columns \code{item@time}), so change over time
#' shows in the item estimates. \emph{Stacking} keeps one column per item
#' and duplicates the persons per time point (rows), so
#' change shows in the person estimates and DIF of items over time can be
#' examined with \code{time} as a within-person factor. The returned
#' \code{id} is the original person identifier and therefore repeats across
#' occasions; \code{row_id} uniquely identifies each person-occasion row.
#'
#' @param data A long data frame with one measurement per row.
#' @param person,time Names of distinct person and time-point columns.
#' @param items Character vector naming the item columns.
#' @return \code{rack_data}: a wide data frame with one row per person and
#'   \code{length(items) * n_times} item columns. \code{stack_data}: a data
#'   frame with one row per person-time, the repeated original \code{id}, a
#'   unique \code{row_id}, the original item columns, and \code{time} as a
#'   factor column for repeated-measures DIF analysis.
#' @examples
#' d <- data.frame(pid = rep(1:100, 2), t = rep(1:2, each = 100),
#'                 Q1 = rbinom(200, 1, 0.6), Q2 = rbinom(200, 1, 0.5))
#' racked <- rack_data(d, person = "pid", time = "t", items = c("Q1", "Q2"))
#' names(racked)
#' stacked <- stack_data(d, person = "pid", time = "t", items = c("Q1", "Q2"))
#' head(stacked)
#' # the follow-up analysis assigns every reshaped column a role: the
#' # repeated person id, time as a within-person factor, and the items
#' fit <- rasch(stacked, id = "id", factors = "time",
#'              items = c("Q1", "Q2"))
#' @export
rack_data <- function(data, person, time, items) {
  .check_column_names(data)
  data <- as.data.frame(data)
  if (!is.character(items) || !length(items) || anyNA(items) ||
      any(!nzchar(trimws(items))))
    stop("`items` must name at least one item column")
  # dereferencing an empty or multiple column name gives a base subscript
  # error rather than a statement of what was wrong
  .check_reshape_column(data, person, "person")
  .check_reshape_column(data, time, "time")
  if (identical(as.character(person), as.character(time)))
    stop("the person and time columns must be distinct")
  pv <- .role_text_values(data[[person]])
  tv <- .role_text_values(data[[time]])
  if (anyNA(pv) || any(!nzchar(pv[!is.na(pv)])))
    stop(sum(is.na(pv) | !nzchar(pv)), " row(s) have a missing or blank ",
         "person identifier; they cannot be aligned across occasions")
  if (anyNA(tv) || any(!nzchar(tv[!is.na(tv)])))
    stop(sum(is.na(tv) | !nzchar(tv)),
         " row(s) have a missing or blank occasion value")
  if (anyDuplicated(items))
    stop("item column(s) named more than once: ",
         paste(unique(items[duplicated(items)]), collapse = ", "))
  overlap <- intersect(items, c(person, time))
  if (length(overlap))
    stop("the person or time column cannot also be an item: ",
         paste(overlap, collapse = ", "))
  for (col in c(person, time)) if (!col %in% names(data))
    stop("column not found: ", col)
  bad <- setdiff(items, names(data))
  if (length(bad)) stop("item column(s) not found: ", paste(bad, collapse = ", "))
  # Use the same canonical values that were validated above. Otherwise
  # visually identical labels such as "T1" and " T1 " pass the blank check
  # but become different occasions, and padded person IDs become different
  # respondents in the reshaped design.
  time_column <- .canonical_role_column(data[[time]])
  times <- sort(unique(time_column))
  time_labels <- as.character(times)
  made <- unlist(lapply(time_labels, function(tt) paste0(items, "@", tt)),
                 use.names = FALSE)
  if (anyDuplicated(c("id", made)))
    stop("generated racked column names are not unique; rename the items or time levels")
  first_id <- !duplicated(pv)
  ids <- pv[first_id]
  id_out <- .canonical_role_column(data[[person]])[first_id]
  out <- data.frame(id = id_out)
  for (j in seq_along(times)) {
    tt <- times[j]
    time_label <- time_labels[j]
    rows_t <- which(time_column == tt)
    d_t <- data[rows_t, , drop = FALSE]
    p_t <- pv[rows_t]
    if (anyDuplicated(p_t))
      stop("more than one row for a person at time ", time_label)
    idx <- match(ids, p_t)
    blk <- d_t[idx, items, drop = FALSE]
    names(blk) <- paste0(items, "@", time_label)
    out <- cbind(out, blk)
  }
  rownames(out) <- NULL
  out
}

#' @rdname rack_data
#' @export
stack_data <- function(data, person, time, items) {
  .check_column_names(data)
  data <- as.data.frame(data)
  if (!is.character(items) || !length(items) || anyNA(items) ||
      any(!nzchar(trimws(items))))
    stop("`items` must name at least one item column")
  # dereferencing an empty or multiple column name gives a base subscript
  # error rather than a statement of what was wrong
  .check_reshape_column(data, person, "person")
  .check_reshape_column(data, time, "time")
  if (identical(as.character(person), as.character(time)))
    stop("the person and time columns must be distinct")
  pv <- .role_text_values(data[[person]])
  tv <- .role_text_values(data[[time]])
  if (anyNA(pv) || any(!nzchar(pv[!is.na(pv)])))
    stop(sum(is.na(pv) | !nzchar(pv)), " row(s) have a missing or blank ",
         "person identifier; they cannot be aligned across occasions")
  if (anyNA(tv) || any(!nzchar(tv[!is.na(tv)])))
    stop(sum(is.na(tv) | !nzchar(tv)),
         " row(s) have a missing or blank occasion value")
  if (anyDuplicated(items))
    stop("item column(s) named more than once: ",
         paste(unique(items[duplicated(items)]), collapse = ", "))
  overlap <- intersect(items, c(person, time))
  if (length(overlap))
    stop("the person or time column cannot also be an item: ",
         paste(overlap, collapse = ", "))
  for (col in c(person, time)) if (!col %in% names(data))
    stop("column not found: ", col)
  bad <- setdiff(items, names(data))
  if (length(bad)) stop("item column(s) not found: ", paste(bad, collapse = ", "))
  reserved <- intersect(items, c("id", "row_id", "time"))
  if (length(reserved))
    stop("item name(s) reserved by the stacked output: ",
         paste(reserved, collapse = ", "), "; rename them before stacking")
  key <- .factor_cells(data.frame(person = pv, time = tv), sep = "\r")
  if (anyDuplicated(key))
    stop("more than one row for a person at the same time point")
  row_id <- .factor_cells(data.frame(person = pv, time = tv), sep = "@")
  time_column <- .canonical_role_column(data[[time]])
  if (!is.factor(time_column)) time_column <- factor(time_column)
  out <- data.frame(id = .canonical_role_column(data[[person]]),
                    row_id = row_id,
                    time = time_column,
                    data[, items, drop = FALSE], check.names = FALSE)
  rownames(out) <- NULL
  out
}
