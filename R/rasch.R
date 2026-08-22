# rasch :: top-level analysis
# ===========================================================================
# rasch() ties the engine together: data preparation (with category
# collapsing and constant-item removal, recorded as notes), pairwise
# conditional ML item estimation, Warm WLE person estimation per
# missing-data pattern, the full test-of-fit suite, and the score-to-measure
# table. An ID variable and any number of person factors carry through to
# the person estimates.
# ===========================================================================

# Map missing-data codes (and any negative score) to NA. Valid item scores
# are non-negative integers from zero; by long-standing convention -1 marks a
# missing response, so any value below zero is read as missing.
.apply_na_codes <- function(v, na_codes) {
  v[v %in% na_codes | (!is.na(v) & v < 0)] <- NA
  v
}

# Prepare the item matrix: integer scores from 0, consecutive observed
# categories, no constant items. Returns the matrix plus human-readable notes.
.prepare_X <- function(X, na_codes = -1) {
  notes <- character(0)
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- sprintf("I%02d", seq_len(ncol(X)))
  if (anyNA(colnames(X)) || any(!nzchar(colnames(X))))
    stop("item column names must be non-missing and non-empty")
  if (anyDuplicated(colnames(X)))
    stop("item column names must be unique: ",
         paste(unique(colnames(X)[duplicated(colnames(X))]), collapse = ", "))
  Xn <- suppressWarnings(apply(X, 2, function(col) as.numeric(as.character(col))))
  Xi <- suppressWarnings(apply(X, 2, function(col) as.integer(as.character(col))))
  dim(Xn) <- dim(X); dim(Xi) <- dim(X); dimnames(Xi) <- dimnames(X)
  # as.integer() TRUNCATES fractional values (1.9 -> 1) without a warning:
  # that silently alters response data, so it must be an error, not a note
  frac <- colSums(!is.na(Xn) & !is.na(Xi) & Xn != Xi) > 0
  if (any(frac))
    stop("non-integer score(s) in: ",
         paste(colnames(X)[frac], collapse = ", "),
         " (e.g. ", format(Xn[!is.na(Xn) & !is.na(Xi) & Xn != Xi][1]),
         "); Rasch categories are integer counts -- round or rescore ",
         "explicitly before analysis")
  bad_num <- colSums(!is.na(X) & is.na(Xi)) > 0
  if (any(bad_num))
    notes <- c(notes, paste0("non-numeric entries set to missing in: ",
                             paste(colnames(X)[bad_num], collapse = ", ")))
  n_na <- sum(!is.na(Xi) & (Xi %in% na_codes | Xi < 0))
  Xi[] <- .apply_na_codes(Xi, na_codes)
  if (n_na > 0) {
    codes <- paste(unique(c(na_codes, "negative")), collapse = ", ")
    notes <- c(notes, sprintf("%d cell(s) with a missing-data code (%s) set to missing",
                              n_na, codes))
  }
  X <- Xi
  const <- apply(X, 2, function(col) length(unique(col[!is.na(col)])) < 2)
  if (any(const)) {
    notes <- c(notes, paste0("dropped constant item(s): ",
                             paste(colnames(X)[const], collapse = ", ")))
    X <- X[, !const, drop = FALSE]
  }
  if (ncol(X) < 2) stop("need at least two non-constant items")
  for (i in seq_len(ncol(X))) {
    v <- X[, i]; obs <- sort(unique(v[!is.na(v)]))
    full <- seq(0L, max(obs))
    if (!identical(obs, full)) {
      X[, i] <- match(v, obs) - 1L
      notes <- c(notes, sprintf("item %s rescored: observed categories [%s] mapped to 0:%d",
                                colnames(X)[i], paste(obs, collapse = ","),
                                length(obs) - 1L))
    }
  }
  list(X = X, notes = notes)
}

#' Fit a Rasch model
#'
#' Fits the partial credit model (PCM) or rating scale model (RSM) by pairwise
#' conditional maximum likelihood. Person locations are Warm weighted
#' likelihood estimates. The fitted object contains item and person fit,
#' targeting, reliability, threshold diagnostics, residuals, and a
#' score-to-measure table.
#'
#' @details
#' For scores \eqn{x=0,\ldots,m_i}, the PCM is
#' \deqn{P(X_{ni}=x)=\frac{\exp\{x\theta_n-\sum_{k=1}^{x}\delta_{ik}\}}
#' {\sum_{y=0}^{m_i}\exp\{y\theta_n-\sum_{k=1}^{y}\delta_{ik}\}}.}
#' The RSM constrains \eqn{\delta_{ik}=\beta_i+\tau_k}, where \eqn{\beta_i}
#' is the item location and \eqn{\tau_k} is common across items. Dichotomous
#' items are the one-threshold case of the PCM.
#'
#' Pairwise conditioning removes \eqn{\theta_n} from the item likelihood.
#' Missing responses are omitted from pairwise contributions, and person
#' measures are estimated within each observed item pattern. The observed
#' item-pair graph must identify a common scale. This covers planned linked
#' designs and ignorable missingness; informative missingness can still bias
#' the estimates.
#'
#' The fit residual is the log-of-mean-square statistic described by Andrich
#' and Marais (2019, ch. 23). It is approximately standard normal under fit;
#' positive values indicate under-discrimination and negative values indicate
#' over-discrimination. The item-trait chi-square and class-interval F tests
#' are large-sample diagnostic approximations and should be considered with
#' the residual statistics, effect sizes, and item content.
#'
#' Multiple-choice responses may be scored from a named item-to-key vector,
#' an item/key table, or an item/option/score table. A slash separates
#' alternative correct options. The third form assigns integer category
#' scores to nominated options and fits the resulting item as polytomous;
#' unlisted options score zero. Raw responses are retained in \code{fit$mc}
#' for distractor analysis.
#'
#' If \code{adjust_N} is supplied, each item-trait chi-square is multiplied by
#' the reference sample size divided by the number of classified persons.
#' The scaling is global: an item answered by a subset retains its
#' proportionally smaller share of the reference sample.
#'
#' @param data Persons-by-items integer score matrix (categories from 0), or a
#'   data frame also containing ID and person-factor columns. Missing values
#'   are allowed subject to the identification and ignorability conditions
#'   described above.
#' @param model Either \code{"PCM"} (partial credit) or \code{"RSM"} (rating
#'   scale).
#' @param id Optional name of an ID column in \code{data}, or a vector of IDs;
#'   carried through to the person estimates.
#' @param factors Optional character vector of person-factor column names in
#'   \code{data} (for DIF analysis), a data frame of factors, or one grouping
#'   vector with one entry per data row.
#' @param items Optional character vector naming the item columns; by default
#'   every column not named in \code{id} or \code{factors}.
#' @param n_groups Number of class intervals for the item-trait chi-square
#'   and ANOVA item fit. The default \code{NULL} applies the rule of Andrich
#'   and Marais (2019, ch. 15): as
#'   many intervals of at least 50 non-extreme persons as the sample allows,
#'   at most 10, at least 2. The resolved value is stored in
#'   \code{fit$n_groups}.
#' @param adjust_N Optional reference sample size used to rescale the
#'   item-trait chi-squares. See Details.
#' @param anchors Optional anchor table for equating: a data frame with
#'   columns \code{item}, \code{k}, and \code{tau}; see \code{\link{pcml}}.
#'   Anchors determine the scale origin.
#' @param na_codes Values to read as missing. Defaults to \code{-1}, the
#'   conventional missing-response code; any negative score is also treated as
#'   missing, since valid category scores start at zero.
#' @param maxit,tol Newton-Raphson iteration cap and convergence
#'   tolerance of the pairwise conditional estimation.
#' @param key Optional multiple-choice key: a named item-to-option vector, an
#'   item/key table, or an item/option/score table. See Details.
#' @param pc_components \code{NULL} (the default) estimates all PCM thresholds
#'   freely. Values from 1 to 4 use the principal-components form in
#'   \code{\link{pcml_pc}}: location, then spread, skewness, and kurtosis.
#'   This can stabilise sparse categories. Component estimates are stored in
#'   the estimation details. Available for PCM fits without anchors.
#' @section Estimated item discrimination:
#' The item summary includes a post-estimation slope \code{disc}. For item
#' \eqn{i}, it maximises that item's response likelihood over \eqn{a_i} while
#' holding the fitted person locations and thresholds fixed:
#' \deqn{\hat a_i=\arg\max_{a_i}
#'   \sum_n\log P(X_{ni}=x_{ni}\mid\hat\theta_n,\hat\delta_i,a_i).}
#' The same slope multiplies every threshold of a polytomous item. It is a
#' descriptive index, not a freely estimated parameter of the Rasch model,
#' and no sampling standard error or hypothesis test is attached to it.
#' @return An object of class \code{"rasch"}. Its principal components are
#'   the item summary, threshold table, person table, score table, residuals,
#'   reliability, targeting, item-trait statistics, threshold diagnostics,
#'   and estimation details. The component \code{summary_stats} contains the
#'   distribution summaries, fit-location correlations, and the cell
#'   degrees-of-freedom factor. The item summary carries a \code{disc}
#'   column described below.
#' @references
#' Rasch, G. (1960). Probabilistic Models for Some Intelligence and
#' Attainment Tests. Copenhagen: Danish Institute for Educational Research.
#' (Expanded edition, 1980, Chicago: University of Chicago Press.)
#'
#' Rasch, G. (1961). On general laws and the meaning of measurement in
#' psychology. In Proceedings of the Fourth Berkeley Symposium on
#' Mathematical Statistics and Probability (Vol. 4, pp. 321--333).
#' Berkeley: University of California Press.
#'
#' Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
#' Rasch model for ordered response categories using principal components.
#' Journal of Applied Measurement, 4(3), 205--221.
#'
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. Psychometrika, 54(3), 427--450.
#' @seealso \code{\link{rasch_mfrm}}, \code{\link{rasch_efrm}},
#'   \code{\link{btl}}, \code{\link{dif_anova}},
#'   \code{\link{test_information}}, and \code{\link{run_app}}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(X, model = "PCM")
#' fit$items
#' fit$psi$PSI
#' @export
rasch <- function(data, model = c("PCM", "RSM"), id = NULL, factors = NULL,
                  items = NULL, n_groups = NULL, adjust_N = NA, anchors = NULL,
                  na_codes = -1, key = NULL, pc_components = NULL,
                  maxit = 60, tol = 1e-8) {
  .check_column_names(data)
  n_groups_requested <- n_groups
  # name for a factors= vector passed by value (not by column name)
  .factors_sym <- substitute(factors)
  .factors_label <- if (is.name(.factors_sym)) as.character(.factors_sym) else "factor"
  model <- match.arg(model)
  # adjust_N rescales the item-trait chi-square by (reference N / classified
  # N); a non-positive reference would zero or negate every statistic
  if (!is.na(adjust_N) && (!is.numeric(adjust_N) || adjust_N <= 0))
    stop("`adjust_N` must be a positive reference sample size")
  if (!is.null(pc_components)) {
    if (model != "PCM")
      stop("pc_components applies to the PCM only")
    if (!is.null(anchors))
      stop("pc_components cannot be combined with anchors")
  }

  # --- split data frame into ID, factors, and item columns ---------------
  id_vec <- NULL; fac_df <- NULL
  # a simulated dataset carries its own person identifier: use it, so the
  # documented bare call rasch(simulate_rasch(...)) keeps person ids
  if (inherits(data, "rasch_sim") && is.null(id) && "id" %in% names(data))
    id <- "id"
  if (is.data.frame(data)) {
    nm <- names(data)
    id_is_col <- is.character(id) && length(id) == 1L
    # a misspelled column name must be an error, never a silent fallback:
    # dropping it quietly produces a valid-looking analysis of the wrong data
    if (id_is_col) {
      if (!id %in% nm)
        stop("id column '", id, "' not found in the data")
      id_vec <- data[[id]]
    } else if (!is.null(id)) {
      # a supplied id vector must line up with the data; a length mismatch
      # (e.g. a stale upstream vector) must error, not be silently dropped
      if (length(id) != nrow(data))
        stop("`id` has ", length(id), " entries but the data has ",
             nrow(data), " rows")
      id_vec <- id
    }
    # Character input is ambiguous: all existing column names means the
    # documented column-name form; a row-length character vector is a
    # grouping vector passed by value. A short non-matching character input
    # remains a misspelled-column error rather than silently changing modes.
    factors_are_cols <- is.character(factors) &&
      (length(factors) == 0L || all(factors %in% nm) ||
         length(factors) != nrow(data))
    factors_by_value <- !is.null(factors) && is.atomic(factors) &&
      !factors_are_cols
    if (factors_are_cols) {
      miss <- setdiff(factors, nm)
      if (length(miss))
        stop("factor column(s) not found in the data: ",
             paste(miss, collapse = ", "))
      fac_df <- data[, factors, drop = FALSE]
    } else if (is.data.frame(factors)) {
      if (nrow(factors) != nrow(data))
        stop("`factors` data frame has ", nrow(factors), " rows but the data ",
             "has ", nrow(data), " rows")
      fac_df <- factors
    } else if (factors_by_value) {
      # a factors= grouping vector passed by VALUE (not by column name):
      # accept it when it lines up, rather than silently ignoring it and
      # leaving any same-named data column to be treated as an item
      if (!is.atomic(factors) || length(factors) != nrow(data))
        stop("`factors` must be column name(s) in the data, a data frame ",
             "with one row per data row, or a vector with one entry per row")
      fac_df <- stats::setNames(data.frame(factors, stringsAsFactors = FALSE),
                                .factors_label)
    } else if (!is.null(factors)) {
      stop("`factors` must be column name(s) in the data, a data frame ",
           "with one row per data row, or a vector with one entry per row")
    }
    # a data column whose values are identical to a by-value factors vector
    # is almost certainly that same variable: exclude it so it is not also
    # scored as a numeric item
    val_factor_cols <- if (factors_by_value)
      nm[vapply(data, function(col)
        length(col) == length(factors) && isTRUE(all.equal(
          as.character(col), as.character(factors))), logical(1))] else NULL
    val_id_cols <- if (!is.null(id) && !id_is_col)
      nm[vapply(data, function(col)
        length(col) == length(id) && isTRUE(all.equal(
          as.character(col), as.character(id))), logical(1))] else NULL
    drop_cols <- c(if (id_is_col) id else val_id_cols,
                   if (factors_are_cols) factors else NULL,
                   # an externally supplied factor data frame whose column
                   # names also appear in `data` almost certainly refers to
                   # those columns: without this they would silently become
                   # numeric ITEMS
                   if (is.data.frame(factors))
                     intersect(names(factors), nm) else NULL,
                   val_factor_cols)
    # identifier-named columns must never be silently SCORED as items: the
    # stacked/racked reshapes emit id/row_id/time columns, and calling
    # rasch(stacked) without id = "id" would otherwise rescore a numeric
    # person identifier as a many-category item with a valid-looking
    # report. Only numeric-convertible columns can be scored, so only they
    # are refused; character identifiers keep the old dropped-with-a-note
    # path (they can never silently enter the item matrix).
    ident_like <- intersect(c("id", "row_id", "time", "person"), nm)
    ident_like <- setdiff(ident_like, c(drop_cols,
                                        if (is.character(items)) items))
    ident_like <- ident_like[vapply(ident_like, function(cn) {
      v <- data[[cn]]
      vn <- suppressWarnings(as.numeric(as.character(v)))
      any(!is.na(vn))
    }, logical(1))]
    if (is.null(items) && length(ident_like))
      stop("the data contain identifier-like column(s) not assigned a role: ",
           paste(ident_like, collapse = ", "),
           " -- pass them via id=/factors= and name the item columns with ",
           "items= (for stack_data output: rasch(stacked, id = \"id\", ",
           "factors = \"time\", items = <the item columns>)), or drop them")
    item_cols <- if (is.null(items)) setdiff(nm, drop_cols)
    else if (is.character(items)) {
      miss <- setdiff(items, nm)
      if (length(miss))
        stop("item column(s) not found in the data: ",
             paste(miss, collapse = ", "))
      items
    } else nm[items]
    # an explicit items= must not silently pull in an id/factor column (a
    # positional items = 1:k over an id-first layout would score the id as
    # an item and drop a real one), nor name the same column twice
    if (!is.null(items)) {
      clash <- intersect(item_cols, drop_cols)
      if (length(clash))
        stop("items= includes id/factor column(s): ",
             paste(clash, collapse = ", "),
             " -- name only item columns, or drop them from id=/factors=")
    }
    dup <- item_cols[duplicated(item_cols)]
    if (length(dup))
      stop("item column(s) named more than once: ",
           paste(unique(dup), collapse = ", "))
    X <- as.matrix(data[, item_cols, drop = FALSE])
  } else {
    X <- as.matrix(data)
    if (!is.null(id)) {
      if (length(id) != nrow(X))
        stop("`id` has ", length(id), " entries but the data has ",
             nrow(X), " rows")
      id_vec <- id
    }
    if (is.data.frame(factors)) {
      if (nrow(factors) != nrow(X))
        stop("`factors` data frame has ", nrow(factors), " rows but the data ",
             "has ", nrow(X), " rows")
      fac_df <- factors
    } else if (!is.null(factors)) {
      if (!is.atomic(factors) || length(factors) != nrow(X))
        stop("`factors` must be a data frame with one row per data row, or ",
             "a vector with one entry per row")
      fac_df <- stats::setNames(data.frame(factors, stringsAsFactors = FALSE),
                                .factors_label)
    }
  }
  if (is.null(id_vec)) id_vec <- seq_len(nrow(X))

  # score multiple-choice items against the key, keeping the raw responses
  mc <- NULL
  if (!is.null(key)) {
    key <- .resolve_key(key)
    sc <- .score_mc(X, key)
    X[, colnames(sc$scored)] <- sc$scored
    mc <- list(key = sc$key, map = sc$map, raw = sc$raw)
  }

  prep <- .prepare_X(X, na_codes = na_codes); X <- prep$X
  if (!is.null(mc)) {
    prep$notes <- c(prep$notes,
                    sprintf("%d item(s) scored 0/1 against the key", ncol(mc$raw)))
    gone <- setdiff(colnames(mc$raw), colnames(X))
    if (length(gone)) mc$raw <- mc$raw[, setdiff(colnames(mc$raw), gone), drop = FALSE]
  }

  if (!is.null(anchors)) {
    a_names <- if (is.character(anchors$item) || is.factor(anchors$item))
      as.character(anchors$item) else colnames(X)[anchors$item]
    gone <- setdiff(a_names, colnames(X))
    if (length(gone))
      stop("anchored item(s) not present after data preparation: ",
           paste(gone, collapse = ", "))
    resc <- grepl("rescored", prep$notes) &
      vapply(prep$notes, function(n) any(vapply(a_names, grepl, TRUE, x = n)), TRUE)
    if (any(resc))
      stop("anchored item(s) were rescored during data preparation; ",
           "anchor values would no longer match the threshold numbering")
    prep$notes <- c(prep$notes,
                    sprintf("%d threshold(s) anchored; scale origin from anchors",
                            nrow(anchors)))
  }

  # --- item estimation ----------------------------------------------------
  est <- if (is.null(pc_components))
    pcml(X, model = model, anchors = anchors, maxit = maxit, tol = tol)
  else pcml_pc(X, n_components = pc_components, maxit = maxit, tol = tol)
  if (!is.null(pc_components))
    prep$notes <- c(prep$notes,
                    sprintf("thresholds estimated through %d principal component(s); see est$components",
                            pc_components))
  if (!isTRUE(est$converged))
    warning("estimation did NOT converge in ", est$iterations,
            " iterations: estimates, standard errors, fit statistics, and ",
            "p-values are unreliable -- increase maxit or check the data ",
            "for unanswerable structure", call. = FALSE)
  fit <- .assemble_fit(model, X, est, id_vec, fac_df, n_groups, adjust_N,
                       c(prep$notes, est$notes))
  fit$mc <- mc
  # Keep the arguments that define the fitted model. Post-fit operations such
  # as drop_items() must not silently change the identification, threshold
  # parameterisation, fit grouping, or optimiser controls when they refit.
  # Anchors are stored by item name so their meaning survives column removal.
  anchors_named <- anchors
  if (!is.null(anchors_named) &&
      !(is.character(anchors_named$item) || is.factor(anchors_named$item)))
    anchors_named$item <- colnames(X)[as.integer(anchors_named$item)]
  key_spec <- NULL
  if (!is.null(mc)) {
    key_spec <- do.call(rbind, lapply(names(mc$map), function(it) {
      data.frame(item = it, option = names(mc$map[[it]]),
                 score = unname(mc$map[[it]]), stringsAsFactors = FALSE)
    }))
    rownames(key_spec) <- NULL
  }
  fit$refit_spec <- list(
    model = model, n_groups = n_groups_requested, adjust_N = adjust_N,
    anchors = anchors_named, na_codes = na_codes, key = key_spec,
    pc_components = pc_components, maxit = maxit, tol = tol)
  fit
}

# Post-estimation pipeline shared by rasch(), rasch_mfrm(), and rasch_efrm():
# person estimation, residuals, the full fit suite, and the assembled tables.
# disc is an optional per-column discrimination (frame unit) vector; with
# unequal discriminations the raw score is no longer sufficient, so person
# estimation switches to the weighted-score routine and the score table is
# replaced by per-unit score curves.
.assemble_fit <- function(model, X, est, id_vec, fac_df, n_groups, adjust_N,
                          notes, disc = NULL) {
  m <- est$m; L <- ncol(X)
  thr <- est$thr
  tau_list <- lapply(seq_len(L), function(i) thr$tau[thr$item == i])
  names(tau_list) <- colnames(X)
  equal_disc <- is.null(disc) || length(unique(disc)) == 1L
  disc_v <- if (is.null(disc)) rep(1, L) else disc

  # --- person estimation and residuals -------------------------------------
  person <- if (equal_disc) .person_estimates(X, tau_list, disc = disc_v[1])
            else .efrm_person_estimates(X, tau_list, disc_v)
  mo <- .moment_arrays(person$theta, tau_list, disc = disc_v)
  Z <- (X - mo$E) / sqrt(mo$V)
  colnames(Z) <- colnames(X)

  # --- fit statistics ------------------------------------------------------
  # an item scored at its floor or ceiling by every non-extreme person has
  # no finite location; exclude it from person fit as extreme persons are
  # excluded from item fit
  m_i <- if (!is.null(est$m)) est$m else apply(X, 2, max, na.rm = TRUE)
  item_extreme <- vapply(seq_len(ncol(X)), function(j) {
    col <- X[!person$extreme, j]
    tot <- sum(col, na.rm = TRUE); nn <- sum(!is.na(col))
    nn == 0 || tot == 0 || tot == nn * m_i[j]
  }, logical(1))
  ifit <- .item_fit(X, Z, mo, disc = if (is.null(disc)) NULL else disc_v,
                    extreme = person$extreme)
  pfit <- .person_fit(X, Z, mo, disc = if (is.null(disc)) NULL else disc_v,
                      item_extreme = item_extreme)
  n_par <- if (is.null(est$n_parameters)) nrow(est$thr) - 1L else est$n_parameters
  rf <- .fitres(Z, mo, person$extreme, n_par)
  ng_req <- n_groups                       # NULL = the automatic rule
  ci <- .class_intervals(person$theta, person$extreme, n_groups)
  n_groups <- attr(ci, "n_groups")
  # class intervals compiled per item when data are missing (the automatic
  # per-item basis; Andrich & Marais 2019, ch. 15), so every item is tested
  # over intervals of
  # its own responders (with the group-count rule applied per item when
  # automatic)
  ci_list <- if (anyNA(X))
    .class_intervals_by_item(X, person$theta, person$extreme, ng_req)
  else NULL
  it <- .item_trait(X, mo, ci, adjust_N = adjust_N, ci_list = ci_list)
  ia <- .item_anova(Z, ci, person$extreme, ci_list = ci_list)
  psi <- .psi(person$theta, person$se)
  psi_noext <- .psi(person$theta, person$se, keep = !person$extreme)
  alpha <- .alpha(X)

  # --- assembled person table ----------------------------------------------
  parts <- list(data.frame(id = id_vec), fac_df, person,
                data.frame(infit_ms = pfit$infit_ms, outfit_ms = pfit$outfit_ms,
                           outfit_z = pfit$outfit_z,
                           fit_resid = rf$persons$fit_resid,
                           natural_resid = rf$persons$natural,
                           df_fit = rf$persons$df, class_interval = ci))
  person <- do.call(cbind, parts[!vapply(parts, is.null, TRUE)])
  rownames(person) <- NULL

  # --- item table ----------------------------------------------------------
  loc <- vapply(tau_list, mean, 0)
  weak_thr <- if (is.null(thr$weak)) rep(FALSE, nrow(thr)) else thr$weak
  se_loc <- vapply(seq_len(L), function(i) {
    rows <- thr$id[thr$item == i]
    # a weakly determined threshold (sparse adjacent category) makes the
    # ridged covariance block spuriously small: report NA, not a number
    if (any(weak_thr[thr$item == i])) return(NA_real_)
    # anchored items have a structurally zero variance that floating-point
    # noise can render as a tiny negative number on some BLAS builds
    sqrt(max(mean(est$cov_tau[rows, rows]), 0))
  }, 0)
  items_df <- data.frame(item = colnames(X), max = m, location = loc,
                         se = se_loc,
                         disc = .item_discrim(person$theta, X, tau_list,
                                              person$extreme),
                         fit_resid = rf$items$fit_resid, df_fit = rf$items$df,
                         natural_resid = rf$items$natural,
                         infit_ms = ifit$infit_ms, outfit_ms = ifit$outfit_ms,
                         infit_z = ifit$infit_z, outfit_z = ifit$outfit_z,
                         chisq = it$chisq, df = it$df, p = it$p,
                         p_adj = it$p_adj, p_bonf = it$p_bonf,
                         F_anova = ia$F_anova, p_anova = ia$p)
  rownames(items_df) <- NULL

  # --- score table (complete responders; raw score is only sufficient when
  # --- discriminations are equal) ---------------------------------------------
  sc <- if (equal_disc) {
    pe_full <- person_wle(tau_list, disc = disc_v[1])
    data.frame(score = 0:sum(m), theta = unname(pe_full$theta),
               se = unname(pe_full$se))
  } else NULL
  if (is.null(sc))
    notes <- c(notes, "person measures use the weighted score; per-group score curves replace the raw-score table (see score_curves)")

  # --- threshold diagnostics --------------------------------------------------
  td <- lapply(seq_len(L), function(i) {
    tau_i <- tau_list[[i]]
    grid <- seq(-8, 8, by = 0.05)
    modal <- unique(vapply(grid, function(th)
      which.max(item_moments(th, tau_i, disc = disc_v[i])$P) - 1L, 1L))
    list(item = colnames(X)[i], thresholds = tau_i,
         ordered = all(diff(tau_i) > 0) || length(tau_i) == 1L,
         reversed_at = which(diff(tau_i) <= 0) + 1L,
         never_modal_categories = setdiff(0:length(tau_i), modal),
         category_counts = as.integer(table(factor(X[, i], levels = 0:length(tau_i)))))
  })
  names(td) <- colnames(X)

  out <- list(model = model, X = X, m = m, items = items_df, thresholds = thr,
              tau_list = tau_list, person = person, score_table = sc,
              residuals = Z, moments = mo, n_groups = n_groups,
              ci_item = ci_list,
              item_trait = it, item_anova = ia,
              psi = psi, psi_noext = psi_noext,
              isi = .psi(items_df$location, items_df$se),
              alpha = alpha,
              targeting = .targeting(person, thr),
              power_of_fit = .fit_power(psi$PSI),
              total_chisq = if (sum(it$df, na.rm = TRUE) > 0)
                sum(it$chisq, na.rm = TRUE) else NA_real_,
              total_df = if (sum(it$df, na.rm = TRUE) > 0)
                sum(it$df, na.rm = TRUE) else NA_integer_,
              total_chisq_p = if (sum(it$df, na.rm = TRUE) > 0)
                pchisq(sum(it$chisq, na.rm = TRUE),
                       sum(it$df, na.rm = TRUE), lower.tail = FALSE)
                else NA_real_,
              item_fit_summary = .dist_stats(rf$items$fit_resid),
              person_fit_summary = .dist_stats(rf$persons$fit_resid),
              summary_stats = list(
                item_location = .dist_stats(items_df$location),
                person_location = .dist_stats(person$theta),
                person_location_noext = .dist_stats(person$theta[!person$extreme]),
                cor_item_fit_location = .safe_cor(items_df$location,
                                                  rf$items$fit_resid),
                cor_person_fit_location = .safe_cor(person$theta,
                                                    rf$persons$fit_resid),
                df_factor = rf$f_cell),
              thresholds_diag = td, est = est, notes = notes,
              factors = fac_df, disc = disc)
  out <- .tag_tables(out)
  class(out) <- "rasch"
  out
}

#' @export
print.rasch <- function(x, ...) {
  cat(sprintf("rasch %s analysis: %d items, %d persons\n",
              x$model, ncol(x$X), nrow(x$X)))
  cat(sprintf("Pairwise conditional ML (Andrich & Luo): %s in %d iterations\n",
              if (x$est$converged) "converged" else "NOT converged",
              x$est$iterations))
  cat(sprintf("PSI %.3f (no extremes %.3f), item SI %.3f, alpha %.3f%s, power of fit: %s\n",
              x$psi$PSI, x$psi_noext$PSI, x$isi$PSI, x$alpha$alpha,
              if (isFALSE(x$alpha$applicable))
                sprintf(" [complete cases only, n = %d]", x$alpha$n) else "",
              x$power_of_fit))
  cat(sprintf("Total item-trait chi-square %.3f on %d df, p = %s\n",
              x$total_chisq, x$total_df, .fmt_p(x$total_chisq_p)))
  if (length(x$notes)) cat(sprintf("Notes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' @export
summary.rasch <- function(object, ...) {
  x <- object
  structural <- inherits(x, c("rasch_mfrm", "rasch_efrm"))
  unit <- if (structural) "Response-cell" else "Item"
  units <- if (structural) "response cells" else "items"
  print(x)
  cat(sprintf("\nTargeting: person mean %.3f (SD %.3f); %sthresholds span %.3f to %.3f\n",
              x$targeting$person_mean, x$targeting$person_sd,
              if (structural) "calibration " else "",
              x$targeting$threshold_range[1], x$targeting$threshold_range[2]))
  cat(sprintf("%s fit residual mean %.3f SD %.3f (skew %.2f, kurt %.2f); person fit residual mean %.3f SD %.3f (skew %.2f, kurt %.2f)\n",
              unit,
              x$item_fit_summary$mean, x$item_fit_summary$sd,
              x$item_fit_summary$skewness, x$item_fit_summary$kurtosis,
              x$person_fit_summary$mean, x$person_fit_summary$sd,
              x$person_fit_summary$skewness, x$person_fit_summary$kurtosis))
  cat(sprintf("Fit residual-location correlation: %s %.3f, persons %.3f; cell df factor %.3f\n",
              units,
              x$summary_stats$cor_item_fit_location,
              x$summary_stats$cor_person_fit_location,
              x$summary_stats$df_factor))
  cat(sprintf("%s with Holm-adjusted chi-square p < 0.05: %d of %d\n\n",
              if (structural) "Response cells" else "Items",
              sum(x$items$p_adj < 0.05, na.rm = TRUE), nrow(x$items)))
  core <- c("item", "max", "location", "se", "fit_resid", "infit_ms",
            "outfit_ms", "chisq", "df", "p_adj")
  print(.fmt_df(x$items[, intersect(core, names(x$items))]), row.names = FALSE)
  cat("(further columns on fit$items: natural and standardised forms,\n",
      " ANOVA fit, Bonferroni probabilities)\n", sep = "")
  dis <- vapply(x$thresholds_diag, function(d) !d$ordered, TRUE) &
    vapply(x$thresholds_diag, function(d) length(d$thresholds) > 1L, TRUE)
  if (any(dis)) cat(sprintf("\nDisordered %sthresholds: %s\n",
                            if (structural) "response-cell " else "",
                            paste(names(dis)[dis], collapse = ", ")))
  invisible(x)
}
