# rasch :: many-facet Rasch model
# ===========================================================================
# The many-facet Rasch model (Linacre 1989) estimated by the same pairwise
# conditional likelihood as the rest of the package. Each combination of an
# item with the levels of the rating facets (for example item x rater) is a
# "virtual item" whose thresholds decompose structurally as
#
#   tau_{(i, r), k} = delta_{ik} + rho_r (+ further facet terms)
#
# so the facet severities enter the design matrix exactly as the rating
# scale structure does, and the person parameter still cancels within every
# pair of virtual items. Identification: item thresholds sum to zero and
# each facet's severities sum to zero. The pairwise conditional likelihood
# is concave in the structural parameters, so Newton-Raphson from zero
# converges globally.
# ===========================================================================

# Pooled fit over a set of columns of the residual matrix (used for facet
# levels and for underlying items across their virtual columns). The fit
# residual is the log-of-mean-square statistic (Andrich & Marais 2019,
# ch. 23) pooled over the group's
# observed cells of non-extreme persons: Y^2 = sum z^2 against the summed
# cell degrees of freedom f = f_cell x cells, symmetrised as in .fitres.
# The published three-facet fit tables report a different margin statistic,
# the mean of the virtual items' fit residuals, which rasch_mfrm() computes
# alongside this pooled form (see its Details); EFRM frame fit keeps the
# pooled form. Infit and outfit mean squares are kept alongside.
.group_col_fit <- function(Z, mo, cols, disc = NULL, extreme = NULL,
                           f_cell = NA_real_) {
  E2 <- .z2_expectation(mo, Z, disc)
  sub <- Z[, cols, drop = FALSE]
  keep <- if (is.null(extreme)) rep(TRUE, nrow(Z)) else !extreme
  ok <- which(!is.na(sub) & keep)
  if (length(ok) < 3)
    return(list(infit_ms = NA_real_, outfit_ms = NA_real_,
                fit_resid = NA_real_, df_fit = NA_real_, n = length(ok)))
  z2 <- sub[ok]^2
  V <- mo$V[, cols, drop = FALSE][ok]
  C4 <- mo$M4[, cols, drop = FALSE][ok]
  e2 <- E2[, cols, drop = FALSE][ok]
  n <- length(ok)
  outfit <- sum(z2) / sum(e2)
  infit <- sum(z2 * V) / sum(e2 * V)
  fr <- df <- NA_real_
  if (is.finite(f_cell) && f_cell > 0) {
    y2 <- sum(z2); v <- sum(C4 / V^2 - 1); f <- f_cell * n
    if (v > 1e-8 && y2 > 0) {
      fr <- f * (log(y2) - log(f)) / sqrt(v)
      df <- f
    }
  }
  if (!is.finite(fr)) {                     # degenerate pooling: fall back
    qi <- sqrt(max(sum(C4 - V^2) / sum(V)^2, 1e-8))
    fr <- .wh(infit, qi)
  }
  list(infit_ms = infit, outfit_ms = outfit, fit_resid = fr, df_fit = df,
       n = n)
}

#' Fit a many-facet Rasch model
#'
#' Fits an additive many-facet Rasch model (Linacre 1989) to scored responses
#' indexed by person, item, and one or more facets such as rater, task, or
#' occasion. Facet severities, item thresholds, person locations, and fit
#' statistics are reported on a common logit scale.
#'
#' @details
#' For person \eqn{n}, item \eqn{i}, and facet levels
#' \eqn{f_1,\ldots,f_Q}, the additive model is
#' \deqn{P(X_{ni\mathbf{f}}=x)=\frac{\exp\{x\theta_n-
#'   \sum_{k=1}^{x}[\delta_{ik}+\sum_{q=1}^{Q}\rho_{qf_q}]\}}
#'   {\sum_{y=0}^{m_i}\exp\{y\theta_n-
#'   \sum_{k=1}^{y}[\delta_{ik}+\sum_{q=1}^{Q}\rho_{qf_q}]\}}.}
#' Positive facet values therefore denote greater severity. The item
#' thresholds have a common sum-zero origin and the levels of each facet sum
#' to zero. If \code{interaction} is requested, an item-by-level term is
#' added with both its item and facet margins constrained to sum to zero.
#'
#' Estimation represents each observed item-by-facet combination as a virtual
#' item and imposes the additive structure in the pairwise conditional
#' likelihood. The person parameter cancels before calibration. The covariance
#' of the structural parameters is the transformed Godambe sandwich covariance.
#'
#' Facet levels must be connected through common persons and items. A facet
#' nested within an item or a person-disjoint block can be confounded with the
#' item location. The function checks the structural rank and response graph
#' before fitting the model.
#'
#' An item-by-facet interaction retains equal discrimination but allows facet
#' differences to vary by item. The omnibus Wald test in
#' \code{interaction_test} is the primary test; cell tests are Holm-adjusted
#' follow-ups. Interaction probabilities require at least
#' \eqn{\max\{30,q+2\}} persons and effective persons in every observed
#' item-by-level cell, where \eqn{q} is the omnibus degrees of freedom.
#' The interaction covariance must also identify the omnibus contrast and
#' leave positive denominator degrees of freedom. Estimates remain descriptive
#' when these conditions are not met.
#'
#' @param data Long-format data frame, or a wide data frame when \code{items}
#'   is supplied.
#' @param person Name of the person identifier column. Person, item, score,
#'   facet, and person-factor columns must define distinct roles.
#' @param item Name of the item column.
#' @param score Name of the integer score column (categories from 0; gaps are
#'   collapsed per item with a note).
#' @param facets Character vector naming one or more facet columns (for
#'   example a rater column).
#' @param n_groups Number of class intervals for the item-trait chi-square;
#'   \code{NULL} (the default) applies the class-interval rule of Andrich and
#'   Marais (2019, ch. 15) (at least 50
#'   non-extreme persons per interval, at most 10 intervals, at least 2).
#' @param na_codes Score values to read as missing (default \code{-1}); any
#'   negative score is also treated as missing.
#' @param items Optional character vector of item score columns for data in
#'   wide format: one row per person-by-facet combination (for example one
#'   row per script per rater) with one column per item or criterion. The
#'   long form (\code{item} + \code{score}) remains available for data
#'   where the facet varies within items.
#' @param interaction Optional name of one facet to interact with the items
#'   (interactive facet mode). See Details.
#' @param factors Optional person factors for DIF analysis: a character
#'   vector naming columns constant within person, or a data frame with one
#'   row per data row or unique person. Facets belong in \code{facets}, not
#'   here.
#' @param maxit,tol Newton-Raphson iteration cap and convergence tolerance.
#' @return An object of classes \code{"rasch_mfrm"} and \code{"rasch"}.
#'   Model-specific components describe the facets, items, thresholds, and
#'   facet specification. Interactive fits also contain an omnibus test and
#'   the corresponding item-by-facet effects. The component \code{fit_resid}
#'   averages virtual-item residuals within a margin. Its response-weighted
#'   counterpart is \code{fit_resid_pooled}; its degrees of freedom are in
#'   \code{df_fit}. A non-converged fit retains estimates and residual patterns
#'   for diagnosis but withholds standard errors and inferential probabilities.
#' @references
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#'
#' Linacre, J. M. (1989). Many-Facet Rasch Measurement. Chicago: MESA Press.
#' @seealso \code{\link{rasch}}, \code{\link{rasch_efrm}},
#'   \code{\link{dif_anova}}, and \code{\link{simulate_mfrm}}.
#' @examples
#' set.seed(1)
#' simP <- function(th, tau) {
#'   x <- 0:length(tau)
#'   p <- exp(x * th - c(0, cumsum(tau)))
#'   p / sum(p)
#' }
#' persons <- sprintf("P%03d", 1:120); raters <- paste0("R", 1:4)
#' th <- setNames(rnorm(120, 0, 1.3), persons)
#' rho <- setNames(c(-0.6, -0.2, 0.2, 0.6), raters)
#' tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
#' d <- expand.grid(person = persons, item = names(tau), rater = raters,
#'                  stringsAsFactors = FALSE)
#' d$score <- mapply(function(p, i, r)
#'   sample(0:2, 1, prob = simP(th[p], tau[[i]] + rho[r])),
#'   d$person, d$item, d$rater)
#' fit <- rasch_mfrm(d, person = "person", item = "item", score = "score",
#'                   facets = "rater")
#' fit$facet_effects$rater
#' @export
rasch_mfrm <- function(data, person, item = NULL, score = NULL, facets,
                       items = NULL, n_groups = NULL,
                       na_codes = -1, interaction = NULL,
                       factors = NULL, maxit = 60, tol = 1e-8) {
  .check_controls(maxit, tol)
  if (!is.data.frame(data))
    stop("`data` must be a data frame in long or wide form", call. = FALSE)
  if (!is.null(n_groups))
    n_groups <- .check_whole(n_groups, "n_groups", 2)
  .check_column_names(data)
  # the person column is dereferenced by BOTH entry forms, so it is resolved
  # to one existing column before either of them runs; item and score are
  # the long form's own and are checked on that path
  .check_reshape_column(data, person, "person")
  if (!is.character(facets) || !is.null(dim(facets)) ||
      !is.null(oldClass(facets)) || !length(facets) || anyNA(facets) ||
      any(!nzchar(facets)))
    stop("`facets` must name at least one facet column")
  if (anyDuplicated(facets))
    stop("facet column(s) named more than once: ",
         paste(unique(facets[duplicated(facets)]), collapse = ", "))
  if (person %in% facets)
    stop("the person column cannot also be a facet")
  if (!is.null(factors) &&
      (!(is.character(factors) && is.null(dim(factors)) &&
         is.null(oldClass(factors))) && !is.data.frame(factors)))
    stop("`factors` must be a plain character vector of column names or a data frame",
         call. = FALSE)
  factor_names <- if (is.character(factors)) factors else
    if (is.data.frame(factors)) names(factors) else character(0)
  if (length(factor_names) &&
      (anyNA(factor_names) || any(!nzchar(factor_names))))
    stop("every person factor needs a non-empty column name")
  if (anyDuplicated(factor_names))
    stop("duplicate factor column name(s): ",
         paste(unique(factor_names[duplicated(factor_names)]), collapse = ", "))
  # wide entry: item score columns are melted to the long form internally
  if (!is.null(items)) {
    if (!is.null(item) || !is.null(score))
      stop("give either `items` (wide: one column per item) or `item` + `score` (long)")
    if (!is.character(items) || !is.null(dim(items)) ||
        !is.null(oldClass(items)) || !length(items) || anyNA(items) ||
        any(!nzchar(items)))
      stop("`items` must name at least one item column")
    if (anyDuplicated(items))
      stop("item column(s) named more than once: ",
           paste(unique(items[duplicated(items)]), collapse = ", "))
    overlap <- intersect(items, c(person, facets))
    if (length(overlap))
      stop("wide item column(s) cannot also be the person or a facet: ",
           paste(overlap, collapse = ", "))
    factor_overlap <- intersect(factor_names, c(person, facets, items))
    if (length(factor_overlap))
      stop("person-factor column(s) cannot also define another model role: ",
           paste(factor_overlap, collapse = ", "))
    miss <- setdiff(c(person, facets, items), names(data))
    if (length(miss)) stop("column(s) not in data: ", paste(miss, collapse = ", "))
    taken <- unique(c(names(data), factor_names))
    temp_name <- function(base) {
      out <- base
      while (out %in% taken) out <- paste0(out, ".")
      taken <<- c(taken, out)
      out
    }
    tmp_person <- temp_name("..person")
    tmp_item <- temp_name("..item")
    tmp_score <- temp_name("..score")
    long <- data.frame(
      rep(.role_text_values(data[[person]]), length(items)),
      rep(items, each = nrow(data)),
      unlist(lapply(items, function(cn) {
        v0 <- as.character(data[[cn]])
        v <- suppressWarnings(as.numeric(v0))
        bad <- !is.na(v0) & is.na(v)
        if (any(bad))
          stop("non-numeric score(s) in item column ", cn, " (e.g. '",
               v0[bad][1], "'); scores must be integer counts", call. = FALSE)
        v
      })), check.names = FALSE,
      stringsAsFactors = FALSE)
    names(long) <- c(tmp_person, tmp_item, tmp_score)
    for (f in facets)
      long[[f]] <- rep(.role_text_values(data[[f]]), length(items))
    # person factors survive the melt: named columns are replicated like
    # facets, a data frame is replicated row-wise to match the long rows
    fac_pass <- factors
    if (is.character(factors)) {
      missf <- setdiff(factors, names(data))
      if (length(missf))
        stop("factor column(s) not found in the data: ",
             paste(missf, collapse = ", "))
      for (cn in factors)
        long[[cn]] <- rep(.role_text_values(data[[cn]]), length(items))
    } else if (is.data.frame(factors)) {
      if (anyDuplicated(names(factors)))
        stop("duplicate factor column name(s): ",
             paste(unique(names(factors)[duplicated(names(factors))]),
                   collapse = ", "))
      persons_row <- .role_text_values(data[[person]])
      pu <- unique(persons_row)
      if (nrow(factors) == nrow(data)) {
        row_idx <- seq_len(nrow(data))
      } else if (nrow(factors) == length(pu)) {
        # documented alternative: one row per unique person -- map each
        # data row to its person's factor row
        row_idx <- match(persons_row, pu)
      } else {
        stop("`factors` data frame needs one row per data row (", nrow(data),
             ") or one per unique person (", length(pu), ")")
      }
      fac_pass <- factors[rep(row_idx, length(items)), , drop = FALSE]
      rownames(fac_pass) <- NULL
    }
    return(rasch_mfrm(long, person = tmp_person, item = tmp_item,
                      score = tmp_score, facets = facets,
                      n_groups = n_groups,
                      na_codes = na_codes, interaction = interaction,
                      factors = fac_pass, maxit = maxit, tol = tol))
  }
  if (is.null(item) || is.null(score))
    stop("give either `items` (wide) or `item` + `score` (long)")
  # the long form dereferences these, so each names one existing column
  .check_reshape_column(data, item, "item")
  .check_reshape_column(data, score, "score")
  roles <- c(person, item, score, facets)
  if (anyDuplicated(roles))
    stop("person, item, score, and facet columns must be distinct; repeated: ",
         paste(unique(roles[duplicated(roles)]), collapse = ", "))
  factor_overlap <- intersect(factor_names, roles)
  if (length(factor_overlap))
    stop("person-factor column(s) cannot also define another model role: ",
         paste(factor_overlap, collapse = ", "))
  if (!is.null(interaction)) {
    if (!is.character(interaction) || !is.null(dim(interaction)) ||
        !is.null(oldClass(interaction)) || length(interaction) != 1L ||
        is.na(interaction))
      stop("'interaction' must name exactly one facet")
    interaction <- as.character(interaction)
    if (!interaction %in% facets)
      stop("'interaction' must name one of the facets")
  }
  need <- c(person, item, score, facets)
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("column(s) not in data: ", paste(miss, collapse = ", "))
  notes <- character(0)

  pid <- .role_text_values(data[[person]])
  itm <- .role_text_values(data[[item]])
  .check_integer_scores(data[[score]], "the score column")
  sc <- suppressWarnings(as.integer(as.character(data[[score]])))
  n_na <- sum(!is.na(sc) & (sc %in% na_codes | sc < 0))
  sc[sc %in% na_codes | (!is.na(sc) & sc < 0)] <- NA
  if (n_na > 0)
    notes <- c(notes, sprintf("%d response(s) with a missing-data code (%s) set to missing",
                              n_na, paste(unique(c(na_codes, "negative")), collapse = ", ")))
  if (all(is.na(sc))) stop("score column has no usable integer values")
  fac <- lapply(facets, function(f) .role_text_values(data[[f]]))
  names(fac) <- facets
  # a whitespace-only label is not an identifier. Unlike NA, which is
  # dropped with a note below, it silently becomes a level of its own: a
  # blank person joins the sample as another respondent, and a blank item
  # or facet level is calibrated alongside the real ones
  cand <- c(stats::setNames(list(pid), person),
            stats::setNames(list(itm), item), fac)
  blanks <- vapply(cand, function(v) any(!is.na(v) & !nzchar(trimws(v))), TRUE)
  if (any(blanks))
    stop("blank identifier(s) in column(s): ",
         paste(names(blanks)[blanks], collapse = ", "),
         "; a whitespace-only label is not a person, item, or facet level")

  # a missing person, item, or facet identifier cannot be attached to any
  # virtual item: paste() would otherwise coerce NA to the literal string
  # "NA", silently pooling unrelated rows under a phantom level
  bad_id <- is.na(pid) | is.na(itm) |
    Reduce(`|`, lapply(fac, is.na), rep(FALSE, length(pid)))
  if (any(bad_id)) {
    notes <- c(notes, sprintf(
      "%d row(s) dropped: missing person, item, or facet identifier",
      sum(bad_id)))
    pid <- pid[!bad_id]; itm <- itm[!bad_id]; sc <- sc[!bad_id]
    fac <- lapply(fac, `[`, !bad_id)
    if (!length(pid)) stop("no rows remain after dropping rows with ",
                           "missing person/item/facet identifiers")
  }

  # rescore each item to consecutive categories from 0
  items_u <- sort(unique(itm))
  item_m <- setNames(integer(length(items_u)), items_u)
  for (it in items_u) {
    sel <- which(itm == it & !is.na(sc))
    obs <- sort(unique(sc[sel]))
    if (length(obs) < 2) stop("item ", it, " has fewer than two observed categories")
    full <- seq(0L, max(obs))
    if (!identical(obs, full)) {
      sc[sel] <- match(sc[sel], obs) - 1L
      notes <- c(notes, sprintf("item %s rescored: observed categories [%s] mapped to 0:%d",
                                it, paste(obs, collapse = ","), length(obs) - 1L))
    }
    item_m[it] <- length(sort(unique(sc[sel]))) - 1L
  }

  # virtual items: item x facet-level combinations present in the data
  fkey <- as.character(.factor_cells(fac, sep = ":"))
  vkey <- as.character(.factor_cells(data.frame(item = itm, cell = fkey),
                                     sep = ":"))
  vlev <- unique(vkey[order(match(itm, items_u), fkey)])
  vmap <- data.frame(vkey = vlev,
                     item = itm[match(vlev, vkey)],
                     stringsAsFactors = FALSE)
  for (f in facets) vmap[[f]] <- fac[[f]][match(vlev, vkey)]

  persons_u <- unique(pid)
  Xv <- matrix(NA_integer_, length(persons_u), length(vlev),
               dimnames = list(NULL, vlev))
  ri <- match(pid, persons_u); cj <- match(vkey, vlev)
  dup <- duplicated(cbind(ri, cj))
  if (any(dup))
    stop(sum(dup), " duplicate person-by-item-by-facet response(s): the ",
         "design has one cell per combination, so duplicates are ",
         "ambiguous (keeping the first would make results depend on row ",
         "order) -- aggregate or de-duplicate explicitly")
  use <- !dup & !is.na(sc)
  Xv[cbind(ri[use], cj[use])] <- sc[use]

  m_v <- item_m[vmap$item]
  thr_v <- threshold_index(m_v)

  # --- structural design matrix --------------------------------------------
  thr_items <- threshold_index(item_m[items_u])     # delta enumeration
  Md <- nrow(thr_items)
  A_delta <- rbind(diag(Md - 1L), rep(-1, Md - 1L))
  flevs <- lapply(facets, function(f) sort(unique(fac[[f]])))
  names(flevs) <- facets
  A_fac <- lapply(flevs, function(lv) {
    Lf <- length(lv)
    if (Lf < 2) stop("a facet needs at least two levels")
    rbind(diag(Lf - 1L), rep(-1, Lf - 1L))
  })
  Li <- length(items_u)
  A_item <- rbind(diag(Li - 1L), rep(-1, Li - 1L))   # item-location margins
  n_gam <- if (is.null(interaction)) 0L else
    (Li - 1L) * (length(flevs[[interaction]]) - 1L)
  P <- (Md - 1L) + sum(vapply(A_fac, ncol, 1L)) + n_gam
  B <- matrix(0, nrow(thr_v), P)
  for (row in seq_len(nrow(thr_v))) {
    v <- thr_v$item[row]; k <- thr_v$k[row]
    iu <- match(vmap$item[v], items_u)
    drow <- which(thr_items$item == iu & thr_items$k == k)
    B[row, seq_len(Md - 1L)] <- A_delta[drow, ]
    cursor <- Md - 1L
    for (f in facets) {
      lev <- match(vmap[[f]][v], flevs[[f]])
      nc <- ncol(A_fac[[f]])
      B[row, cursor + seq_len(nc)] <- A_fac[[f]][lev, ]
      cursor <- cursor + nc
    }
    if (n_gam > 0L) {
      lev <- match(vmap[[interaction]][v], flevs[[interaction]])
      B[row, cursor + seq_len(n_gam)] <-
        as.vector(outer(A_item[iu, ], A_fac[[interaction]][lev, ]))
    }
  }

  # concave likelihood: Newton-Raphson from zero with step halving
  # the structural design must have full column rank: a facet level that
  # never shares items (or persons) with the others is confounded with
  # the item parameters, and the ridged solve would return a
  # valid-looking but unidentified decomposition
  qrB <- qr(B)
  if (qrB$rank < ncol(B))
    stop("the facet design is structurally unidentified (rank ",
         qrB$rank, " for ", ncol(B), " parameters): some facet level(s) ",
         "are confounded with the items or with each other -- every ",
         "facet level needs items (and persons) in common with the rest")
  # algebraic rank of B is necessary but not sufficient: the pairwise
  # conditional likelihood carries information only within blocks of
  # virtual items that share persons, so a relative shift between
  # person-disjoint blocks is a flat direction of the likelihood whenever
  # the structural map can express it -- the solve would then land
  # wherever the ridge sends it while reporting convergence
  # informative co-observation only: a person at the extreme total of a
  # pair (both responses 0, or both at their maxima) has one feasible
  # conditional allocation and links nothing
  obs_m <- !is.na(Xv)
  zeros <- obs_m & !is.na(Xv) & Xv == 0
  maxs <- obs_m & sweep(Xv, 2, m_v, `==`); maxs[is.na(maxs)] <- FALSE
  zeros[is.na(zeros)] <- FALSE
  co_obs <- (crossprod(obs_m) - crossprod(zeros) - crossprod(maxs)) > 0
  edges_v <- which(co_obs & upper.tri(co_obs), arr.ind = TRUE)
  comp_v <- .btlef_components(ncol(Xv), edges_v)
  if (length(unique(comp_v)) > 1L) {
    comps_u <- sort(unique(comp_v))
    U <- vapply(comps_u, function(cc) as.numeric(comp_v[thr_v$item] == cc),
                numeric(nrow(thr_v)))
    if (qr(cbind(B, U))$rank < qrB$rank + qr(U)$rank) {
      blocks <- vapply(comps_u, function(cc)
        paste(sort(unique(vmap$item[comp_v == cc])), collapse = ", "), "")
      stop("the response design is disconnected and the facet structure ",
           "does not bridge it: no person links the blocks {",
           paste(blocks, collapse = "} and {"), "}, so their relative ",
           "locations are unidentified -- link the blocks through common ",
           "persons")
    }
  }
  sol <- .pcml_solve(Xv, thr_v, m_v, B, rep(0, P), maxit = maxit, tol = tol)
  if (!isTRUE(sol$converged))
    warning("MFRM estimation did NOT converge in ", sol$iterations,
            " iterations; increase maxit or inspect the design",
            call. = FALSE)

  thr_v$tau <- sol$tau; thr_v$se <- sol$se_tau; thr_v$anchored <- FALSE
  # a virtual item threshold resting on a near-empty category is a boundary
  # artefact: flag it and report its SE as NA, the same honesty rasch()/
  # pcml() apply -- the facet decomposition does not exempt it
  weak <- .pcml_weak_thresholds(Xv, m_v, thr_v, colnames(Xv))
  thr_v$weak <- weak$flag
  thr_v$se[weak$flag] <- NA_real_
  if (length(weak$notes)) notes <- c(notes, weak$notes)
  est <- list(model = "MFRM", thr = thr_v, cov_tau = sol$cov_tau,
              loglik = sol$loglik, iterations = sol$iterations,
              converged = sol$converged, m = m_v, anchors = NULL,
              n_parameters = P)

  # person factors for DIF: columns of `data` (constant within person) or a
  # data frame keyed to the unique persons, carried through so dif_anova()
  # works on an MFRM fit the way it does on any rasch fit. Facets are NOT
  # person factors -- facet DIF is an item-by-facet interaction and belongs
  # to `interaction=`.
  fac_df <- NULL
  if (!is.null(factors)) {
    if (is.character(factors)) {
      miss <- setdiff(factors, names(data))
      if (length(miss))
        stop("factor column(s) not found in the data: ",
             paste(miss, collapse = ", "))
      fac_df <- as.data.frame(lapply(factors, function(cn) {
        v <- .role_text_values(data[[cn]])[!bad_id]
        nvar <- tapply(v, pid, function(x) length(unique(x[!is.na(x)])))
        if (any(nvar > 1L, na.rm = TRUE))
          stop("factor '", cn, "' varies within person(s) ",
               paste(names(nvar)[which(nvar > 1L)], collapse = ", "),
               ": person factors must be constant per person (a facet ",
               "is not a person factor; see `interaction=`)")
        vv <- tapply(v, pid, function(x) {
          z <- x[!is.na(x)]
          if (length(z)) z[1L] else NA_character_
        })
        unname(vv[match(persons_u, names(vv))])
      }), col.names = factors, stringsAsFactors = FALSE)
      names(fac_df) <- factors
    } else {
      fac_df <- as.data.frame(factors, stringsAsFactors = FALSE)
      # a row-per-data-row frame is collapsed to one row per person: a
      # column that varies within a person has no person-level value, and
      # keeping the first row's would silently pick one occasion
      .check_person_constant <- function(df, key) {
        for (cn in names(df)) {
          v <- .role_text_values(df[[cn]])
          nvar <- tapply(v, key, function(x) length(unique(x[!is.na(x)])))
          if (any(nvar > 1L, na.rm = TRUE))
            stop("factor '", cn, "' varies within person(s) ",
                 paste(names(nvar)[which(nvar > 1L)], collapse = ", "),
                 ": person factors must be constant per person (a facet ",
                 "is not a person factor; see `interaction=`)")
        }
      }
      if (nrow(fac_df) == length(bad_id)) {
        # one row per ORIGINAL data row: rows dropped for missing
        # identifiers drop from the factors too, keeping them aligned
        fac_df <- fac_df[!bad_id, , drop = FALSE]
        .check_person_constant(fac_df, pid)
        fac_df <- fac_df[match(persons_u, pid), , drop = FALSE]
      } else if (nrow(fac_df) == length(pid)) {
        .check_person_constant(fac_df, pid)
        fac_df <- fac_df[match(persons_u, pid), , drop = FALSE]
      } else if (nrow(fac_df) != length(persons_u))
        stop("`factors` needs one row per data row or one per unique ",
             "person (", length(persons_u), ")")
      rownames(fac_df) <- NULL
    }
  }
  .check_factor_frame(fac_df)
  fit <- .assemble_fit("MFRM", Xv, est, persons_u, fac_df, n_groups, notes)
  # When an item is represented by several facet cells, the expanded
  # columns are not one administered item set. Alpha and a universal
  # raw-score conversion over those columns have no test-level interpretation.
  # Retain both for the one-cell-per-item reduction, which is an ordinary
  # administered matrix despite having been fitted through this interface.
  expanded_cells <- any(table(vmap$item) > 1L)
  if (expanded_cells) {
    fit$alpha <- list(
      alpha = NA_real_, n = NA_integer_, applicable = FALSE,
      design_applicable = FALSE,
      reason = "not applicable when an item has several facet response cells")
    fit$score_table <- NULL
    fit$notes <- unique(c(fit$notes, paste(
      "a universal raw-score conversion is not defined across the expanded",
      "facet response cells; use the design-specific information curves")))
  } else fit$alpha$design_applicable <- TRUE

  # --- structural effects -----------------------------------------------------
  covb <- sol$cov_beta
  d_idx <- seq_len(Md - 1L)
  delta <- drop(A_delta %*% sol$beta[d_idx])
  cov_d <- A_delta %*% covb[d_idx, d_idx, drop = FALSE] %*% t(A_delta)
  fit$item_thresholds <- data.frame(item = items_u[thr_items$item],
                                    k = thr_items$k, tau = delta,
                                    se = sqrt(pmax(diag(cov_d), 0)))
  item_fit <- lapply(items_u, function(it)
    .group_col_fit(fit$residuals, fit$moments, which(vmap$item == it),
                    extreme = fit$person$extreme,
                    f_cell = fit$summary_stats$df_factor))
  # The facet-margin fit residual of the published three-facet tables
  # (Andrich & Marais 2019, ch. 26 and app. C) is the MEAN of the
  # constituent virtual items' fit residuals; the pooled log-residual over
  # the margin's cells is kept alongside with its degrees of freedom.
  vmean <- function(sel) {
    z <- fit$items$fit_resid[sel]
    z <- z[is.finite(z)]
    if (length(z)) mean(z) else NA_real_
  }
  fit$item_effects <- data.frame(
    item = items_u,
    location = vapply(seq_along(items_u), function(i)
      mean(delta[thr_items$item == i]), 0),
    se = vapply(seq_along(items_u), function(i) {
      rows <- which(thr_items$item == i)
      sqrt(max(mean(cov_d[rows, rows]), 0))
    }, 0),
    n = vapply(item_fit, `[[`, 0, "n"),
    infit_ms = vapply(item_fit, `[[`, 0, "infit_ms"),
    outfit_ms = vapply(item_fit, `[[`, 0, "outfit_ms"),
    fit_resid = vapply(items_u, function(it) vmean(vmap$item == it), 0),
    fit_resid_pooled = vapply(item_fit, `[[`, 0, "fit_resid"),
    df_fit = vapply(item_fit, `[[`, 0, "df_fit"))

  cursor <- Md - 1L
  fit$facet_effects <- list()
  for (f in facets) {
    nc <- ncol(A_fac[[f]]); idx <- cursor + seq_len(nc); cursor <- cursor + nc
    rho <- drop(A_fac[[f]] %*% sol$beta[idx])
    cov_r <- A_fac[[f]] %*% covb[idx, idx, drop = FALSE] %*% t(A_fac[[f]])
    lev_fit <- lapply(flevs[[f]], function(lv)
      .group_col_fit(fit$residuals, fit$moments, which(vmap[[f]] == lv),
                      extreme = fit$person$extreme,
                      f_cell = fit$summary_stats$df_factor))
    fit$facet_effects[[f]] <- data.frame(
      level = flevs[[f]], severity = rho,
      se = sqrt(pmax(diag(cov_r), 0)),
      n = vapply(lev_fit, `[[`, 0, "n"),
      infit_ms = vapply(lev_fit, `[[`, 0, "infit_ms"),
      outfit_ms = vapply(lev_fit, `[[`, 0, "outfit_ms"),
      fit_resid = vapply(flevs[[f]], function(lv) vmean(vmap[[f]] == lv), 0),
      fit_resid_pooled = vapply(lev_fit, `[[`, 0, "fit_resid"),
      df_fit = vapply(lev_fit, `[[`, 0, "df_fit"))
  }
  if (n_gam > 0L) {
    idx <- cursor + seq_len(n_gam)
    R0 <- length(flevs[[interaction]])
    K <- A_fac[[interaction]] %x% A_item        # vec(gamma) = K beta_gamma
    gvec <- drop(K %*% sol$beta[idx])
    cov_g <- K %*% covb[idx, idx, drop = FALSE] %*% t(K)
    fit$interaction_effects <- data.frame(
      item = rep(items_u, R0),
      level = rep(flevs[[interaction]], each = Li),
      gamma = gvec, se = sqrt(pmax(diag(cov_g), 0)))
    # Inferential support is set by the least-observed item-by-level cell, not
    # by the facet-level or total calibration sample. A sparse interaction
    # cell cannot borrow denominator degrees of freedom from observations of
    # the same facet level on other items. Response counts supply Kish
    # weights, so highly unequal coverage also reduces the effective number
    # of persons.
    cell_support <- lapply(flevs[[interaction]], function(lv)
      lapply(items_u, function(it) {
        cc <- which(vmap[[interaction]] == lv & vmap$item == it)
        nr <- if (length(cc)) rowSums(!is.na(fit$X[, cc, drop = FALSE])) else
          rep(0, nrow(fit$X))
        nr[fit$person$extreme] <- 0
        ww <- nr[nr > 0]
        data.frame(item = it, level = lv, n_persons = length(ww),
                   effective_persons = if (length(ww))
                     sum(ww)^2 / sum(ww^2) else 0,
                   stringsAsFactors = FALSE)
      }))
    fit$interaction_support <- do.call(rbind, unlist(cell_support,
                                                      recursive = FALSE))
    q_int <- length(sol$beta[idx])
    min_required <- max(30L, q_int + 2L)
    fit$interaction_support$minimum_required <- min_required
    support_ok <- all(fit$interaction_support$n_persons >= min_required &
      fit$interaction_support$effective_persons >=
        min_required - sqrt(.Machine$double.eps))
    n_units <- floor(min(fit$interaction_support$effective_persons))

    # inferential reference: the sandwich covariance is ESTIMATED from the
    # persons' score contributions, so a chi-square reference for the
    # multi-degree-of-freedom Wald is anticonservative in realistic samples
    # (the Hotelling effect: with n persons and q parameters the statistic
    # behaves as a scaled F, not chi-square; a null simulation at n = 50
    # showed ~13% rejection at nominal 5% under the chi-square reference).
    # Use the T-squared-style F reference with persons as the units, and a
    # t reference for the per-cell follow-ups.
    fit$interaction_effects$z <- .wald_ratio(
      fit$interaction_effects$gamma, fit$interaction_effects$se)
    fit$interaction_effects$p <- if (support_ok)
      2 * stats::pt(-abs(fit$interaction_effects$z),
                    df = max(n_units - 1L, 1L)) else NA_real_
    fit$interaction_effects$p_adj <- .p_adjust_family(
      fit$interaction_effects$p, method = "holm")
    fit$interaction_effects$significant <-
      fit$interaction_effects$p_adj < 0.05
    bg <- sol$beta[idx]
    Vg <- covb[idx, idx, drop = FALSE]
    Wg <- if (.covariance_is_psd(Vg))
      tryCatch(drop(t(bg) %*% solve(Vg) %*% bg),
               error = function(e) NA_real_) else NA_real_
    if (is.finite(Wg) && Wg < 0) Wg <- NA_real_
    q_int <- length(bg)
    test_ok <- support_ok && is.finite(Wg) && n_units > q_int + 1L
    if (test_ok) {
      Fg <- Wg * (n_units - q_int) / (q_int * (n_units - 1L))
      pg <- stats::pf(Fg, q_int, n_units - q_int, lower.tail = FALSE)
    } else { Fg <- NA_real_; pg <- NA_real_ }
    fit$interaction_test <- data.frame(
      facet = interaction, df = q_int, wald = Wg,
      f = Fg, df2 = if (test_ok) n_units - q_int else NA_real_,
      p = pg, min_effective_persons = n_units,
      minimum_required = min_required,
      inference_available = test_ok)
    if (!support_ok) fit$notes <- unique(c(fit$notes, sprintf(
      "the %s interaction estimates are descriptive because at least one item-by-level cell has fewer than %d persons or effective persons; probabilities are withheld",
      interaction, min_required)))
    if (support_ok && !test_ok) fit$notes <- unique(c(fit$notes, sprintf(
      "the %s interaction omnibus is unavailable because its estimated covariance is singular or does not leave positive denominator degrees of freedom",
      interaction)))
    fit$interaction <- interaction
  }
  fit$facet_spec <- facets
  fit$virtual_map <- vmap
  # Retain the public fitting contract needed to reproduce this structural
  # calibration from a generated virtual response matrix.  The response
  # matrix itself and the person factors already live on `fit`; keeping the
  # roles and controls here avoids trying to infer them from display tables.
  fit$refit_spec <- list(
    facets = facets, interaction = interaction,
    n_groups = n_groups, na_codes = na_codes,
    maxit = maxit, tol = tol)
  if (!isTRUE(sol$converged)) {
    for (f in names(fit$facet_effects)) fit$facet_effects[[f]]$se[] <- NA_real_
    if (!is.null(fit$interaction_effects)) {
      for (nm in intersect(c("se", "z", "p", "p_adj"),
                           names(fit$interaction_effects)))
        fit$interaction_effects[[nm]][] <- NA_real_
      fit$interaction_effects$significant[] <- NA
    }
    if (!is.null(fit$interaction_test)) {
      for (nm in intersect(c("wald", "f", "df2", "p"),
                           names(fit$interaction_test)))
        fit$interaction_test[[nm]][] <- NA_real_
      fit$interaction_test$inference_available[] <- FALSE
    }
  }
  fit <- .tag_tables(fit)
  class(fit) <- c("rasch_mfrm", "rasch")
  fit
}

#' @export
print.rasch_mfrm <- function(x, ...) {
  separation_quality <- x$separation_quality %||% x$power_of_fit %||%
    .separation_quality(x$psi$PSI)
  cat(sprintf("rasch multiple ratings analysis: %d items x %s = %d response cells, %d persons\n",
              nrow(x$item_effects),
              paste(vapply(x$facet_spec, function(f)
                sprintf("%d %s level(s)", nrow(x$facet_effects[[f]]), f), ""),
                collapse = " x "),
              ncol(x$X), nrow(x$X)))
  cat(sprintf("Pairwise conditional ML: %s in %d iterations\n",
              if (x$est$converged) "converged" else "NOT converged",
              x$est$iterations))
  cat(sprintf("PSI %.3f, separation quality: %s\n", x$psi$PSI,
              separation_quality))
  for (f in x$facet_spec) {
    fe <- x$facet_effects[[f]]
    core <- c("level", "severity", "se", "n", "fit_resid")
    cat(sprintf("\nFacet '%s' severities (logits):\n", f))
    print(.fmt_df(fe[, intersect(core, names(fe))]), row.names = FALSE)
  }
  cat("(pooled fit residuals and their df on fit$facet_effects)\n")
  if (!is.null(x$interaction)) {
    it <- x$interaction_test
    cat(sprintf("\nItem-by-%s omnibus test: Wald %.3f -> F(%d, %d) = %.3f, p = %s\n",
                x$interaction, it$wald, it$df, it$df2,
                if (is.finite(it$f)) it$f else NA, .fmt_p(it$p)))
    big <- x$interaction_effects[x$interaction_effects$significant %in% TRUE,
                                 , drop = FALSE]
    cat(sprintf("Holm-adjusted exploratory cells: %d significant of %d\n",
                nrow(big), nrow(x$interaction_effects)))
    if (nrow(big)) print(big, digits = 3, row.names = FALSE)
  }
  if (length(x$notes)) cat("\nNotes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}

#' Plot facet severities
#'
#' Caterpillar plot of the severity of each level of a facet from a
#' many-facet analysis, with 95 per cent error bars; levels with pooled fit
#' residuals beyond the band are highlighted.
#'
#' @param fit A fitted object from \code{\link{rasch_mfrm}}.
#' @param facet Facet name; defaults to the first facet.
#' @param band Fit residual band beyond which a level is highlighted.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' set.seed(1)
#' simP <- function(th, tau) {
#'   x <- 0:length(tau)
#'   p <- exp(x * th - c(0, cumsum(tau)))
#'   p / sum(p)
#' }
#' persons <- sprintf("P%03d", 1:120); raters <- paste0("R", 1:4)
#' th <- setNames(rnorm(120, 0, 1.3), persons)
#' rho <- setNames(c(-0.6, -0.2, 0.2, 0.6), raters)
#' tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
#' d <- expand.grid(person = persons, item = names(tau), rater = raters,
#'                  stringsAsFactors = FALSE)
#' d$score <- mapply(function(p, i, r)
#'   sample(0:2, 1, prob = simP(th[p], tau[[i]] + rho[r])),
#'   d$person, d$item, d$rater)
#' plot_facets(rasch_mfrm(d, "person", "item", "score", facets = "rater"))
#' }
#' @export
plot_facets <- function(fit, facet = NULL, band = 2.5) {
  if (!inherits(fit, "rasch_mfrm")) stop("plot_facets needs a rasch_mfrm fit")
  .check_response_display_fit(fit, "facet plots")
  .check_band(band)
  if (!is.null(facet) && (!is.character(facet) || !is.null(dim(facet)) ||
                          !is.null(oldClass(facet)) || length(facet) != 1L ||
                          is.na(facet) || !nzchar(trimws(facet))))
    stop("`facet` must be one non-empty facet name")
  if (is.null(facet)) facet <- fit$facet_spec[1]
  fe <- fit$facet_effects[[facet]]
  if (is.null(fe)) stop("no such facet: ", facet)
  fe <- fe[order(fe$severity), ]
  lo <- fe$severity - 1.96 * fe$se; hi <- fe$severity + 1.96 * fe$se
  point_ok <- is.finite(fe$severity)
  if (!any(point_ok))
    .refuse("this facet has no finite severity estimates to display")
  ci_ok <- point_ok & is.finite(lo) & is.finite(hi)
  n <- nrow(fe)
  op <- par(mar = c(4.2, 7.5, 3.2, 1.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  xr <- range(c(fe$severity[point_ok], lo[ci_ok], hi[ci_ok], 0))
  plot(NA, xlim = xr + c(-0.2, 0.2), ylim = c(0.5, n + 0.5),
       xlab = "Severity (logits)", ylab = "", axes = FALSE, main = "")
  title(main = facet, adj = 0, line = 1.4)
  abline(h = seq_len(n), col = .rr$grid, lwd = 0.8)
  abline(v = 0, lty = 2, col = .rr$soft)
  .rr_axis(1)
  axis(2, at = seq_len(n), labels = fe$level, cex.axis = 0.8,
       col = .rr$grid, col.ticks = NA)
  misfit <- !is.na(fe$fit_resid) & abs(fe$fit_resid) > band
  segments(lo[ci_ok], seq_len(n)[ci_ok], hi[ci_ok], seq_len(n)[ci_ok],
           lwd = 2.2, col = ifelse(misfit[ci_ok], .rr$red, .rr$soft))
  points(fe$severity[point_ok], seq_len(n)[point_ok], pch = 21, cex = 1.5,
         bg = ifelse(misfit[point_ok], .rr$red, .rr$blue),
         col = "white", lwd = 1.2)
  if (any(misfit))
    mtext(sprintf("%d level(s) with |fit residual| > %.1f", sum(misfit), band),
          side = 3, line = 0.2, adj = 0, cex = 0.8, col = .rr$red)
  invisible(NULL)
}
