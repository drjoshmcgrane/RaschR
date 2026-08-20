# rasch :: differential item functioning
# ===========================================================================
# DIF by analysis of variance of the standardised residuals (Hagquist &
# Andrich 2017). For each item, residuals are analysed by person factor and
# trait class interval: a factor main effect indicates uniform DIF and a
# factor-by-interval interaction indicates non-uniform DIF. With several
# person factors they are modelled jointly by dif_anova (main effects by
# default, factor-by-factor interactions optional), with Tukey HSD
# comparisons on the significant group terms and the convention that a
# significant interaction supersedes the main effects of the factors
# involved. Multiplicity across items is handled by Benjamini-Hochberg
# false-discovery-rate adjustment.
# ===========================================================================

.dif_factors <- function(fit, factors) {
  if (is.null(factors)) factors <- fit$factors
  if (is.null(factors)) stop("no person factors supplied or stored in the fit")
  if (is.character(factors) && !is.null(fit$factors) &&
      all(factors %in% names(fit$factors)))
    factors <- fit$factors[, factors, drop = FALSE]
  if (!is.data.frame(factors)) factors <- data.frame(group = factors)
  factors
}

# Class intervals for a DIF analysis are set from the cells the analysis
# actually uses: the residual ANOVA crosses trait intervals with group
# levels (or with the factor-combination cells in the factorial), so the
# interval count is chosen to keep the smallest group's expected cell size
# adequate -- independently of the interval count of the overall fit.
.dif_n_groups <- function(fit, grp, cell_min = 30L, id = NULL) {
  ok <- !is.na(grp) & !is.na(fit$person$theta)
  if (!any(ok)) return(2L)
  # with repeated ids the cells are counted in PERSONS, not rows: stacked
  # or duplicated observations must not widen the interval rule
  n_min <- if (!is.null(id))
    min(tapply(as.character(id)[ok], droplevels(factor(grp[ok])),
               function(v) length(unique(v))))
  else min(table(droplevels(factor(grp[ok]))))
  max(2L, min(10L, as.integer(n_min) %/% as.integer(cell_min)))
}

.dif_class_intervals <- function(fit, n_groups) {
  ci <- fit$person$class_interval
  if (is.null(ci) || !identical(n_groups, fit$n_groups))
    ci <- .class_intervals(fit$person$theta, fit$person$extreme, n_groups)
  factor(ci)
}

# Person-level class intervals for a within-subjects analysis: each person
# gets one interval from their mean location, so the interval is a clean
# whole-plot factor. Returned aligned to the rows of the fit.
.dif_person_ci <- function(fit, id, n_groups) {
  th <- fit$person$theta; ex <- fit$person$extreme
  pth <- tapply(th, id, mean, na.rm = TRUE)
  pex <- tapply(ex, id, function(v) all(v, na.rm = TRUE))
  pci <- .class_intervals(as.numeric(pth), as.logical(pex), n_groups)
  factor(pci[match(as.character(id), names(pth))])
}

# variables of an ANOVA term label, e.g. "g1:ci" -> c("g1", "ci")
.term_vars <- function(term) strsplit(term, ":", fixed = TRUE)[[1]]

# ---------------------------------------------------------------------------
# Order-invariant person-level DIF tests. Between-person terms get Type II
# sums of squares on person-level residual means: each term is adjusted for
# every term NOT containing it (so the class interval is always adjusted
# out of every group test), with F against the full model's between-person
# residual -- sequential (Type I) tests let a group factor absorb trait or
# correlated-factor variance in unbalanced designs, flipping which factor
# flags with entry order. Within-person terms are tested on person-by-cell
# means through orthonormal contrasts with the Greenhouse-Geisser epsilon
# correction (Maxwell & Delaney 2004): classical split-plot strata assume
# sphericity, and a nonspherical 4-level null rejected at ~9% nominal 5%.
# ---------------------------------------------------------------------------
.dif_type2 <- function(d, term_labels, resp = "z") {
  mk <- function(tl) stats::as.formula(paste(
    resp, "~", if (length(tl)) paste(tl, collapse = " + ") else "1"))
  full <- tryCatch(stats::lm(mk(term_labels), data = d),
                   error = function(e) NULL)
  if (is.null(full)) return(NULL)
  rss_full <- sum(stats::resid(full)^2)
  df_res <- stats::df.residual(full)
  if (df_res < 1 || rss_full <= 0) return(NULL)
  mse <- rss_full / df_res
  out <- list()
  for (tt in term_labels) {
    tv <- .term_vars(tt)
    not_cont <- term_labels[!vapply(term_labels, function(u)
      all(tv %in% .term_vars(u)), TRUE)]
    m0 <- stats::lm(mk(not_cont), data = d)
    m1 <- stats::lm(mk(c(not_cont, tt)), data = d)
    df_t <- stats::df.residual(m0) - stats::df.residual(m1)
    if (df_t < 1) next
    ss_t <- max(sum(stats::resid(m0)^2) - sum(stats::resid(m1)^2), 0)
    Fv <- (ss_t / df_t) / mse
    out[[length(out) + 1L]] <- data.frame(
      term = tt, df = df_t, df_denom = df_res, gg_epsilon = NA_real_,
      sum_sq = ss_t, mean_sq = ss_t / df_t,
      F_value = Fv, p = stats::pf(Fv, df_t, df_res, lower.tail = FALSE),
      resid_ss = rss_full, stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  rbind(do.call(rbind, out),
        data.frame(term = "Residuals", df = df_res, df_denom = NA_real_,
                   gg_epsilon = NA_real_, sum_sq = rss_full,
                   mean_sq = mse, F_value = NA_real_, p = NA_real_,
                   resid_ss = NA_real_, stringsAsFactors = FALSE))
}

# Within-stratum tests on the person-by-within-cell mean matrix Y (complete
# cases over cells). For a term pairing the within subspace (Kronecker
# contrast matrix over the within factors) with a between portion, the
# contrast scores S = Y C are tested by Type II model comparison of each
# score column on the between design, pooled over columns, with df scaled
# by the Greenhouse-Geisser epsilon of the residual score covariance.
.dif_within_tests <- function(Y, pdat, wname, wlv, within_terms,
                              bterms_all) {
  n <- nrow(Y)
  contr_of <- function(k) {
    C <- stats::contr.helmert(k)
    sweep(C, 2, sqrt(colSums(C^2)), "/")
  }
  meanvec_of <- function(k) matrix(1 / sqrt(k), k, 1)
  mk <- function(tl, resp) stats::as.formula(paste(
    resp, "~", if (length(tl)) paste(tl, collapse = " + ") else "1"))
  out <- list(); resid_pool <- 0; resid_df <- 0
  na_row <- function(tt) data.frame(
    term = tt, df = NA_real_, df_denom = NA_real_, gg_epsilon = NA_real_,
    sum_sq = NA_real_, mean_sq = NA_real_, F_value = NA_real_, p = NA_real_,
    resid_ss = NA_real_, stringsAsFactors = FALSE)
  # between design for the scores: only terms whose factors survive the
  # complete-panel filtering with at least two levels (a group observed at
  # a single occasion pattern can lose every complete panel; its
  # interactions are then non-estimable and are reported NA rather than
  # crashing lm with a one-level factor)
  pdat <- droplevels(pdat)
  ok_var <- vapply(names(pdat), function(cn)
    !is.factor(pdat[[cn]]) || nlevels(pdat[[cn]]) >= 2L, TRUE)
  bad_vars <- names(pdat)[!ok_var]
  bt_full <- bterms_all[!vapply(bterms_all, function(u)
    any(.term_vars(u) %in% bad_vars), TRUE)]
  for (tt in within_terms) {
    tv <- .term_vars(tt)
    w_t <- intersect(tv, names(wlv))
    b_t <- setdiff(tv, w_t)
    if (any(b_t %in% bad_vars)) {         # non-estimable after filtering
      out[[length(out) + 1L]] <- na_row(tt)
      next
    }
    Cm <- matrix(1, 1, 1)
    for (wf in names(wlv)) {
      k <- wlv[[wf]]
      Cm <- Cm %x% (if (wf %in% w_t) contr_of(k) else meanvec_of(k))
    }
    S <- Y %*% Cm
    m <- ncol(S)
    fits_j <- lapply(seq_len(m), function(j) {
      dd <- pdat; dd$s_ <- S[, j]
      stats::lm(mk(bt_full, "s_"), data = dd)
    })
    rss_f <- sum(vapply(fits_j, function(f) sum(stats::resid(f)^2), 0))
    dfr1 <- stats::df.residual(fits_j[[1]])
    df_err <- m * dfr1
    if (df_err < 1 || rss_f <= 0) next
    # Greenhouse-Geisser epsilon from the residual score covariance
    E <- vapply(fits_j, stats::resid, numeric(n))
    Sg <- crossprod(as.matrix(E)) / dfr1
    lam <- eigen(Sg, symmetric = TRUE, only.values = TRUE)$values
    lam <- pmax(lam, 0)
    eps <- if (m == 1L || sum(lam^2) <= 0) 1 else
      max(min(sum(lam)^2 / (m * sum(lam^2)), 1), 1 / m)
    if (!length(b_t)) {
      # the within main effect is the grand mean of the contrast scores,
      # adjusted for the between design. Removing the intercept from a
      # FORMULA does nothing when factors are present (R re-parameterises
      # them to absorb the constant), so the design is built explicitly
      # with sum-to-zero factor coding, where the intercept column is the
      # balanced grand mean and genuinely separable.
      Xb <- if (length(bt_full)) {
        fml <- stats::as.formula(paste("~", paste(bt_full, collapse = " + ")))
        used <- unique(unlist(lapply(bt_full, .term_vars)))
        fac_cols <- intersect(
          names(pdat)[vapply(pdat, is.factor, TRUE)], used)
        ctr <- stats::setNames(
          rep(list("contr.sum"), length(fac_cols)), fac_cols)
        stats::model.matrix(fml, data = pdat, contrasts.arg = ctr)
      } else matrix(1, n, 1)
      ss_t <- 0
      for (j in seq_len(m)) {
        r_full <- stats::lm.fit(Xb, S[, j])$residuals
        r_red <- stats::lm.fit(Xb[, -1, drop = FALSE], S[, j])$residuals
        ss_t <- ss_t + max(sum(r_red^2) - sum(r_full^2), 0)
      }
      df_t <- m
    } else {
      not_cont <- bt_full[!vapply(bt_full, function(u)
        all(b_t %in% .term_vars(u)), TRUE)]
      rss0 <- rss1 <- 0; df_t1 <- NA_integer_
      for (j in seq_len(m)) {
        dd <- pdat; dd$s_ <- S[, j]
        f0 <- stats::lm(mk(not_cont, "s_"), data = dd)
        f1 <- stats::lm(mk(c(not_cont, paste(b_t, collapse = ":")), "s_"),
                        data = dd)
        rss0 <- rss0 + sum(stats::resid(f0)^2)
        rss1 <- rss1 + sum(stats::resid(f1)^2)
        df_t1 <- stats::df.residual(f0) - stats::df.residual(f1)
      }
      if (is.na(df_t1) || df_t1 < 1) next
      ss_t <- max(rss0 - rss1, 0)
      df_t <- m * df_t1
    }
    Fv <- (ss_t / df_t) / (rss_f / df_err)
    # the p-value is computed at the Greenhouse-Geisser corrected degrees
    # of freedom (eps * df, eps * df_denom); the nominal df, the
    # denominator df, and epsilon are all returned so the test is
    # reproducible and reportable
    out[[length(out) + 1L]] <- data.frame(
      term = tt, df = df_t, df_denom = df_err, gg_epsilon = eps,
      sum_sq = ss_t, mean_sq = ss_t / df_t,
      F_value = Fv,
      p = stats::pf(Fv, eps * df_t, eps * df_err, lower.tail = FALSE),
      resid_ss = rss_f, stringsAsFactors = FALSE)
    resid_pool <- rss_f; resid_df <- df_err
  }
  if (!length(out)) return(NULL)
  rbind(do.call(rbind, out),
        data.frame(term = "Residuals", df = resid_df, df_denom = NA_real_,
                   gg_epsilon = NA_real_, sum_sq = resid_pool,
                   mean_sq = if (resid_df > 0) resid_pool / resid_df else
                     NA_real_, F_value = NA_real_, p = NA_real_,
                   resid_ss = NA_real_, stringsAsFactors = FALSE))
}

# Flatten an aov (single- or multi-stratum) into one row per term, carrying
# the residual sum of squares of the term's own stratum so a partial
# eta-squared can be formed against the right error in a mixed design.
.aov_terms_flat <- function(a) {
  sm <- summary(a)
  strata <- if (inherits(a, "aovlist")) sm else list(sm)
  rows <- list()
  for (st in strata) {
    tab <- st[[1L]]; rn <- trimws(rownames(tab))
    nm <- names(tab)
    rss <- if ("Residuals" %in% rn) tab[["Sum Sq"]][rn == "Residuals"] else NA_real_
    for (k in seq_len(nrow(tab)))
      rows[[length(rows) + 1L]] <- data.frame(
        term = rn[k], df = tab[["Df"]][k], sum_sq = tab[["Sum Sq"]][k],
        mean_sq = tab[["Mean Sq"]][k],
        F_value = if ("F value" %in% nm) tab[["F value"]][k] else NA_real_,
        p = if ("Pr(>F)" %in% nm) tab[["Pr(>F)"]][k] else NA_real_,
        resid_ss = rss, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  # missing responses unbalance a mixed design, and aov then projects a
  # within term onto an earlier (between-person) error stratum as well as
  # its own; only the deepest stratum tests the term against its proper
  # error, so keep each term's last occurrence (Residuals rows stay
  # per-stratum for the partial eta-squared)
  dup <- out$term != "Residuals" & duplicated(out$term, fromLast = TRUE)
  out[!dup, , drop = FALSE]
}

#' Differential item functioning by residual analysis of variance
#'
#' Tests uniform and non-uniform DIF by analysing each item's standardised
#' residuals over person factors and trait class intervals (Andrich and Marais
#' 2019, ch. 16). Several person factors are fitted jointly. The function also
#' supports designs containing both between-person and within-person factors.
#'
#' @details
#' With one factor \eqn{G} and class interval \eqn{C}, the residual model is
#' \deqn{z=\mu+G+C+G\mathbin{:}C+\varepsilon.}
#' The factor term tests uniform DIF and its interaction with class interval
#' tests non-uniform DIF. With several factors, \code{effects = "main"} fits
#' \code{(f1 + f2 + ...) * ci}; \code{effects = "factorial"} also includes
#' factor-by-factor interactions. Type II sums of squares are used, and
#' probabilities are adjusted across items separately within each term.
#'
#' When identifiers repeat, the person is the unit of analysis. Between-person
#' terms use person means and the between-person error stratum. Within-person
#' terms use orthonormal contrasts of person-by-cell means. A
#' Greenhouse--Geisser correction is applied to within-person factors with
#' more than two levels. Persons missing a required cell are excluded from the
#' corresponding within-person test. In incomplete mixed designs, within-cell
#' effects are removed before the between-person analysis.
#'
#' A significant higher-order factor term supersedes its component terms in
#' the summary. For EFRM fits, frame-defining factors are excluded because
#' they define the model rather than a separate DIF contrast; testing such a
#' factor means stepping outside the model, which is what
#' \code{\link{frame_invariance}} does. MFRM residuals
#' are pooled to underlying items unless \code{pool_facets = FALSE}.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param factors A vector (one factor), a data frame of person factors, or a
#'   character vector naming factor columns nominated in the fit. Defaults to
#'   every factor stored in the fit.
#' @param n_groups Number of trait class intervals. The default uses the
#'   smallest joint factor cell to retain about 30 expected responses per
#'   interval and cell, with between 2 and 10 intervals. The selected value is
#'   returned in \code{n_groups}.
#' @param p_adjust Multiplicity adjustment across items within each term;
#'   default \code{"BH"}.
#' @param alpha Significance level applied to the adjusted probabilities.
#' @param effects \code{"main"} (default) models several factors additively
#'   (each factor's main effect and its class-interval interaction, but no
#'   factor-by-factor terms); \code{"factorial"} also crosses the factors
#'   with each other. Immaterial with a single factor.
#' @param id Person identifier for stacked or repeated-measures data. It may
#'   be a column name stored in the fit or a vector with one value per row;
#'   by default the identifier carried by the fit is used.
#' @param within Names of within-person factors. With repeated identifiers,
#'   varying factors are detected automatically when this is omitted. See
#'   Details for the mixed-design analysis.
#' @param pool_facets For MFRM fits: pool residuals to the underlying
#'   items (the default), so DIF is tested per item rather than per
#'   item-by-facet cell; \code{FALSE} tests each cell as its own item.
#'   Ignored for other fits.
#' @param sizes If \code{TRUE}, refit each flagged item-term and calculate
#'   pairwise DIF differences in logits using \code{\link{dif_size}}.
#' @return A list with:
#' \describe{
#'   \item{\code{summary}}{One row per item and group term, containing the
#'   uniform and non-uniform tests, partial eta-squared, adjusted
#'   probabilities, DIF flags, and supersession flag.}
#'   \item{\code{terms}}{The complete item-wise analysis-of-variance tables.}
#'   \item{\code{tukey}}{Tukey comparisons for significant,
#'   non-superseded terms with more than two levels.}
#'   \item{\code{sizes}}{When requested, pairwise logit differences for the
#'   significant, non-superseded item-terms.}
#'   \item{\code{posthoc}}{When \code{sizes = TRUE}, marginal pairwise
#'   differences for main effects and difference-in-differences magnitudes
#'   for interactions, calculated by \code{\link{dif_posthoc}}.}
#' }
#' The remaining components record the factors, class intervals, adjustment,
#' significance level, and design settings.
#' @references
#' Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing. Journal of
#' the Royal Statistical Society: Series B, 57(1), 289--300.
#'
#' Hagquist, C. and Andrich, D. (2017). Recent advances in analysis of
#' differential item functioning in health research using the Rasch model.
#' Health and Quality of Life Outcomes, 15, 181.
#'
#' Maxwell, S. E. and Delaney, H. D. (2004). Designing Experiments and
#' Analyzing Data: A Model Comparison Perspective (2nd ed.). Lawrence Erlbaum.
#' @seealso \code{\link{dif_size}}, \code{\link{dif_contrasts}}, and
#'   \code{\link{resolve_dif}}; and \code{\link{frame_invariance}} for the
#'   frame-defining factor this function excludes.
#' @examples
#' set.seed(1); n <- 800
#' d <- seq(-1.5, 1.5, length.out = 6)
#' g1 <- rep(c("a", "b"), each = n / 2)
#' g2 <- rep(c("x", "y"), times = n / 2)
#' sh <- matrix(0, n, 6); sh[g1 == "b", 2] <- 0.8
#' X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
#' colnames(X) <- paste0("I", 1:6)
#' fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))
#' dif_anova(fit)$summary
#'
#' \donttest{
#' # Mixed design: group is between persons and occasion is within persons.
#' N <- 320; theta <- rnorm(N); group <- rep(c("A", "B"), each = N / 2)
#' make_wave <- function(occasion_shift) {
#'   shift <- matrix(0, N, 6)
#'   shift[group == "B", 2] <- 0.9
#'   shift[, 5] <- occasion_shift
#'   matrix(rbinom(N * 6, 1,
#'          plogis(outer(theta, d, "-") - shift)), N, 6)
#' }
#' Xm <- rbind(make_wave(0), make_wave(1.0))
#' colnames(Xm) <- paste0("I", 1:6)
#' repeated <- data.frame(Xm, group = rep(group, 2),
#'                        occasion = rep(c("T1", "T2"), each = N))
#' mixed_fit <- rasch(repeated, id = rep(seq_len(N), 2),
#'                    factors = c("group", "occasion"))
#' mixed_dif <- dif_anova(mixed_fit, within = "occasion")
#' subset(mixed_dif$summary, uniform_DIF | nonuniform_DIF)
#' }
#' @export
dif_anova <- function(fit, factors = NULL, n_groups = NULL,
                                p_adjust = "BH", alpha = 0.05,
                                effects = c("main", "factorial"),
                                sizes = FALSE, id = NULL, within = NULL,
                                pool_facets = TRUE) {
  effects <- match.arg(effects)
  Z <- fit$residuals; L <- ncol(Z)
  # MFRM residuals pool to the UNDERLYING items by default: users ask
  # whether item A shows DIF, not whether the virtual cell A:R2 does; the
  # per-virtual-item tests remain available with pool_facets = FALSE
  pooled_note <- NULL
  if (inherits(fit, "rasch_mfrm") && isTRUE(pool_facets) &&
      !is.null(fit$virtual_map)) {
    vm <- fit$virtual_map
    items_u <- unique(vm$item)
    Zp <- vapply(items_u, function(it) {
      zz <- Z[, vm$vkey[vm$item == it], drop = FALSE]
      nn <- rowSums(is.finite(zz))
      out <- rowSums(zz, na.rm = TRUE) / sqrt(pmax(nn, 1L))
      out[nn == 0L] <- NA_real_
      out
    }, numeric(nrow(Z)))
    Zp[!is.finite(Zp)] <- NA_real_
    colnames(Zp) <- items_u
    Z <- Zp; L <- ncol(Z)
    pooled_note <- paste(
      "MFRM residuals pooled to the underlying items (standardised sum over",
      "each item's observed facet cells, so rows with different facet",
      "coverage retain comparable null variance); pool_facets = FALSE tests",
      "each item-by-facet cell as its own item")
  }
  factors <- .dif_factors(fit, factors)
  # the EFRM frame group IS the frame structure: each frame has its own
  # virtual items, so the group factor has a single level among any
  # virtual item's responders and cannot be tested as DIF
  drop_frame_note <- NULL
  if (!is.null(fit$frame_group) && any(names(factors) %in% fit$frame_group)) {
    hit <- intersect(names(factors), fit$frame_group)
    if (length(hit) == length(factors))
      .refuse("'", paste(hit, collapse = "', '"),
           "' define(s) the EFRM frame structure itself: the factor is ",
           "constant among the persons responding within any frame, so ",
           "no within-frame comparison remains to test as DIF; nominate ",
           "other person factors")
    factors <- factors[!names(factors) %in% fit$frame_group]
    drop_frame_note <- paste0("frame factor(s) '",
                              paste(hit, collapse = "', '"),
                              "' excluded: they are the frame structure, ",
                              "not testable DIF factors")
  }

  # within-subject factors (levels repeating within a person) turn the
  # analysis into a mixed (split-plot) one: the class interval is taken at
  # the person level so it is a clean whole-plot factor, and the within
  # factors carry a person error stratum.
  if (is.character(id) && length(id) == 1L) {
    if (is.null(fit$factors) || !id %in% names(fit$factors))
      stop("id column '", id, "' not found among the fit's factors")
    id <- fit$factors[[id]]
  }
  if (is.null(id) && !is.null(fit$person$id)) id <- as.character(fit$person$id)
  if (!is.null(id)) id <- as.character(id)
  repeated <- !is.null(id) && anyDuplicated(id) > 0L
  if (!is.null(within)) {
    unknown <- setdiff(within, names(factors))
    if (length(unknown))
      stop("within-subject factor(s) not among the nominated factors: ",
           paste(unknown, collapse = ", "))
    if (!repeated)
      stop("within-subject factors need repeated person ids (each id ",
           "observed more than once); no id repeats here")
    varies <- vapply(within, function(fn)
      any(tapply(as.character(factors[[fn]]), id, function(v)
        length(unique(v[!is.na(v)])) > 1L), na.rm = TRUE), TRUE)
    if (any(!varies))
      stop("factor(s) declared within-subject never vary within any id: ",
           paste(within[!varies], collapse = ", "))
    # the converse is equally ill-defined: a factor that varies within
    # persons has no person-level value, so it cannot be treated as
    # between-subjects (the old row-level treatment pseudo-replicated)
    other <- setdiff(names(factors), within)
    ovaries <- vapply(other, function(fn)
      any(tapply(as.character(factors[[fn]]), id, function(v)
        length(unique(v[!is.na(v)])) > 1L), na.rm = TRUE), TRUE)
    if (length(other) && any(ovaries))
      stop("factor(s) vary within persons but are not declared in ",
           "`within`: ", paste(other[ovaries], collapse = ", "),
           "; declare them within-subject (or aggregate the data to one ",
           "row per person) -- treating repeated observations as ",
           "independent between-person rows manufactures information")
  }
  if (is.null(within) && repeated) {
    within <- names(factors)[vapply(names(factors), function(fn)
      any(tapply(as.character(factors[[fn]]), id, function(v)
        length(unique(v[!is.na(v)])) > 1L), na.rm = TRUE), TRUE)]
  }
  if (is.null(within)) within <- character(0)
  within <- intersect(within, names(factors))
  mixed <- length(within) > 0L

  if (is.null(n_groups)) {
    cells <- interaction(factors, drop = TRUE)
    n_groups <- .dif_n_groups(fit, cells,
                              id = if (repeated) id else NULL)
  }
  # PERSONS are the units whenever ids repeat (stacked or duplicated rows):
  # observation-level tests would let copied observations manufacture
  # information. The class interval is taken at the person level so it is a
  # clean whole-plot covariate.
  ci <- if (repeated) .dif_person_ci(fit, id, n_groups) else
    .dif_class_intervals(fit, n_groups)

  fnames <- names(factors)
  safe <- paste0("f", seq_along(fnames))           # syntactic stand-ins
  wsafe <- safe[match(within, fnames)]
  bsafe <- setdiff(safe, wsafe)
  op <- if (effects == "factorial") " * " else " + "
  form_all <- stats::as.formula(
    paste("z ~ (", paste(safe, collapse = op), ") * ci"))
  all_terms <- attr(stats::terms(form_all), "term.labels")
  bterms <- all_terms[!vapply(all_terms, function(tt)
    any(.term_vars(tt) %in% wsafe), TRUE)]
  wterms <- setdiff(all_terms, bterms)

  fits <- vector("list", L); rows <- list()
  incomplete_note <- 0L
  for (i in seq_len(L)) {
    d <- data.frame(z = Z[, i], ci = ci)
    d$pid <- if (is.null(id)) sprintf("p%06d", seq_len(nrow(d))) else id
    for (j in seq_along(fnames)) d[[safe[j]]] <- factor(factors[[fnames[j]]])
    d <- d[stats::complete.cases(d), ]
    if (nrow(d) < 10 || any(vapply(safe, function(s)
      length(unique(d[[s]])) < 2, TRUE))) next

    # aggregate to one mean residual per person per within-cell (persons
    # without repeats aggregate to themselves). Within cells are ordered
    # with the LAST within factor varying fastest, matching the Kronecker
    # construction of the contrast matrices -- interaction()'s
    # first-fastest order silently rotated multi-within effects into the
    # wrong subspaces (a pure w1 effect flagged w1:w2).
    if (mixed) {
      wl0 <- lapply(wsafe, function(sn)
        sort(unique(as.character(d[[sn]]))))
      names(wl0) <- wsafe
      grid <- expand.grid(rev(wl0), stringsAsFactors = FALSE,
                          KEEP.OUT.ATTRS = FALSE)
      grid <- grid[, rev(seq_along(wl0)), drop = FALSE]
      names(grid) <- wsafe
      cell_levels <- do.call(paste, c(grid, list(sep = "\r")))
      wcell <- factor(do.call(paste, c(lapply(wsafe, function(sn)
        as.character(d[[sn]])), list(sep = "\r"))), levels = cell_levels)
    } else wcell <- factor(rep("all", nrow(d)))
    key <- interaction(d$pid, wcell, drop = TRUE)
    agz <- tapply(d$z, key, mean)
    firsts <- which(!duplicated(key))
    ag <- d[firsts[match(levels(key), as.character(key[firsts]))],
            c("pid", "ci", safe), drop = FALSE]
    ag$z <- as.numeric(agz)
    ag$wcell <- wcell[firsts[match(levels(key),
                                   as.character(key[firsts]))]]

    # person-level frame for the between stratum: one row per person, the
    # mean over that person's OCCASION-ADJUSTED within-cell values. With
    # differentially incomplete within panels, raw person means are not
    # comparable between groups (a common occasion effect plus one group
    # missing an occasion masqueraded as group DIF at F = 37.6); centring
    # each within cell at its all-person mean removes the common within
    # effects from the between comparison.
    # centring must remove within effects that VARY BY TRAIT LEVEL too: a
    # common occasion-by-class-interval structure plus differential
    # missingness otherwise masquerades as non-uniform group DIF (observed
    # F = 214.7 on grp:ci with no group effect). Centre each within cell
    # within each class interval; empty combinations fall back to the
    # cell's overall mean.
    cellci <- interaction(ag$wcell, ag$ci, drop = TRUE)
    m_cellci <- tapply(ag$z, cellci, mean)
    m_cell <- tapply(ag$z, ag$wcell, mean)
    ctr <- as.numeric(m_cellci[as.character(cellci)])
    miss_ctr <- is.na(ctr)
    if (any(miss_ctr))
      ctr[miss_ctr] <- as.numeric(m_cell[as.character(ag$wcell)[miss_ctr]])
    zc <- ag$z - ctr
    pkey <- factor(ag$pid)
    pz <- tapply(zc, pkey, mean)
    pfirst <- which(!duplicated(pkey))
    pdat <- ag[pfirst[match(levels(pkey), as.character(pkey[pfirst]))],
               c("pid", "ci", bsafe), drop = FALSE]
    pdat$z <- as.numeric(pz)

    ft_b <- .dif_type2(pdat, bterms)
    ft_w <- NULL
    if (mixed && length(wterms)) {
      # complete within-cell matrix per person; incomplete persons are
      # dropped explicitly (multi-stratum projections on unbalanced data
      # are exactly the murky territory this engine replaces)
      wl <- lapply(wsafe, function(sn) sort(unique(as.character(ag[[sn]]))))
      names(wl) <- wsafe
      Ywide <- tapply(ag$z, list(factor(ag$pid), ag$wcell), mean)
      complete <- rowSums(is.na(Ywide)) == 0L
      incomplete_note <- incomplete_note + sum(!complete)
      if (sum(complete) >= 6L) {
        Y <- Ywide[complete, , drop = FALSE]
        pd2 <- pdat[match(rownames(Y), as.character(pdat$pid)), ,
                    drop = FALSE]
        wlv <- lapply(wl, length)
        ft_w <- .dif_within_tests(Y, pd2, paste(wsafe, collapse = ":"),
                                  wlv, wterms, bterms)
      }
    }
    ft <- rbind(ft_b, ft_w)
    if (is.null(ft)) next
    # a person-level aov (class interval first, BETWEEN factors only, one
    # row per person) retained solely for the Tukey HSD follow-ups: cell
    # means and the between-person error MS are what Tukey needs. Within
    # factors are excluded -- ordinary Tukey on repeated cells would treat
    # them as independent, so within-term follow-ups are not offered.
    fits[i] <- list(if (length(bsafe)) tryCatch(stats::aov(
      stats::as.formula(paste("z ~ ci + (",
                              paste(bsafe, collapse = op), ")")),
      data = pdat), error = function(e) NULL) else NULL)
    rows[[length(rows) + 1L]] <- data.frame(item = colnames(Z)[i], ft)
  }
  if (!length(rows)) stop("no item yielded an estimable factorial ANOVA")
  terms <- do.call(rbind, rows)
  rownames(terms) <- NULL

  # partial eta-squared per term against its own stratum residual
  terms$eta2_partial <- ifelse(
    terms$term == "Residuals" | is.na(terms$resid_ss), NA_real_,
    terms$sum_sq / (terms$sum_sq + terms$resid_ss))
  terms$resid_ss <- NULL

  # adjust across items within each term (the residual rows carry no test)
  terms$p_adj <- NA_real_
  for (tt in setdiff(unique(terms$term), "Residuals")) {
    sel <- terms$term == tt
    terms$p_adj[sel] <- p.adjust(terms$p[sel], method = p_adjust)
  }
  terms$significant <- !is.na(terms$p_adj) & terms$p_adj < alpha

  # a significant higher-order GROUP interaction supersedes the lower-order
  # group terms built from a subset of its factors (within the same item).
  # Terms crossing the class interval are excluded from the pass: a term's
  # own ci interaction is reported WITH it as non-uniform DIF, so it must
  # not supersede it (an item significant on both would otherwise drop out
  # of the follow-ups entirely).
  terms$superseded <- FALSE
  is_group_t <- !vapply(terms$term, function(t)
    "ci" %in% .term_vars(t), TRUE)
  for (it in unique(terms$item)) {
    sel <- which(terms$item == it & terms$significant & is_group_t)
    if (length(sel) < 2) next
    vlist <- lapply(terms$term[sel], .term_vars)
    for (a_i in seq_along(sel)) for (b_i in seq_along(sel)) {
      if (a_i == b_i) next
      if (all(vlist[[a_i]] %in% vlist[[b_i]]) &&
          length(vlist[[a_i]]) < length(vlist[[b_i]]))
        terms$superseded[sel[a_i]] <- TRUE
    }
  }

  # Tukey HSD for significant, non-superseded terms that do not involve the
  # class interval (the group structure itself)
  tk <- list(); tukey_note <- NULL
  for (i in seq_len(L)) {
    a <- fits[[i]]; if (is.null(a)) next
    it <- colnames(Z)[i]
    cand <- terms[terms$item == it & terms$significant & !terms$superseded, ]
    # group terms only; and no comparisons for a two-level main effect,
    # where the F test is already the only contrast
    keep_t <- !vapply(cand$term, function(tt) "ci" %in% .term_vars(tt), TRUE) &
      !(cand$df == 1L & !grepl(":", cand$term, fixed = TRUE)) &
      # within-subject terms have no place in an ordinary Tukey HSD (the
      # follow-up aov is person-level and between-factors only)
      !vapply(cand$term, function(tt)
        any(.term_vars(tt) %in% wsafe), TRUE)
    cand <- cand$term[keep_t]
    if (!length(cand)) next
    # TukeyHSD has no method for the multi-stratum aov of a mixed design;
    # say so rather than return an empty table silently
    if (inherits(a, "aovlist")) {
      tukey_note <- paste("Tukey comparisons are unavailable for",
                          "within-subject (mixed) designs; use the",
                          "resolved DIF magnitudes (sizes) instead")
      next
    }
    th <- tryCatch(stats::TukeyHSD(a, which = cand), error = function(e) NULL)
    if (is.null(th)) next
    for (tt in names(th)) {
      tb <- as.data.frame(th[[tt]])
      tk[[length(tk) + 1L]] <- data.frame(
        item = it, term = tt, comparison = rownames(tb),
        difference = tb$diff, lower = tb$lwr, upper = tb$upr,
        p_tukey = tb$`p adj`, row.names = NULL)
    }
  }
  tukey <- if (length(tk)) do.call(rbind, tk) else
    data.frame(item = character(), term = character(),
               comparison = character(), difference = numeric(),
               lower = numeric(), upper = numeric(), p_tukey = numeric())

  # map a term's syntactic stand-ins (f1..fk) back to the nominated factor
  # names by exact whole-token match, so a factor named like a stand-in
  # ("f1") or like the class interval ("ci") cannot be re-substituted or
  # collide. Applied only for display, after all classification is done on
  # the stand-ins.
  relabel <- function(x) vapply(x, function(t) {
    toks <- strsplit(t, ":", fixed = TRUE)[[1]]
    i <- match(toks, safe); toks[!is.na(i)] <- fnames[i[!is.na(i)]]
    if ("ci" %in% fnames)
      toks[is.na(i) & toks == "ci"] <- "(class interval)"
    paste(toks, collapse = ":")
  }, character(1), USE.NAMES = FALSE)

  # DIF magnitudes in logits for the significant, non-superseded group
  # terms (interaction terms resolved by their cells)
  size_tab <- posthoc_tab <- NULL
  if (isTRUE(sizes)) {
    sz <- ph <- list(); size_fail <- posthoc_fail <- character(0)
    cand <- terms[terms$significant & !terms$superseded &
                  !vapply(terms$term, function(tt)
                    "ci" %in% .term_vars(tt), TRUE), , drop = FALSE]
    for (r in seq_len(nrow(cand))) {
      it <- cand$item[r]; tt <- cand$term[r]
      # dif_size takes the nominated factor names, not the stand-ins
      by_user <- fnames[match(.term_vars(tt), safe)]
      ds <- tryCatch(dif_size(fit, it, by = by_user),
                     error = function(e) e)
      if (inherits(ds, "error")) {
        size_fail <- c(size_fail, conditionMessage(ds))
        next
      }
      p <- ds$pairs
      sz[[length(sz) + 1L]] <- data.frame(item = it, term = tt, p,
                                          row.names = NULL)
      dp <- tryCatch(dif_posthoc(
        fit, it, term = paste(by_user, collapse = ":"), factors = fnames,
        within = fnames[match(wsafe, safe)], id = id),
        error = function(e) e)
      if (inherits(dp, "error"))
        posthoc_fail <- c(posthoc_fail, conditionMessage(dp))
      else
        ph[[length(ph) + 1L]] <- data.frame(
          item = it, term = tt, dp$table, row.names = NULL)
    }
    size_tab <- if (length(sz)) do.call(rbind, sz) else
      data.frame(item = character(), term = character())
    size_note <- if (length(size_fail)) paste0(
      "DIF magnitudes unavailable for some flagged term(s): ",
      paste(unique(size_fail), collapse = "; ")) else NULL
    posthoc_tab <- if (length(ph)) do.call(rbind, ph) else
      data.frame(item = character(), term = character())
    posthoc_note <- if (length(posthoc_fail)) paste0(
      "DIF post-hoc comparisons unavailable for some flagged term(s): ",
      paste(unique(posthoc_fail), collapse = "; ")) else NULL
  } else size_note <- posthoc_note <- NULL

  # compact reading: one row per item and group term, its own effect being
  # uniform DIF and its crossing with the class interval non-uniform DIF
  gterms <- setdiff(unique(terms$term), "Residuals")
  gterms <- gterms[!vapply(gterms, function(tt)
    "ci" %in% .term_vars(tt), TRUE)]
  srows <- list()
  for (it in unique(terms$item)) for (tt in gterms) {
    u <- terms[terms$item == it & terms$term == tt, , drop = FALSE]
    nu <- terms[terms$item == it & terms$term == paste0(tt, ":ci"), ,
                drop = FALSE]
    if (!nrow(u)) next
    # .aov_terms_flat deduplicates strata-split terms; a term reaching here
    # twice would silently break every isTRUE() below, so fail loudly
    if (nrow(u) > 1 || nrow(nu) > 1)
      stop("internal error: term '", tt, "' appears in several error strata")
    srows[[length(srows) + 1L]] <- data.frame(
      item = it, term = tt,
      F_uniform = u$F_value, p_uniform = u$p,
      p_uniform_adj = u$p_adj, eta2_uniform = u$eta2_partial,
      uniform_DIF = isTRUE(u$significant),
      F_nonuniform = if (nrow(nu)) nu$F_value else NA_real_,
      p_nonuniform = if (nrow(nu)) nu$p else NA_real_,
      p_nonuniform_adj = if (nrow(nu)) nu$p_adj else NA_real_,
      eta2_nonuniform = if (nrow(nu)) nu$eta2_partial else NA_real_,
      nonuniform_DIF = nrow(nu) > 0 && isTRUE(nu$significant),
      superseded = isTRUE(u$superseded))
  }
  summary_tab <- do.call(rbind, srows)
  rownames(summary_tab) <- NULL

  # relabel the stand-ins to the nominated names for display, now that all
  # classification is done
  terms$term <- relabel(terms$term)
  tukey$term <- relabel(tukey$term)
  summary_tab$term <- relabel(summary_tab$term)
  if (!is.null(size_tab) && nrow(size_tab))
    size_tab$term <- relabel(size_tab$term)
  if (!is.null(posthoc_tab) && nrow(posthoc_tab))
    posthoc_tab$term <- relabel(posthoc_tab$term)

  notes <- character(0)
  if (!is.null(pooled_note)) notes <- c(notes, pooled_note)
  if (!is.null(drop_frame_note)) notes <- c(notes, drop_frame_note)
  if (!is.null(size_note)) notes <- c(notes, size_note)
  if (!is.null(posthoc_note)) notes <- c(notes, posthoc_note)
  if (incomplete_note > 0L)
    notes <- c(notes, sprintf(
      "%d person-by-item panel(s) missing a within-subject cell were dropped from the within-person tests (their between-person information is retained)",
      incomplete_note))
  if (mixed && any(is.na(terms$F_value) & terms$term != "Residuals"))
    notes <- c(notes, paste(
      "term(s) reported NA were non-estimable after complete-panel",
      "filtering (a between level lost all its complete within panels)"))
  out <- list(summary = summary_tab, terms = terms, tukey = tukey,
              n_groups = nlevels(as.factor(ci)), within = within,
              factor_names = fnames,
              effects = effects, alpha = alpha, p_adjust = p_adjust,
              notes = notes)
  if (!is.null(tukey_note)) out$tukey_note <- tukey_note
  if (isTRUE(sizes)) {
    out$sizes <- size_tab
    out$posthoc <- posthoc_tab
  }
  out <- .tag_tables(out)
  class(out) <- "rasch_dif"
  out
}

#' @export
print.rasch_dif <- function(x, ...) {
  s <- x$summary
  cat(sprintf("DIF by residual analysis of variance (%s; %d class intervals%s)\n",
              if (length(unique(s$term)) > length(unique(s$item)) ||
                  x$effects == "factorial")
                sprintf("%d terms, %s effects", length(unique(s$term)),
                        x$effects) else "one-way",
              x$n_groups[1],
              if (length(x$within))
                sprintf("; within-subject: %s", paste(x$within, collapse = ", "))
              else ""))
  show <- s[, c("item", "term", "F_uniform", "p_uniform_adj", "uniform_DIF",
                "F_nonuniform", "p_nonuniform_adj", "nonuniform_DIF")]
  print(.fmt_df(show), row.names = FALSE)
  cat(sprintf("%d uniform, %d non-uniform DIF flag(s) after %s adjustment.\n",
              sum(s$uniform_DIF, na.rm = TRUE),
              sum(s$nonuniform_DIF, na.rm = TRUE), x$p_adjust))
  invisible(x)
}

#' DIF differences between factor levels
#'
#' Resolves an item into one parameter per factor level, refits the model, and
#' reports pairwise differences between the resolved locations. Wald tests use
#' the full sandwich covariance and are adjusted over the family of level
#' comparisons. Several names in \code{by} define factor-combination cells for
#' following up an interaction in \code{\link{dif_anova}}.
#'
#' @details
#' Let \eqn{\delta_i} contain the resolved locations of item
#' \eqn{i}, and let \eqn{\mathbf{c}_{ab}} place 1 on level \eqn{a}, -1 on
#' level \eqn{b}, and zero elsewhere. The reported difference and its
#' standard error are
#' \deqn{\Delta_{i,ab}=\mathbf{c}_{ab}^{\mathsf T}
#' \delta_i,}
#' \deqn{\operatorname{SE}(\Delta_{i,ab})=
#' \sqrt{\mathbf{c}_{ab}^{\mathsf T}\mathbf{V}_i
#' \mathbf{c}_{ab}},}
#' where \eqn{\mathbf{V}_i} is the full covariance of the resolved
#' locations. Wald probabilities are adjusted over the pairwise family.
#'
#' With repeated person identifiers, the row-level calibration covariance does
#' not represent within-person sampling dependence. Logit differences and
#' practical flags are retained, but their standard errors and Wald tests are
#' withheld. Use \code{\link{dif_contrasts}} for person-level inference in a
#' repeated-measures design.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or index.
#' @param by One or more person-factor names nominated in the fit (several
#'   names give interaction cells), or a grouping vector/data frame with
#'   one entry per person.
#' @param p_adjust Familywise adjustment over the pairwise comparisons;
#'   default \code{"holm"}.
#' @param alpha Significance level for the adjusted probabilities.
#' @param flag_logits Absolute difference flagged as practically
#'   significant.
#' @param min_n Levels with fewer responders to the item are dropped (their
#'   resolved locations would be too unstable to compare), with a note.
#' @return A list of class \code{"rasch_dif_size"}: \code{levels} (resolved
#'   location and SE per level, with its n), \code{pairs} (per comparison:
#'   difference in logits, SE, z, raw and adjusted p, 95 per cent interval,
#'   \code{significant}, \code{practical}, \code{ets}), the settings, and
#'   any notes. Sampling-uncertainty fields are \code{NA} when person IDs
#'   repeat.
#'
#' @section The ETS categories:
#' \code{pairs$ets} reports the A, B and C categories used at ETS, signed
#' for direction as they are there. They are defined from the
#' Mantel-Haenszel common odds ratio on the delta scale, where
#' \eqn{\mathrm{MH\ D\text{-}DIF} = -2.35\log\hat\alpha_{MH}}. Under the
#' Rasch model that log odds ratio and the difference in item location
#' estimate the same quantity, both being conditional on the total score, so
#' the delta thresholds convert exactly at 2.35 delta units to the logit: A
#' is not significant or below 0.426 logits, C is at or above 0.638 logits
#' and significantly beyond 0.426, and B is the remainder. The C rule tests
#' against a non-zero null, so it uses the standard error rather than the
#' probability alone.
#'
#' Read the letter beside the magnitude rather than instead of it, because
#' the categories were built to triage items for an operational bank and not
#' to answer whether an item is invariant. A tops out at 0.426 logits, which
#' is a difference of 10.6 percentage points in success at the item's own
#' location -- plainly not invariance. What justifies calling it negligible
#' is its effect on the score rather than on the item: in simulation at
#' 3,000 persons, one item sitting at that ceiling moved the comparison
#' between the two groups by 0.026 logits on a ten-item test and 0.018 on a
#' twenty-item one, while three such items moved it by 0.086 and 0.050. So
#' the letter answers "does this item distort the total score", and
#' \code{practical} against \code{flag_logits} answers "is this item
#' behaving the same way in both groups". They are different questions.
#'
#' The column is \code{NA} for a polytomous item. ETS classifies those from
#' a standardised mean difference in the observed-score metric, which is a
#' different statistic rather than a rescaling of this one.
#' @references
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' Scandinavian Journal of Statistics, 6(2), 65--70.
#' @seealso \code{\link{dif_anova}} and \code{\link{dif_contrasts}}.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' dif_size(fit, "I3", by = "grp")
#' @export
dif_size <- function(fit, item, by, p_adjust = "holm", alpha = 0.05,
                     flag_logits = 0.5, min_n = 20) {
  mfrm_item <- inherits(fit, "rasch_mfrm") && !is.null(fit$virtual_map) &&
    !(item %in% colnames(fit$X)) && item %in% fit$virtual_map$item
  if (!mfrm_item) {
    i <- .item_idx(fit, item)
    if (is.na(i)) stop("no such item")
    item <- fit$items$item[i]
  }
  if (is.character(by) && length(by) < nrow(fit$X)) {
    bad <- if (is.null(fit$factors)) by else setdiff(by, names(fit$factors))
    if (length(bad))
      stop("not a person factor nominated in the fit: ",
           paste(bad, collapse = ", "))
  }
  factors <- .dif_factors(fit, by)
  grp <- if (ncol(factors) == 1L) factor(factors[[1]])
         else interaction(factors, sep = ":", drop = TRUE)
  notes <- character(0)
  # Repeated persons couple the resolved locations. The calibration sandwich
  # treats rows as independent, so it is not a sampling covariance for this
  # design and cannot support Wald inference.
  repeated_person <- FALSE
  if (!is.null(fit$person$id)) {
    idv <- as.character(fit$person$id)
    seen <- !is.na(grp)
    repeated_person <- anyDuplicated(idv[seen]) > 0L
    if (repeated_person)
      notes <- c(notes, paste(
        "person identifiers repeat across response rows: resolved point",
        "differences remain descriptive, but sampling SEs, confidence",
        "intervals and Wald tests are withheld; use dif_contrasts for",
        "person-level inference or a whole-person bootstrap"))
  }

  # drop levels too thin on this item to resolve
  obs_i <- if (mfrm_item)
    rowSums(!is.na(fit$X[, fit$virtual_map$vkey[
      fit$virtual_map$item == item], drop = FALSE])) > 0L
  else !is.na(fit$X[, i])
  n_lev <- table(grp[obs_i & !is.na(grp)])
  thin <- names(n_lev)[n_lev < min_n]
  if (length(thin)) {
    notes <- c(notes, sprintf("level(s) dropped with fewer than %d responders: %s",
                              min_n, paste(thin, collapse = ", ")))
    grp <- factor(ifelse(as.character(grp) %in% thin, NA, as.character(grp)))
  }
  if (nlevels(droplevels(grp)) < 2)
    stop("fewer than two usable levels for item ", item)
  grp <- droplevels(grp)

  if (mfrm_item) {
    # underlying MFRM item: pooled virtual-level resolution (one joint
    # unstructured refit of the virtual matrix; see .dif_resolve)
    rs <- .dif_resolve(fit, item, grp, min_n)
    if (is.null(rs))
      stop("could not resolve item ", item, " (too little data per level)")
    levs <- rs$levs; loc <- rs$loc; vloc <- rs$vloc; weak_lev <- rs$weak
    notes <- c(notes, rs$notes)
  } else {
    refit <- split_items(fit, item, by = grp)
    levs <- levels(grp)
    split_names <- paste0(item, " (", levs, ")")
    idx <- match(split_names, refit$items$item)
    if (anyNA(idx))
      stop("resolved item(s) missing after re-analysis (too little data): ",
           paste(split_names[is.na(idx)], collapse = ", "))

    # location covariance from the sandwich: var(mean of a threshold block)
    thr <- refit$thresholds; cv <- refit$est$cov_tau
    block <- lapply(idx, function(k) thr$id[thr$item == k])
    loc <- refit$items$location[idx]
    vloc <- matrix(NA_real_, length(levs), length(levs))
    for (a in seq_along(levs)) for (b in seq_along(levs))
      vloc[a, b] <- mean(cv[block[[a]], block[[b]], drop = FALSE])
    weak_lev <- .dif_weak_levels(refit, as.list(idx)); names(weak_lev) <- levs
    if (any(weak_lev))
      notes <- c(notes, sprintf(
        "location(s) for level(s) %s rest on a near-empty category and are weakly identified; their DIF magnitude and significance are withheld",
        paste(levs[weak_lev], collapse = ", ")))
  }
  n_item <- as.integer(table(grp[obs_i & !is.na(grp)])[levs])
  lev_se <- sqrt(pmax(diag(vloc), 0)); lev_se[weak_lev] <- NA_real_
  if (repeated_person) lev_se[] <- NA_real_
  levels_df <- data.frame(level = levs, location = loc,
                          se = lev_se, weak = unname(weak_lev), n = n_item)

  pr <- t(utils::combn(seq_along(levs), 2))
  pair_weak <- weak_lev[pr[, 1]] | weak_lev[pr[, 2]]
  pairs <- data.frame(
    level_a = levs[pr[, 1]], level_b = levs[pr[, 2]],
    difference = loc[pr[, 1]] - loc[pr[, 2]],
    se = sqrt(pmax(diag(vloc)[pr[, 1]] + diag(vloc)[pr[, 2]] -
                   2 * vloc[cbind(pr[, 1], pr[, 2])], 1e-12)))
  # a pair touching a weakly-identified level carries no trustworthy
  # magnitude: withhold its SE and every SE-derived verdict
  pairs$se[pair_weak] <- NA_real_
  if (repeated_person) pairs$se[] <- NA_real_
  pairs$z <- pairs$difference / pairs$se
  pairs$p <- 2 * pnorm(-abs(pairs$z))
  pairs$p_adj <- p.adjust(pairs$p, method = p_adjust)
  pairs$lower <- pairs$difference - qnorm(0.975) * pairs$se
  pairs$upper <- pairs$difference + qnorm(0.975) * pairs$se
  pairs$significant <- ifelse(pair_weak | repeated_person, NA,
                              pairs$p_adj < alpha)
  pairs$practical <- ifelse(pair_weak, NA, abs(pairs$difference) >= flag_logits)
  pairs$ets <- .ets_category(pairs$difference, pairs$se, pairs$p_adj, alpha,
                             dichotomous = max(fit$m[i]) == 1L)

  out <- list(item = item, by = paste(names(factors), collapse = ":"),
              levels = levels_df, pairs = pairs, alpha = alpha,
              p_adjust = p_adjust, flag_logits = flag_logits, notes = notes)
  out <- .tag_tables(out)
  class(out) <- "rasch_dif_size"
  out
}

#' @export
print.rasch_dif_size <- function(x, ...) {
  cat(sprintf("DIF size for %s by %s (resolved locations, logits)\n",
              x$item, x$by))
  lv <- x$levels; lv[-1] <- lapply(lv[-1], round, 3)
  print(lv, row.names = FALSE)
  pr <- x$pairs
  num <- vapply(pr, is.numeric, TRUE)
  pr[num] <- lapply(pr[num], round, 3)
  pr$significant <- ifelse(pr$significant, "*", "")
  pr$practical <- ifelse(pr$practical, sprintf(">= %.2f", x$flag_logits), "")
  print(pr, row.names = FALSE)
  cat(sprintf("p adjusted by %s over %d pairwise comparison(s); practical criterion %.2f logits\n",
              x$p_adjust, nrow(pr), x$flag_logits))
  if (length(x$notes)) cat("notes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Planned contrasts: the confirmatory alternative to exhaustive post-hoc
# pairwise comparison. The family of questions is derived from the structure
# of the nominated factors (or supplied), estimated on the logit scale from
# resolved item locations, and -- when persons repeat across rows of a
# stacked design -- tested from person-level contrast scores so that
# within-subject dependence is respected.
# ---------------------------------------------------------------------------

# Resolve one item over grouping cells: locations and sandwich covariance.
# A resolved level is weakly identified when its split copy rests on a
# near-empty category: split_items() already flags such thresholds
# (weak = TRUE, se = NA) and leaves the item-location SE NA. A location
# built on such a threshold is a boundary artefact, so its DIF magnitude
# and significance must be withheld rather than recomputed from the ridged
# covariance -- otherwise dif_size()/dif_contrasts() report a fabricated
# finite SE and a spurious 'significant'/'practical' verdict. item_rows is
# a list, one entry per level, of the refit$items row-index(es) whose
# thresholds back that level's location.
.dif_weak_levels <- function(refit, item_rows) {
  wk <- refit$thresholds$weak
  se <- refit$items$se
  vapply(item_rows, function(ks) any(vapply(ks, function(k) {
    (!is.null(wk) && isTRUE(any(wk[refit$thresholds$item == k], na.rm = TRUE))) ||
      (!is.null(se) && k <= length(se) && is.na(se[k]))
  }, logical(1))), logical(1))
}

.dif_resolve <- function(fit, item, grp, min_n) {
  # an UNDERLYING MFRM item resolves at the virtual level: every one of
  # its facet cells is split by the groups in one joint unstructured
  # refit of the virtual matrix (the facet decomposition is not
  # reimposed), and the per-level locations are precision-weighted means
  # over the item's cells with the full covariance carried
  if (inherits(fit, "rasch_mfrm") && !is.null(fit$virtual_map) &&
      !(item %in% colnames(fit$X)) && item %in% fit$virtual_map$item) {
    vm <- fit$virtual_map
    cols <- vm$vkey[vm$item == item]
    notes <- paste0(item, ": resolved at the virtual-item level, pooled ",
                    "over its facet cells (facet structure not reimposed)")
    obs <- rowSums(!is.na(fit$X[, cols, drop = FALSE])) > 0L
    n_lev <- table(grp[obs & !is.na(grp)])
    thin <- names(n_lev)[n_lev < min_n]
    if (length(thin)) {
      notes <- c(notes, sprintf(
        "%s: level(s) dropped with fewer than %d responders: %s",
        item, min_n, paste(thin, collapse = ", ")))
      grp <- factor(ifelse(as.character(grp) %in% thin, NA,
                           as.character(grp)))
    }
    grp <- droplevels(grp)
    if (nlevels(grp) < 2) return(NULL)
    vfit <- fit; class(vfit) <- "rasch"
    vfit$model <- "PCM"   # unstructured virtual thresholds refit as PCM
    refit <- tryCatch(split_items(vfit, cols, by = grp),
                      error = function(e) NULL)
    if (is.null(refit)) return(NULL)
    levs <- levels(grp)
    thr <- refit$thresholds; cv <- refit$est$cov_tau
    # COMMON cells with COMMON weights: every facet cell used must be
    # resolved for EVERY level, and each cell gets one weight shared by
    # all levels, so the cell's facet severity cancels exactly from every
    # level contrast. Group-specific precision weights let severity leak
    # into the DIF magnitude when groups have different facet exposure
    # (a no-DIF design with sex-linked rater allocation read -1.75
    # logits, z = -6.5).
    idx_m <- sapply(levs, function(l)
      match(paste0(cols, " (", l, ")"), refit$items$item))
    if (is.null(dim(idx_m))) idx_m <- matrix(idx_m, nrow = length(cols))
    common <- rowSums(is.na(idx_m)) == 0L
    if (sum(common) < 1L) return(NULL)
    if (any(!common))
      notes <- c(notes, sprintf(
        "%s: facet cell(s) dropped from the magnitude (not resolvable for every level): %s",
        item, paste(cols[!common], collapse = ", ")))
    idx_m <- idx_m[common, , drop = FALSE]
    blocks <- lapply(seq_len(nrow(idx_m)), function(ci)
      lapply(idx_m[ci, ], function(k) thr$id[thr$item == k]))
    # one weight per cell: inverse of the level-averaged location variance
    vr_c <- vapply(blocks, function(bl)
      mean(vapply(bl, function(rws) mean(cv[rws, rws]), 0)), 0)
    w_c <- 1 / pmax(vr_c, 1e-10); w_c <- w_c / sum(w_c)
    loc <- vapply(seq_along(levs), function(a)
      sum(w_c * vapply(seq_along(blocks), function(ci)
        refit$items$location[idx_m[ci, a]], 0)), 0)
    vloc <- matrix(NA_real_, length(levs), length(levs))
    for (a in seq_along(levs)) for (b in seq_along(levs)) {
      acc <- 0
      for (ca in seq_along(blocks)) for (cb in seq_along(blocks))
        acc <- acc + w_c[ca] * w_c[cb] *
          mean(cv[blocks[[ca]][[a]], blocks[[cb]][[b]], drop = FALSE])
      vloc[a, b] <- acc
    }
    weak_lev <- .dif_weak_levels(refit, lapply(seq_along(levs),
                                               function(a) idx_m[, a]))
    names(weak_lev) <- levs
    if (any(weak_lev))
      notes <- c(notes, sprintf(
        "%s: location(s) for level(s) %s rest on a near-empty category and are weakly identified; their DIF magnitude and significance are withheld",
        item, paste(levs[weak_lev], collapse = ", ")))
    return(list(levs = levs, loc = loc, vloc = vloc, weak = weak_lev,
                notes = notes))
  }
  i <- .item_idx(fit, item)
  item <- fit$items$item[i]
  notes <- character(0)
  n_lev <- table(grp[!is.na(fit$X[, i]) & !is.na(grp)])
  thin <- names(n_lev)[n_lev < min_n]
  if (length(thin)) {
    notes <- c(notes, sprintf(
      "%s: level(s) dropped with fewer than %d responders: %s",
      item, min_n, paste(thin, collapse = ", ")))
    grp <- factor(ifelse(as.character(grp) %in% thin, NA, as.character(grp)))
  }
  grp <- droplevels(grp)
  if (nlevels(grp) < 2) return(NULL)
  refit <- split_items(fit, item, by = grp)
  levs <- levels(grp)
  idx <- match(paste0(item, " (", levs, ")"), refit$items$item)
  if (anyNA(idx)) return(NULL)
  thr <- refit$thresholds; cv <- refit$est$cov_tau
  block <- lapply(idx, function(k) thr$id[thr$item == k])
  loc <- refit$items$location[idx]
  vloc <- matrix(NA_real_, length(levs), length(levs))
  for (a in seq_along(levs)) for (b in seq_along(levs))
    vloc[a, b] <- mean(cv[block[[a]], block[[b]], drop = FALSE])
  weak_lev <- .dif_weak_levels(refit, as.list(idx)); names(weak_lev) <- levs
  if (any(weak_lev))
    notes <- c(notes, sprintf(
      "%s: location(s) for level(s) %s rest on a near-empty category and are weakly identified; their DIF magnitude and significance are withheld",
      item, paste(levs[weak_lev], collapse = ", ")))
  list(levs = levs, loc = loc, vloc = vloc, weak = weak_lev, notes = notes)
}

# A factor is treated as ordered when declared ordered or when its levels
# parse as numbers (ages, waves, doses).
.dif_is_ordered <- function(f)
  is.ordered(f) || !any(is.na(suppressWarnings(as.numeric(levels(f)))))

# The leading contrast of a factor: the difference for two levels, the
# linear trend for an ordered factor, none for a nominal many-level factor.
.dif_leading <- function(f) {
  K <- nlevels(f)
  if (K == 2L) {
    w <- c(-1, 1); names(w) <- levels(f)
    list(weights = w,
         label = sprintf("%s - %s", levels(f)[2], levels(f)[1]))
  } else if (.dif_is_ordered(f)) {
    sc <- suppressWarnings(as.numeric(levels(f)))
    cp <- if (!any(is.na(sc))) stats::contr.poly(K, scores = sc)
          else stats::contr.poly(K)
    w <- cp[, 1]; names(w) <- levels(f)
    list(weights = w, label = "linear")
  } else NULL
}

# The planned questions a single factor admits.
.dif_factor_contrasts <- function(f, fname) {
  K <- nlevels(f); out <- list()
  if (K == 2L) {
    lead <- .dif_leading(f)
    out[[sprintf("%s: %s", fname, lead$label)]] <- lead$weights
  } else if (.dif_is_ordered(f)) {
    sc <- suppressWarnings(as.numeric(levels(f)))
    cp <- if (!any(is.na(sc))) stats::contr.poly(K, scores = sc)
          else stats::contr.poly(K)
    w1 <- cp[, 1]; names(w1) <- levels(f)
    out[[sprintf("%s: linear", fname)]] <- w1
    w2 <- cp[, 2]; names(w2) <- levels(f)
    out[[sprintf("%s: quadratic", fname)]] <- w2
  } else if (K <= 4L) {
    pr <- utils::combn(levels(f), 2)
    for (j in seq_len(ncol(pr))) {
      w <- stats::setNames(numeric(K), levels(f))
      w[pr[2, j]] <- 1; w[pr[1, j]] <- -1
      out[[sprintf("%s: %s - %s", fname, pr[2, j], pr[1, j])]] <- w
    }
  } else {
    for (l in levels(f)) {
      w <- stats::setNames(rep(-1 / (K - 1), K), levels(f)); w[l] <- 1
      out[[sprintf("%s: %s - others", fname, l)]] <- w
    }
  }
  out
}

# Scale cell weights so the positive and negative parts each sum to one:
# every contrast then reads as a difference between two weighted averages,
# in logits, comparable across contrasts and against the practical flag.
.dif_norm <- function(w) {
  w[is.na(w)] <- 0
  s <- sum(abs(w))
  if (s < 1e-10) return(NULL)
  w * 2 / s
}

# Spread factor-level weights over the design cells (unweighted marginal
# means: each level's weight is shared equally by the cells carrying it).
.dif_cell_weights <- function(cellmap, fw, fname) {
  lev <- as.character(cellmap[[fname]])
  keep <- lev %in% names(fw)
  w <- stats::setNames(numeric(nrow(cellmap)), cellmap$cell)
  w[keep] <- fw[lev[keep]] / as.numeric(table(lev[keep])[lev[keep]])
  w
}

# Derive the planned family from the factor structure.
.dif_contrast_family <- function(factors, cellmap, within_names) {
  fam <- list(); meta <- list()
  for (fname in names(factors)) {
    fc <- .dif_factor_contrasts(factors[[fname]], fname)
    for (nm in names(fc)) {
      w <- .dif_norm(.dif_cell_weights(cellmap, fc[[nm]], fname))
      if (is.null(w)) next
      fam[[nm]] <- w
      meta[[nm]] <- list(factors = fname, fweights = fc[nm],
                         within = fname %in% within_names)
    }
  }
  fns <- names(factors)
  if (length(fns) >= 2) for (a in seq_len(length(fns) - 1))
    for (b in seq(a + 1, length(fns))) {
      la <- .dif_leading(factors[[fns[a]]])
      lb <- .dif_leading(factors[[fns[b]]])
      if (is.null(la) || is.null(lb)) next
      wa <- la$weights[as.character(cellmap[[fns[a]]])]
      wb <- lb$weights[as.character(cellmap[[fns[b]]])]
      key <- paste(cellmap[[fns[a]]], cellmap[[fns[b]]])
      w <- .dif_norm(stats::setNames(
        wa * wb / as.numeric(table(key)[key]), cellmap$cell))
      if (is.null(w)) next
      nm <- sprintf("%s(%s) x %s(%s)", fns[a], la$label, fns[b], lb$label)
      fam[[nm]] <- w
      meta[[nm]] <- list(factors = c(fns[a], fns[b]),
                         fweights = list(la$weights, lb$weights),
                         within = any(c(fns[a], fns[b]) %in% within_names))
    }
  list(family = fam, meta = meta)
}

# Welch test of a linear combination of independent group means.
.welch_contrast <- function(vals, g, w) {
  ok <- !is.na(vals) & !is.na(g)
  vals <- vals[ok]; g <- droplevels(factor(g[ok]))
  w <- w[levels(g)]
  if (any(is.na(w)) || sum(abs(w)) < 1e-10) return(NULL)
  m <- tapply(vals, g, mean); v <- tapply(vals, g, stats::var)
  n <- tapply(vals, g, length)
  if (any(n < 2)) return(NULL)
  vv <- sum(w^2 * v / n)
  if (!is.finite(vv) || vv <= 0) return(NULL)
  df <- vv^2 / sum((w^2 * v / n)^2 / (n - 1))
  t <- sum(w * m) / sqrt(vv)
  list(stat = t, df = df, p = 2 * stats::pt(-abs(t), df))
}

# Test any resolved-cell contrast in a stacked design without treating rows
# from the same person as independent. The supplied cell weights already
# encode the desired marginal comparison. Within each between-person cell we
# first form one weighted residual score per person, then combine the
# independent cell means with a Welch--Satterthwaite reference.
.dif_paired_cell_contrast <- function(z, factors, grp, id, within,
                                      cellmap, weights) {
  id <- as.character(id)
  between <- setdiff(names(factors), within)
  bkey <- if (length(between))
    do.call(paste, c(lapply(factors[between], as.character), list(sep = "\r")))
  else rep("all", nrow(factors))
  map_bkey <- if (length(between))
    do.call(paste, c(lapply(cellmap[between], as.character), list(sep = "\r")))
  else rep("all", nrow(cellmap))
  names(weights) <- cellmap$cell

  people <- split(seq_along(id), id)
  score <- group <- rep(NA_character_, length(people))
  score_num <- rep(NA_real_, length(people))
  for (j in seq_along(people)) {
    r <- people[[j]]
    ok <- is.finite(z[r]) & !is.na(grp[r]) & !is.na(bkey[r])
    if (!any(ok)) next
    r <- r[ok]
    bg <- unique(bkey[r])
    if (length(bg) != 1L) next
    use_cells <- cellmap$cell[map_bkey == bg & weights != 0]
    if (!length(use_cells)) next
    means <- tapply(z[r], as.character(grp[r]), mean)
    if (anyNA(match(use_cells, names(means)))) next
    score_num[j] <- sum(weights[use_cells] * means[use_cells])
    group[j] <- bg
  }
  ok <- is.finite(score_num) & !is.na(group)
  if (sum(ok) < 3L) return(NULL)
  sp <- split(score_num[ok], group[ok])
  sp <- sp[lengths(sp) >= 2L]
  if (!length(sp)) return(NULL)
  means <- vapply(sp, mean, 0)
  vars <- vapply(sp, stats::var, 0)
  ns <- lengths(sp)
  parts <- vars / ns
  vv <- sum(parts)
  if (!is.finite(vv) || vv <= 0) return(NULL)
  df <- vv^2 / sum(parts^2 / (ns - 1))
  stat <- sum(means) / sqrt(vv)
  list(stat = stat, df = df, p = 2 * stats::pt(-abs(stat), df))
}

# Pairwise marginal comparisons for one term. For a main effect these are
# differences between factor levels, averaged equally over complete cells of
# the remaining factors. For an interaction they are tensor products of the
# level differences: difference-in-differences for two factors and the direct
# higher-order analogue beyond two.
.dif_posthoc_family <- function(factors, cellmap, term, within) {
  target <- .term_vars(term)
  bad <- setdiff(target, names(factors))
  if (length(bad))
    stop("term factor(s) not found: ", paste(bad, collapse = ", "))
  pairs <- lapply(target, function(fn) {
    lv <- levels(factors[[fn]])
    if (length(lv) < 2L) return(list())
    pr <- utils::combn(lv, 2)
    lapply(seq_len(ncol(pr)), function(j) {
      w <- stats::setNames(numeric(length(lv)), lv)
      w[pr[1, j]] <- -1; w[pr[2, j]] <- 1
      list(weights = w, label = sprintf("%s - %s", pr[2, j], pr[1, j]))
    })
  })
  if (any(!lengths(pairs))) stop("every term factor needs at least two levels")
  grid <- expand.grid(lapply(pairs, seq_along), KEEP.OUT.ATTRS = FALSE)
  nuisance <- setdiff(names(factors), target)
  nkey <- if (length(nuisance))
    do.call(paste, c(lapply(cellmap[nuisance], as.character), list(sep = "\r")))
  else rep("all", nrow(cellmap))
  family <- meta <- list()
  for (r in seq_len(nrow(grid))) {
    chosen <- lapply(seq_along(target), function(j) pairs[[j]][[grid[r, j]]])
    raw <- rep(1, nrow(cellmap))
    active <- rep(TRUE, nrow(cellmap))
    for (j in seq_along(target)) {
      fw <- chosen[[j]]$weights
      cw <- unname(fw[as.character(cellmap[[target[j]]])])
      cw[is.na(cw)] <- 0
      raw <- raw * cw
      active <- active & cw != 0
    }
    # Marginalise only over nuisance strata containing the complete target
    # contrast. This avoids changing the estimand when an unbalanced design
    # has a structurally absent target cell.
    complete_n <- names(which(tapply(active, nkey, sum) == 2^length(target)))
    raw[!active | !nkey %in% complete_n] <- 0
    if (!length(complete_n)) next
    w <- stats::setNames(raw / length(complete_n), cellmap$cell)
    label <- paste(vapply(chosen, `[[`, "", "label"), collapse = " x ")
    family[[label]] <- w
    meta[[label]] <- list(
      factors = target,
      fweights = lapply(chosen, `[[`, "weights"),
      within = any(target %in% within))
  }
  if (!length(family))
    stop("no complete cells support post-hoc comparisons for term '", term, "'")
  list(family = family, meta = meta)
}

#' Planned DIF contrasts
#'
#' Tests a specified family of one-degree-of-freedom DIF contrasts. By default,
#' contrasts are derived from the factor structure: differences for two-level
#' factors, polynomial trends for ordered factors, and pairwise or
#' level-against-rest comparisons for nominal factors. Leading contrasts are
#' crossed to form two-factor interactions. User-supplied cell weights are
#' also accepted.
#'
#' @details
#' Each logit contrast is calculated from resolved item locations. Weights are
#' scaled so their positive and negative parts each sum to one. With repeated
#' persons, inference uses person-level residual contrast scores: within-person
#' contrasts are tested against zero, between-person contrasts use person
#' means, and mixed interactions compare within-person contrast scores between
#' groups. The resolved logit estimate is retained, but its calibration-based
#' standard error is withheld because it does not include repeated-person
#' dependence.
#'
#' For independent rows, a contrast with weights \eqn{\mathbf{c}} is
#' \deqn{\Delta_i=\mathbf{c}^{\mathsf T}\delta_i,\qquad
#' \operatorname{SE}(\Delta_i)=
#' \sqrt{\mathbf{c}^{\mathsf T}\mathbf{V}_i\mathbf{c}}.}
#' In a repeated-measures design, a within-person contrast is formed from the
#' standardised residuals,
#' \deqn{s_p=\sum_l c_l z_{pl},}
#' and tested over persons. Between-person contrasts use Welch tests of person
#' means; mixed interactions use Welch tests of the within-person contrast
#' scores. The sign of each residual test is aligned with the resolved logit
#' contrast.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param factors A data frame of person factors, a character vector naming
#'   factors nominated in the fit, or a single grouping vector. Defaults to
#'   every factor stored in the fit.
#' @param items Item names or indices to test; all items by default.
#' @param within Names of factors that vary within person (for example
#'   time). Detected automatically when \code{id} is supplied and a factor
#'   varies within an id.
#' @param id Person identifier with one entry per row, or the name of a
#'   nominated factor holding it; required for stacked designs where the
#'   same person occupies several rows.
#' @param contrasts \code{"auto"} (derive the family from the factor
#'   structure) or a named list of numeric cell-weight vectors, each named
#'   by the design-cell labels (factor levels joined by \code{":"}).
#'   Weights are rescaled so the positive and negative parts each sum to
#'   one.
#' @param p_adjust Familywise adjustment across items and contrasts. The
#'   default is \code{"holm"}.
#' @param alpha Significance level for the adjusted probabilities.
#' @param flag_logits Absolute estimate flagged as practically significant.
#' @param min_n Cells with fewer responders to an item are dropped from that
#'   item's resolution, with a note.
#' @return A list of class \code{"rasch_dif_contrasts"}: \code{table} (one row
#'   per item and contrast: estimate in logits, SE, statistic, df where a t
#'   test was used, raw and adjusted p, 95 per cent interval,
#'   \code{significant}, \code{practical}, \code{within}), \code{family}
#'   (the derived questions with their cell weights), the settings, and any
#'   \code{notes}.
#' @references
#' Maxwell, S. E. and Delaney, H. D. (2004). Designing Experiments and
#' Analyzing Data (2nd ed.). Erlbaum.
#'
#' Andrich, D. and Hagquist, C. (2015). Real and artificial differential item
#' functioning in polytomous items. Educational and Psychological Measurement,
#' 75(2), 185--207.
#'
#' Hagquist, C. and Andrich, D. (2017). Recent advances in analysis of
#' differential item functioning in health research using the Rasch model.
#' Health and Quality of Life Outcomes, 15, 181.
#' @seealso \code{\link{dif_anova}} and \code{\link{dif_size}}.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' dif_contrasts(fit, items = c("I3", "I5"))
#' @export
dif_contrasts <- function(fit, factors = NULL, items = NULL, within = NULL,
                          id = NULL, contrasts = "auto", p_adjust = "holm",
                          alpha = 0.05, flag_logits = 0.5, min_n = 20) {
  factors <- .dif_factors(fit, factors)
  factors <- as.data.frame(lapply(factors, function(v) {
    f <- droplevels(if (is.ordered(v)) v else factor(v))
    f
  }), check.names = FALSE, stringsAsFactors = FALSE)
  grp <- interaction(factors, sep = ":", drop = TRUE)
  cellmap <- unique(data.frame(cell = as.character(grp), factors,
                               check.names = FALSE))
  cellmap <- cellmap[match(levels(grp), cellmap$cell), , drop = FALSE]

  if (is.character(id) && length(id) == 1L && !is.null(fit$factors) &&
      id %in% names(fit$factors)) id <- fit$factors[[id]]
  if (is.null(within) && !is.null(id) && anyDuplicated(id)) {
    within <- names(factors)[vapply(names(factors), function(fn)
      any(tapply(as.character(factors[[fn]]), id, function(v)
        length(unique(v[!is.na(v)])) > 1L), na.rm = TRUE), TRUE)]
  }
  if (is.null(within)) within <- character(0)
  within <- intersect(within, names(factors))
  if (length(within) && is.null(id))
    stop("`within` requires `id` to pair a person's rows")
  paired <- !is.null(id) && anyDuplicated(id) > 0

  if (identical(contrasts, "auto")) {
    fam <- .dif_contrast_family(factors, cellmap, within)
  } else {
    if (!is.list(contrasts) || is.null(names(contrasts)))
      stop("`contrasts` must be \"auto\" or a named list of cell weights")
    supplied_meta <- attr(contrasts, "dif_meta", exact = TRUE)
    fam <- list(family = list(), meta = list())
    for (nm in names(contrasts)) {
      w <- contrasts[[nm]]
      if (is.null(names(w)) || !all(names(w) %in% cellmap$cell))
        stop("weights of contrast '", nm, "' must be named by design cells: ",
             paste(cellmap$cell, collapse = ", "))
      full <- stats::setNames(numeric(nrow(cellmap)), cellmap$cell)
      full[names(w)] <- w
      preserve_scale <- isTRUE(supplied_meta[[nm]]$preserve_scale)
      w <- if (preserve_scale) full else .dif_norm(full)
      if (is.null(w) || !any(w > 0) || !any(w < 0) ||
          abs(sum(w)) > 1e-8)
        stop("contrast '", nm, "' needs positive and negative weights that sum to zero")
      fam$family[[nm]] <- w
      fam$meta[[nm]] <- if (!is.null(supplied_meta[[nm]]))
        supplied_meta[[nm]]
      else list(factors = names(factors), fweights = NULL, within = FALSE)
    }
  }
  if (!length(fam$family)) stop("no contrasts could be formed")

  its <- if (is.null(items)) fit$items$item else
    fit$items$item[vapply(items, function(x) .item_idx(fit, x), 1L)]
  Z <- fit$residuals
  notes <- character(0)
  rows <- list()

  for (item in its) {
    i <- .item_idx(fit, item)
    rs <- .dif_resolve(fit, item, grp, min_n)
    if (!is.null(rs)) notes <- c(notes, rs$notes)
    for (nm in names(fam$family)) {
      w_full <- fam$family[[nm]]
      mt <- fam$meta[[nm]]
      est <- se <- stat <- df <- p <- NA_real_
      if (!is.null(rs)) {
        w <- w_full[rs$levs]
        w[is.na(w)] <- 0
        # a contrast placing weight on a weakly-identified level rests on a
        # boundary-artefact location: withhold its estimate and SE
        touches_weak <- !is.null(rs$weak) && any(w != 0 & rs$weak)
        preserve_scale <- isTRUE(mt$preserve_scale)
        complete_support <- all(names(w_full)[w_full != 0] %in% rs$levs)
        valid_weights <- sum(w > 0) > 0 && sum(w < 0) > 0 &&
          abs(sum(w)) < 1e-8 &&
          (if (preserve_scale) complete_support
           else abs(sum(abs(w)) - 2) < 0.5)
        if (valid_weights && !touches_weak) {
          if (!preserve_scale) w <- .dif_norm(w)
          est <- sum(w * rs$loc)
          se <- sqrt(max(drop(t(w) %*% rs$vloc %*% w), 1e-12))
        }
      }
      if (!paired) {
        if (is.finite(est)) {
          stat <- est / se
          p <- 2 * stats::pnorm(-abs(stat))
        }
      } else if (isTRUE(mt$within) && length(mt$factors) == 1L &&
                 mt$factors %in% within && !is.null(mt$fweights)) {
        # A single within-factor question is the ordinary paired contrast:
        # one complete contrast score per person, tested against zero.
        fw <- .dif_norm(mt$fweights[[1]])
        lev <- as.character(factors[[mt$factors]])
        zi <- Z[, i]
        ps <- split(seq_along(zi), id)
        psi <- vapply(ps, function(rws) {
          l <- lev[rws]
          if (anyDuplicated(l) || !all(names(fw) %in% l)) return(NA_real_)
          sum(fw * zi[rws][match(names(fw), l)])
        }, 0)
        psi <- psi[!is.na(psi)]
        if (length(psi) >= 10) {
          stat <- -mean(psi) / (stats::sd(psi) / sqrt(length(psi)))
          df <- length(psi) - 1
          p <- 2 * stats::pt(-abs(stat), df)
        }
      } else if (isTRUE(mt$within) && length(mt$factors) == 2L &&
                 sum(mt$factors %in% within) == 1L &&
                 !is.null(mt$fweights)) {
        # A two-factor mixed contrast compares one within-person contrast
        # score across the between-person groups.
        wf <- intersect(mt$factors, within)[1]
        bf <- setdiff(mt$factors, wf)[1]
        fw <- .dif_norm(mt$fweights[[match(wf, mt$factors)]])
        bw <- mt$fweights[[match(bf, mt$factors)]]
        lev <- as.character(factors[[wf]])
        blev <- as.character(factors[[bf]])
        zi <- Z[, i]
        ps <- split(seq_along(zi), id)
        psi <- vapply(ps, function(rws) {
          l <- lev[rws]
          if (anyDuplicated(l) || !all(names(fw) %in% l)) return(NA_real_)
          sum(fw * zi[rws][match(names(fw), l)])
        }, 0)
        pb <- vapply(ps, function(rws) blev[rws][1], "")
        ok <- !is.na(psi)
        bw_test <- if (isTRUE(mt$preserve_scale)) bw else bw / 2
        wc <- .welch_contrast(psi[ok], pb[ok], bw_test)
        if (!is.null(wc)) { stat <- -wc$stat; df <- wc$df; p <- wc$p }
      } else {
        # One generic person-level calculation covers between, within and
        # mixed contrasts, including interactions of several within factors.
        wc <- .dif_paired_cell_contrast(
          Z[, i], factors, grp, id, within, cellmap, w_full)
        if (!is.null(wc)) { stat <- -wc$stat; df <- wc$df; p <- wc$p }
      }
      if (paired) se <- NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        item = item, contrast = nm, within = isTRUE(mt$within),
        estimate = est, se = se, statistic = stat, df = df, p = p)
    }
  }
  tab <- do.call(rbind, rows)
  tab$p_adj <- stats::p.adjust(tab$p, method = p_adjust)
  tab$lower <- tab$estimate - stats::qnorm(0.975) * tab$se
  tab$upper <- tab$estimate + stats::qnorm(0.975) * tab$se
  tab$significant <- !is.na(tab$p_adj) & tab$p_adj < alpha
  tab$practical <- !is.na(tab$estimate) & abs(tab$estimate) >= flag_logits
  rownames(tab) <- NULL

  fam_df <- data.frame(
    contrast = names(fam$family),
    within = vapply(fam$meta, function(m) isTRUE(m$within), TRUE),
    cells = vapply(fam$family, function(w)
      paste(sprintf("%s %+0.2f", names(w)[w != 0], w[w != 0]),
            collapse = ", "), ""))
  rownames(fam_df) <- NULL

  out <- list(table = tab, family = fam_df, within = within,
              paired = paired, alpha = alpha, p_adjust = p_adjust,
              flag_logits = flag_logits, notes = unique(notes))
  out <- .tag_tables(out)
  class(out) <- "rasch_dif_contrasts"
  out
}

#' Pairwise follow-up comparisons for a DIF term
#'
#' Resolves one item's locations over the complete person-factor design and
#' follows up a selected main effect or interaction. Main effects are pairwise
#' marginal differences. Interactions are differences between those
#' differences, providing a logit-scale magnitude for the interaction itself.
#'
#' @details
#' For levels \eqn{a,b} of one factor, the comparison is
#' \deqn{\Delta_{ba}=\bar\delta_b-\bar\delta_a,}
#' where the bars average equally over complete cells of the other nominated
#' factors. For a two-factor interaction, levels \eqn{a,b} and \eqn{c,d} give
#' \deqn{\Delta_{ba\mathbin{:}dc}=
#' (\delta_{bd}-\delta_{ad})-(\delta_{bc}-\delta_{ac}).}
#' Higher-order interactions use the corresponding tensor-product contrast.
#' Standard errors use the full covariance of the resolved locations.
#'
#' This is the preferred follow-up to a significant DIF term with more than
#' two levels. It reports effects in Rasch logits, adjusts the chosen family
#' of comparisons, and uses person-level scores for repeated-measures designs.
#' Tukey's HSD in \code{\link{dif_anova}} instead compares residual means and
#' is limited to between-person analysis-of-variance terms.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or index.
#' @param term A factor name for a main effect, or factor names joined by
#'   \code{":"} for an interaction.
#' @param factors The complete person-factor design, specified as for
#'   \code{\link{dif_contrasts}}. Other factors are retained when calculating
#'   marginal comparisons.
#' @param within,id Within-person factor names and person identifiers, specified
#'   as for \code{\link{dif_contrasts}}.
#' @param p_adjust Familywise adjustment over this post-hoc family; default
#'   \code{"holm"}.
#' @param alpha Significance level for adjusted probabilities.
#' @param flag_logits Absolute logit magnitude flagged as practically
#'   important.
#' @param min_n Minimum responders required in a resolved design cell.
#' @return An object of class \code{"rasch_dif_posthoc"}, extending the
#'   \code{\link{dif_contrasts}} result. Its \code{table} contains the pairwise
#'   marginal differences or interaction contrasts, with logit estimates,
#'   standard errors where available, confidence intervals, raw and adjusted
#'   probabilities, and statistical and practical flags.
#' @references Holm, S. (1979). A simple sequentially rejective multiple test
#'   procedure. Scandinavian Journal of Statistics, 6(2), 65--70.
#' @seealso \code{\link{dif_anova}}, \code{\link{dif_size}}, and
#'   \code{\link{dif_contrasts}}.
#' @examples
#' set.seed(1); n <- 800
#' g <- factor(rep(c("A", "B", "C", "D"), each = n / 4))
#' sex <- factor(rep(c("female", "male"), length.out = n))
#' d <- seq(-1.5, 1.5, length.out = 6)
#' sh <- matrix(0, n, 6); sh[g == "D", 2] <- 0.8
#' X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
#' colnames(X) <- paste0("I", 1:6)
#' fit <- rasch(data.frame(X, group = g, sex = sex),
#'              factors = c("group", "sex"))
#' dif_posthoc(fit, "I2", term = "group")
#' @export
dif_posthoc <- function(fit, item, term, factors = NULL, within = NULL,
                        id = NULL, p_adjust = "holm", alpha = 0.05,
                        flag_logits = 0.5, min_n = 20) {
  if (length(term) != 1L || is.na(term) || !nzchar(term))
    stop("`term` must be one main-effect or interaction term")
  if (length(item) != 1L)
    stop("`item` must name one item; run dif_posthoc() per item so the ",
         "multiplicity adjustment covers one post-hoc family at a time")
  item_names <- if (!is.null(fit$items$item)) fit$items$item else colnames(fit$X)
  ok_item <- if (is.character(item)) item %in% item_names
             else is.finite(item) && item >= 1 && item <= length(item_names)
  if (!isTRUE(ok_item))
    stop("item '", item, "' not found in the fit (items: ",
         paste(utils::head(item_names, 8), collapse = ", "),
         if (length(item_names) > 8) ", ..." else "", ")")
  factors <- .dif_factors(fit, factors)
  factors <- as.data.frame(lapply(factors, function(v)
    droplevels(if (is.ordered(v)) v else factor(v))),
    check.names = FALSE, stringsAsFactors = FALSE)
  if (is.character(id) && length(id) == 1L && !is.null(fit$factors) &&
      id %in% names(fit$factors)) id <- fit$factors[[id]]
  if (is.null(id) && !is.null(fit$person$id)) id <- fit$person$id
  if (is.null(within) && !is.null(id) && anyDuplicated(id)) {
    within <- names(factors)[vapply(names(factors), function(fn)
      any(tapply(as.character(factors[[fn]]), id, function(v)
        length(unique(v[!is.na(v)])) > 1L), na.rm = TRUE), TRUE)]
  }
  if (is.null(within)) within <- character(0)
  unknown_within <- setdiff(within, names(factors))
  if (length(unknown_within))
    stop("within-subject factor(s) not found: ",
         paste(unknown_within, collapse = ", "))

  grp <- interaction(factors, sep = ":", drop = TRUE)
  cellmap <- unique(data.frame(cell = as.character(grp), factors,
                               check.names = FALSE))
  cellmap <- cellmap[match(levels(grp), cellmap$cell), , drop = FALSE]
  fam <- .dif_posthoc_family(factors, cellmap, term, within)
  contrasts <- fam$family
  fam$meta <- lapply(fam$meta, function(x) {
    x$preserve_scale <- TRUE
    x
  })
  attr(contrasts, "dif_meta") <- fam$meta
  out <- dif_contrasts(
    fit, factors = factors, items = item, within = within, id = id,
    contrasts = contrasts, p_adjust = p_adjust, alpha = alpha,
    flag_logits = flag_logits, min_n = min_n)
  out$term <- term
  out$type <- if (grepl(":", term, fixed = TRUE))
    "interaction magnitude" else "pairwise marginal difference"
  if (nrow(out$table) && all(!is.finite(out$table$estimate)))
    stop("no contrast in the '", term, "' family is estimable for item '",
         item, "': every design cell fell below min_n = ", min_n,
         " responders, or the resolved refits were not identified; lower ",
         "min_n, pool sparse levels, or check the factor coding")
  class(out) <- c("rasch_dif_posthoc", class(out))
  out
}

#' @export
print.rasch_dif_posthoc <- function(x, ...) {
  cat(sprintf("DIF follow-up for %s (%s; %s)\n",
              x$term, x$type, x$p_adjust))
  show <- x$table[, c("item", "contrast", "estimate", "se", "statistic",
                      "p_adj", "significant", "practical")]
  print(.fmt_df(show), row.names = FALSE)
  if (length(x$notes)) cat("\n", paste(x$notes, collapse = "\n"), "\n", sep = "")
  invisible(x)
}

#' @export
print.rasch_dif_contrasts <- function(x, ...) {
  cat("Planned DIF contrasts (", nrow(x$family), " questions x ",
      length(unique(x$table$item)), " items; ", x$p_adjust,
      " over the family)\n", sep = "")
  for (r in seq_len(nrow(x$family)))
    cat(sprintf("  %s%s\n", x$family$contrast[r],
                if (x$family$within[r]) "  [within subjects]" else ""))
  if (x$paired)
    cat("Stacked design: tests use person-level residual scores;",
        "logit SEs and intervals are withheld.\n")
  cat("\n")
  tab <- x$table
  show <- tab[, c("item", "contrast", "estimate", "se", "statistic", "p_adj",
                  "significant", "practical")]
  print(.fmt_df(show), row.names = FALSE)
  if (length(x$notes)) cat("\n", paste(x$notes, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
