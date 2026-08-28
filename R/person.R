# rasch :: person estimation
# ===========================================================================
# Warm (1989) weighted likelihood person estimates. The WLE is finite at the
# extreme (zero and maximum) raw scores, so extreme persons receive usable
# locations; they are still flagged so reliability and test-of-fit statistics
# can exclude them. Persons with missing responses are
# estimated on their own observed item subset, grouped by missing-data
# pattern for speed.
# ===========================================================================

#' Category-score moments for a polytomous item
#'
#' @param theta Person location, in logits.
#' @param tau_i Numeric vector of the item's threshold parameters.
#' @param disc Discrimination (frame unit) multiplier on the exponent; 1 for
#'   the ordinary Rasch model.
#' @return A list with category probabilities \code{P}, expected score \code{E},
#'   variance \code{V}, and third and fourth central moments \code{mu3},
#'   \code{mu4}.
#' @examples
#' item_moments(0.5, c(-1, 0, 1))
#' @export
item_moments <- function(theta, tau_i, disc = 1) {
  if (length(theta) != 1L || !is.numeric(theta) || !is.finite(theta))
    stop("`theta` must be one finite location")
  if (!is.numeric(tau_i) || !length(tau_i) || any(!is.finite(tau_i)))
    stop("`tau_i` must be a non-empty vector of finite thresholds")
  if (length(disc) != 1L || !is.numeric(disc) || !is.finite(disc) || disc <= 0)
    stop("`disc` must be one positive finite discrimination")
  m <- length(tau_i); x <- 0:m
  lp <- disc * (x * theta - c(0, cumsum(tau_i)))
  num <- exp(lp - max(lp)); P <- num / sum(num)   # log-sum-exp: no overflow
  E <- sum(x * P); d <- x - E
  list(P = P, E = E, V = sum(d^2 * P), mu3 = sum(d^3 * P), mu4 = sum(d^4 * P))
}

#' Warm's weighted likelihood estimates by raw score
#'
#' Computes the weighted likelihood estimate (WLE) of person location for every
#' possible raw score on a set of items, with standard errors. WLE estimates
#' are finite at the extreme (zero and maximum) scores, unlike the maximum
#' likelihood estimate.
#'
#' @details
#' For raw score \eqn{R}, let \eqn{E(\theta)}, \eqn{V(\theta)}, and
#' \eqn{\mu_3(\theta)} be the sums of the item expected scores, variances, and
#' third central moments. The estimate solves Warm's weighted score equation
#' \deqn{R-E(\theta)+\frac{\mu_3(\theta)}{2V(\theta)}=0.}
#' With common discrimination \eqn{d}, its explicit multiplier cancels from
#' this equation, although the moments are evaluated under \eqn{d}. The
#' reported standard error is
#' \deqn{\operatorname{SE}(\hat{\theta})=
#' \{d^2V(\hat{\theta})\}^{-1/2}.}
#'
#' @param tau_list List of per-item threshold vectors.
#' @param disc Common discrimination (frame unit) of the items; with a
#'   constant discrimination the raw score remains sufficient.
#' @return A list with \code{theta} and \code{se}, each named by raw score.
#' @references
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. Psychometrika, 54(3), 427--450.
#' @seealso \code{\link{score_table}} and \code{\link{person_extrapolated}}.
#' @examples
#' person_wle(list(c(-1, 0), c(-0.5, 0.5), c(0, 1)))
#' @export
person_wle <- function(tau_list, disc = 1) {
  if (!is.list(tau_list) || !length(tau_list) ||
      !all(vapply(tau_list, function(t)
        is.numeric(t) && length(t) >= 1L && all(is.finite(t)), TRUE)))
    stop("`tau_list` must be a list of non-empty finite threshold vectors")
  if (length(disc) != 1L || !is.numeric(disc) || !is.finite(disc) ||
      disc <= 0)
    stop("`disc` must be one positive finite discrimination: the raw-score ",
         "WLE requires a common discrimination")
  Smax <- sum(vapply(tau_list, length, 1L))
  theta <- se <- setNames(rep(NA_real_, Smax + 1L), as.character(0:Smax))
  for (R in 0:Smax) {
    g <- function(th) {
      mo <- lapply(tau_list, item_moments, theta = th, disc = disc)
      E  <- sum(vapply(mo, `[[`, 0, "E"));  V <- sum(vapply(mo, `[[`, 0, "V"))
      m3 <- sum(vapply(mo, `[[`, 0, "mu3"))
      # Warm's weighted score is disc*(R - E) + disc^3 mu3 / (2 disc^2 V):
      # the discrimination cancels throughout, so the correction carries NO
      # disc factor (a stray disc here biased WLEs by up to 0.26 logits at
      # disc = 0.5; the vector-disc EFRM path was already correct)
      (R - E) + m3 / (2 * V)
    }
    root <- tryCatch(uniroot(g, c(-30, 30), tol = 1e-9)$root, error = function(e) NA_real_)
    theta[as.character(R)] <- root
    if (!is.na(root)) {
      V <- sum(vapply(lapply(tau_list, item_moments, theta = root, disc = disc),
                      `[[`, 0, "V"))
      se[as.character(R)] <- 1 / sqrt(disc^2 * V)      # WLE SE ~ 1/sqrt(information)
    }
  }
  list(theta = theta, se = se)
}

#' Person estimates with externally imposed weights
#'
#' Calculates a second set of person measures after assigning relative
#' weights to items or item sets. The fitted calibration is not changed.
#' In particular, these estimates do not replace the ordinary Rasch person
#' measures used for fit, reliability, targeting or DIF.
#'
#' @details
#' Let \eqn{q_i} be the external weight and \eqn{a_i} the model unit for
#' response cell \eqn{i}. Write
#' \eqn{H(\theta)=\sum_i q_i a_i^2V_i(\theta)} and
#' \eqn{J(\theta)=\sum_i q_i^2a_i^2V_i(\theta)}. The estimate solves the
#' externally weighted Warm score equation
#' \deqn{\sum_i q_i a_i\{x_i-E_i(\theta)\}+
#' \frac{J(\theta)\sum_i q_i a_i^3\mu_{3i}(\theta)}
#' {2H(\theta)^2}=0.}
#' Its standard error is the sandwich form
#' \deqn{\operatorname{SE}(\hat\theta)=
#' \frac{\{\sum_iq_i^2a_i^2V_i(\hat\theta)\}^{1/2}}
#' {\sum_iq_ia_i^2V_i(\hat\theta)}.}
#' This matters because an external weight changes the estimating equation;
#' it does not create independent replications of an item. With equal weights
#' the equation and standard error reduce to the person estimates in
#' \code{fit$person}. Weights are normalised to mean one over the fitted
#' response cells, so only their relative values matter. A zero weight omits
#' that item or set from the secondary measure.
#'
#' For an MFRM fit, an item weight applies to all response cells belonging to
#' that item. For an EFRM fit, \code{by = "set"} uses the fitted item-set map
#' unless \code{sets} is supplied.
#'
#' @param fit A fitted object from \code{\link{rasch}},
#'   \code{\link{rasch_mfrm}} or \code{\link{rasch_efrm}}.
#' @param weights A named numeric vector of non-negative finite relative
#'   weights. Names must identify every fitted item when \code{by = "item"},
#'   or every item set when \code{by = "set"}.
#' @param by Whether \code{weights} names individual \code{"item"}s or
#'   \code{"set"}s.
#' @param sets For \code{by = "set"}, a named character vector mapping items
#'   to sets, or a named list whose elements contain item names. It can be
#'   omitted for an EFRM fit, which already contains this map.
#' @return A data frame with the person identifiers and factors, observed
#'   item count, raw and externally weighted scores, maximum scores, weighted
#'   location, sandwich standard error and extreme-score flag. The resolved
#'   item weights are retained in the \code{"weighting"} attribute.
#' @references
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. Psychometrika, 54(3), 427--450.
#' @examples
#' set.seed(1)
#' d <- simulate_rasch(n_persons = 200, n_items = 6)
#' fit <- rasch(d)
#' weighted_person_estimates(
#'   fit, setNames(c(2, 2, 1, 1, 0.5, 0.5), colnames(fit$X)))
#' @export
weighted_person_estimates <- function(fit, weights,
                                      by = c("item", "set"), sets = NULL) {
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    stop("`fit` must be an ordinary, multiple-ratings or extended-frames Rasch fit",
         call. = FALSE)
  by <- match.arg(by)
  if (!is.numeric(weights) || !length(weights) || is.null(names(weights)) ||
      any(is.na(names(weights))) || any(!nzchar(names(weights))) ||
      anyDuplicated(names(weights)) || any(!is.finite(weights)) ||
      any(weights < 0))
    stop("`weights` must be a non-empty named numeric vector of non-negative finite values",
         call. = FALSE)
  if (!any(weights > 0))
    stop("at least one external weight must be positive", call. = FALSE)

  X <- fit$X
  if (!is.matrix(X)) X <- as.matrix(X)
  cells <- colnames(X)
  if (is.null(cells) || !length(cells) || length(fit$tau_list) != ncol(X))
    stop("the fitted response cells and thresholds are incomplete",
         call. = FALSE)
  if (inherits(fit, c("rasch_mfrm", "rasch_efrm"))) {
    vm <- fit$virtual_map
    if (is.null(vm) || !all(c("vkey", "item") %in% names(vm)))
      stop("the fitted response-cell map is unavailable", call. = FALSE)
    ii <- match(cells, as.character(vm$vkey))
    if (anyNA(ii))
      stop("the fitted response-cell map does not cover every response cell",
           call. = FALSE)
    source_item <- as.character(vm$item[ii])
  } else {
    source_item <- cells
  }
  items <- unique(source_item)

  set_of <- NULL
  if (by == "item") {
    if (!is.null(sets))
      stop("`sets` is used only when `by = \"set\"`", call. = FALSE)
    missing_w <- setdiff(items, names(weights))
    extra_w <- setdiff(names(weights), items)
    if (length(missing_w) || length(extra_w))
      stop("item weights must name every fitted item exactly",
           if (length(missing_w)) paste0("; missing: ", paste(missing_w, collapse = ", ")) else "",
           if (length(extra_w)) paste0("; unknown: ", paste(extra_w, collapse = ", ")) else "",
           call. = FALSE)
    q <- unname(weights[source_item])
  } else {
    if (is.null(sets) && inherits(fit, "rasch_efrm")) sets <- fit$set_of
    if (is.null(sets))
      stop("`sets` must map every fitted item when `by = \"set\"`",
           call. = FALSE)
    if (is.list(sets)) {
      if (!length(sets) || is.null(names(sets)) ||
          any(is.na(names(sets))) || any(!nzchar(names(sets))) ||
          anyDuplicated(names(sets)) ||
          !all(vapply(sets, function(z)
            (is.character(z) || is.factor(z)) && length(z), logical(1))))
        stop("a set list needs unique non-empty names and item-name elements",
             call. = FALSE)
      set_items <- unlist(lapply(sets, as.character), use.names = FALSE)
      set_names <- rep(names(sets), lengths(sets))
      if (anyDuplicated(set_items))
        stop("an item cannot belong to more than one externally weighted set",
             call. = FALSE)
      set_of <- stats::setNames(set_names, set_items)
    } else {
      if (!(is.character(sets) || is.factor(sets)) || !length(sets) ||
          is.null(names(sets)) || any(is.na(names(sets))) ||
          any(!nzchar(names(sets))) || anyDuplicated(names(sets)) ||
          anyNA(sets) || any(!nzchar(as.character(sets))))
        stop("`sets` must be a named item-to-set vector or a named list",
             call. = FALSE)
      set_of <- stats::setNames(as.character(sets), names(sets))
    }
    missing_set <- setdiff(items, names(set_of))
    extra_set <- setdiff(names(set_of), items)
    if (length(missing_set) || length(extra_set))
      stop("the set map must name every fitted item exactly",
           if (length(missing_set)) paste0("; missing: ", paste(missing_set, collapse = ", ")) else "",
           if (length(extra_set)) paste0("; unknown: ", paste(extra_set, collapse = ", ")) else "",
           call. = FALSE)
    set_names <- unique(unname(set_of[items]))
    missing_w <- setdiff(set_names, names(weights))
    extra_w <- setdiff(names(weights), set_names)
    if (length(missing_w) || length(extra_w))
      stop("set weights must name every fitted set exactly",
           if (length(missing_w)) paste0("; missing: ", paste(missing_w, collapse = ", ")) else "",
           if (length(extra_w)) paste0("; unknown: ", paste(extra_w, collapse = ", ")) else "",
           call. = FALSE)
    q <- unname(weights[set_of[source_item]])
  }
  q <- q / mean(q)

  disc <- fit$disc
  if (is.null(disc)) disc <- rep(1, ncol(X))
  if (length(disc) == 1L) disc <- rep(disc, ncol(X))
  if (!is.numeric(disc) || length(disc) != ncol(X) ||
      any(!is.finite(disc)) || any(disc <= 0))
    stop("the fitted model units are unavailable", call. = FALSE)

  obs <- !is.na(X) & rep(q > 0, each = nrow(X))
  m <- vapply(fit$tau_list, length, integer(1))
  theta <- se <- rep(NA_real_, nrow(X))
  weighted_score <- rowSums(sweep(X, 2L, q * disc, `*`), na.rm = TRUE)
  weighted_score[rowSums(obs) == 0L] <- NA_real_
  max_weighted_score <- as.numeric(obs %*% (q * disc * m))
  raw <- rowSums(replace(X, !obs, NA_real_), na.rm = TRUE)
  raw[rowSums(obs) == 0L] <- NA_real_
  max_raw <- as.numeric(obs %*% m)
  extreme <- !is.na(weighted_score) &
    (weighted_score <= 1e-12 |
       weighted_score >= max_weighted_score - 1e-12)
  pat <- apply(obs, 1L, function(z) paste(which(z), collapse = ","))

  for (key in unique(pat)) {
    cols <- as.integer(strsplit(key, ",", fixed = TRUE)[[1L]])
    if (!length(cols)) next
    who_pat <- which(pat == key)
    for (Wu in unique(signif(weighted_score[who_pat], 12L))) {
      who <- who_pat[signif(weighted_score[who_pat], 12L) == Wu]
      score <- function(th) {
        mo <- lapply(cols, function(j)
          item_moments(th, fit$tau_list[[j]], disc = disc[j]))
        E <- vapply(mo, `[[`, 0, "E")
        V <- vapply(mo, `[[`, 0, "V")
        m3 <- vapply(mo, `[[`, 0, "mu3")
        H <- sum(q[cols] * disc[cols]^2 * V)
        J <- sum(q[cols]^2 * disc[cols]^2 * V)
        sum(q[cols] * disc[cols] * (X[who[1L], cols] - E)) +
          J * sum(q[cols] * disc[cols]^3 * m3) / (2 * H^2)
      }
      root <- tryCatch(stats::uniroot(score, c(-30, 30), tol = 1e-9)$root,
                       error = function(e) NA_real_)
      theta[who] <- root
      if (is.finite(root)) {
        V <- vapply(cols, function(j)
          item_moments(root, fit$tau_list[[j]], disc = disc[j])$V, 0)
        H <- sum(q[cols] * disc[cols]^2 * V)
        J <- sum(q[cols]^2 * disc[cols]^2 * V)
        se[who] <- sqrt(J) / H
      }
    }
  }

  identity_cols <- intersect(unique(c("id", names(fit$factors))),
                             names(fit$person))
  out <- cbind(fit$person[, identity_cols, drop = FALSE],
               data.frame(n_items = rowSums(obs), raw = raw,
                          max_raw = max_raw,
                          weighted_score = weighted_score,
                          max_weighted_score = max_weighted_score,
                          theta = theta, se = se, extreme = extreme,
                          check.names = FALSE))
  weighting <- data.frame(response_cell = cells, item = source_item,
                          set = if (is.null(set_of)) NA_character_
                            else unname(set_of[source_item]),
                          supplied_weight = if (by == "item")
                            unname(weights[source_item])
                          else unname(weights[set_of[source_item]]),
                          normalised_weight = q,
                          stringsAsFactors = FALSE)
  attr(out, "weighting") <- weighting
  attr(out, "by") <- by
  out
}

# Person locations for an arbitrary response matrix, grouped by
# missing-data pattern so each WLE score table is solved once. disc is a
# single common discrimination (frame unit); raw scores stay sufficient.
.person_estimates <- function(X, tau_list, disc = 1) {
  N <- nrow(X)
  obs <- !is.na(X)
  m <- vapply(tau_list, length, 1L)
  pat <- apply(obs, 1, function(z) paste(which(z), collapse = ","))
  theta <- se <- rep(NA_real_, N)
  raw <- rowSums(X, na.rm = TRUE); raw[rowSums(obs) == 0L] <- NA
  max_raw <- as.numeric(obs %*% m)
  for (key in unique(pat)) {
    cols <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    if (!length(cols)) next
    sel <- which(pat == key)
    pe <- person_wle(tau_list[cols], disc = disc)
    r <- rowSums(X[sel, cols, drop = FALSE])
    theta[sel] <- pe$theta[as.character(r)]
    se[sel]    <- pe$se[as.character(r)]
  }
  data.frame(n_items = rowSums(obs), raw = raw, max_raw = max_raw,
             theta = theta, se = se,
             extreme = !is.na(raw) & (raw == 0L | raw == max_raw))
}

# Model moments evaluated at each person's location, observed cells only.
# Persons sharing a location (same pattern and raw score) share a row of the
# unique-theta moment tables. disc may be a per-column vector (frame units).
.moment_arrays <- function(theta, tau_list, disc = NULL) {
  L <- length(tau_list)
  if (is.null(disc)) disc <- rep(1, L)
  if (length(disc) == 1L) disc <- rep(disc, L)
  ut <- sort(unique(theta[!is.na(theta)]))
  E <- V <- M3 <- M4 <- matrix(NA_real_, length(ut), L)
  for (u in seq_along(ut)) {
    mo <- lapply(seq_len(L), function(i)
      item_moments(ut[u], tau_list[[i]], disc = disc[i]))
    E[u, ]  <- vapply(mo, `[[`, 0, "E")
    V[u, ]  <- vapply(mo, `[[`, 0, "V")
    M3[u, ] <- vapply(mo, `[[`, 0, "mu3")
    M4[u, ] <- vapply(mo, `[[`, 0, "mu4")
  }
  idx <- match(theta, ut)
  list(E = E[idx, , drop = FALSE], V = V[idx, , drop = FALSE],
       M3 = M3[idx, , drop = FALSE], M4 = M4[idx, , drop = FALSE])
}

#' Raw score to measure conversion table
#'
#' The score-to-logit conversion for complete responders: every possible raw
#' score with its location, standard error, and the frequency and cumulative
#' percentage of complete responders at that score (the complete-data
#' estimates table of Andrich and Marais 2019, ch. 10).
#'
#' Two estimators are available. \code{"wle"} (the default) is Warm's
#' weighted likelihood estimate, finite at the extreme scores. \code{"mle"}
#' is the plain maximum likelihood estimate, infinite at the
#' extremes. \code{extremes = "extrapolated"} replaces the extreme-score
#' entries by the geometric extrapolation described in Andrich and Marais
#' (2019, ch. 10): successive score-to-score
#' differences grow towards the extremes, so the last difference is
#' continued geometrically -- the extrapolated top difference \eqn{d} solves
#' \eqn{b = \sqrt{a d}} where \eqn{a, b} are the two preceding differences
#' (equivalently \eqn{d = b^2/a}), and symmetrically at zero. The standard
#' error at an extrapolated location is \eqn{1/\sqrt{I(\theta)}} evaluated
#' there. With \code{method = "wle"} the extrapolation replaces the finite
#' Warm estimates at the extremes, giving the extrapolated form of the
#' conversion table from a WLE analysis.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param method \code{"wle"} (Warm, default) or \code{"mle"}.
#' @param extremes Treatment of the extreme scores. \code{"model"} keeps the
#'   estimator's own values; these are \code{NA} for MLE.
#'   \code{"extrapolated"} applies the geometric extrapolation.
#' @return A data frame with \code{score}, \code{theta}, \code{se},
#'   \code{freq}, \code{cum_pct} (omitted when no complete responders
#'   exist), and \code{extrapolated}; \code{NULL} when the fitted items do not
#'   share one discrimination or an item is represented by several MFRM or
#'   EFRM response cells.
#' @references
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. Psychometrika, 54(3), 427--450.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' score_table(rasch(X), method = "mle", extremes = "extrapolated")
#' @export
score_table <- function(fit, method = c("wle", "mle"),
                        extremes = c("model", "extrapolated")) {
  method <- match.arg(method); extremes <- match.arg(extremes)
  if (is.null(fit$score_table)) return(NULL)
  tab <- fit$score_table[, c("score", "theta", "se")]
  M <- max(tab$score)
  disc <- if (is.null(fit$disc)) 1 else fit$disc[1]
  info <- function(th) sum(vapply(fit$tau_list, function(tt)
    disc^2 * item_moments(th, tt, disc = disc)$V, 0))
  if (method == "mle") {
    for (r in seq_len(M - 1)) {
      # the ML equation is sum(disc*(x_i - E_i)) = 0: the common
      # discrimination cancels, so the root solves r = sum(E), not
      # r = sum(disc*E)
      tab$theta[tab$score == r] <- uniroot(function(th)
        r - sum(vapply(fit$tau_list, function(tt)
          item_moments(th, tt, disc = disc)$E, 0)),
        c(-30, 30), tol = 1e-9)$root
    }
    tab$theta[c(1, M + 1)] <- NA_real_
    tab$se <- 1 / sqrt(vapply(tab$theta, function(th)
      if (is.na(th)) NA_real_ else info(th), 0))
  }
  tab$extrapolated <- FALSE
  if (extremes == "extrapolated") {
    if (M < 4) stop("extrapolation needs at least three interior scores")
    th <- tab$theta
    lo <- th[2] - (th[3] - th[2])^2 / (th[4] - th[3])
    hi <- th[M] + (th[M] - th[M - 1])^2 / (th[M - 1] - th[M - 2])
    tab$theta[c(1, M + 1)] <- c(lo, hi)
    tab$se[c(1, M + 1)] <- 1 / sqrt(c(info(lo), info(hi)))
    tab$extrapolated[c(1, M + 1)] <- TRUE
  }
  raw <- rowSums(fit$X)
  freq <- as.integer(table(factor(raw[stats::complete.cases(fit$X)],
                                  levels = 0:M)))
  if (sum(freq) > 0) {
    tab$freq <- freq
    tab$cum_pct <- 100 * cumsum(freq) / sum(freq)
  }
  tab
}

#' Person measures with extrapolated extreme scores
#'
#' Extreme persons (zero or maximum raw score on their observed items) are
#' excluded from calibration, but they cannot be left out of group
#' comparisons; Andrich and Marais (2019, ch. 10) therefore describe an
#' extrapolated measure for them,
#' continuing the growth of the score-to-score differences so the last
#' difference is the geometric mean of its neighbours (see
#' \code{\link{score_table}}). This helper applies the same rule to the
#' person table: for each missing-data pattern with extreme persons, the
#' score-to-measure conversion over that pattern's items is extrapolated at
#' its ends, and the extreme persons receive the extrapolated location with
#' the standard error \eqn{1/\sqrt{I(\theta)}} evaluated there. Non-extreme
#' persons keep their estimates unchanged. The extrapolation continues the
#' Warm (weighted likelihood) conversion, matching the package's person
#' estimates.
#'
#' @param fit A fitted object from \code{\link{rasch}} (equal
#'   discriminations; not EFRM).
#' @return The fit's person table with two added columns,
#'   \code{theta_extrapolated} and \code{se_extrapolated}: equal to
#'   \code{theta} and \code{se} for non-extreme persons, extrapolated for
#'   extreme persons. Patterns with fewer than three interior scores cannot
#'   be extrapolated and keep their Warm values.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(300 * 8, 1, plogis(outer(rnorm(300, 0, 2), d, "-"))), 300, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(X)
#' pe <- person_extrapolated(fit)
#' head(pe[pe$extreme, c("theta", "theta_extrapolated", "se", "se_extrapolated")])
#' @export
person_extrapolated <- function(fit) {
  if (!is.null(fit$disc) && length(unique(fit$disc)) > 1L)
    stop("extrapolation over a common raw score needs equal discriminations; ",
         "EFRM fits report per-group score curves instead")
  p <- fit$person
  p$theta_extrapolated <- p$theta
  p$se_extrapolated <- p$se
  ext <- which(p$extreme)
  if (!length(ext)) return(p)
  X <- fit$X
  pat <- apply(!is.na(X), 1, paste, collapse = "")
  for (pt in unique(pat[ext])) {
    rows <- intersect(which(pat == pt), ext)
    obs <- which(!is.na(X[rows[1], ]))
    if (length(obs) < 2) next
    tl <- fit$tau_list[obs]
    pe <- person_wle(tl)
    th <- unname(pe$theta); M <- length(th) - 1L
    if (M < 4) next                        # too few interior scores
    lo <- th[2] - (th[3] - th[2])^2 / (th[4] - th[3])
    hi <- th[M] + (th[M] - th[M - 1])^2 / (th[M - 1] - th[M - 2])
    info <- function(t) sum(vapply(tl, function(tt)
      item_moments(t, tt)$V, 0))
    for (r in rows) {
      raw <- sum(X[r, obs])
      if (raw == 0L) {
        p$theta_extrapolated[r] <- lo
        p$se_extrapolated[r] <- 1 / sqrt(info(lo))
      } else {
        p$theta_extrapolated[r] <- hi
        p$se_extrapolated[r] <- 1 / sqrt(info(hi))
      }
    }
  }
  p
}
