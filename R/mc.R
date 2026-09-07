# rasch :: multiple choice
# ===========================================================================
# Scoring of raw multiple-choice responses against a key, with the raw
# responses retained for distractor analysis: per option, the count,
# proportion, mean person location, and point-biserial with the person
# measure. A distractor chosen by abler persons than the keyed option is the
# classic signature of a miskeyed item and is flagged. Option curves show
# the proportion choosing each option across trait class intervals.
# Three key forms are supported:
# a single correct option per item (scored 0/1), double keying (several
# options separated by "/" all scoring 1), and full polytomous option
# scoring (a data frame item/option/score assigning partial credit to
# informative distractors; Andrich & Styles 2011). distractor_rescore()
# proposes such a scoring from the rest-measure evidence.
# ===========================================================================

# Resolve a key into a per-item scoring map (named integer vector
# option -> score). Accepted forms: a named vector or data frame
# (item, key) where key is the correct option, possibly several separated
# by "/" (double keying, all scoring 1); or a data frame
# (item, option, score) assigning an integer score to every credited
# option (unlisted options score 0). Matching is case-insensitive.
.resolve_key <- function(key) {
  if (is.data.frame(key)) .check_column_names(key)
  if (is.data.frame(key) && all(c("item", "option", "score") %in% names(key))) {
    # Item values are column selectors. Preserve an exact, deliberately
    # spaced column name; only response options are canonicalised for the
    # case-insensitive matching used by .score_mc().
    item_missing <- is.na(key$item)
    key$item <- as.character(key$item)
    key$item[item_missing | is.na(key$item)] <- NA_character_
    # a row with no item names no item: split() would file it under a
    # phantom group and the real item would be left unscored
    blank_item <- is.na(key$item) | !nzchar(trimws(key$item))
    if (any(blank_item))
      stop("missing or blank item name in the scoring table (row(s) ",
           paste(which(blank_item), collapse = ", "), ")")
    key$option <- toupper(.role_text_values(key$option))
    # an option nobody can have chosen scores its item zero for everyone,
    # and the item is then dropped as constant under a misleading message
    blank_opt <- is.na(key$option) | !nzchar(key$option)
    if (any(blank_opt))
      stop("missing or blank option in the scoring table for item(s): ",
           paste(unique(as.character(key$item)[blank_opt]), collapse = ", "))
    sc <- suppressWarnings(as.numeric(as.character(key$score)))
    if (anyNA(sc) || any(!is.finite(sc)) || any(sc != floor(sc)) ||
        any(sc < 0) || any(sc > .Machine$integer.max))
      stop("option scores must be non-negative integers")
    key$score <- as.integer(sc)
    map <- lapply(split(key, key$item), function(d) {
      if (anyDuplicated(d$option))
        stop("duplicate option row(s) for item ", d$item[1])
      if (max(d$score) < 1)
        stop("item ", d$item[1], " credits no option")
      setNames(d$score, d$option)
    })
    return(map)
  }
  if (is.data.frame(key)) {
    if (!all(c("item", "key") %in% names(key)))
      stop("a key data frame needs columns item, key ",
           "(or item, option, score for polytomous option scoring)")
    item_missing <- is.na(key$item)
    item_names <- as.character(key$item)
    item_names[item_missing | is.na(item_names)] <- NA_character_
    key <- setNames(.role_text_values(key$key), item_names)
  }
  if (is.null(names(key))) stop("the key must be named by item")
  # a key entry that names no item scores no item: the real item would be
  # left unscored and dropped from the analysis without a word
  blank_nm <- is.na(names(key)) | !nzchar(trimws(names(key)))
  if (any(blank_nm))
    stop("missing or blank item name in the key (entr(ies) ",
         paste(which(blank_nm), collapse = ", "), ")")
  if (anyDuplicated(names(key)))
    stop("duplicate key entr(ies) for item(s): ",
         paste(unique(names(key)[duplicated(names(key))]), collapse = ", "),
         "; give one key per item (use '/' for double keying)")
  # a missing key value cannot score its item: without this it would
  # zero-score every response and the item would be dropped as constant
  # under a misleading message
  if (anyNA(key))
    stop("missing (NA) key value for item(s): ",
         paste(names(key)[is.na(key)], collapse = ", "))
  blank_key <- !nzchar(trimws(as.character(key)))
  if (any(blank_key))
    stop("blank key value for item(s): ",
         paste(names(key)[blank_key], collapse = ", "))
  lapply(setNames(trimws(toupper(as.character(key))), names(key)),
         function(k) {
           opts <- trimws(strsplit(k, "/", fixed = TRUE)[[1]])
           if (startsWith(k, "/") || endsWith(k, "/") ||
               any(!nzchar(opts)))
             stop("empty credited option within an item key")
           if (anyDuplicated(opts))
             stop("each credited option must be named once within an item key")
           setNames(rep(1L, length(opts)), opts)
         })
}

# One-line display form of an item's scoring map (e.g. "C", "A/C", or
# "C=2, B=1").
.key_label <- function(m) {
  m <- sort(m[m > 0], decreasing = TRUE)
  if (all(m == 1L)) paste(names(m), collapse = "/")
  else paste(sprintf("%s=%d", names(m), m), collapse = ", ")
}

# Score raw responses against the scoring maps. Observed options absent
# from an item's map score 0; blank, NA, and missing-data codes become NA.
.score_mc <- function(X, map, na_codes = -1) {
  .check_na_codes(na_codes)
  supplied <- names(map)
  idx <- match(supplied, colnames(X))
  fallback <- is.na(idx)
  if (any(fallback))
    idx[fallback] <- match(.role_text_values(supplied[fallback]), colnames(X))
  if (all(is.na(idx))) stop("no key item matches an item column")
  # A key entry naming no item column is almost always a typo: surface it
  # rather than drop it silently. Exact names win before surrounding
  # whitespace is treated as selector syntax.
  if (anyNA(idx))
    stop("key item(s) with no matching data column: ",
         paste(supplied[is.na(idx)], collapse = ", "), call. = FALSE)
  keyed <- colnames(X)[idx]
  if (anyDuplicated(keyed))
    stop("key entries resolve to the same item column: ",
         paste(unique(keyed[duplicated(keyed)]), collapse = ", "),
         call. = FALSE)
  names(map) <- keyed
  raw <- matrix(trimws(toupper(as.character(X[, keyed]))), nrow(X),
                length(keyed), dimnames = list(NULL, keyed))
  codes <- trimws(toupper(as.character(na_codes)))
  codes <- codes[!is.na(codes)]
  raw_num <- suppressWarnings(as.numeric(raw))
  code_num <- suppressWarnings(as.numeric(codes))
  missing <- is.na(raw) | raw %in% c("", "NA") | raw %in% codes |
    (!is.na(raw_num) & raw_num < 0)
  if (any(is.finite(code_num)))
    missing <- missing | (!is.na(raw_num) & raw_num %in% code_num[is.finite(code_num)])
  raw[missing] <- NA_character_
  scored <- matrix(NA_integer_, nrow(raw), ncol(raw),
                   dimnames = dimnames(raw))
  for (j in seq_along(keyed)) {
    m <- map[[keyed[j]]]
    s <- unname(m[raw[, j]])
    s[is.na(s) & !is.na(raw[, j])] <- 0L   # observed but uncredited option
    scored[, j] <- s
  }
  list(scored = scored, raw = raw, map = map[keyed],
       key = vapply(map[keyed], .key_label, ""))
}

#' Distractor analysis for multiple-choice items
#'
#' For every keyed item and response option: the count and proportion
#' choosing it among respondents with a non-extreme rest measure, their mean
#' location, and the point-biserial correlation between choosing the option
#' and the person measure. These summaries use the rest measure (the person
#' estimate from the other items), so the analysed item cannot credit its own
#' takers. The keyed
#' option should attract able persons and usually carry a positive
#' point-biserial; a distractor whose takers are abler than the pooled takers
#' of the full-credit option or options (with at least \code{min_n} takers) is
#' flagged as a possible
#' miskey. The analysis requires one response row per person; repeated rows
#' do not supply independent taker counts or rest measures.
#'
#' @param fit A fitted object from \code{\link{rasch}} run with a \code{key}.
#' @param items Optional subset of item names; defaults to every keyed item.
#' @param min_n Minimum takers for an option to be eligible for the miskey
#'   flag.
#' @return A data frame with one row per item-option: \code{item},
#'   \code{option}, its assigned \code{score}, \code{keyed} (full credit),
#'   \code{n}, \code{prop}, \code{mean_location}, \code{point_biserial},
#'   and \code{flag}.
#' @examples
#' set.seed(1); Np <- 400
#' th <- rnorm(Np)
#' raw <- sapply(seq(-1, 1, length.out = 6), function(d) {
#'   ok <- rbinom(Np, 1, plogis(th - d))
#'   ifelse(ok == 1, "A", sample(c("B", "C", "D"), Np, replace = TRUE))
#' })
#' colnames(raw) <- paste0("M", 1:6)
#' fit <- rasch(raw, key = setNames(rep("A", 6), colnames(raw)))
#' head(distractor_analysis(fit))
#' @export
distractor_analysis <- function(fit, items = NULL, min_n = 10) {
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    stop("`fit` must be a fitted response-data Rasch model", call. = FALSE)
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; rest-measure distractor analysis is unavailable",
         call. = FALSE)
  if (!.efrm_link_converged(fit))
    stop("the fitted set-unit link did not converge; rest-measure distractor analysis is unavailable",
         call. = FALSE)
  if (.has_repeated_residual_units(fit))
    .refuse("distractor analysis needs one response row per person; repeated ",
            "identifiers would inflate option counts and treat dependent ",
            "rest measures as independent")
  if (is.null(fit$mc)) stop("the fit has no key: run rasch(..., key = )")
  min_n <- .check_whole(min_n, "min_n", 1)
  raw <- fit$mc$raw; map <- fit$mc$map
  if (is.null(items)) items <- colnames(raw)
  if (!is.atomic(items) || !is.null(dim(items)) || !length(items))
    stop("`items` must be an ordinary vector naming at least one item")
  supplied <- as.character(items)
  if (anyNA(supplied) || any(!nzchar(trimws(supplied))))
    stop("`items` must contain non-missing, non-empty item names")
  idx <- match(supplied, colnames(raw))
  fallback <- is.na(idx)
  if (any(fallback))
    idx[fallback] <- match(.role_text_values(supplied[fallback]), colnames(raw))
  if (anyNA(idx))
    stop("item(s) without raw multiple-choice responses: ",
         paste(supplied[is.na(idx)], collapse = ", "))
  items <- colnames(raw)[idx]
  if (anyDuplicated(items))
    stop("`items` must not name the same keyed item more than once")
  out <- list()
  for (it in items) {
    r <- raw[, it]
    idx <- match(it, colnames(fit$X))
    rp <- .person_estimates(fit$X[, -idx, drop = FALSE],
                            fit$tau_list[-idx])         # rest measure
    th <- rp$theta
    ok <- !is.na(r) & is.finite(th) & !rp$extreme
    # Retain every observed option even when all of its takers have an
    # extreme or unavailable rest score. Its rest-measure summaries are then
    # honestly unavailable rather than the option silently disappearing.
    opts <- sort(unique(r[!is.na(r)]))
    m <- map[[it]]
    sc <- unname(m[opts]); sc[is.na(sc)] <- 0L
    rows <- data.frame(item = it, option = opts, score = sc,
                       keyed = sc == max(m),
                       n = NA_integer_, prop = NA_real_,
                       mean_location = NA_real_, point_biserial = NA_real_)
    for (j in seq_along(opts)) {
      sel <- ok & r == opts[j]
      rows$n[j] <- sum(sel)
      rows$prop[j] <- if (any(ok)) sum(sel) / sum(ok) else NA_real_
      rows$mean_location[j] <- if (any(sel)) mean(th[sel]) else NA_real_
      ind <- as.integer(r[ok] == opts[j])
      # var() is NA for a single usable rest measure.  That is insufficient
      # to define a correlation, but it is not an analysis fault.
      rows$point_biserial[j] <- if (length(ind) > 1L &&
                                      isTRUE(stats::var(ind) > 0) &&
                                      isTRUE(stats::var(th[ok]) > 0))
        cor(ind, th[ok]) else NA_real_
    }
    # Several options can share full credit. They form one scored response
    # category, so compare a distractor with their pooled takers rather than
    # selecting whichever keyed option happened to have the highest sample
    # mean.
    keyed_options <- rows$option[rows$keyed]
    key_values <- th[ok & r %in% keyed_options]
    key_mean <- if (length(key_values)) mean(key_values) else NA_real_
    rows$flag <- if (is.finite(key_mean))
      !rows$keyed & rows$n >= min_n & rows$mean_location > key_mean
    else rep(FALSE, nrow(rows))
    out[[it]] <- rows
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

#' Plot multiple-choice option curves
#'
#' The proportion choosing each response option across class intervals of the
#' rest measure (the person estimate from the other items), with the keyed
#' option drawn solid and bold. The curves are descriptive. Under ordered
#' option scoring, higher-scored options should tend to occur at higher rest
#' measures; an intermediate-credit option may peak in the middle.
#'
#' @param fit A fitted object from \code{\link{rasch}} run with a \code{key}.
#' @param item Keyed item name.
#' @param n_groups Number of class intervals. By default, use the fit's
#'   count, with at least two requested intervals. Tied locations may
#'   produce fewer intervals.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1); Np <- 400
#' th <- rnorm(Np)
#' raw <- sapply(seq(-1, 1, length.out = 6), function(d) {
#'   ok <- rbinom(Np, 1, plogis(th - d))
#'   ifelse(ok == 1, "A", sample(c("B", "C", "D"), Np, replace = TRUE))
#' })
#' colnames(raw) <- paste0("M", 1:6)
#' fit <- rasch(raw, key = setNames(rep("A", 6), colnames(raw)))
#' plot_distractors(fit, "M3")
#' @export
plot_distractors <- function(fit, item, n_groups = NULL) {
  .check_response_display_fit(fit, "rest-measure distractor plots")
  if (.has_repeated_residual_units(fit))
    .refuse("distractor plots need one response row per person; repeated ",
            "identifiers would give dependent rows equal weight")
  if (!is.atomic(item) || !is.null(dim(item)) ||
      length(item) != 1L || is.na(item))
    stop("`item` must name exactly one item")
  item <- as.character(item)
  n_groups <- if (is.null(n_groups)) max(2L, fit$n_groups)
              else .check_whole(n_groups, "n_groups", 2)
  if (is.null(fit$mc)) stop("the fit has no key: run rasch(..., key = )")
  if (!item %in% colnames(fit$mc$raw)) stop("no such keyed item: ", item)
  r <- fit$mc$raw[, item]
  idx <- match(item, colnames(fit$X))
  rp <- .person_estimates(fit$X[, -idx, drop = FALSE],
                          fit$tau_list[-idx])          # rest measure
  th <- rp$theta
  candidate <- !is.na(r) & !is.na(th) & !rp$extreme
  if (sum(candidate) < 4L)
    stop("fewer than 4 non-extreme rest measures are available for this item")
  ng <- min(n_groups, max(2, floor(sum(candidate) / 25)))
  ci <- .class_intervals(ifelse(is.na(r), NA_real_, th), rp$extreme, ng)
  ok <- !is.na(ci)
  mid <- tapply(th[ok], ci[ok], mean)
  opts <- sort(unique(r[ok]))
  m <- fit$mc$map[[item]]
  sc <- unname(m[opts]); sc[is.na(sc)] <- 0L
  keyed_v <- sc == max(m)
  op <- .rr_canvas(range(mid) + c(-0.2, 0.2), c(0, 1),
                   "Person location (logits)", "Proportion choosing option",
                   paste0(item, "  (key: ", fit$mc$key[item], ")"))
  on.exit(par(op))
  for (j in seq_along(opts)) {
    pr <- tapply(r[ok] == opts[j], ci[ok], mean)
    colr <- .rr$pal[(j - 1L) %% length(.rr$pal) + 1L]
    lines(mid, pr, lwd = if (keyed_v[j]) 3.2 else if (sc[j] > 0) 2.4 else 1.8,
          lty = if (keyed_v[j]) 1 else if (sc[j] > 0) 2 else 5, col = colr)
    points(mid, pr, pch = 21, bg = colr, col = "white",
           cex = if (keyed_v[j]) 1.5 else 1.1, lwd = 1)
  }
  labs <- paste0(opts,
                 ifelse(keyed_v, " (key)",
                        ifelse(sc > 0, sprintf(" (credit %d)", sc), "")))
  .rr_legend("right", labs,
             lwd = ifelse(keyed_v, 3.2, ifelse(sc > 0, 2.4, 1.8)),
             lty = ifelse(keyed_v, 1, ifelse(sc > 0, 2, 5)),
             col = .rr$pal[(seq_along(opts) - 1L) %% length(.rr$pal) + 1L])
  invisible(NULL)
}

#' Propose polytomous option scores from the distractor evidence
#'
#' Multiple-choice items can be rescored polytomously so
#' that a distractor carrying information about the trait receives partial
#' credit (Andrich and Styles 2011). This function proposes such a scoring
#' from the rest-measure distractor analysis: within each keyed item, a
#' distractor qualifies for credit when it attracts at least \code{min_n}
#' takers, its takers' mean rest location exceeds that of the pooled remaining
#' distractors by more than \code{z} standard errors of the difference,
#' and it remains below the keyed option. Qualifying distractors are
#' ranked by mean location and scored 1, 2, ... below the keyed option's
#' top score. The result is a proposal for substantive review, not an
#' automatic decision: inspect \code{\link{plot_distractors}} and the item
#' content, edit as needed, then refit with
#' \code{rasch(raw_data, key = proposal$option_scores)}.
#'
#' @param fit A fitted object from \code{\link{rasch}} run with a \code{key}.
#' @param items Optional subset of keyed item names.
#' @param min_n Minimum takers for a distractor to be considered.
#' @param z Required separation, in standard errors, between a credited
#'   distractor and the pooled remaining distractors.
#' @return A list of class \code{"rasch_rescore"}: \code{option_scores}, a
#'   data frame (\code{item}, \code{option}, \code{score}) ready for
#'   \code{rasch(key = )}, covering every observed option of the examined
#'   items and retaining the existing scoring of other keyed items; and
#'   \code{evidence}, the distractor analysis with the proposed scores and
#'   the separation z per option.
#' @references Andrich, D. and Styles, I. (2011). Distractors with
#'   information in multiple choice items: A rationale based on the Rasch
#'   model. Journal of Applied Measurement, 12, 67-95.
#' @examples
#' set.seed(1); Np <- 600
#' th <- rnorm(Np)
#' raw <- sapply(seq(-0.5, 0.5, length.out = 4), function(d) {
#'   x <- vapply(th, function(b) sample(0:2, 1,
#'     prob = item_moments(b, c(d - 0.7, d + 0.7))$P), 0L)
#'   c("D", "B", "A")[x + 1]   # B is an informative distractor
#' })
#' colnames(raw) <- paste0("M", 1:4)
#' fit <- rasch(raw, key = setNames(rep("A", 4), colnames(raw)))
#' pr <- distractor_rescore(fit)
#' pr$option_scores
#' @export
distractor_rescore <- function(fit, items = NULL, min_n = 20, z = 1.96) {
  min_n <- .check_whole(min_n, "min_n", 1)
  if (length(z) != 1L || !is.numeric(z) || is.complex(z) ||
      !is.null(dim(z)) || !is.null(oldClass(z)) || !is.finite(z) || z <= 0)
    stop("`z` must be one positive finite separation threshold")
  if (!is.null(items) && !length(items))
    stop("`items` must name at least one keyed item")
  da <- distractor_analysis(fit, items = items, min_n = min_n)
  raw <- fit$mc$raw
  ev <- list(); os <- list()
  for (it in unique(da$item)) {
    d <- da[da$item == it, ]
    if (!any(d$n > 0L))
      stop("item ", it, " has no non-extreme rest measures; a rescoring ",
           "proposal cannot be estimated", call. = FALSE)
    usable_key <- which(d$keyed %in% TRUE & d$n > 0L &
                          is.finite(d$mean_location))
    if (!length(usable_key))
      stop("item ", it, " has no observed full-credit option with a ",
           "non-extreme rest measure; a rescoring proposal cannot locate ",
           "the keyed response on the rest-measure scale", call. = FALSE)
    idx <- match(it, colnames(fit$X))
    rp <- .person_estimates(fit$X[, -idx, drop = FALSE],
                            fit$tau_list[-idx])
    th <- rp$theta
    r <- raw[, it]
    ok <- !is.na(r) & is.finite(th) & !rp$extreme
    # per-option SE of the mean rest location
    d$se_location <- vapply(d$option, function(o) {
      x <- th[ok & r == o]
      if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
    }, 0)
    key_options <- d$option[d$keyed]
    key_values <- th[ok & r %in% key_options]
    key_mean <- if (length(key_values)) mean(key_values) else NA_real_
    cand <- which(!d$keyed & d$n >= min_n)
    # the uncredited baseline: distractors not currently under consideration
    d$z_sep <- NA_real_
    credited <- integer(0)
    for (j in cand) {
      # min_n governs whether the nominated option is stable enough to
      # receive credit. The comparison group is the documented pool of all
      # other non-keyed options; requiring each of those options separately
      # to meet min_n can leave a single eligible candidate with no baseline.
      others <- setdiff(which(!d$keyed), j)
      if (!length(others)) next
      base <- th[ok & r %in% d$option[others]]
      xj <- th[ok & r == d$option[j]]
      if (length(base) < 2 || length(xj) < 2) next
      den <- sqrt(var(xj) / length(xj) + var(base) / length(base))
      # zero spread in both groups (every chooser at the same location)
      # gives a 0/0 separation z: leave it NA rather than NaN
      if (!is.finite(den) || den <= 0) next
      sep <- (mean(xj) - mean(base)) / den
      d$z_sep[j] <- sep
      if (sep > z && mean(xj) < key_mean) credited <- c(credited, j)
    }
    credited <- credited[order(d$mean_location[credited])]
    d$proposed <- 0L
    if (length(credited)) d$proposed[credited] <- seq_along(credited)
    d$proposed[d$keyed] <- length(credited) + 1L
    ev[[it]] <- d
    os[[it]] <- data.frame(item = it, option = d$option, score = d$proposed)
  }
  # A subset proposal must still be usable with the original raw matrix.
  # Leaving the other character-valued item columns out of the key would
  # leave them unscored and make the documented refit fail. Preserve their
  # existing maps unchanged.
  untouched <- setdiff(colnames(raw), names(os))
  for (it in untouched) {
    m <- fit$mc$map[[it]]
    os[[it]] <- data.frame(item = it, option = names(m),
                           score = as.integer(m))
  }
  item_order <- intersect(colnames(raw), names(os))
  out <- list(option_scores = do.call(rbind, os[item_order]),
              evidence = do.call(rbind, ev))
  rownames(out$option_scores) <- rownames(out$evidence) <- NULL
  out <- .tag_tables(out)
  class(out) <- "rasch_rescore"
  out
}

#' @export
print.rasch_rescore <- function(x, ...) {
  n_credit <- sum(x$option_scores$score > 0 &
                  x$option_scores$score < ave(x$option_scores$score,
                                              x$option_scores$item,
                                              FUN = max))
  cat(sprintf("Polytomous option-scoring proposal (Andrich & Styles 2011): %d distractor(s) credited\n",
              n_credit))
  ev <- x$evidence
  num <- vapply(ev, is.numeric, TRUE)
  ev[num] <- lapply(ev[num], round, 3)
  print(ev, row.names = FALSE)
  cat("review substantively, edit, then refit: rasch(raw, key = proposal$option_scores)\n")
  invisible(x)
}
