# rasch :: a null distribution for the item fit statistics
# ===========================================================================
# Every item fit statistic in this package is computed at ESTIMATED person
# locations, and referred to a distribution derived as though those locations
# were known. Each pays for that differently, and all of them pay more as the
# sample grows.
#
# The class-interval chi-square forms its intervals by selecting on the
# estimates, so within an interval the true locations regress toward the mean
# relative to the estimates that selected them and the expected interval mean
# is biased. The bias does not shrink as the sample grows while the interval
# sizes do, so the statistic grows about linearly in N against a fixed
# reference: at 8 items a correctly fitting item is rejected 9% of the time
# at 250 persons and 100% at 4,000, and the whole-test total is rejecting
# every correctly fitting instrument by 1,000.
#
# The fit residual and the Wilson-Hilferty standardisations fail in a
# different direction. Under a true model the fit residual's null SD runs
# 0.71 at 250 persons to 1.00 at 4,000 with a mean drifting from -0.12 to
# -0.37, so the conventional +/-2.5 cut is far too lenient at small samples
# and only reaches its nominal meaning at a few thousand persons. The
# standardised statistics are worse: infit z beyond 1.96 flags 12% of
# correctly fitting items at 250 persons and 69% at 4,000. Only the mean
# squares themselves are stable, sitting at 1.06 at every sample size -- it
# is the standardisation that drifts, because the spread it divides by
# shrinks as N grows.
#
# Rescaling the chi-square to a reference sample size, the remedy RUMM2030
# offers and this package carried until 1.12.1, is one constant applied
# across items: it can move how many items are flagged but not which, and at
# a reference of 500 it silenced a planted 2.5-slope item in every one of 120
# replicates. Analysing a subsample instead inherits whatever miscalibration
# holds at the subsample's size. Neither addresses the reference
# distribution, which is the thing that is wrong.
#
# The parametric bootstrap replaces it, once, for all of them. Data are
# generated from the fitted model and each replicate travels the same road as
# the observed data -- estimate the item side, estimate person locations from
# it, form class intervals on those estimates, take residuals against those
# same estimates -- so the same bias enters the null as entered the observed
# statistic and cancels in the comparison. One set of replicates serves every
# statistic, because a replicate has to refit the whole model anyway.
#
# What the correction cannot supply is power the statistic never had. A
# flatter-than-Rasch item varies less across the class intervals and so
# carries LESS of the selection bias than a fitting item does: its chi-square
# comes out smaller than a fitting item's, and no reference distribution can
# make a small statistic significant. The fit residual sees that departure
# plainly. They have complementary blind spots, which is the argument for
# calibrating both rather than choosing between them.
# ===========================================================================

# Everything the item fit statistics need from a response matrix and nothing
# else: the person tables, score curves and targeting that rasch() also builds
# are not used here and would multiply the cost by the number of replicates.
# The class-interval count is imposed rather than re-derived per replicate, so
# every replicate computes the same statistic as the observed fit rather than
# one on a neighbouring number of intervals.
.fit_refit <- function(X, model, n_groups, anchors, maxit, tol) {
  est <- tryCatch(pcml(X, model = model, anchors = anchors,
                       maxit = maxit, tol = tol),
                  error = function(e) NULL)
  if (is.null(est)) return(NULL)
  L <- ncol(X)
  tau_list <- lapply(seq_len(L), function(i) est$thr$tau[est$thr$item == i])
  if (!all(vapply(tau_list, function(t) length(t) && all(is.finite(t)), TRUE)))
    return(NULL)
  person <- .person_estimates(X, tau_list, disc = 1)
  if (all(is.na(person$theta))) return(NULL)
  mo <- .moment_arrays(person$theta, tau_list, disc = rep(1, L))
  Z <- (X - mo$E) / sqrt(mo$V)
  colnames(Z) <- colnames(X)
  ci <- .class_intervals(person$theta, person$extreme, n_groups)
  ci_list <- if (anyNA(X))
    .class_intervals_by_item(X, person$theta, person$extreme, n_groups) else NULL
  it <- .item_trait(X, mo, ci, ci_list = ci_list)
  ifit <- .item_fit(X, Z, mo, extreme = person$extreme)
  n_par <- if (is.null(est$n_parameters)) nrow(est$thr) - 1L else est$n_parameters
  rf <- .fitres(Z, mo, person$extreme, n_par)
  list(chisq = it$chisq, fit_resid = rf$items$fit_resid,
       infit_ms = ifit$infit_ms, outfit_ms = ifit$outfit_ms,
       infit_z = ifit$infit_z, outfit_z = ifit$outfit_z)
}

# The statistics this bootstrap calibrates, and how each is read. The
# chi-square is a discrepancy and only its upper tail is misfit; the others
# depart in two directions -- above for a flatter item than the model
# predicts, below for a steeper one -- and both readings are misfit, so their
# p-values are equal-tailed.
.FIT_STATS <- c(chisq = "upper", fit_resid = "two", infit_ms = "two",
                outfit_ms = "two", infit_z = "two", outfit_z = "two")

# Bootstrap p-value for one statistic against its own replicated null.
# (1 + r)/(1 + B) never returns zero; the two-sided form doubles the smaller
# tail and is capped, so its resolution is half the one-sided form's.
.boot_p <- function(obs, null, side) {
  n <- sum(!is.na(null))
  if (!is.finite(obs) || n < 1L) return(NA_real_)
  up <- (1 + sum(null >= obs, na.rm = TRUE)) / (1 + n)
  if (side == "upper") return(up)
  lo <- (1 + sum(null <= obs, na.rm = TRUE)) / (1 + n)
  min(1, 2 * min(up, lo))
}

# One replicate conditional on the observed scores: for each person, a
# response pattern drawn from the Rasch conditional distribution given their
# raw score over their own observed items. Sufficiency cancels the person
# parameter, so no ability has to be drawn at all -- which is what makes
# this the default. Generating at sampled abilities either inflates the
# person variance (resampled estimates carry their error variance, and the
# standardised fit statistics feel it as the sample grows) or breaks the tie
# between who a person is and what they answered (a booklet design's
# missingness is informative, and an independently drawn ability severs it).
# Conditioning keeps each person's score and missingness glued together and
# reproduces the observed score margins exactly.
#
# Sampling is sequential with suffix elementary-symmetric functions: for
# items j = 1..k and remaining score r, P(x_j = x) is proportional to
# w_j(x) G_{j+1}(r - x), where w_j(x) = exp(-cumulative tau_j(x)) and
# G_{j+1} is the ESF of the items after j. Persons sharing a missingness
# pattern share the suffix tables, so booklets cost one recursion each.
.fit_gen_conditional <- function(X, tau_list, na_mask) {
  N <- nrow(X); L <- length(tau_list)
  obs_mask <- if (is.null(na_mask)) matrix(TRUE, N, L) else !na_mask
  # per-item category weights, scaled per item (constants cancel in the
  # conditional distribution; scaling guards the exp against extreme taus)
  W <- lapply(tau_list, function(tau) {
    w <- exp(-c(0, cumsum(tau)))
    w / max(w)
  })
  out <- matrix(NA_integer_, N, L, dimnames = dimnames(X))
  key <- apply(obs_mask, 1, function(z) paste(which(z), collapse = ","))
  for (kk in unique(key)) {
    rows <- which(key == kk)
    items <- as.integer(strsplit(kk, ",", fixed = TRUE)[[1]])
    k <- length(items)
    if (!k) next
    m <- vapply(items, function(j) length(W[[j]]) - 1L, 0L)
    # suffix ESFs: G[[t]][s + 1] carries items t..k at score s, scaled at
    # each step so long tests cannot underflow
    G <- vector("list", k + 1L)
    G[[k + 1L]] <- 1
    for (t in k:1) {
      prev <- G[[t + 1L]]
      cur <- numeric(m[t] + length(prev))
      w <- W[[items[t]]]
      for (x in 0:m[t])
        cur[(x + 1):(x + length(prev))] <-
          cur[(x + 1):(x + length(prev))] + w[x + 1] * prev
      G[[t]] <- cur / max(cur)
    }
    for (i in rows) {
      r <- sum(X[i, items])
      for (t in seq_len(k)) {
        w <- W[[items[t]]]
        xs <- 0:m[t]
        rest <- G[[t + 1L]]
        ok <- (r - xs) >= 0 & (r - xs) <= (length(rest) - 1L)
        pr <- numeric(length(xs))
        pr[ok] <- w[xs[ok] + 1] * rest[r - xs[ok] + 1L]
        x <- if (t == k) r else sample.int(length(xs), 1L, prob = pr) - 1L
        out[i, items[t]] <- x
        r <- r - x
      }
    }
  }
  out
}

# One replicate of the fitted model: responses for the given locations under
# the estimated thresholds, carrying the observed missing-data pattern so
# that items are tested on the same number of responders they were fitted
# on.
.fit_gen <- function(theta, tau_list, na_mask, item_names) {
  X <- vapply(seq_along(tau_list),
              function(i) .sim_item(theta, tau_list[[i]]),
              integer(length(theta)))
  dim(X) <- c(length(theta), length(tau_list))
  colnames(X) <- item_names
  if (!is.null(na_mask)) X[na_mask] <- NA_integer_
  X
}

# Person locations for the non-default, ability-sampling schemes. Each has a
# documented cost, which is why the score-conditional generator above is the
# default. "resample" generates at the spread of the ESTIMATES -- the true
# variance plus their error variance -- and the standardised fit statistics
# feel that inflation as the sample grows (infit z reached 14% at 4,000
# persons); it also draws abilities independently of each person's
# missingness pattern, which severs a linked-booklet design's informative
# missingness. "normal" corrects the variance but leaves the extreme
# categories of easy and hard items unvisited often enough to lose a
# non-random share of replicates (a third, on a four-category rating scale
# at 400 persons). "fixed" keeps each person's own estimate and mask paired,
# at the estimates' inflated spread. All three remain for comparison and for
# designs where conditioning is not wanted.
.fit_theta <- function(scheme, theta, psi) {
  n <- length(theta)
  switch(scheme,
    normal = {
      vt <- psi$var_theta; mse <- psi$mean_error_var
      sd_c <- if (is.finite(vt) && is.finite(mse)) sqrt(max(vt - mse, 0)) else NA_real_
      if (!is.finite(sd_c) || sd_c <= 0) sample(theta, n, replace = TRUE)
      else stats::rnorm(n, mean(theta), sd_c)
    },
    resample = sample(theta, n, replace = TRUE),
    fixed = theta)
}

#' Bootstrap null distribution for the item fit statistics
#'
#' Every item fit statistic \code{\link{rasch}} reports is computed at
#' estimated person locations and referred to a distribution derived as
#' though those locations were known, and each is miscalibrated by an amount
#' that grows with the sample. The item-trait chi-square, whose class
#' intervals are formed by selecting on the estimates, grows about linearly
#' in N against a fixed reference. The fit residual's null SD runs from about
#' 0.71 at 250 persons to 1.00 at 4,000, so the conventional \eqn{\pm}2.5 cut
#' means different things at different sizes; infit z beyond 1.96 flags 12%
#' of correctly fitting items at 250 persons and 69% at 4,000. This function
#' replaces the reference distribution for all of them at once.
#'
#' Each replicate generates responses from the fitted thresholds at person
#' locations drawn under \code{theta}, then re-estimates the item side,
#' re-estimates person locations, re-forms class intervals and re-takes
#' residuals exactly as the observed fit did, so the same bias enters the
#' null as entered the observed statistics. One set of replicates serves
#' every statistic, since a replicate must refit the whole model in any case.
#' The class-interval count is held at the fitted value, and observed missing
#' data are carried into every replicate.
#'
#' The chi-square is read in its upper tail alone. The fit residual, the mean
#' squares and the standardised statistics depart in both directions -- above
#' for an item flatter than the model predicts, below for a steeper one --
#' and both are misfit, so their p-values are equal-tailed. A bootstrap
#' p-value is \code{(1 + r) / (1 + B)}, so it is never zero and its
#' resolution is \code{1 / (1 + B)} one-sided and \code{2 / (1 + B)}
#' two-sided. Familywise flagging multiplies that floor by the item count:
#' Holm across L items cannot reach .05 below \code{B = 20 L - 1} for the
#' chi-square and \code{B = 40 L - 1} for the two-sided statistics, so the
#' default \code{B = 200} resolves adjusted two-sided tests only to five
#' items and serious familywise use wants \code{B} in the hundreds to
#' thousands.
#'
#' Calibration is not power. A flatter-than-Rasch item carries less of the
#' class-interval selection bias than a fitting item does, so its chi-square
#' comes out \emph{smaller} than a fitting item's and no reference
#' distribution can make it significant; the fit residual detects that
#' departure readily. Calibrating both is what covers the two blind spots.
#' A null estimated from data that contain a misfitting item is also mildly
#' contaminated by it, which leaves the remaining items flagging somewhat
#' above nominal.
#'
#' @param fit A fitted object from \code{\link{rasch}}. Extended-frame,
#'   many-facet and explanatory fits are not supported.
#' @param B Number of bootstrap replicates.
#' @param theta How each replicate is generated. \code{"conditional"} (the
#'   default) draws each person's responses from the Rasch conditional
#'   distribution given their observed raw score over their own observed
#'   items: sufficiency cancels the person parameter, so no ability is drawn
#'   at all, the observed score margins are reproduced exactly, and any tie
#'   between who answers and what is missing (a linked-booklet design,
#'   informative missingness) is preserved. The ability-sampling schemes
#'   remain: \code{"resample"} resamples the person estimates --- whose
#'   spread carries their estimation error, which the standardised fit
#'   statistics feel as anticonservatism at several thousand persons ---
#'   \code{"fixed"} reuses them as they stand, and \code{"normal"} draws
#'   from a normal with the error-corrected variance, at the price of losing
#'   replicates whose extreme categories go unvisited.
#' @param workers Number of parallel bootstrap workers. The default is four,
#'   reduced when fewer physical cores are available or the R process has a
#'   lower system limit. Starting them costs about half a second, which even
#'   the smallest useful run earns back; a larger count pays at large samples
#'   or large \code{B}, where eight workers run about five times faster than
#'   one. Person locations and the per-replicate seeds are drawn before
#'   distribution, so a fixed seed gives the same result for any worker
#'   count.
#' @param seed Optional seed. The caller's random stream is restored on
#'   exit.
#' @return A list with \code{items}, one row per item carrying each observed
#'   statistic beside its bootstrap probability and that probability's Holm
#'   adjustment across items; \code{total} for the whole-test chi-square and
#'   for the mean and SD of the item fit residuals, each against its own
#'   bootstrap null; \code{replicates}, the raw replicated statistics as one
#'   matrix per statistic; \code{B} requested and \code{B_used} (replicates
#'   that estimated); and the \code{theta} scheme.
#' @references Andrich, D. and Marais, I. (2019) \emph{A Course in Rasch
#'   Measurement Theory}. Springer.
#' @seealso \code{\link{chisq_detail}} for the class-interval breakdown
#'   behind one item's chi-square.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' # an exploratory run, kept small for speed: raw probabilities are usable,
#' # and the warning says what Holm-adjusted flagging at .05 would need
#' bs <- suppressWarnings(fit_bootstrap(rasch(X), B = 49, seed = 1))
#' bs$items[c("item", "chisq", "chisq_p_boot", "fit_resid", "fit_resid_p_boot")]
#' @export
fit_bootstrap <- function(fit, B = 200,
                          theta = c("conditional", "resample", "fixed", "normal"),
                          workers = 4L, seed = NULL) {
  if (!inherits(fit, "rasch"))
    stop("`fit` must be a fitted model from rasch()")
  if (inherits(fit, c("rasch_efrm", "rasch_mfrm", "rasch_explanatory")))
    .refuse("the item fit bootstrap generates from a single-facet Rasch ",
            "model; an extended-frame, many-facet or explanatory fit has a ",
            "generating structure this function does not reproduce")
  if (!is.null(fit$disc) && length(unique(fit$disc)) > 1L)
    .refuse("the item fit bootstrap generates under equal discriminations; ",
            "this fit carries frame units that differ across items")
  if (!is.null((fit$refit_spec %||% list())$pc_components))
    .refuse("the item fit bootstrap re-estimates each replicate the way the ",
            "fit was estimated; thresholds estimated through principal ",
            "components are not reproduced")
  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) || B != floor(B) ||
      B < 1 || B > .Machine$integer.max)
    stop("`B` must be one whole positive number of replicates")
  theta <- match.arg(theta)
  if (length(workers) != 1L || !is.numeric(workers) || !is.finite(workers) ||
      workers != floor(workers) || workers < 1 ||
      workers > .Machine$integer.max)
    stop("`workers` must be one whole positive number of workers")
  B <- as.integer(B)
  workers <- min(as.integer(workers), .rasch_available_workers())

  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed))
      stop("`seed` must be one finite number")
    old <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old), add = TRUE)
    set.seed(seed)
  }

  X <- fit$X
  L <- ncol(X)
  # the resolution floor: a two-sided bootstrap probability cannot fall
  # below 2/(B + 1), so Holm across the items cannot reach .05 unless
  # B is at least 40 times the item count (20 times for the one-sided
  # chi-square). Warn now rather than let a whole run produce probabilities
  # that could never have flagged anything.
  if (2 * L >= 0.05 * (B + 1))
    warning("B = ", B, " cannot reach Holm-adjusted significance at .05 ",
            "for ", L, " two-sided tests (floor 2L/(B+1) = ",
            signif(2 * L / (B + 1), 2), "); use B >= ", 40 * L,
            " for adjusted two-sided inference", call. = FALSE)
  spec <- fit$refit_spec %||% list()
  obs <- list(chisq = fit$items$chisq, fit_resid = fit$items$fit_resid,
              infit_ms = fit$items$infit_ms, outfit_ms = fit$items$outfit_ms,
              infit_z = fit$items$infit_z, outfit_z = fit$items$outfit_z)
  th <- fit$person$theta
  th <- th[is.finite(th)]
  if (length(th) < 2L)
    .refuse("the fit has fewer than two finite person locations to generate ",
            "from; the item fit bootstrap is unavailable")
  na_mask <- if (anyNA(X)) is.na(X) else NULL
  # a replicate must have the same number of persons as the fit, so a
  # location distribution shortened by non-finite estimates is resampled
  # back up to the fitted N
  if (length(th) < nrow(X)) th <- sample(th, nrow(X), replace = TRUE)

  # every random draw a replicate depends on is made here, in the parent:
  # its person locations (for the ability-sampling schemes), and the seed its
  # responses are generated under. A worker therefore reproduces its
  # replicate from what it is handed, and the worker count cannot move the
  # result.
  th_b <- if (theta == "conditional") NULL
          else lapply(seq_len(B), function(b) .fit_theta(theta, th, fit$psi))
  seeds <- sample.int(.Machine$integer.max, B)
  tau_list <- fit$tau_list; item_names <- colnames(X)
  model <- fit$model; ng <- fit$n_groups; anchors <- spec$anchors
  maxit <- spec$maxit %||% 60L; tol <- spec$tol %||% 1e-8
  one <- function(b) {
    old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old_stream), add = TRUE)
    set.seed(seeds[b])
    Xb <- if (is.null(th_b)) .fit_gen_conditional(X, tau_list, na_mask)
          else .fit_gen(th_b[[b]], tau_list, na_mask, item_names)
    .fit_refit(Xb, model, ng, anchors, maxit, tol)
  }
  reps <- .rasch_boot_apply(B, one, workers = workers,
                            label = "item fit bootstrap")

  keep <- !vapply(reps, is.null, TRUE)
  B_used <- sum(keep)
  if (B_used < 1L)
    .refuse("no bootstrap replicate could be estimated; the fitted model ",
            "generates data this estimator cannot fit")
  reps <- reps[keep]
  # replicates are not lost at random -- an unvisited extreme category is
  # likelier in the less dispersed ones -- so a materially thinned null is
  # worth saying out loud rather than leaving to be read off B_used
  if (B_used < 0.9 * B)
    warning(B - B_used, " of ", B, " bootstrap replicates could not be ",
            "estimated and the null is formed from the remaining ", B_used,
            "; replicates are not lost at random, so read the result with ",
            "that in mind", call. = FALSE)

  M <- lapply(names(.FIT_STATS), function(st)
    do.call(rbind, lapply(reps, `[[`, st)))
  names(M) <- names(.FIT_STATS)

  items <- data.frame(item = fit$items$item, stringsAsFactors = FALSE)
  for (st in names(.FIT_STATS)) {
    side <- .FIT_STATS[[st]]
    p <- vapply(seq_len(L), function(i)
      .boot_p(obs[[st]][i], M[[st]][, i], side), 0)
    usable <- is.finite(p)
    adj <- rep(NA_real_, L)
    adj[usable] <- stats::p.adjust(p[usable], method = "holm")
    items[[st]] <- obs[[st]]
    items[[paste0(st, "_p_boot")]] <- p
    items[[paste0(st, "_p_boot_adj")]] <- adj
  }
  items$n_boot <- colSums(!is.na(M$chisq))
  rownames(items) <- NULL

  # whole-test readings: the summed chi-square, and the mean and SD of the
  # item fit residuals that the fit summary reports. The SD is the reason
  # this one matters -- the convention reads it against 1, and under a true
  # model it is nowhere near 1 until the sample runs to a few thousand.
  ok <- !is.na(obs$chisq)
  tot_rep <- rowSums(M$chisq[, ok, drop = FALSE])
  fr_mean <- rowMeans(M$fit_resid, na.rm = TRUE)
  fr_sd <- apply(M$fit_resid, 1, stats::sd, na.rm = TRUE)
  total <- list(
    chisq = fit$total_chisq, df = fit$total_df, p = fit$total_chisq_p,
    chisq_p_boot = .boot_p(fit$total_chisq, tot_rep, "upper"),
    fit_resid_mean = fit$item_fit_summary$mean,
    fit_resid_mean_p_boot = .boot_p(fit$item_fit_summary$mean, fr_mean, "two"),
    fit_resid_sd = fit$item_fit_summary$sd,
    fit_resid_sd_p_boot = .boot_p(fit$item_fit_summary$sd, fr_sd, "two"),
    n_boot = sum(!is.na(tot_rep)))

  .tag_tables(list(items = items, total = total, replicates = M,
                   B = B, B_used = B_used, theta = theta))
}
