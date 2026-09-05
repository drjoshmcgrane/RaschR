# rasch :: model comparison
# ===========================================================================
# Side-by-side comparison of fitted models. Three kinds of evidence are
# reported. (1) The pairwise conditional log-likelihood with the number of
# structural parameters and, for fits of the SAME response data (identical
# items, categories, and persons, hence identical conditional information),
# twice the log-likelihood difference from the reference fit. Because the
# likelihood is a composite (pairwise) one, the difference is descriptive
# and is not chi-square calibrated; it is most meaningful for nested
# structures (for example the rating scale model inside the partial credit
# model, or equal units inside the extended frame of reference model).
# (2) Composite-likelihood information criteria that ARE calibrated for the
# pairwise over-counting: CL-AIC (Varin & Vidoni 2005) and CL-BIC (Gao &
# Song 2010) penalise -2 cl with the effective parameter count
# tr(H^-1 J) from the Godambe matrices -- the same quantity whose
# eigenvalues calibrate lr_test() -- instead of the nominal count, with
# the CL-BIC log(n) taken over the independent units (persons; judges for
# paired comparisons). Smaller is better, valid across models of the same
# data whether or not they nest.
# (3) Likelihood-free fit descriptors retained as descriptive context when
# the data preparation changes: the total trait chi-square against its degrees
# of freedom, the spread of the calibration and person fit residuals (ideal
# SD 1), the person separation index, and coefficient alpha where it is
# defined for an administered item matrix. They are not a formal selection
# criterion across different responses.
# ===========================================================================

# Composite-likelihood information criteria for one fit: effective parameter
# count tr(H^-1 J), CL-AIC, CL-BIC. NA when the fit does not carry its
# Godambe matrices (MFRM and EFRM assemble their own estimation structures).
.cl_ic <- function(f) {
  if (inherits(f, "rasch_btl")) {
    if (is.null(f$cl) || !isTRUE(f$cl$inference_available))
      return(c(eff = NA_real_, aic = NA_real_, bic = NA_real_))
    eff <- f$cl$eff_params; n <- f$cl$n_units; ll <- f$loglik
  } else {
    est <- f$est
    if (is.null(est$cov_beta) || is.null(est$H_beta))
      return(c(eff = NA_real_, aic = NA_real_, bic = NA_real_))
    p <- nrow(est$cov_beta)
    sensitivity <- -est$H_beta
    if (!.covariance_supports_wald(est$cov_beta, p) ||
        !.covariance_supports_wald(sensitivity, p))
      return(c(eff = NA_real_, aic = NA_real_, bic = NA_real_))
    # The sensitivity must identify every fitted beta direction. A singular
    # positive-semidefinite matrix cannot define an effective parameter count,
    # even though it passes the covariance gate used for identified contrasts.
    ev <- eigen((sensitivity + t(sensitivity)) / 2,
                symmetric = TRUE, only.values = TRUE)$values
    scale <- max(abs(ev))
    if (!is.finite(scale) || scale <= 0 ||
        min(ev) <= sqrt(.Machine$double.eps) * scale)
      return(c(eff = NA_real_, aic = NA_real_, bic = NA_real_))
    eff <- sum(diag(est$cov_beta %*% sensitivity))
    # Independent units are persons contributing at least one informative
    # pair, not response occasions. A missing identifier denotes its own
    # person, consistently with the clustered calibration sandwich.
    eligible <- !f$person$extreme & rowSums(!is.na(f$X)) >= 2L
    if (is.null(f$person$id)) {
      n <- sum(eligible)
    } else {
      id <- .role_text_values(f$person$id[eligible])
      missing <- is.na(id) | !nzchar(id)
      n <- length(unique(id[!missing])) + sum(missing)
    }
    ll <- est$loglik
  }
  if (!is.finite(eff) || eff < 0 || !is.finite(ll) ||
      !is.finite(n) || n < 1L)
    return(c(eff = NA_real_, aic = NA_real_, bic = NA_real_))
  c(eff = eff, aic = -2 * ll + 2 * eff, bic = -2 * ll + log(n) * eff)
}

.compare_person_count <- function(f) {
  if (is.null(f$person$id)) return(nrow(f$X))
  length(unique(.dif_ids(f$person$id)))
}

#' Compare fitted Rasch models
#'
#' Builds a comparison table for two or more fits from \code{\link{rasch}},
#' \code{\link{rasch_mfrm}}, \code{\link{rasch_efrm}}, or (all together)
#' \code{\link{btl}}. For fits of the
#' same response data (identical item columns, maximum scores, and number of
#' persons) the pairwise conditional log-likelihoods share their conditional
#' information, and twice the difference from the reference fit is reported
#' with the difference in parameter counts; this is descriptive (composite
#' likelihood), and most meaningful for nested structures such as RSM inside
#' PCM.
#'
#' The calibrated comparison is carried by the composite-likelihood
#' information criteria \code{cl_aic} (Varin and Vidoni 2005) and
#' \code{cl_bic} (Gao and Song 2010): \eqn{-2\,cl + c \cdot tr(H^{-1}J)},
#' with \eqn{c = 2} or \eqn{\log n}. Because every response enters every
#' pair its item forms, the pairwise log-likelihood over-counts the data;
#' the effective parameter count \eqn{tr(H^{-1}J)} from the Godambe
#' matrices -- the same quantity whose eigenvalues calibrate
#' \code{\link{lr_test}} -- absorbs exactly that over-counting, where the
#' nominal parameter count would not. \eqn{n} counts independent units:
#' persons contributing at least one informative pair, or judges for
#' paired-comparison fits (count-weighted comparisons when unclustered).
#' Smaller is better; the criteria are valid across models of the same data
#' whether or not they nest, and are \code{NA} (with the reason in the
#' printed note) for MFRM and EFRM fits, which do not carry their Godambe
#' matrices.
#'
#' Across different data preparations (subtests, splits, facet or frame
#' structures), or different allocations of response rows to persons, the
#' likelihood-based criteria are not comparable and are withheld. The table retains
#' descriptive context: total trait chi-square per degree of freedom,
#' calibration and person fit-residual SDs (ideal 1), PSI, and alpha where
#' applicable (OSI for paired comparisons). Alpha is \code{NA} when an MFRM
#' or EFRM item is represented by several response cells. These columns do
#' not provide a formal selection test across different response data.
#'
#' @param ... Two or more fitted objects, preferably given unique names. Supply either all
#'   Rasch-family fits or all \code{btl} fits. For \code{btl}, fits of the same
#'   comparison data (same objects, comparisons, and judges) support the
#'   likelihood columns -- e.g. free versus principal-component thresholds,
#'   with and without a position effect or within-judge dependence.
#'   Row order, arbitrary person or judge labels, and expansion versus count
#'   compression do not change data identity; allocation to independent
#'   persons or judges does.
#' @param reference Index or name of the reference fit for the
#'   log-likelihood difference; defaults to the first.
#' @return A data frame with one row per fit: label, model, persons, items
#'   (judges, objects, comparisons for \code{btl}),
#'   parameters, log-likelihood, \code{eff_params}, \code{cl_aic},
#'   \code{cl_bic}, comparability with the reference,
#'   \code{two_delta_ll} and \code{delta_parameters} (same-data fits only),
#'   chi-square per df, fit residual SDs, PSI, and alpha (OSI for
#'   \code{btl}).
#' @examples
#' set.seed(1)
#' simP <- function(th, tau) {
#'   x <- 0:length(tau)
#'   p <- exp(x * th - c(0, cumsum(tau)))
#'   p / sum(p)
#' }
#' th <- rnorm(400)
#' X <- sapply(seq(-1, 1, length.out = 6), function(b)
#'   sapply(th, function(t)
#'     sample(0:3, 1, prob = simP(t, b + c(-0.8, 0, 0.8)))))
#' colnames(X) <- paste0("R", 1:6)
#' compare_fits(PCM = rasch(X, model = "PCM"),
#'              RSM = rasch(X, model = "RSM"))
#' @export
compare_fits <- function(..., reference = 1) {
  if (is.character(reference)) {
    if (length(reference) != 1L || !is.null(dim(reference)) ||
        !is.null(oldClass(reference)) || is.na(reference) ||
        !nzchar(trimws(reference)))
      stop("`reference` must be one fit name or one whole number")
  } else reference <- .check_whole(reference, "reference", 1)
  fits <- list(...)
  if (length(fits) < 2) stop("supply at least two fits to compare")
  is_btl <- vapply(fits, inherits, TRUE, what = "rasch_btl")
  bad <- !vapply(fits, inherits, TRUE, what = "rasch") & !is_btl
  if (any(bad)) stop("argument(s) ", paste(which(bad), collapse = ", "),
                     " are not rasch or btl fits")
  if (any(is_btl) && !all(is_btl))
    stop("compare either Rasch-family fits or btl fits, not a mixture: ",
         "their likelihoods are over different data")
  labs <- names(fits)
  if (is.null(labs)) labs <- rep("", length(fits))
  if (anyNA(labs)) stop("fit names must not be missing")
  labs[labs == ""] <- paste0("fit", seq_along(fits))[labs == ""]
  if (anyDuplicated(labs))
    stop("fit names must be unique: ",
         paste(unique(labs[duplicated(labs)]), collapse = ", "))
  if (is.character(reference)) reference <- match(reference, labs)
  if (is.na(reference) || reference < 1 || reference > length(fits))
    stop("no such reference fit")

  if (all(is_btl)) {
    # Compare the actual comparison records exactly after canonicalising
    # orientation and row order. A scalar checksum can collide, and a
    # row-position-weighted checksum treats a harmless reorder as new data.
    # Weight and judge allocation are part of the composite data too: they
    # determine the likelihood contribution and independent clusters.
    seq_names <- function(f) intersect(c("exposure", "carry_over", "order"),
                                      names(f$comparisons))
    use_presented <- any(vapply(fits, function(f)
      !is.null(f$dependence) ||
        length(intersect(c("exposure", "carry_over", "position"),
                         names(f$comparisons))) > 0L, TRUE))
    # only sequence columns every compared fit carries can distinguish them
    seq_cols <- Reduce(intersect, lapply(fits, seq_names))
    if (is.null(seq_cols)) seq_cols <- character(0)
    sig <- function(f) {
      cmp <- f$comparisons
      if (is.null(cmp)) return(NULL)
      ca <- as.character(cmp$object_a); cb <- as.character(cmp$object_b)
      swap <- ca > cb
      lo <- ifelse(swap, cb, ca); hi <- ifelse(swap, ca, cb)
      resp <- as.numeric(cmp$response)
      resp[swap] <- max(f$m) - resp[swap]
      z <- data.frame(object_a = lo, object_b = hi, response = resp,
                      weight = as.numeric(cmp$weight),
                      stringsAsFactors = FALSE)
      # Orientation is arbitrary in a plain BTL, so the canonical form above
      # is the right fingerprint. It stops being arbitrary once any compared
      # fit models position or sequence: there the presented order enters
      # the likelihood, and two presentation designs would otherwise compare
      # as the same data. The decision is made across ALL the fits, so a
      # plain fit and a position-effect fit of the SAME comparisons still
      # match; the derived position column is a fitted covariate, not a
      # datum, so it never enters.
      if (use_presented) { z$presented_a <- ca; z$presented_b <- cb }
      for (cn in seq_cols) z[[cn]] <- as.character(cmp[[cn]])
      row_keys <- function(d, omit = character(0)) {
        x <- d[, setdiff(names(d), omit), drop = FALSE]
        vapply(seq_len(nrow(x)), function(i) {
          one <- x[i, , drop = FALSE]
          rownames(one) <- NULL
          paste(format(serialize(one, NULL, version = 3L)), collapse = "")
        }, character(1))
      }
      canonical_rows <- function(d) {
        d <- d[do.call(order, c(d, list(na.last = TRUE))), , drop = FALSE]
        rownames(d) <- NULL
        d
      }
      aggregate_rows <- function(d) {
        # A count-compressed row and the corresponding repeated rows are the
        # same comparison data.  Sum weights only over rows identical in
        # every field that enters the likelihood or its design.
        key <- row_keys(d, "weight")
        first <- !duplicated(key)
        out <- d[first, , drop = FALSE]
        out$weight <- as.numeric(rowsum(d$weight, match(key, unique(key)),
                                        reorder = FALSE))
        canonical_rows(out)
      }
      block_key <- function(d)
        paste(format(serialize(canonical_rows(d), NULL, version = 3L)),
              collapse = "")
      judge <- as.character(cmp$judge)
      if (all(is.na(judge))) {
        ia <- f$independent_allocation
        if (is.data.frame(ia) && nrow(ia) == nrow(z) &&
            all(c("unit", "replicates") %in% names(ia))) {
          # Express each independent unit once, plus its replication count.
          # This makes an expanded dataset identical to its count-compressed
          # form.  Half-scored ties remain a two-row atomic unit and therefore
          # cannot be confused with two independent opposing judgements.
          blocks <- lapply(split(seq_len(nrow(z)), as.character(ia$unit)),
                           function(ii) {
            nr <- unique(as.numeric(ia$replicates[ii]))
            if (length(nr) != 1L || !is.finite(nr) || nr <= 0)
              return(NULL)
            d <- z[ii, , drop = FALSE]
            d$weight <- d$weight / nr
            list(data = canonical_rows(d), n = nr)
          })
          if (length(blocks) && !any(vapply(blocks, is.null, logical(1)))) {
            keys <- vapply(blocks, function(b) block_key(b$data), "")
            uk <- unique(keys)
            allocation <- lapply(uk, function(k) list(
              data = blocks[[which(keys == k)[1L]]]$data,
              n = sum(vapply(blocks[keys == k], `[[`, numeric(1), "n"))))
            allocation <- unname(allocation[order(vapply(
              allocation, function(b) block_key(b$data), ""))])
          } else allocation <- canonical_rows(z)
        } else {
          # Compatibility with fits made before the allocation was retained.
          allocation <- canonical_rows(z)
        }
      } else {
        # Judge labels are arbitrary; the allocation of rows to judges is not.
        # Count compression inside a judge is representational, so identical
        # rows are combined before anonymous judge blocks are ordered.
        blocks <- lapply(split(seq_len(nrow(z)), judge), function(ii)
          aggregate_rows(z[ii, , drop = FALSE]))
        key <- vapply(blocks, block_key, "")
        allocation <- unname(blocks[order(key)])
      }
      list(objects = sort(f$objects$object), comparisons = allocation)
    }
    ref_sig <- sig(fits[[reference]])
    rows <- lapply(seq_along(fits), function(i) {
      f <- fits[[i]]
      conv <- isTRUE(f$converged)
      ic <- if (conv) .cl_ic(f)
            else c(eff = NA_real_, aic = NA_real_, bic = NA_real_)
      dep <- if (is.null(f$dependence)) "" else
        paste0(" + ", paste(f$dependence$effect, collapse = " + "))
      model_label <- if (inherits(f, "rasch_btl_efrm"))
        "BTL with frame-dependent units"
      else if (inherits(f, "rasch_btl_explanatory"))
        paste0("Explanatory BTL (", if (max(f$m) > 1L)
          "polytomous" else "dichotomous", ")")
      else paste0("BTL (", if (max(f$m) > 1L)
        paste0("polytomous, ", f$thr_structure, " thresholds") else
          "dichotomous", dep, ")")
      data.frame(
        label = labs[i], converged = conv,
        model = model_label,
        judges = if (is.null(f$judges)) NA_integer_ else nrow(f$judges),
        objects = nrow(f$objects), comparisons = f$n_comparisons,
        parameters = if (!is.null(f$n_parameters)) f$n_parameters
          else if (is.null(f$cl)) NA_integer_ else f$cl$n_parameters,
        loglik = if (conv) f$loglik else NA_real_,
        eff_params = unname(ic["eff"]), cl_aic = unname(ic["aic"]),
        cl_bic = unname(ic["bic"]),
        same_data = identical(sig(f), ref_sig),
        two_delta_ll = NA_real_, delta_parameters = NA_integer_,
        chisq_per_df = if (is.finite(f$total_df) && f$total_df > 0)
          f$total_chisq / f$total_df else NA_real_,
        OSI = f$osi$PSI)
    })
  } else {
    # same_data must compare the ACTUAL responses, not just the item names,
    # maximum scores, and person count: two different datasets sharing those
    # margins would otherwise be declared the same data and get a spurious
    # two_delta_ll. The full response matrix is the exact fingerprint.
    sig <- function(f) {
      # Row and column order are presentation details, not different response
      # data.  Canonicalise the items by name, then represent each independent
      # person by the (sorted) multiset of response rows belonging to that
      # person.  Sorting those person blocks also makes an arbitrary relabeling
      # of the IDs immaterial while retaining repeated-person allocation.
      item_order <- order(colnames(f$X))
      X <- unname(as.matrix(f$X[, item_order, drop = FALSE]))
      row_key <- apply(X, 1L, function(z)
        paste(ifelse(is.na(z), "NA", format(z, scientific = FALSE,
                                             trim = TRUE)), collapse = ","))
      person <- if (is.null(f$person$id)) seq_len(nrow(X)) else
        .dif_ids(f$person$id)
      blocks <- unname(vapply(split(row_key, person), function(z)
        paste(sort(z), collapse = "|"), character(1)))
      list(items = colnames(f$X)[item_order],
           m = unname(f$m[item_order]), n = nrow(X),
           person_blocks = sort(blocks))
    }
    # underlying item count: MFRM/EFRM columns of X are virtual item-by-facet
    # or item-by-group cells, not the real items
    n_items <- function(f)
      if (!is.null(f$item_effects)) nrow(f$item_effects)
      else if (!is.null(f$item_arbitrary)) nrow(f$item_arbitrary)
      else ncol(f$X)
    ref_sig <- sig(fits[[reference]])
    rows <- lapply(seq_along(fits), function(i) {
      f <- fits[[i]]
      conv <- isTRUE(f$est$converged)
      # an unconverged fit has no trustworthy log-likelihood or information:
      # withhold its information criteria rather than rank on them
      ic <- if (conv) .cl_ic(f)
            else c(eff = NA_real_, aic = NA_real_, bic = NA_real_)
      data.frame(
        label = labs[i],
        model = if (inherits(f, "rasch_explanatory"))
          f$explanatory_model else f$model,
        converged = conv,
        persons = .compare_person_count(f), items = n_items(f),
        parameters = if (is.null(f$est$n_parameters)) NA_integer_
                     else f$est$n_parameters,
        loglik = if (conv) f$est$loglik else NA_real_,
        eff_params = unname(ic["eff"]), cl_aic = unname(ic["aic"]),
        cl_bic = unname(ic["bic"]),
        same_data = identical(sig(f), ref_sig),
        two_delta_ll = NA_real_, delta_parameters = NA_integer_,
        chisq_per_df = if (is.finite(f$total_df) && f$total_df > 0)
          f$total_chisq / f$total_df else NA_real_,
        item_fit_sd = f$item_fit_summary$sd,
        person_fit_sd = f$person_fit_summary$sd,
        PSI = f$psi$PSI, alpha = f$alpha$alpha)
    })
  }
  out <- do.call(rbind, rows)
  ref <- out[reference, ]
  # Composite information criteria compare models only when the responses
  # and independent-unit allocation match the reference. Retaining numeric
  # criteria for a different dataset invites a comparison with no common
  # target likelihood.
  out$cl_aic[!out$same_data] <- NA_real_
  out$cl_bic[!out$same_data] <- NA_real_
  # the descriptive two_delta_ll needs the same data AND two trustworthy
  # (converged) log-likelihoods
  cmp <- out$same_data & out$converged & isTRUE(ref$converged) &
    seq_len(nrow(out)) != reference
  out$two_delta_ll[cmp] <- 2 * (out$loglik[cmp] - ref$loglik)
  out$delta_parameters[cmp] <- out$parameters[cmp] - ref$parameters
  rownames(out) <- NULL
  attr(out, "reference") <- labs[reference]
  attr(out, "note") <- paste0(
    "cl_aic and cl_bic are composite-likelihood information criteria ",
    "(Varin & Vidoni 2005; Gao & Song 2010): -2 cl penalised by the ",
    "effective parameter count tr(H^-1 J), which absorbs the pairwise ",
    "over-counting that the nominal count would not; smaller is better, ",
    "valid across models of the same data",
    if (any(!vapply(fits, function(f) {
      converged <- if (inherits(f, "rasch_btl")) isTRUE(f$converged) else
        isTRUE(f$est$converged)
      converged && is.finite(.cl_ic(f)["eff"])
    }, TRUE)))
      paste0(" (NA for MFRM/EFRM fits without the required Godambe ",
             "matrices, non-converged fits, and judge-clustered BTL fits ",
             "with too few independent clusters)")
    else "",
    ". Information criteria and two_delta_ll are withheld when the response ",
    "data or independent-unit allocation differ from the reference. ",
    "two_delta_ll is the raw composite difference against the reference, ",
    "descriptive only. Across different data preparations, chisq_per_df, ",
    "the fit residual SDs and separation/reliability columns provide ",
    "descriptive context rather than a formal selection criterion.")
  class(out) <- c("rasch_compare", "data.frame")
  out
}

#' @export
print.rasch_compare <- function(x, ...) {
  cat(sprintf("Model comparison (reference: %s)\n\n", attr(x, "reference")))
  core <- c("label", "model", "persons", "items", "judges", "objects",
            "eff_params", "cl_aic", "cl_bic", "two_delta_ll",
            "chisq_per_df", "item_fit_sd", "person_fit_sd", "PSI", "alpha",
            "OSI")
  y <- as.data.frame(x)[, intersect(core, names(x)), drop = FALSE]
  print(.fmt_df(y), row.names = FALSE)
  cat("(further columns on the object: loglik, parameters, same_data)\n")
  cat("\n", attr(x, "note"), "\n", sep = "")
  invisible(x)
}

#' Compare the partial credit and rating scale models
#'
#' Compares a fitted partial credit model with the rating scale
#' reparameterisation of the same data. Both raw and composite-likelihood
#' adjusted statistics are returned.
#'
#' @details
#' The pairwise conditional likelihood is a composite likelihood: each
#' response contributes to every item pair in which it appears. Consequently,
#' the raw statistic \eqn{W=2(cl_{PCM}-cl_{RSM})} does not have an ordinary
#' chi-square reference distribution. Its limiting distribution is
#' \eqn{\sum_j\lambda_j\chi^2_1} (Kent 1982; Varin, Reid and Firth 2011), where
#' the \eqn{\lambda_j} are obtained from the sensitivity matrix \eqn{H},
#' variability matrix \eqn{J}, and the constraints defining the RSM. The
#' mean-matched statistic is
#' \deqn{W_{adj}=rW/\sum_j\lambda_j,}
#' with \eqn{r} degrees of freedom.
#'
#' Use \code{p_adj} for inference. The unadjusted \code{p} is retained for
#' descriptive comparison with conventional displays. The adjustment is a
#' first-order approximation and can be mildly anti-conservative in small
#' samples with long polytomous tests. Interpret values near the nominal
#' level cautiously in such designs.
#'
#' @param fit An unrestricted, unanchored \code{"PCM"} fit from
#'   \code{\link{rasch}} with equal maximum scores across items (the rating
#'   parameterisation requires them).
#' @param maxit,tol Passed to the rating-scale refit.
#' @return A list of class \code{"rasch_lr"}: raw \code{chisq}, \code{df},
#'   \code{p} (the conventional display); adjusted \code{chisq_adj}, \code{p_adj},
#'   and the eigenvalues \code{lambda}; the two log-likelihoods; and the
#'   rating-scale refit (\code{fit_rsm}).
#' @references Kent, J. T. (1982). Robust properties of likelihood ratio
#'   tests. Biometrika, 69, 19-27. Varin, C., Reid, N. and Firth, D. (2011).
#'   An overview of composite likelihood methods. Statistica Sinica, 21,
#'   5-42.
#' @examples
#' set.seed(1)
#' tau <- c(-0.7, 0.7)
#' X <- sapply(seq(-1, 1, length.out = 6), function(d) vapply(rnorm(300),
#'   function(b) sample(0:2, 1, prob = item_moments(b, tau + d)$P), 0L))
#' colnames(X) <- paste0("Q", 1:6)
#' lr_test(rasch(X, model = "PCM"))
#' @export
lr_test <- function(fit, maxit = 60, tol = 1e-8) {
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    stop("`fit` must be a fitted Rasch-family object from rasch()")
  if (!identical(fit$model, "PCM"))
    stop("lr_test() compares an unrestricted (PCM) fit with its rating ",
         "re-parameterisation; supply a PCM fit")
  if (length(unique(fit$m)) != 1L)
    stop("the rating parameterisation requires equal maximum scores across items")
  if (max(fit$m) < 2L)
    stop("with dichotomous items the two parameterisations coincide")
  if (!isTRUE(fit$est$converged))
    stop("the PCM fit did not converge; its likelihood cannot support a model comparison")
  if (inherits(fit, "rasch_explanatory"))
    stop("lr_test() compares an unrestricted PCM fit with its rating ",
         "re-parameterisation; an explanatory fit already restricts the ",
         "thresholds to its design, and the rating refit would drop that ",
         "restriction, so the two models are not nested")
  spec <- fit$refit_spec
  if (!is.null(spec$anchors) && nrow(spec$anchors))
    stop("lr_test() requires an unrestricted PCM fit; fixed threshold anchors change the null constraints")
  if (!is.null(spec$pc_components))
    stop("lr_test() requires an unrestricted PCM fit; principal-component threshold constraints are already a restricted model")
  # the refit is meant to sit beside the PCM fit, so it must be grouped the
  # same way: dropping the person factors makes the returned fit useless for
  # the follow-up analyses
  rsm <- rasch(fit$X, model = "RSM", n_groups = fit$n_groups,
               id = fit$person$id, factors = fit$factors,
               maxit = maxit, tol = tol)
  if (!isTRUE(rsm$est$converged))
    stop("the rating-scale refit did not converge; the model comparison is unavailable")
  chisq <- max(0, 2 * (fit$est$loglik - rsm$est$loglik))
  df <- fit$est$n_parameters - rsm$est$n_parameters

  # composite-likelihood calibration: eigenvalues of the Godambe ratio over
  # the constrained directions (see Details)
  chisq_adj <- p_adj <- NA_real_; lambda <- NULL
  Bp <- fit$est$B; Br <- rsm$est$B
  if (!is.null(Bp) && !is.null(Br) && !is.null(fit$est$H_beta)) {
    M <- nrow(Bp)
    S <- cbind(Br, rep(1, M))            # rating subspace + the null shift
    ss <- svd(S)
    rs <- sum(ss$d > max(1e-10, max(ss$d) * 1e-8))
    U <- ss$u[, seq_len(rs), drop = FALSE]
    Portho <- diag(M) - tcrossprod(U)
    A <- Portho %*% Bp                   # constraint: A beta = 0
    qa <- qr(t(A))
    r <- qa$rank
    if (r > 0) {
      C <- qr.Q(qa)[, seq_len(r), drop = FALSE]
      Hinv <- tryCatch(solve(-fit$est$H_beta), error = function(e) NULL)
      kc <- if (is.null(Hinv))
        list(chisq = NA_real_, p = NA_real_, lambda = numeric(0)) else
        .kent_calibration(chisq, C, fit$est$cov_beta, Hinv)
      chisq_adj <- kc$chisq
      p_adj <- kc$p
      lambda <- kc$lambda
      df <- r
    }
  }
  out <- list(chisq = chisq, df = df,
              p = pchisq(chisq, df, lower.tail = FALSE),
              chisq_adj = chisq_adj, p_adj = p_adj, lambda = lambda,
              loglik_pcm = fit$est$loglik, loglik_rsm = rsm$est$loglik,
              fit_rsm = rsm)
  out <- .tag_tables(out)
  class(out) <- "rasch_lr"
  out
}

#' @export
print.rasch_lr <- function(x, ...) {
  cat("Likelihood-ratio test: partial credit vs rating parameterisation\n")
  cat(sprintf("  Raw composite chi-square %.3f on %d df, p = %s (conventional display; anticonservative)\n",
              x$chisq, x$df, .fmt_p(x$p)))
  if (is.finite(x$chisq_adj))
    cat(sprintf("  Adjusted chi-square %.3f on %d df, p = %s (Kent 1982 first-order calibration)\n",
                x$chisq_adj, x$df, .fmt_p(x$p_adj)))
  cat(sprintf("  log-likelihood (pairwise composite): PCM %.3f, RSM %.3f\n",
              x$loglik_pcm, x$loglik_rsm))
  invisible(x)
}
