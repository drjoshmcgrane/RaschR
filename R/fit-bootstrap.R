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

.fit_boot_failure <- function(status = c("error", "nonconverged")) {
  status <- match.arg(status)
  structure(list(status = status), class = "rasch_fit_boot_failure")
}

.fit_boot_refuse <- function(..., B, B_used, B_nonconverged, B_errors) {
  stop(structure(
    class = c("rasch_fit_bootstrap_refusal", "rasch_refusal", "error",
              "condition"),
    list(message = paste0(...), call = NULL, B = B, B_used = B_used,
         B_nonconverged = B_nonconverged, B_errors = B_errors)))
}

.fit_boot_status <- function(x) {
  if (inherits(x, "rasch_fit_boot_failure")) return(x$status)
  if (is.null(x)) return("error")
  "ok"
}

.fit_boot_md5 <- function(x) {
  path <- tempfile("rasch-fit-signature-", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  dput(x, path, control = c("keepNA", "keepInteger", "showAttributes",
                            "hexNumeric"))
  unname(tools::md5sum(path))
}

.fit_boot_signature <- function(fit) {
  if (inherits(fit, "rasch_btl")) {
    return(list(
      kind = "btl", m = fit$m, n_comparisons = fit$n_comparisons,
      objects = as.character(fit$objects$object),
      locations = unname(as.numeric(fit$objects$location)),
      pairs = .btl_pair_key(fit$pairs$object_a, fit$pairs$object_b),
      total_chisq = unname(as.numeric(fit$total_chisq)),
      fingerprint = .fit_boot_md5(list(
        comparisons = fit$comparisons, refit_spec = fit$refit_spec,
        objects = fit$objects, total_chisq = fit$total_chisq))))
  }
  list(
    kind = "rasch", model = fit$model, dimensions = dim(fit$X),
    items = as.character(fit$items$item),
    item_chisq = unname(as.numeric(fit$items$chisq)),
    raw_scores = unname(as.numeric(fit$person$raw)),
    fingerprint = .fit_boot_md5(list(
      X = fit$X, model = fit$model, thresholds = fit$thresholds,
      n_groups = fit$n_groups, refit_spec = fit$refit_spec)))
}

.new_fit_bootstrap <- function(x, fit, kind) {
  x$model_kind <- kind
  x$fit_signature <- .fit_boot_signature(fit)
  x <- .tag_tables(x)
  class(x) <- c("rasch_fit_bootstrap", "list")
  x
}

.validate_fit_bootstrap <- function(bootstrap, fit) {
  if (is.null(bootstrap)) return(invisible(NULL))
  if (!inherits(bootstrap, "rasch_fit_bootstrap"))
    stop("`bootstrap` must be a current fit_bootstrap() result")
  expected <- if (inherits(fit, "rasch_btl")) "btl" else "rasch"
  if (!identical(bootstrap$model_kind, expected))
    stop("`bootstrap` was computed for a different model family")
  required <- if (expected == "btl")
    c("pairs", "objects", "total", "B", "B_used") else
    c("items", "total", "B", "B_used")
  if (!all(required %in% names(bootstrap)))
    stop("`bootstrap` is incomplete and cannot be exported")
  if (!identical(bootstrap$fit_signature, .fit_boot_signature(fit)))
    stop("`bootstrap` was computed from a different fitted model")
  invisible(NULL)
}

# Everything the item fit statistics need from a response matrix and nothing
# else: the person tables, score curves and targeting that rasch() also builds
# are not used here and would multiply the cost by the number of replicates.
# The class-interval count is imposed rather than re-derived per replicate, so
# every replicate computes the same statistic as the observed fit rather than
# one on a neighbouring number of intervals.
.fit_refit <- function(X, model, n_groups, anchors, maxit, tol) {
  est <- tryCatch(pcml(X, model = model, anchors = anchors,
                       maxit = maxit, tol = tol),
                  error = function(e) .fit_boot_failure("error"))
  if (inherits(est, "rasch_fit_boot_failure")) return(est)
  if (!isTRUE(est$converged)) return(.fit_boot_failure("nonconverged"))
  L <- ncol(X)
  tau_list <- lapply(seq_len(L), function(i) est$thr$tau[est$thr$item == i])
  if (!all(vapply(tau_list, function(t) length(t) && all(is.finite(t)), TRUE)))
    return(.fit_boot_failure("error"))
  person <- .person_estimates(X, tau_list, disc = 1)
  if (all(is.na(person$theta))) return(.fit_boot_failure("error"))
  mo <- .moment_arrays(person$theta, tau_list, disc = rep(1, L))
  Z <- (X - mo$E) / sqrt(mo$V)
  colnames(Z) <- colnames(X)
  m_i <- vapply(tau_list, length, 1L)
  item_extreme <- vapply(seq_len(L), function(j) {
    col <- X[!person$extreme, j]
    tot <- sum(col, na.rm = TRUE); nn <- sum(!is.na(col))
    nn == 0L || tot == 0 || tot == nn * m_i[j]
  }, logical(1))
  ci <- .class_intervals(person$theta, person$extreme, n_groups)
  ci_list <- if (anyNA(X))
    .class_intervals_by_item(X, person$theta, person$extreme, n_groups) else NULL
  it <- .item_trait(X, mo, ci, ci_list = ci_list)
  ifit <- .item_fit(X, Z, mo, extreme = person$extreme)
  pfit <- .person_fit(X, Z, mo, item_extreme = item_extreme)
  n_par <- if (is.null(est$n_parameters)) nrow(est$thr) - 1L else est$n_parameters
  rf <- .fitres(Z, mo, person$extreme, n_par)
  list(chisq = it$chisq, fit_resid = rf$items$fit_resid,
       infit_ms = ifit$infit_ms, outfit_ms = ifit$outfit_ms,
       infit_z = ifit$infit_z, outfit_z = ifit$outfit_z,
       persons = list(fit_resid = rf$persons$fit_resid,
                      infit_ms = pfit$infit_ms,
                      outfit_ms = pfit$outfit_ms,
                      infit_z = pfit$infit_z,
                      outfit_z = pfit$outfit_z))
}

# The statistics this bootstrap calibrates, and how each is read. The
# chi-square is a discrepancy and only its upper tail is misfit; the others
# depart in two directions -- above for a flatter item than the model
# predicts, below for a steeper one -- and both readings are misfit, so their
# p-values are equal-tailed.
.FIT_STATS <- c(chisq = "upper", fit_resid = "two", infit_ms = "two",
                outfit_ms = "two", infit_z = "two", outfit_z = "two")
.PERSON_FIT_STATS <- c(fit_resid = "two", infit_ms = "two",
                       outfit_ms = "two", infit_z = "two", outfit_z = "two")

# Bootstrap p-value for one statistic against its own replicated null.
# (1 + r)/(1 + B) never returns zero; the two-sided form doubles the smaller
# tail and is capped, so its resolution is half the one-sided form's.
.boot_p <- function(obs, null, side) {
  null <- null[is.finite(null)]
  n <- length(null)
  if (!is.finite(obs) || n < 1L) return(NA_real_)
  up <- (1 + sum(null >= obs)) / (1 + n)
  if (side == "upper") return(up)
  lo <- (1 + sum(null <= obs)) / (1 + n)
  min(1, 2 * min(up, lo))
}

# A requested exploratory run below 30 draws may still be useful at its very
# coarse Monte Carlo resolution, but it must retain a majority.  At 30 or
# above, use the package-wide floor of 30 successful draws as well as a
# majority.  This is a p-value null, not a covariance-rank condition.
.fit_min_boot_success <- function(B) {
  as.integer(max(if (B >= 30L) 30L else 1L, floor(B / 2) + 1L))
}

# Single-step maximum-statistic adjustment from the joint bootstrap null. The
# centring and scale are statistic-specific, while each replicate retains the
# dependence among the quantities. Unlike an empirical minimum-p reference,
# maxT does not collapse to the endpoint when the family is larger than B.
.boot_maxt <- function(obs, null, side, min_success,
                       mode = c("studentised", "centred", "raw"),
                       transform = identity) {
  mode <- match.arg(mode)
  if (!is.matrix(null)) null <- as.matrix(null)
  n <- ncol(null)
  n_stat <- colSums(is.finite(null))
  raw <- vapply(seq_len(n), function(i) {
    if (n_stat[i] < min_success) return(NA_real_)
    .boot_p(obs[i], null[, i], side)
  }, 0)
  family <- is.finite(obs)
  adj <- rep(NA_real_, n)
  family_n <- sum(family)
  if (!family_n)
    return(list(p = raw, p_adj = adj, n_boot = n_stat,
                family_n = 0L, family_boot = 0L, adjusted_n = 0L))

  V <- transform(null[, family, drop = FALSE])
  obs_e <- transform(obs[family])
  complete <- stats::complete.cases(V)
  family_boot <- sum(complete)
  # The family is declared by the finite observed statistics, not by which
  # null columns happen to survive. A partial maxT family would make every
  # remaining adjusted probability depend on an estimation failure.
  if (all(is.finite(raw[family])) && all(is.finite(obs_e)) &&
      family_boot >= min_success) {
    V <- V[complete, , drop = FALSE]
    centre <- if (mode == "raw") rep(0, ncol(V)) else colMeans(V)
    scale <- if (mode == "studentised") apply(V, 2L, stats::sd) else
      rep(1, ncol(V))
    spread <- apply(V, 2L, stats::sd)
    stable <- is.finite(spread) & spread > sqrt(.Machine$double.eps) &
      is.finite(scale) & scale > sqrt(.Machine$double.eps)
    if (all(stable)) {
      Z <- sweep(sweep(V, 2L, centre, "-"), 2L, scale, "/")
      z_obs <- (obs_e - centre) / scale
      if (side == "two") {
        Z <- abs(Z)
        z_obs <- abs(z_obs)
      }
      max_null <- apply(Z, 1L, max)
      a <- vapply(z_obs, function(z)
        (1 + sum(max_null >= z)) / (1 + length(max_null)), 0)
      pos <- which(family)
      adj[pos] <- pmax(raw[pos], a)
    }
  }
  list(p = raw, p_adj = adj, n_boot = n_stat,
       family_n = family_n, family_boot = family_boot,
       adjusted_n = sum(is.finite(adj)))
}

.btl_pair_key <- function(a, b) {
  lo <- ifelse(a <= b, a, b)
  hi <- ifelse(a <= b, b, a)
  paste0(nchar(lo, type = "bytes"), ":", lo,
         nchar(hi, type = "bytes"), ":", hi)
}

.btl_boot_data <- function(fit) {
  cmp <- fit$comparisons
  spec <- fit$refit_spec %||% list()
  if (is.null(cmp) || !all(c("object_a", "object_b", "weight") %in% names(cmp)))
    .refuse("the paired-comparison fit does not carry the comparison design ",
            "needed for a fit bootstrap; refit it with the current package version")
  if (any(!is.finite(cmp$weight)) || any(cmp$weight <= 0) ||
      any(cmp$weight != floor(cmp$weight)))
    .refuse("the paired-comparison fit bootstrap requires whole positive ",
            "comparison counts; a fit using half-weighted ties cannot be generated ",
            "as independent comparisons. Code ties as an ordered middle category")
  if (any(grepl("half a win", fit$notes %||% character(0), fixed = TRUE)))
    .refuse("the paired-comparison fit used half-weighted ties, which the fit ",
            "bootstrap cannot generate as independent comparisons. Code ties ",
            "as an ordered middle category")

  requested_effects <- intersect(c("exposure", "carry_over", "position"),
                                 names(cmp))
  fitted_effects <- if (is.null(fit$dependence)) character(0) else
    as.character(fit$dependence$effect)
  dropped_effects <- setdiff(requested_effects, fitted_effects)
  if (length(dropped_effects))
    .refuse("the paired-comparison fit dropped the following requested ",
            "history or position effect(s): ",
            paste(dropped_effects, collapse = ", "),
            ". The bootstrap cannot refit a reduced dependence model through ",
            "the public btl() specification")

  has_order <- isTRUE(spec$has_order)
  if (has_order && (any(cmp$weight != 1) || !"order" %in% names(cmp)))
    .refuse("an ordered paired-comparison bootstrap needs one row per judgement ",
            "with a unique sequence value; count-weighted ordered rows are not a ",
            "generative history")

  if (!has_order) {
    P <- fit$fitted_prob
    if (!is.matrix(P) || nrow(P) != nrow(cmp) || ncol(P) != fit$m + 1L ||
        any(!is.finite(P)) || any(abs(rowSums(P) - 1) > 1e-8))
      .refuse("the paired-comparison fit does not carry valid fitted category ",
              "probabilities; refit it with the current package version")
    rows <- vector("list", nrow(cmp))
    for (r in seq_len(nrow(cmp))) {
      z <- as.vector(stats::rmultinom(1L, size = cmp$weight[r], prob = P[r, ]))
      k <- which(z > 0L)
      rows[[r]] <- data.frame(
        object_a = cmp$object_a[r], object_b = cmp$object_b[r],
        response = k - 1L, count = z[k],
        judge = cmp$judge[r], stringsAsFactors = FALSE)
    }
    out <- do.call(rbind, rows)
    # Keep the fitted response scale even when a replicate happens not to
    # visit an extreme category. Numeric scores would make btl() infer a
    # shorter scale from that replicate and fit a different model.
    out$response <- ordered(out$response, levels = 0:fit$m)
    return(out)
  }

  # History-dependent effects must be generated in sequence. Holding the
  # observed carry-over covariate fixed would condition on the outcome being
  # tested and produce the wrong null.
  beta <- stats::setNames(fit$objects$location, fit$objects$object)
  if (anyNA(beta[cmp$object_a]) || anyNA(beta[cmp$object_b]))
    .refuse("the ordered comparison design names an object without a fitted location")
  tau <- if (is.null(fit$thresholds)) numeric(fit$m) else fit$thresholds$tau
  dep <- if (is.null(fit$dependence)) numeric(0) else
    stats::setNames(fit$dependence$estimate, fit$dependence$effect)
  out <- cmp[, c("object_a", "object_b", "judge", "order"), drop = FALSE]
  out$response <- NA_integer_
  out$count <- 1L
  cnt <- new.env(hash = TRUE, parent = emptyenv())
  tot <- new.env(hash = TRUE, parent = emptyenv())
  gets <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
  key_a <- .factor_keys(data.frame(judge = cmp$judge,
                                   object = cmp$object_a,
                                   stringsAsFactors = FALSE))
  key_b <- .factor_keys(data.frame(judge = cmp$judge,
                                   object = cmp$object_b,
                                   stringsAsFactors = FALSE))
  for (r in order(cmp$judge, cmp$order)) {
    ka <- key_a[r]; kb <- key_b[r]
    na <- gets(cnt, ka); nb <- gets(cnt, kb)
    exposure <- as.numeric(na > 0) - as.numeric(nb > 0)
    carry <- (if (na > 0) gets(tot, ka) / na else 0) -
      (if (nb > 0) gets(tot, kb) / nb else 0)
    d <- beta[[cmp$object_a[r]]] - beta[[cmp$object_b[r]]]
    if ("exposure" %in% names(dep)) d <- d + dep[["exposure"]] * exposure
    if ("carry_over" %in% names(dep)) d <- d + dep[["carry_over"]] * carry
    if ("position" %in% names(dep)) d <- d + dep[["position"]]
    eta <- d * (0:fit$m) - c(0, cumsum(tau))
    eta <- eta - max(eta)
    pr <- exp(eta) / sum(exp(eta))
    x <- sample.int(fit$m + 1L, 1L, prob = pr) - 1L
    out$response[r] <- x
    cnt[[ka]] <- na + 1; cnt[[kb]] <- nb + 1
    tot[[ka]] <- gets(tot, ka) + (2 * x / fit$m - 1)
    tot[[kb]] <- gets(tot, kb) + (2 * (fit$m - x) / fit$m - 1)
  }
  out$response <- ordered(out$response, levels = 0:fit$m)
  out
}

.btl_boot_refit <- function(fit) {
  d <- .btl_boot_data(fit)
  spec <- fit$refit_spec %||% list()
  args <- list(data = d, object_a = "object_a", object_b = "object_b",
               response = "response", count = "count",
               thresholds = spec$thresholds %||% "free",
               position = isTRUE(spec$position),
               anchors = spec$anchors,
               maxit = spec$maxit %||% 60L, tol = spec$tol %||% 1e-8,
               .object_design = spec$object_design)
  if (any(!is.na(d$judge))) args$judge <- "judge"
  if (isTRUE(spec$has_order)) args$order <- "order"
  nonconv <- FALSE
  z <- tryCatch(withCallingHandlers(
    do.call(btl, args), warning = function(w) {
      # A non-convergence warning makes that replicate unusable. Other
      # warnings are expected occasionally at sparse bootstrap boundaries;
      # the returned fit decides whether the replicate is estimable.
      if (grepl("did NOT converge", conditionMessage(w), fixed = TRUE))
        nonconv <<- TRUE
      invokeRestart("muffleWarning")
    }), error = function(e) .fit_boot_failure("error"))
  if (inherits(z, "rasch_fit_boot_failure")) return(z)
  if (nonconv || !isTRUE(z$converged))
    return(.fit_boot_failure("nonconverged"))
  wanted_effects <- if (is.null(fit$dependence)) character(0) else
    sort(as.character(fit$dependence$effect))
  fitted_effects <- if (is.null(z$dependence)) character(0) else
    sort(as.character(z$dependence$effect))
  if (!identical(wanted_effects, fitted_effects))
    return(.fit_boot_failure("error"))
  active <- if (is.null(fit$objects$extreme))
    rep(TRUE, nrow(fit$objects)) else !fit$objects$extreme
  active_objects <- fit$objects$object[active]
  if (!setequal(as.character(z$objects$object), as.character(active_objects)) ||
      any(z$objects$extreme %||% FALSE)) return(.fit_boot_failure("error"))

  align <- function(tab, key, target, stats) {
    if (is.null(tab)) return(NULL)
    j <- match(target, key)
    ans <- lapply(stats, function(st) {
      v <- rep(NA_real_, length(target)); ok <- !is.na(j)
      v[ok] <- tab[[st]][j[ok]]; v
    })
    names(ans) <- stats
    ans
  }
  pair_target <- .btl_pair_key(fit$pairs$object_a, fit$pairs$object_b)
  pair_key <- .btl_pair_key(z$pairs$object_a, z$pairs$object_b)
  list(
    total_chisq = z$total_chisq,
    pairs = align(z$pairs, pair_key, pair_target, c("chisq")),
    objects = align(z$objects, z$objects$object, fit$objects$object,
                    c("fit_resid", "infit_ms", "outfit_ms")),
    judges = if (is.null(fit$judges)) NULL else
      align(z$judges, z$judges$judge, fit$judges$judge,
            c("fit_resid", "infit_ms", "outfit_ms")))
}

.btl_boot_table <- function(base, reps, stats, side, id, min_success) {
  out <- base
  meta <- list()
  matrices <- list()
  for (st in stats) {
    M <- do.call(rbind, lapply(reps, function(z) z[[st]]))
    mode <- if (st == "chisq") "raw" else if (st == "fit_resid")
      "centred" else "studentised"
    trans <- if (st %in% c("infit_ms", "outfit_ms")) log else identity
    z <- .boot_maxt(base[[st]], M, side[[st]], min_success,
                    mode = mode, transform = trans)
    out[[paste0(st, "_p_boot")]] <- z$p
    out[[paste0(st, "_p_boot_adj")]] <- z$p_adj
    out[[paste0("n_boot_", st)]] <- z$n_boot
    matrices[[st]] <- M
    meta[[st]] <- list(family_n = z$family_n,
                       adjusted_n = z$adjusted_n,
                       joint_replicates = z$family_boot)
  }
  list(table = out, adjustment = meta, replicates = matrices)
}

.btl_fit_bootstrap <- function(fit, B, workers, seed) {
  if (inherits(fit, "rasch_btl_efrm"))
    .refuse("the paired-comparison fit bootstrap does not yet generate the ",
            "frame-dependent unit structure; use it with the single-frame ",
            "paired-comparison fit")
  if (!isTRUE(fit$converged))
    .refuse("the observed paired-comparison fit did not converge; refit it ",
            "successfully before bootstrapping fit")
  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) || B != floor(B) ||
      B < 1 || B > .Machine$integer.max)
    stop("`B` must be one whole positive number of replicates")
  if (length(workers) != 1L || !is.numeric(workers) || !is.finite(workers) ||
      workers != floor(workers) || workers < 1 || workers > .Machine$integer.max)
    stop("`workers` must be one whole positive number of workers")
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
       seed < 0 || seed != floor(seed) || seed > .Machine$integer.max))
    stop("`seed` must be one non-negative whole number within the integer range")
  B <- as.integer(B)
  workers <- min(as.integer(workers), .rasch_available_workers())
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    old <- .sim_seed_capture(); on.exit(.sim_seed_restore(old), add = TRUE)
    set.seed(seed)
  }
  # Check the design before starting the workers and before spending B refits.
  # The check generates one candidate dataset; restore the stream so that
  # validation cannot change the requested bootstrap draws.
  validation_stream <- .sim_seed_capture()
  tryCatch(invisible(.btl_boot_data(fit)),
           finally = .sim_seed_restore(validation_stream))
  seeds <- sample.int(.Machine$integer.max, B)
  one <- function(b) {
    old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old_stream), add = TRUE)
    set.seed(seeds[b])
    .btl_boot_refit(fit)
  }
  raw <- .rasch_boot_apply(B, one, workers = workers,
                           label = "paired-comparison fit bootstrap")
  status <- vapply(raw, .fit_boot_status, "")
  keep <- status == "ok"
  B_used <- sum(keep); min_success <- .fit_min_boot_success(B)
  if (B_used < min_success)
    .fit_boot_refuse(
      "only ", B_used, " of ", B, " paired-comparison bootstrap replicates ",
      "were usable (", sum(status == "nonconverged"),
      " did not converge; ", sum(status == "error"),
      " otherwise failed); at least ", min_success,
      " are required for the bootstrap null",
      B = B, B_used = B_used,
      B_nonconverged = sum(status == "nonconverged"),
      B_errors = sum(status == "error"))
  reps <- raw[keep]
  if (B_used < .9 * B)
    warning(B - B_used, " of ", B, " paired-comparison bootstrap replicates ",
            "were unusable (", sum(status == "nonconverged"),
            " did not converge; ", sum(status == "error"),
            " otherwise failed); the null uses the remaining ", B_used,
            call. = FALSE)

  pair_reps <- lapply(reps, `[[`, "pairs")
  object_reps <- lapply(reps, `[[`, "objects")
  judge_reps <- if (is.null(fit$judges)) NULL else lapply(reps, `[[`, "judges")
  pairs <- .btl_boot_table(fit$pairs, pair_reps, "chisq",
                           c(chisq = "upper"), "pair", min_success)
  objects <- .btl_boot_table(
    fit$objects, object_reps, c("fit_resid", "infit_ms", "outfit_ms"),
    c(fit_resid = "two", infit_ms = "two", outfit_ms = "two"),
    "object", min_success)
  judges <- if (is.null(fit$judges)) NULL else .btl_boot_table(
    fit$judges, judge_reps, c("fit_resid", "infit_ms", "outfit_ms"),
    c(fit_resid = "two", infit_ms = "two", outfit_ms = "two"),
    "judge", min_success)
  tot <- vapply(reps, `[[`, 0, "total_chisq")
  total <- list(chisq = fit$total_chisq, df = fit$total_df, p = fit$total_p,
                chisq_p_boot = if (sum(is.finite(tot)) >= min_success)
                  .boot_p(fit$total_chisq, tot, "upper") else NA_real_,
                n_boot = sum(is.finite(tot)))
  adjustment <- list(
    method = paste("single-step maximum-statistic bootstrap within each",
                   "statistic and parameter family"),
    pairs = pairs$adjustment, objects = objects$adjustment,
    judges = if (is.null(judges)) NULL else judges$adjustment)
  flat_meta <- c(pairs$adjustment, objects$adjustment,
                 if (is.null(judges)) list() else judges$adjustment)
  incomplete <- vapply(flat_meta, function(z)
    z$family_n > z$adjusted_n ||
      (z$family_n > 0L && z$joint_replicates < min_success), logical(1))
  if (any(incomplete))
    warning("at least one paired-comparison fit family had too few joint ",
            "replicates or no replicated variation for every requested ",
            "adjusted probability; affected values are NA", call. = FALSE)
  .new_fit_bootstrap(list(
    model = "paired comparisons", pairs = pairs$table,
    objects = objects$table,
    judges = if (is.null(judges)) NULL else judges$table,
    total = total, adjustment = adjustment,
    replicates = list(total_chisq = tot, pairs = pairs$replicates,
                      objects = objects$replicates,
                      judges = if (is.null(judges)) NULL else judges$replicates),
    B = B, B_used = B_used, B_failed = B - B_used,
    B_nonconverged = sum(status == "nonconverged"),
    B_errors = sum(status == "error"),
    minimum_usable = min_success, theta = NULL), fit, "btl")
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
.fit_theta_source <- function(theta, X, scheme) {
  if (length(theta) != nrow(X))
    .refuse("the fit does not contain one person location per response row; ",
            "the item fit bootstrap is unavailable")
  finite <- theta[is.finite(theta)]
  if (scheme == "conditional") return(NULL)
  if (scheme %in% c("normal", "resample")) {
    if (length(finite) < 2L)
      .refuse("the fit has fewer than two finite person locations to generate ",
              "from; the item fit bootstrap is unavailable")
    return(finite)
  }
  observed <- rowSums(!is.na(X)) > 0L
  if (any(!is.finite(theta) & observed))
    .refuse("`theta = \"fixed\"` requires a finite person location for ",
            "every row with an observed response")
  theta[!is.finite(theta)] <- if (length(finite)) mean(finite) else 0
  theta
}

.fit_theta <- function(scheme, theta, psi, n = length(theta)) {
  switch(scheme,
    normal = {
      vt <- psi$var_theta; mse <- psi$mean_error_var
      sd_c <- if (is.finite(vt) && is.finite(mse)) sqrt(max(vt - mse, 0)) else NA_real_
      if (!is.finite(sd_c) || sd_c <= 0) sample(theta, n, replace = TRUE)
      else stats::rnorm(n, mean(theta), sd_c)
    },
    resample = sample(theta, n, replace = TRUE),
    fixed = {
      if (length(theta) != n)
        stop("fixed person locations must have one value per response row")
      theta
    })
}

#' Bootstrap fit statistics
#'
#' Refers fit statistics to replicated datasets fitted in the same way as the
#' observed data. For a person-by-item Rasch model, the default generator
#' conditions on each person's observed raw score and missingness pattern.
#' The person parameter then cancels by sufficiency. Item parameters and
#' person locations are re-estimated in every replicate.
#'
#' Item chi-squares use the upper tail. Fit residuals, infit and outfit use
#' equal-tailed probabilities. Holm adjustment is applied separately to each
#' predeclared item-statistic family; an unavailable item remains in the
#' multiplicity count. Under the conditional generator, the same replicates
#' also give person-specific null distributions. Person probabilities are
#' adjusted with a single-step maximum-statistic distribution across persons
#' for each statistic. A maximum-statistic adjustment is withheld for the
#' complete family if any testable member lacks a usable joint null.
#'
#' The adjustments describe the fitted global null. They do not guarantee
#' familywise error control among otherwise fitting items, persons, objects or
#' judges when another member departs from the model. Each fit statistic is a
#' separate family. The marginal probabilities are retained.
#'
#' For a \code{\link{btl}} fit, outcomes are generated on the fitted comparison
#' design and the model is refitted. The result covers the total pairwise
#' chi-square and pair, object and judge fit. Ordered response thresholds,
#' judge allocation, counts, anchors, position effects and explanatory object
#' restrictions are retained. History-dependent effects are generated in
#' sequence. Half-weighted ties and frame-dependent paired-comparison fits are
#' refused because the present generator does not reproduce those models.
#'
#' Bootstrap probabilities use \eqn{(1+r)/(1+B)}. Their one-sided resolution
#' is therefore \eqn{1/(1+B)} and their equal-tailed resolution is
#' \eqn{2/(1+B)}. For \eqn{L} items, Holm-adjusted inference at .05 requires
#' at least \eqn{20L} replicates for the chi-square and \eqn{40L} for the
#' two-sided statistics.
#'
#' @param fit A fitted object from \code{\link{rasch}} or \code{\link{btl}}.
#'   Extended-frame, many-facet and explanatory person-by-item fits are not
#'   supported. Explanatory paired-comparison fits are supported.
#' @param B Positive whole number of bootstrap replicates.
#' @param theta Generator for a person-by-item fit. \code{"conditional"}
#'   retains each observed raw score and missingness pattern. \code{"resample"}
#'   resamples estimated person locations, \code{"fixed"} reuses each
#'   location with the same response row and missingness pattern, and
#'   \code{"normal"} draws from a normal distribution with error-corrected
#'   variance. Fixed generation requires a finite location for each row with
#'   an observed response. This argument does not apply to paired comparisons.
#' @param workers Number of parallel bootstrap workers. The default is four,
#'   reduced when fewer physical cores are available or the R process has a
#'   lower system limit. Per-replicate seeds are fixed before distribution,
#'   so results do not depend on the worker count.
#' @param seed Optional non-negative whole-number seed within the integer
#'   range. The caller's random stream is restored on exit.
#' @return An object of class \code{rasch_fit_bootstrap}. For a person-by-item
#'   fit, it contains \code{items},
#'   \code{persons}, \code{total}, \code{replicates}, adjustment metadata and
#'   replicate counts, including separate non-convergence and other-failure
#'   counts. For a paired-comparison fit, the corresponding tables are
#'   \code{pairs}, \code{objects}, \code{judges} and \code{total}.
#' @references Andrich, D. and Marais, I. (2019) \emph{A Course in Rasch
#'   Measurement Theory}. Springer.
#' Molenaar, I. W. and Hoijtink, H. (1996). Person-fit test statistics for the
#'   Rasch model. Applied Measurement in Education, 9, 87--106.
#'
#' Westfall, P. H. and Young, S. S. (1993). \emph{Resampling-Based Multiple
#'   Testing}. Wiley.
#' @seealso \code{\link{chisq_detail}} and \code{\link{btl}}.
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
  if (inherits(fit, "rasch_btl")) {
    if (!missing(theta))
      stop("`theta` applies to person-by-item Rasch fits, not paired comparisons")
    return(.btl_fit_bootstrap(fit, B = B, workers = workers, seed = seed))
  }
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
  if (!isTRUE(fit$est$converged))
    .refuse("the observed Rasch fit did not converge; refit it successfully ",
            "before bootstrapping fit")
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
    if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
        seed < 0 || seed != floor(seed) || seed > .Machine$integer.max)
      stop("`seed` must be one non-negative whole number within the integer range")
    seed <- as.integer(seed)
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
  th_full <- fit$person$theta
  na_mask <- if (anyNA(X)) is.na(X) else NULL
  # A fixed-location replicate must keep row n paired with row n's observed
  # missingness pattern. Non-finite locations are harmless only for entirely
  # unobserved rows, whose generated values are all masked again.
  th <- .fit_theta_source(th_full, X, theta)

  # every random draw a replicate depends on is made here, in the parent:
  # its person locations (for the ability-sampling schemes), and the seed its
  # responses are generated under. A worker therefore reproduces its
  # replicate from what it is handed, and the worker count cannot move the
  # result.
  th_b <- if (theta == "conditional") NULL
          else lapply(seq_len(B), function(b)
            .fit_theta(theta, th, fit$psi, n = nrow(X)))
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

  status <- vapply(reps, .fit_boot_status, "")
  keep <- status == "ok"
  B_used <- sum(keep)
  min_success <- .fit_min_boot_success(B)
  if (B_used < min_success)
    .fit_boot_refuse(
      "only ", B_used, " of ", B, " bootstrap replicates were usable (",
      sum(status == "nonconverged"), " did not converge; ",
      sum(status == "error"), " otherwise failed); at least ", min_success,
      " are required for the bootstrap null. The fitted model generates data ",
      "this estimator cannot fit reliably",
      B = B, B_used = B_used,
      B_nonconverged = sum(status == "nonconverged"),
      B_errors = sum(status == "error"))
  reps <- reps[keep]
  # replicates are not lost at random -- an unvisited extreme category is
  # likelier in the less dispersed ones -- so a materially thinned null is
  # worth saying out loud rather than leaving to be read off B_used
  if (B_used < 0.9 * B)
    warning(B - B_used, " of ", B, " bootstrap replicates were unusable (",
            sum(status == "nonconverged"), " did not converge; ",
            sum(status == "error"), " otherwise failed); the null is formed ",
            "from the remaining ", B_used,
            "; replicates are not lost at random, so read the result with ",
            "that in mind", call. = FALSE)

  M <- lapply(names(.FIT_STATS), function(st)
    do.call(rbind, lapply(reps, `[[`, st)))
  names(M) <- names(.FIT_STATS)

  items <- data.frame(item = fit$items$item, stringsAsFactors = FALSE)
  resolution_lost <- character(0)
  insufficient <- character(0)
  for (st in names(.FIT_STATS)) {
    side <- .FIT_STATS[[st]]
    n_stat <- colSums(is.finite(M[[st]]))
    p <- vapply(seq_len(L), function(i) {
      if (n_stat[i] < min_success) return(NA_real_)
      .boot_p(obs[[st]][i], M[[st]][, i], side)
    }, 0)
    usable <- is.finite(p)
    # Keep the predeclared family intact when one item's bootstrap null is
    # unavailable. Dropping that item from Holm's family would make the
    # remaining adjusted probabilities less conservative because a failed
    # test happened to be unreportable.
    family_n <- L
    if (any(is.finite(obs[[st]]) & n_stat < min_success))
      insufficient <- c(insufficient, st)
    adj <- rep(NA_real_, L)
    adj[usable] <- stats::p.adjust(p[usable], method = "holm", n = family_n)
    items[[st]] <- obs[[st]]
    items[[paste0(st, "_p_boot")]] <- p
    items[[paste0(st, "_p_boot_adj")]] <- adj
    items[[paste0("n_boot_", st)]] <- n_stat

    # The pre-run warning uses requested B.  Recheck the resolution from the
    # usable null for each statistic, since failures or statistic-specific NAs
    # can make an initially attainable Holm threshold unattainable.
    mult <- if (side == "upper") 1 else 2
    requested_possible <- mult * family_n / (B + 1) < 0.05
    actual_possible <- any(usable) &&
      family_n * min(mult / (n_stat[usable] + 1)) < 0.05
    if (requested_possible && !actual_possible)
      resolution_lost <- c(resolution_lost, st)
  }
  # Backward-compatible alias for the original chi-square count.
  items$n_boot <- items$n_boot_chisq
  rownames(items) <- NULL
  if (length(insufficient))
    warning("too few usable replicated statistics for ",
            paste(unique(insufficient), collapse = ", "),
            " in at least one item; the affected bootstrap probabilities are NA",
            call. = FALSE)
  if (length(resolution_lost))
    warning("after unusable replicates were removed, Holm-adjusted significance ",
            "at .05 is no longer resolvable for ",
            paste(unique(resolution_lost), collapse = ", "),
            "; increase `B`", call. = FALSE)

  # Person fit has a different nuisance structure from item fit. Under the
  # score-conditional generator, row n keeps its observed raw score and
  # missingness pattern in every replicate, so its replicated statistics are
  # the relevant null for that same response pattern. A joint maximum-
  # statistic reference retains the dependence among people in every replicate.
  persons <- NULL
  person_adjustment <- NULL
  if (theta == "conditional") {
    keep_cols <- intersect(c("id", "raw", "theta", "se"), names(fit$person))
    persons <- fit$person[, keep_cols, drop = FALSE]
    person_adjustment <- list(method = "single-step maximum-statistic bootstrap",
                              statistics = list())
    person_warnings <- character(0)
    for (st in names(.PERSON_FIT_STATS)) {
      Mp <- do.call(rbind, lapply(reps, function(z) z$persons[[st]]))
      mode <- if (st %in% c("fit_resid", "infit_z", "outfit_z"))
        "centred" else "studentised"
      trans <- if (st %in% c("infit_ms", "outfit_ms")) log else identity
      z <- .boot_maxt(fit$person[[st]], Mp, .PERSON_FIT_STATS[[st]],
                      min_success, mode = mode, transform = trans)
      persons[[st]] <- fit$person[[st]]
      persons[[paste0(st, "_p_boot")]] <- z$p
      persons[[paste0(st, "_p_boot_adj")]] <- z$p_adj
      persons[[paste0("n_boot_", st)]] <- z$n_boot
      person_adjustment$statistics[[st]] <- list(
        family_n = z$family_n, adjusted_n = z$adjusted_n,
        joint_replicates = z$family_boot)
      if (z$family_n > z$adjusted_n ||
          (z$family_n > 0L && z$family_boot < min_success))
        person_warnings <- c(person_warnings, st)
    }
    rownames(persons) <- NULL
    if (length(person_warnings))
      warning("too few joint person-fit bootstrap replicates, or no replicated ",
              "variation, for the maxT adjustment of ",
              paste(unique(person_warnings), collapse = ", "),
              "; marginal probabilities remain available but affected ",
              "adjusted probabilities are NA", call. = FALSE)
  }

  # whole-test readings: the summed chi-square, and the mean and SD of the
  # item fit residuals that the fit summary reports. The SD is the reason
  # this one matters -- the convention reads it against 1, and under a true
  # model it is nowhere near 1 until the sample runs to a few thousand.
  ok <- is.finite(obs$chisq)
  tot_rep <- rep(NA_real_, B_used)
  if (any(ok)) {
    complete_chisq <- rowSums(is.finite(M$chisq[, ok, drop = FALSE])) == sum(ok)
    tot_rep[complete_chisq] <- rowSums(
      M$chisq[complete_chisq, ok, drop = FALSE])
  }
  # Use the same observed item set in every replicate. Allowing na.rm to
  # change the set from row to row compares means and SDs with different
  # composition and therefore does not form one bootstrap null statistic.
  fr_ok <- is.finite(obs$fit_resid)
  fr_mean <- fr_sd <- rep(NA_real_, B_used)
  if (any(fr_ok)) {
    complete_fr <- rowSums(
      is.finite(M$fit_resid[, fr_ok, drop = FALSE])) == sum(fr_ok)
    fr_mean[complete_fr] <- rowMeans(
      M$fit_resid[complete_fr, fr_ok, drop = FALSE])
    if (sum(fr_ok) >= 2L)
      fr_sd[complete_fr] <- apply(
        M$fit_resid[complete_fr, fr_ok, drop = FALSE], 1, stats::sd)
  }
  total_p <- function(obs_value, null, side) {
    if (sum(is.finite(null)) < min_success) return(NA_real_)
    .boot_p(obs_value, null, side)
  }
  total <- list(
    chisq = fit$total_chisq, df = fit$total_df, p = fit$total_chisq_p,
    chisq_p_boot = total_p(fit$total_chisq, tot_rep, "upper"),
    fit_resid_mean = fit$item_fit_summary$mean,
    fit_resid_mean_p_boot = total_p(
      fit$item_fit_summary$mean, fr_mean, "two"),
    fit_resid_sd = fit$item_fit_summary$sd,
    fit_resid_sd_p_boot = total_p(fit$item_fit_summary$sd, fr_sd, "two"),
    n_boot = sum(is.finite(tot_rep)),
    n_boot_fit_resid_mean = sum(is.finite(fr_mean)),
    n_boot_fit_resid_sd = sum(is.finite(fr_sd)))
  total_counts <- c(chisq = total$n_boot,
                    fit_resid_mean = total$n_boot_fit_resid_mean,
                    fit_resid_sd = total$n_boot_fit_resid_sd)
  if (any(total_counts < min_success))
    warning("too few usable replicates for whole-test ",
            paste(names(total_counts)[total_counts < min_success],
                  collapse = ", "),
            "; the affected bootstrap probabilities are NA", call. = FALSE)

  .new_fit_bootstrap(list(
    items = items, total = total, replicates = M,
    persons = persons, person_adjustment = person_adjustment,
    B = B, B_used = B_used, B_failed = B - B_used,
    B_nonconverged = sum(status == "nonconverged"),
    B_errors = sum(status == "error"),
    minimum_usable = min_success, theta = theta), fit, "rasch")
}
