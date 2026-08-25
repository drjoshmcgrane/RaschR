# rasch :: explanatory item and threshold models
# =============================================================================
# The LLTM and LPCM retain the Rasch response function and replace freely
# estimated thresholds by linear functions of observed item or threshold
# characteristics. Estimation uses the same pairwise conditional likelihood as
# pcml(); only its threshold design matrix changes. Fixed departures are added
# as nominated design columns, so relaxation remains a conditional Rasch model.
# =============================================================================

.explanatory_projector <- function(m, thr) {
  L <- length(m); M <- nrow(thr)
  a <- 1 / (L * m[thr$item])
  diag(M) - matrix(1, M, 1L) %*% t(a)
}

.explanatory_ordinal_contrasts <- function(x) {
  k <- nlevels(x)
  if (k < 2L) stop("an ordinal predictor needs at least two observed levels")
  out <- outer(seq_len(k), seq_len(k - 1L), `>`) * 1
  colnames(out) <- paste0(make.names(levels(x)[-1L]), "_vs_",
                          make.names(levels(x)[-k]))
  out
}

.explanatory_metadata <- function(predictors, formula, X,
                                  level = c("item", "threshold")) {
  level <- match.arg(level)
  if (!is.data.frame(predictors))
    stop("`predictors` must be a data frame")
  if (!"item" %in% names(predictors))
    stop("`predictors` needs an `item` column")
  if (!inherits(formula, "formula") || length(formula) != 2L)
    stop("`formula` must be one-sided, for example ~ format + threshold")

  items <- colnames(X); m <- apply(X, 2L, max, na.rm = TRUE)
  thr <- threshold_index(m)
  index <- data.frame(item = items[thr$item],
                      threshold_number = thr$k,
                      stringsAsFactors = FALSE)
  predictors$item <- as.character(predictors$item)
  unknown <- setdiff(unique(predictors$item), items)
  if (length(unknown))
    stop("predictor item(s) are not present in the fitted response data: ",
         paste(unknown, collapse = ", "))

  if (level == "item") {
    if (anyDuplicated(predictors$item))
      stop("item-level predictors need exactly one row per item")
    missing <- setdiff(items, predictors$item)
    if (length(missing))
      stop("item-level predictors are missing: ", paste(missing, collapse = ", "))
    meta <- predictors[match(index$item, predictors$item), , drop = FALSE]
    meta$threshold_number <- index$threshold_number
  } else {
    if (!"threshold" %in% names(predictors) &&
        !"threshold_number" %in% names(predictors))
      stop("threshold-level predictors need a `threshold` column")
    kn <- if ("threshold" %in% names(predictors))
      predictors$threshold else predictors$threshold_number
    kn_num <- suppressWarnings(as.numeric(as.character(kn)))
    if (anyNA(kn_num) || any(!is.finite(kn_num)) ||
        any(kn_num != floor(kn_num)) || any(kn_num < 1))
      stop("threshold numbers must be positive integers")
    kn <- as.integer(kn_num)
    key <- paste(predictors$item, kn, sep = "\r")
    if (anyDuplicated(key))
      stop("threshold-level predictors need one row per item and threshold")
    wanted <- paste(index$item, index$threshold_number, sep = "\r")
    missing <- which(!wanted %in% key)
    extra <- setdiff(key, wanted)
    if (length(missing))
      stop("threshold-level predictors are missing: ",
           paste(paste0(index$item[missing], " threshold ",
                        index$threshold_number[missing]), collapse = ", "))
    if (length(extra))
      stop("predictor rows do not correspond to an observed item threshold: ",
           paste(gsub("\r", " threshold ", extra, fixed = TRUE),
                 collapse = ", "))
    meta <- predictors[match(wanted, key), , drop = FALSE]
    meta$threshold_number <- index$threshold_number
  }
  meta$item <- index$item
  meta$threshold <- factor(index$threshold_number,
                           levels = sort(unique(index$threshold_number)))
  rownames(meta) <- NULL

  reserved <- c("item", "threshold", "threshold_number")
  for (nm in setdiff(names(meta), reserved)) {
    if (is.character(meta[[nm]]) || is.logical(meta[[nm]]))
      meta[[nm]] <- factor(meta[[nm]])
  }
  model_meta <- meta
  for (nm in setdiff(names(model_meta), reserved))
    if (is.ordered(model_meta[[nm]]))
      contrasts(model_meta[[nm]]) <-
        .explanatory_ordinal_contrasts(model_meta[[nm]])
  mf <- tryCatch(stats::model.frame(formula, data = model_meta,
                                    na.action = stats::na.fail),
                 error = function(e) stop("cannot construct the explanatory ",
                   "model: ", conditionMessage(e), call. = FALSE))
  mm <- tryCatch(stats::model.matrix(formula, data = mf),
                 error = function(e) stop("cannot construct the explanatory ",
                   "model matrix: ", conditionMessage(e), call. = FALSE))
  if (!ncol(mm)) stop("the explanatory formula produced no predictors")
  if (any(!is.finite(mm)))
    stop("the explanatory predictors produce non-finite model-matrix values")

  A <- .explanatory_projector(m, thr)
  B0 <- A %*% mm
  norms <- sqrt(colSums(B0^2))
  zero <- norms < 1e-10
  if (any(zero & colnames(mm) != "(Intercept)"))
    stop("predictor(s) have no estimable variation after fixing the scale ",
         "origin: ", paste(colnames(mm)[zero &
           colnames(mm) != "(Intercept)"], collapse = ", "))
  B <- B0[, !zero, drop = FALSE]
  mm_keep <- mm[, !zero, drop = FALSE]
  if (!ncol(B)) stop("the explanatory formula contains only an intercept")
  qb <- qr(B, tol = 1e-10)
  if (qb$rank < ncol(B)) {
    aliased <- colnames(B)[qb$pivot[seq.int(qb$rank + 1L, ncol(B))]]
    stop("the explanatory design is not identified; aliased term(s): ",
         paste(aliased, collapse = ", "))
  }
  if (ncol(B) > nrow(B) - 1L)
    stop("the explanatory design has more parameters than the free calibration")
  colnames(B) <- colnames(mm_keep)
  list(B = B, matrix = mm_keep, metadata = meta,
       source_predictors = predictors, threshold_index = thr,
       m = m, level = level, formula = formula, projector = A)
}

.pcml_design <- function(X, B, parameter_names = colnames(B), maxit = 60,
                         tol = 1e-8) {
  X <- as.matrix(X); .check_integer_scores(X, "the score matrix")
  storage.mode(X) <- "integer"
  m <- apply(X, 2L, max, na.rm = TRUE); L <- ncol(X)
  thr <- threshold_index(m); M <- nrow(thr)
  inames <- colnames(X) %||% paste0("V", seq_len(L))
  if (!is.matrix(B) || nrow(B) != M)
    stop("the explanatory design must have one row per fitted threshold")
  if (!ncol(B) || qr(B, tol = 1e-10)$rank < ncol(B))
    stop("the explanatory design matrix is not full column rank")
  pairs <- .pair_counts(X, m)
  .pcml_check_connected(pairs, L, inames)
  weak <- .pcml_weak_thresholds(X, m, thr, inames)
  st <- .start_tau(X, thr)
  beta0 <- tryCatch(qr.solve(B, st, tol = 1e-10), error = function(e)
    rep(0, ncol(B)))
  beta0[!is.finite(beta0)] <- 0
  sol <- .pcml_solve(X, thr, m, B, beta0, maxit = maxit, tol = tol,
                     pairs = pairs)
  names(sol$beta) <- parameter_names
  dimnames(sol$cov_beta) <- list(parameter_names, parameter_names)
  thr$tau <- sol$tau
  thr$se <- sol$se_tau
  thr$anchored <- FALSE
  thr$weak <- weak$flag
  thr$se[thr$weak] <- NA_real_
  se <- sqrt(pmax(diag(sol$cov_beta), 0))
  z <- sol$beta / se
  coef <- data.frame(term = parameter_names, estimate = sol$beta, se = se,
                     z = z, p = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
                     stringsAsFactors = FALSE)
  coef$p_adj <- stats::p.adjust(coef$p, method = "holm")
  rownames(coef) <- coef$term
  list(model = "explanatory", thr = thr, cov_tau = sol$cov_tau,
       loglik = sol$loglik, iterations = sol$iterations,
       converged = sol$converged, m = m, anchors = NULL,
       n_parameters = ncol(B), B = B, beta = sol$beta,
       coefficients = coef, cov_beta = sol$cov_beta,
       H_beta = sol$H_beta, notes = weak$notes)
}

.pcml_nested_test <- function(full, restricted) {
  if (!isTRUE(full$converged) || !isTRUE(restricted$converged))
    stop("both conditional calibrations must converge before comparison")
  Bf <- full$B; Br <- restricted$B
  if (nrow(Bf) != nrow(Br))
    stop("the compared models do not describe the same thresholds")
  M <- nrow(Bf)
  S <- cbind(Br, rep(1, M))
  ss <- svd(S)
  rs <- sum(ss$d > max(1e-10, max(ss$d) * 1e-8))
  U <- ss$u[, seq_len(rs), drop = FALSE]
  Portho <- diag(M) - tcrossprod(U)
  A <- Portho %*% Bf
  sa <- svd(t(A))
  r <- sum(sa$d > max(1e-8, max(sa$d) * 1e-8))
  W <- max(0, 2 * (full$loglik - restricted$loglik))
  if (!r) return(list(chisq = W, df = 0L, p = NA_real_,
                      chisq_kent = NA_real_, p_kent = NA_real_,
                      lambda = numeric(0)))
  C <- sa$u[, seq_len(r), drop = FALSE]
  Hinv <- solve(-full$H_beta)
  num <- crossprod(C, full$cov_beta %*% C)
  den <- crossprod(C, Hinv %*% C)
  lambda <- Re(eigen(solve(den, num), only.values = TRUE)$values)
  kent <- W * r / sum(lambda)
  list(chisq = W, df = r, p = stats::pchisq(W, r, lower.tail = FALSE),
       chisq_kent = kent,
       p_kent = stats::pchisq(kent, r, lower.tail = FALSE),
       lambda = lambda)
}

.btl_explanatory_design <- function(predictors, formula, objects) {
  if (!is.data.frame(predictors) || !"object" %in% names(predictors))
    stop("`predictors` must be a data frame with an `object` column")
  if (!inherits(formula, "formula") || length(formula) != 2L)
    stop("`formula` must be one-sided, for example ~ domain + format")
  predictors$object <- as.character(predictors$object)
  if (anyDuplicated(predictors$object))
    stop("object predictors need exactly one row per object")
  missing <- setdiff(objects, predictors$object)
  if (length(missing))
    stop("object predictors are missing: ", paste(missing, collapse = ", "))
  extra <- setdiff(predictors$object, objects)
  if (length(extra))
    stop("predictor row(s) for object(s) not present in the comparisons: ",
         paste(extra, collapse = ", "))
  meta <- predictors[match(objects, predictors$object), , drop = FALSE]
  for (nm in setdiff(names(meta), "object"))
    if (is.character(meta[[nm]]) || is.logical(meta[[nm]]))
      meta[[nm]] <- factor(meta[[nm]])
  model_meta <- meta
  for (nm in setdiff(names(model_meta), "object"))
    if (is.ordered(model_meta[[nm]]))
      contrasts(model_meta[[nm]]) <-
        .explanatory_ordinal_contrasts(model_meta[[nm]])
  mf <- tryCatch(stats::model.frame(formula, data = model_meta,
                                    na.action = stats::na.fail),
    error = function(e) stop("cannot construct the explanatory object model: ",
                             conditionMessage(e), call. = FALSE))
  mm <- tryCatch(stats::model.matrix(formula, data = mf),
    error = function(e) stop("cannot construct the explanatory object model ",
                             "matrix: ", conditionMessage(e), call. = FALSE))
  if (any(!is.finite(mm)))
    stop("the explanatory object predictors produce non-finite model-matrix values")
  C <- diag(length(objects)) - 1 / length(objects)
  B0 <- C %*% mm
  keep <- sqrt(colSums(B0^2)) >= 1e-10
  bad <- colnames(mm)[!keep & colnames(mm) != "(Intercept)"]
  if (length(bad))
    stop("predictor(s) have no estimable variation after fixing the scale origin: ",
         paste(bad, collapse = ", "))
  B <- B0[, keep, drop = FALSE]
  if (!ncol(B)) stop("the explanatory formula contains only an intercept")
  if (qr(B, tol = 1e-10)$rank < ncol(B))
    stop("the explanatory object design is not identified")
  rownames(B) <- objects
  list(B = B, offset = stats::setNames(numeric(length(objects)), objects),
       parameter_names = colnames(B), metadata = meta,
       source_predictors = predictors, formula = formula)
}

.btl_explanatory_nested_test <- function(full, restricted) {
  if (!isTRUE(full$converged) || !isTRUE(restricted$converged))
    stop("both comparative judgement models must converge before comparison")
  if (!identical(as.character(full$objects$object),
                 as.character(restricted$objects$object)))
    stop("the compared models do not contain the same objects")
  Bf <- full$location_design; Br <- restricted$location_design
  K <- nrow(Bf)
  S <- cbind(Br, rep(1, K))
  ss <- svd(S); rs <- sum(ss$d > max(1e-10, max(ss$d) * 1e-8))
  P <- diag(K) - tcrossprod(ss$u[, seq_len(rs), drop = FALSE])
  A <- P %*% Bf
  sa <- svd(t(A)); r <- sum(sa$d > max(1e-8, max(sa$d) * 1e-8))
  W <- max(0, 2 * (full$loglik - restricted$loglik))
  if (!r) return(list(chisq = W, df = 0L, p = NA_real_,
                      chisq_kent = NA_real_, p_kent = NA_real_,
                      lambda = numeric(0)))
  Cobj <- sa$u[, seq_len(r), drop = FALSE]
  C <- rbind(Cobj,
             matrix(0, nrow(full$sensitivity) - nrow(Cobj), r))
  Hinv <- solve(full$sensitivity)
  num <- crossprod(C, full$cov_parameters %*% C)
  den <- crossprod(C, Hinv %*% C)
  lambda <- Re(eigen(solve(den, num), only.values = TRUE)$values)
  kent <- W * r / sum(lambda)
  available <- isTRUE(full$cl$inference_available)
  list(chisq = W, df = r,
       p = stats::pchisq(W, r, lower.tail = FALSE),
       chisq_kent = if (available) kent else NA_real_,
       p_kent = if (available)
         stats::pchisq(kent, r, lower.tail = FALSE) else NA_real_,
       lambda = if (available) lambda else rep(NA_real_, length(lambda)))
}

.btl_explanatory_refit <- function(fit, B, relaxations) {
  spec <- fit$explanatory$refit_spec
  design <- list(B = B,
                 offset = stats::setNames(numeric(nrow(B)), rownames(B)),
                 parameter_names = colnames(B))
  out <- do.call(btl, c(list(data = spec$data), spec$args,
                        list(.object_design = design)))
  out$reference_fit <- fit$reference_fit
  out$explanatory <- fit$explanatory
  out$explanatory$active_B <- B
  out$explanatory$relaxations <- relaxations
  class(out) <- c("rasch_btl_explanatory", "rasch_btl")
  out
}

#' Fit an explanatory comparative judgement model
#'
#' Constrains Bradley--Terry--Luce object locations to linear functions of
#' observed object characteristics. The formulation applies to dichotomous
#' and ordered comparative judgements; ordered-response thresholds retain the
#' structure selected in \code{thresholds}.
#'
#' @details For objects \eqn{a} and \eqn{b},
#' \deqn{\log\{P(a \succ b)/P(b \succ a)\}=\beta_a-\beta_b,\qquad
#' \beta_i=\mathbf z_i^{T}\boldsymbol\gamma.}
#' The scale origin is fixed at mean object location zero. Numeric predictors
#' are continuous, unordered factors are categorical, and ordered factors use
#' successive contrasts between adjacent levels. Character predictors are
#' converted to unordered factors. Selected
#' interactions may be included in \code{formula}. A free calibration is
#' retained for \code{explanatory_test()}. Standard errors use the same
#' sandwich covariance as \code{btl()}; when judges are identified, coefficient
#' tests use the judge-clustered covariance and a \eqn{t} reference with
#' judge-cluster degrees of freedom. Holm adjustment covers the coefficient
#' family.
#'
#' @inheritParams btl
#' @param predictors Data frame with one row per object, an \code{object}
#'   column, and the predictors named in \code{formula}.
#' @param formula One-sided explanatory formula, including selected
#'   interactions if required.
#' @return An object of class \code{"rasch_btl_explanatory"}, inheriting from
#'   \code{"rasch_btl"}.
#' @references Bradley, R. A. and Terry, M. E. (1952). Rank analysis of
#' incomplete block designs: I. The method of paired comparisons. Biometrika,
#' 39, 324--345.
#'
#' Fischer, G. H. (1973). The linear logistic test model as an instrument in
#' educational research. Acta Psychologica, 37, 359--374.
#' @examples
#' set.seed(1)
#' q <- data.frame(object = LETTERS[1:6],
#'                 domain = rep(0:1, each = 3))
#' beta <- setNames(0.8 * q$domain, q$object)
#' pr <- t(combn(q$object, 2))
#' d <- data.frame(a = rep(pr[, 1], each = 20),
#'                 b = rep(pr[, 2], each = 20))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$winner <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' fit <- btl_explanatory(d, q, ~ domain, "a", "b", winner = "winner")
#' fit$object_coefficients
#' explanatory_test(fit)
#' @seealso \code{\link{btl}}, \code{\link{explanatory_test}},
#'   \code{\link{explanatory_diagnostics}}, and
#'   \code{\link{relax_btl_explanatory}}.
#' @export
btl_explanatory <- function(data, predictors, formula, object_a, object_b,
                            winner = NULL, response = NULL, margin = NULL,
                            judge = NULL, count = NULL, order = NULL,
                            position = FALSE,
                            ties = c("drop", "half", "error"),
                            thresholds = c("free", "pc"), maxit = 60,
                            tol = 1e-8) {
  ties <- match.arg(ties); thresholds <- match.arg(thresholds)
  if (!is.data.frame(data) ||
      !all(c(object_a, object_b) %in% names(data)))
    stop("`object_a` and `object_b` must name columns in `data`")
  if (!is.data.frame(predictors) || !"object" %in% names(predictors))
    stop("`predictors` must be a data frame with an `object` column")
  observed_objects <- unique(c(trimws(as.character(data[[object_a]])),
                               trimws(as.character(data[[object_b]]))))
  observed_objects <- observed_objects[!is.na(observed_objects) &
                                         nzchar(observed_objects)]
  unknown <- setdiff(as.character(predictors$object), observed_objects)
  if (length(unknown))
    stop("predictor object(s) are not present in the comparison data: ",
         paste(unknown, collapse = ", "))
  args <- list(object_a = object_a, object_b = object_b, winner = winner,
               response = response, margin = margin, judge = judge,
               count = count, order = order, position = position,
               ties = ties, thresholds = thresholds, maxit = maxit, tol = tol)
  reference <- do.call(btl, c(list(data = data), args))
  design <- .btl_explanatory_design(predictors, formula,
                                    as.character(reference$objects$object))
  fit <- do.call(btl, c(list(data = data), args,
                        list(.object_design = design)))
  fit$reference_fit <- reference
  fit$explanatory <- list(
    formula = formula, formula_text = paste(deparse(formula), collapse = " "),
    metadata = design$metadata, source_predictors = design$source_predictors,
    base_B = design$B, active_B = fit$location_design,
    relaxations = data.frame(), refit_spec = list(data = data, args = args))
  class(fit) <- c("rasch_btl_explanatory", "rasch_btl")
  fit
}

#' Add a fixed object departure to an explanatory comparative judgement model
#'
#' @param fit A fitted object from \code{\link{btl_explanatory}}.
#' @param object Object name.
#' @return A refitted explanatory comparative judgement model.
#' @export
relax_btl_explanatory <- function(fit, object) {
  if (!inherits(fit, "rasch_btl_explanatory"))
    stop("relax_btl_explanatory() needs an explanatory comparative judgement fit")
  obj_tab <- fit$objects
  if ("extreme" %in% names(obj_tab)) {
    at <- match(object, obj_tab$object)
    if (!is.na(at) && isTRUE(obj_tab$extreme[at]))
      .refuse(object, " was set aside at a response boundary; its location ",
              "is an extrapolation for display and cannot be relaxed")
    obj_tab <- obj_tab[!(obj_tab$extreme %in% TRUE), ]
  }
  objects <- as.character(obj_tab$object)
  j <- match(object, objects)
  if (is.na(j)) stop("object not found in the explanatory fit: ", object)
  D <- diag(length(objects))[, j, drop = FALSE]
  D <- (diag(length(objects)) - 1 / length(objects)) %*% D
  rownames(D) <- objects; colnames(D) <- paste0("departure[", object, "]")
  B <- fit$explanatory$active_B
  if (qr(cbind(B, D), tol = 1e-10)$rank == qr(B, tol = 1e-10)$rank)
    stop("that object departure is already represented by the active model")
  rel <- fit$explanatory$relaxations
  rel <- rbind(rel, data.frame(order = nrow(rel) + 1L, object = object,
                               component = "Object location",
                               parameters_added = 1L,
                               stringsAsFactors = FALSE))
  .btl_explanatory_refit(fit, cbind(B, D), rel)
}

#' @export
print.rasch_btl_explanatory <- function(x, ...) {
  cat("Explanatory comparative judgement model\n")
  cat("Formula: ", x$explanatory$formula_text, "\n", sep = "")
  print.rasch_btl(x, ...)
}

.explanatory_attach <- function(out, reference, design, formula, level,
                                relaxations, n_groups_requested, adjust_N,
                                maxit, tol) {
  out$reference_fit <- reference
  out$mc <- reference$mc
  out$explanatory <- list(
    formula = formula, formula_text = paste(deparse(formula), collapse = " "),
    level = level, metadata = design$metadata,
    model_matrix = design$matrix, base_B = design$B,
    active_B = out$est$B, threshold_index = design$threshold_index,
    source_predictors = design$source_predictors,
    relaxations = relaxations)
  out$refit_spec <- list(model = "PCM", n_groups = n_groups_requested,
    adjust_N = adjust_N, anchors = NULL,
    na_codes = reference$refit_spec$na_codes %||% -1,
    key = reference$refit_spec$key,
    pc_components = NULL, maxit = maxit, tol = tol,
    explanatory = TRUE)
  class(out) <- c("rasch_explanatory", "rasch")
  out
}

#' Fit an explanatory Rasch model
#'
#' Fits the linear logistic test model (LLTM) for dichotomous responses or
#' the linear partial credit model (LPCM) for polytomous responses. Item or
#' threshold locations are linear functions of observed predictors. The
#' response model remains Rasch and is estimated by pairwise conditional
#' maximum likelihood.
#'
#' @details
#' For threshold \eqn{k} of item \eqn{i},
#' \deqn{\delta_{ik}=z_{ik}^{T}\gamma.}
#' The adjacent-category log odds are
#' \deqn{\log\{P(X_{ni}=k)/P(X_{ni}=k-1)\}=\theta_n-\delta_{ik}.}
#' The threshold origin is fixed to the same mean-item-location zero used by
#' \code{\link{rasch}}. An intercept therefore sets the arbitrary origin and
#' is not separately estimated. Numeric predictors are continuous, unordered
#' factors are categorical, and ordered factors use successive contrasts
#' between adjacent levels. Character predictors are converted to unordered
#' factors. The reserved factor
#' \code{threshold} identifies the within-item threshold number;
#' \code{threshold_number} supplies its integer value.
#'
#' A free PCM reference is fitted to the same prepared responses and retained
#' on the object. \code{\link{explanatory_test}} applies the first-order Kent
#' calibration required for the pairwise composite likelihood.
#'
#' @param data,items,id,factors,n_groups,adjust_N,na_codes,key,maxit,tol As in
#'   \code{\link{rasch}}.
#' @param predictors Data frame containing an \code{item} column and the
#'   predictors named in \code{formula}. With \code{level = "threshold"}, it
#'   must also contain \code{threshold}, with one row for every fitted item
#'   threshold.
#' @param formula One-sided explanatory formula. For example,
#'   \code{~ format + operation + format:operation}. The reserved
#'   \code{threshold} factor permits threshold-specific effects.
#' @param level Whether \code{predictors} contains one row per \code{"item"}
#'   or per \code{"threshold"}. Item rows are expanded over their thresholds.
#' @return An object of class \code{"rasch_explanatory"} inheriting from
#'   \code{"rasch"}. Standard item, person, fit and diagnostic components use
#'   the explanatory thresholds. The \code{explanatory} component contains the
#'   formula, metadata and design matrices; \code{reference_fit} is the free
#'   PCM calibration.
#' @references
#' Fischer, G. H. (1973). The linear logistic test model as an instrument in
#' educational research. Acta Psychologica, 37, 359--374.
#'
#' Fischer, G. H. and Ponocny, I. (1994). An extension of the partial credit
#' model with an application to the measurement of change. Psychometrika, 59,
#' 177--192.
#' @examples
#' set.seed(1)
#' q <- data.frame(item = paste0("I", 1:8),
#'                 operation = rep(0:1, each = 4),
#'                 format = rep(c("A", "B"), 4))
#' difficulty <- -1 + 0.7 * q$operation + 0.4 * (q$format == "B")
#' X <- matrix(rbinom(500 * 8, 1,
#'   plogis(outer(rnorm(500), difficulty, "-"))), 500, 8)
#' colnames(X) <- q$item
#' fit <- rasch_explanatory(X, predictors = q,
#'                          formula = ~ operation + format)
#' fit$est$coefficients
#' explanatory_test(fit)
#' @seealso \code{\link{explanatory_test}},
#'   \code{\link{explanatory_diagnostics}}, and
#'   \code{\link{relax_explanatory}}.
#' @export
rasch_explanatory <- function(data, predictors, formula, items = NULL,
                              level = c("item", "threshold"), id = NULL,
                              factors = NULL, n_groups = NULL, adjust_N = NA,
                              na_codes = -1, key = NULL, maxit = 60,
                              tol = 1e-8) {
  level <- match.arg(level)
  reference <- rasch(data, model = "PCM", id = id, factors = factors,
                     items = items, n_groups = n_groups, adjust_N = adjust_N,
                     na_codes = na_codes, key = key, maxit = maxit, tol = tol)
  design <- .explanatory_metadata(predictors, formula, reference$X, level)
  est <- .pcml_design(reference$X, design$B,
                      parameter_names = colnames(design$B),
                      maxit = maxit, tol = tol)
  if (!isTRUE(est$converged))
    warning("the explanatory calibration did not converge; estimates, ",
            "diagnostics and probabilities are unreliable", call. = FALSE)
  kind <- if (all(reference$m == 1L)) "LLTM" else "LPCM"
  notes <- c(reference$notes,
             sprintf("%s explanatory calibration: %s", kind,
                     paste(deparse(formula), collapse = " ")),
             est$notes)
  out <- .assemble_fit("PCM", reference$X, est, reference$person$id,
                       reference$factors, n_groups, adjust_N, notes)
  out$explanatory_model <- kind
  .explanatory_attach(out, reference, design, formula, level,
                      relaxations = data.frame(),
                      n_groups_requested = n_groups, adjust_N = adjust_N,
                      maxit = maxit, tol = tol)
}

#' Compare an explanatory model with its free calibration
#'
#' Tests explanatory item, threshold or object restrictions against the
#' corresponding free calibration of the same responses. The inferential
#' result uses the first-order Kent calibration for the fitted likelihood and
#' sandwich covariance. The calibration coefficient of determination is
#' \deqn{R^2_{cal}=1-\frac{\sum_j(\hat\eta^{free}_j-
#' \hat\eta^{expl}_j-\bar d)^2}{\sum_j(\hat\eta^{free}_j-
#' \bar\eta^{free})^2},}
#' where \eqn{\bar d} removes the arbitrary scale origin. It describes the
#' proportion of variation in the well-determined free threshold calibration
#' (Rasch models) or free object calibration (comparative judgement)
#' reproduced by the explanatory model. It is at most one and may be negative.
#' It is not adjusted for the number of predictors, so with few calibrated
#' parameters it reads above zero even for an uninformative design.
#' \code{r_squared_adj} divides the unexplained proportion by its share of
#' the degrees of freedom, \eqn{1-(1-R^2_{cal})(n-1)/\mathit{df}}, where
#' \eqn{n} counts the calibrated parameters compared and \eqn{\mathit{df}}
#' is \eqn{n} minus the rank of the retained explanatory design with its
#' origin, so exclusions that remove a level's only support reduce it. The
#' correction is exact for independent homoskedastic estimates fitted by
#' least squares, which these calibrations are not, so read it as a
#' descriptive optimism adjustment. Read either beside the test rather than
#' in place of it.
#'
#' @param fit A fitted explanatory Rasch or comparative judgement model.
#' @return A one-row data frame containing the raw and Kent-calibrated
#'   statistics, degrees of freedom and parameter counts. The primary
#'   \code{p} and the retained \code{p_kent} are the Kent-calibrated
#'   probability. \code{p_naive} is the unscaled composite-likelihood
#'   probability and is provided for methodological inspection, not
#'   inference. \code{r_squared} is the calibration coefficient of
#'   determination, \code{r_squared_adj} its degrees-of-freedom-adjusted
#'   counterpart, and \code{r2_basis} names the calibrated parameters used.
#' @export
explanatory_test <- function(fit) {
  if (inherits(fit, "rasch_btl_explanatory")) {
    z <- .btl_explanatory_nested_test(fit$reference_fit, fit)
    free <- fit$reference_fit$objects
    active <- fit$objects
    # an extrapolated boundary row is not a calibrated location; the design
    # rows span the calibrated objects only
    if ("extreme" %in% names(active)) active <- active[!(active$extreme %in% TRUE), ]
    if ("extreme" %in% names(free)) free <- free[!(free$extreme %in% TRUE), ]
    f <- free$location[match(active$object, free$object)]
    a <- active$location
    ok <- is.finite(f) & is.finite(a)
    d <- f[ok] - a[ok]
    den <- sum((f[ok] - mean(f[ok]))^2)
    r2 <- if (sum(ok) > 1L && is.finite(den) && den > 0)
      1 - sum((d - mean(d))^2) / den else NA_real_
    df_sub <- if (sum(ok) > 1L) sum(ok) -
      qr(cbind(1, fit$location_design[ok, , drop = FALSE]))$rank else 0L
    r2_adj <- if (is.finite(r2) && df_sub > 0)
      1 - (1 - r2) * (sum(ok) - 1) / df_sub else NA_real_
    out <- data.frame(
      model = if (nrow(fit$explanatory$relaxations))
        "Partially relaxed explanatory CJ model" else "Explanatory CJ",
      parameters = ncol(fit$location_design),
      free_parameters = ncol(fit$reference_fit$location_design),
      r_squared = r2, r_squared_adj = r2_adj,
      r2_basis = "object calibration",
      chisq = z$chisq, df = z$df, p_naive = z$p,
      chisq_kent = z$chisq_kent, p = z$p_kent, p_kent = z$p_kent,
      stringsAsFactors = FALSE)
    attr(out, "lambda") <- z$lambda
    return(.tag_tables(out))
  }
  if (!inherits(fit, "rasch_explanatory"))
    stop("explanatory_test() needs an explanatory Rasch fit")
  z <- .pcml_nested_test(fit$reference_fit$est, fit$est)
  free <- fit$reference_fit$est$thr$tau
  active <- fit$est$thr$tau
  ok <- is.finite(free) & is.finite(active)
  free_weak <- fit$reference_fit$est$thr$weak %||% rep(FALSE, length(free))
  active_weak <- fit$est$thr$weak %||% rep(FALSE, length(active))
  ok <- ok & !free_weak & !active_weak
  d <- free[ok] - active[ok]
  den <- sum((free[ok] - mean(free[ok]))^2)
  r2 <- if (sum(ok) > 1L && is.finite(den) && den > 0)
    1 - sum((d - mean(d))^2) / den else NA_real_
  df_sub <- if (sum(ok) > 1L)
    sum(ok) - qr(cbind(1, fit$est$B[ok, , drop = FALSE]))$rank else 0L
  r2_adj <- if (is.finite(r2) && df_sub > 0)
    1 - (1 - r2) * (sum(ok) - 1) / df_sub else NA_real_
  out <- data.frame(
    model = if (nrow(fit$explanatory$relaxations))
      "Partially relaxed explanatory model" else fit$explanatory_model,
    parameters = fit$est$n_parameters,
    free_parameters = fit$reference_fit$est$n_parameters,
    r_squared = r2, r_squared_adj = r2_adj,
    r2_basis = "threshold calibration",
    chisq = z$chisq, df = z$df, p_naive = z$p,
    chisq_kent = z$chisq_kent, p = z$p_kent, p_kent = z$p_kent,
    stringsAsFactors = FALSE)
  attr(out, "lambda") <- z$lambda
  .tag_tables(out)
}

.explanatory_candidate <- function(fit, item,
                                   component = c("location", "thresholds")) {
  component <- match.arg(component)
  j <- match(item, colnames(fit$X))
  if (is.na(j)) stop("item not found in the explanatory fit: ", item)
  thr <- fit$explanatory$threshold_index %||%
    threshold_index(fit$m)
  rows <- which(thr$item == j); mi <- length(rows); M <- nrow(thr)
  if (component == "location") {
    D <- matrix(0, M, 1L); D[rows, 1L] <- 1
    colnames(D) <- paste0("departure_location[", item, "]")
  } else {
    if (mi < 2L) stop("a dichotomous item has no threshold-structure block")
    D <- matrix(0, M, mi - 1L)
    D[rows, ] <- rbind(diag(mi - 1L), rep(-1, mi - 1L))
    colnames(D) <- paste0("departure_threshold[", item, ",",
                          seq_len(mi - 1L), "]")
  }
  A <- .explanatory_projector(fit$m, thr)
  A %*% D
}

.explanatory_addable <- function(B, D) {
  qr(cbind(B, D), tol = 1e-10)$rank - qr(B, tol = 1e-10)$rank
}

#' Diagnose fixed departures from an explanatory model
#'
#' Fits each available item-location, polytomous threshold-structure or
#' comparative-judgement object departure separately from the active model.
#' Probabilities use Kent calibration and Holm adjustment over the complete
#' candidate family.
#'
#' @param fit A fitted explanatory Rasch or comparative judgement model.
#' @param p_adjust Multiplicity adjustment over the candidate departures.
#' @return A data frame ordered by adjusted probability. For item fits, a
#'   \code{weak} column marks items whose thresholds the calibration flags
#'   as weakly identified; their probabilities are withheld, since the
#'   departure test rests on the same sparse categories, and a note on the
#'   table records the withholding.
#' @export
explanatory_diagnostics <- function(fit, p_adjust = "holm") {
  if (inherits(fit, "rasch_btl_explanatory")) {
    if (!p_adjust %in% stats::p.adjust.methods)
      stop("p_adjust must name a method in stats::p.adjust.methods")
    obj_tab <- fit$objects
    if ("extreme" %in% names(obj_tab)) obj_tab <- obj_tab[!(obj_tab$extreme %in% TRUE), ]
    objects <- as.character(obj_tab$object)
    B <- fit$explanatory$active_B
    rows <- list()
    for (j in seq_along(objects)) {
      D <- diag(length(objects))[, j, drop = FALSE]
      D <- (diag(length(objects)) - 1 / length(objects)) %*% D
      rownames(D) <- objects
      colnames(D) <- paste0("departure[", objects[j], "]")
      if (qr(cbind(B, D), tol = 1e-10)$rank == qr(B, tol = 1e-10)$rank)
        next
      cand <- .btl_explanatory_refit(fit, cbind(B, D),
                                     fit$explanatory$relaxations)
      tst <- .btl_explanatory_nested_test(cand, fit)
      rows[[length(rows) + 1L]] <- data.frame(
        object = objects[j], component = "Object location",
        parameters_added = 1L,
        departure = utils::tail(cand$object_coefficients$estimate, 1L),
        deviance_reduction = tst$chisq, df = tst$df, p = tst$p_kent,
        stringsAsFactors = FALSE)
    }
    if (!length(rows)) return(.tag_tables(data.frame(
      object = character(0), component = character(0),
      parameters_added = integer(0), departure = numeric(0),
      deviance_reduction = numeric(0), df = integer(0), p = numeric(0),
      p_adj = numeric(0))))
    out <- do.call(rbind, rows)
    out$p_adj <- stats::p.adjust(out$p, method = p_adjust)
    out <- out[order(out$p_adj, -out$deviance_reduction), , drop = FALSE]
    rownames(out) <- NULL
    attr(out, "p_adjust") <- p_adjust
    return(.tag_tables(out))
  }
  if (!inherits(fit, "rasch_explanatory"))
    stop("explanatory_diagnostics() needs an explanatory Rasch fit")
  if (!p_adjust %in% stats::p.adjust.methods)
    stop("p_adjust must name a method in stats::p.adjust.methods")
  B <- fit$est$B; rows <- list()
  spec <- fit$refit_spec
  for (item in colnames(fit$X)) for (component in c("location", "thresholds")) {
    if (component == "thresholds" && fit$m[match(item, colnames(fit$X))] < 2L)
      next
    D <- .explanatory_candidate(fit, item, component)
    add <- .explanatory_addable(B, D)
    if (!add) next
    candB <- cbind(B, D)
    est <- .pcml_design(fit$X, candB, colnames(candB),
                        maxit = spec$maxit, tol = spec$tol)
    tst <- .pcml_nested_test(est, fit$est)
    b <- utils::tail(est$beta, ncol(D))
    departure <- if (component == "location") unname(b[1L]) else
      max(abs(drop(D %*% b)))
    # a departure test rests on the same sparse categories that made the
    # item's thresholds weak; the probability is withheld there, as the
    # threshold standard errors already are, and the departure stays
    # descriptive
    ii <- match(item, colnames(fit$X))
    weak_item <- isTRUE(any(fit$thresholds$weak[fit$thresholds$item == ii]))
    rows[[length(rows) + 1L]] <- data.frame(
      item = item, component = if (component == "location")
        "Item location" else "Threshold structure",
      parameters_added = add, departure = departure,
      deviance_reduction = tst$chisq, df = tst$df,
      p = if (weak_item) NA_real_ else tst$p_kent,
      weak = weak_item, stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(.tag_tables(data.frame(
    item = character(0), component = character(0), parameters_added = integer(0),
    departure = numeric(0), deviance_reduction = numeric(0), df = integer(0),
    p = numeric(0), weak = logical(0), p_adj = numeric(0))))
  out <- do.call(rbind, rows)
  usable <- is.finite(out$p)
  out$p_adj <- NA_real_
  out$p_adj[usable] <- stats::p.adjust(out$p[usable], method = p_adjust)
  if (any(out$weak))
    attr(out, "note") <- paste("departure probabilities are withheld for",
      "item(s) with weak thresholds:",
      paste(unique(out$item[out$weak]), collapse = ", "))
  out <- out[order(out$p_adj, -out$deviance_reduction), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "p_adjust") <- p_adjust
  .tag_tables(out)
}

#' Relax a nominated explanatory restriction
#'
#' Adds either one fixed item-location departure or a fixed block describing
#' an item's threshold structure, then repeats the complete conditional
#' calibration and downstream Rasch analysis. The departure is fixed rather
#' than random; raw-score sufficiency and the common discrimination remain.
#'
#' @param fit A fitted explanatory Rasch model.
#' @param item Item name.
#' @param component Either \code{"location"} or \code{"thresholds"}.
#' @return A partially relaxed \code{"rasch_explanatory"} fit.
#' @export
relax_explanatory <- function(fit, item,
                              component = c("location", "thresholds")) {
  if (!inherits(fit, "rasch_explanatory"))
    stop("relax_explanatory() needs an explanatory Rasch fit")
  component <- match.arg(component)
  D <- .explanatory_candidate(fit, item, component)
  add <- .explanatory_addable(fit$est$B, D)
  if (!add)
    stop("that departure is already represented by the active explanatory model")
  B <- cbind(fit$est$B, D)
  spec <- fit$refit_spec
  est <- .pcml_design(fit$X, B, colnames(B), maxit = spec$maxit,
                      tol = spec$tol)
  if (!isTRUE(est$converged))
    stop("the relaxed explanatory calibration did not converge")
  rel <- fit$explanatory$relaxations
  new <- data.frame(order = nrow(rel) + 1L, item = item,
                    component = if (component == "location")
                      "Item location" else "Threshold structure",
                    parameters_added = add, stringsAsFactors = FALSE)
  rel <- rbind(rel, new)
  out <- .assemble_fit("PCM", fit$X, est, fit$person$id, fit$factors,
                       fit$n_groups, spec$adjust_N,
                       c(fit$notes, sprintf("fixed explanatory departure: %s, %s",
                         item, tolower(new$component))))
  out$explanatory_model <- fit$explanatory_model
  design <- list(B = fit$explanatory$base_B,
                 matrix = fit$explanatory$model_matrix,
                 metadata = fit$explanatory$metadata,
                 source_predictors = fit$explanatory$source_predictors,
                 threshold_index = fit$explanatory$threshold_index)
  out <- .explanatory_attach(out, fit$reference_fit, design,
                      fit$explanatory$formula, fit$explanatory$level, rel,
                      n_groups_requested = fit$n_groups,
                      adjust_N = spec$adjust_N, maxit = spec$maxit,
                      tol = spec$tol)
  out$mc <- fit$mc
  out
}

.explanatory_inherit_mc <- function(fit, source, inherit,
                                    exclude = character(0)) {
  if (is.null(fit$mc) || is.null(fit$mc$raw)) return(NULL)
  items <- colnames(source)
  old <- unname(inherit[items])
  keep <- items[old %in% colnames(fit$mc$raw) & !items %in% exclude]
  if (!length(keep)) return(NULL)
  raw <- vapply(keep, function(it) {
    value <- fit$mc$raw[, unname(inherit[it])]
    value[is.na(source[, it])] <- NA_character_
    value
  }, character(nrow(source)))
  if (!is.matrix(raw)) raw <- matrix(raw, ncol = length(keep))
  colnames(raw) <- keep
  map <- lapply(keep, function(it) fit$mc$map[[unname(inherit[it])]])
  names(map) <- keep
  list(raw = raw, map = map,
       key = vapply(map, .key_label, character(1)))
}

.explanatory_refit_modified <- function(fit, source, inherit = NULL,
                                        location_relaxed = character(0),
                                        fully_relaxed = character(0)) {
  source <- as.matrix(source)
  new_items <- colnames(source)
  old_items <- colnames(fit$X)
  if (is.null(inherit))
    inherit <- stats::setNames(intersect(new_items, old_items),
                               intersect(new_items, old_items))
  if (is.null(names(inherit)) || any(!new_items %in% names(inherit)))
    stop("the explanatory refit needs an inherited predictor source for every item")
  if (any(!unname(inherit[new_items]) %in% old_items))
    stop("an inherited explanatory item is absent from the original fit")
  src <- fit$explanatory$source_predictors
  level <- fit$explanatory$level
  pred <- list()
  m_new <- apply(source, 2L, max, na.rm = TRUE)
  for (it in new_items) {
    old <- unname(inherit[it])
    z <- src[as.character(src$item) == old, , drop = FALSE]
    if (!nrow(z)) stop("predictors are unavailable for inherited item ", old)
    if (level == "item") {
      z <- z[1L, , drop = FALSE]; z$item <- it
    } else {
      kn <- if ("threshold" %in% names(z))
        as.integer(as.character(z$threshold)) else
          as.integer(as.character(z$threshold_number))
      take <- vapply(seq_len(m_new[it]), function(k) {
        hit <- which(kn == k)
        if (length(hit)) hit[1L] else which.max(kn)
      }, integer(1))
      z <- z[take, , drop = FALSE]; z$item <- it
      if ("threshold" %in% names(z)) z$threshold <- seq_len(m_new[it])
      if ("threshold_number" %in% names(z))
        z$threshold_number <- seq_len(m_new[it])
    }
    pred[[it]] <- z
  }
  pred <- do.call(rbind, pred); rownames(pred) <- NULL
  spec <- fit$refit_spec
  out <- rasch_explanatory(source, predictors = pred,
    formula = fit$explanatory$formula, level = level,
    id = fit$person$id, factors = fit$factors,
    n_groups = spec$n_groups %||% fit$n_groups,
    adjust_N = spec$adjust_N %||% NA_real_, maxit = spec$maxit %||% 60,
    tol = spec$tol %||% 1e-8)

  # Preserve prior analyst-approved departures on items that survive or are
  # replaced by inherited copies. A split copy receives the source item's
  # departure before its own location is freed.
  rel <- fit$explanatory$relaxations
  if (nrow(rel)) for (r in seq_len(nrow(rel))) {
    targets <- names(inherit)[unname(inherit) == rel$item[r]]
    for (it in targets) {
      component <- if (rel$component[r] == "Item location")
        "location" else "thresholds"
      D <- tryCatch(.explanatory_candidate(out, it, component),
                    error = function(e) NULL)
      if (!is.null(D) && .explanatory_addable(out$est$B, D) > 0L)
        out <- relax_explanatory(out, it, component)
    }
  }
  for (it in unique(location_relaxed)) {
    D <- .explanatory_candidate(out, it, "location")
    if (.explanatory_addable(out$est$B, D) > 0L)
      out <- relax_explanatory(out, it, "location")
  }
  for (it in unique(fully_relaxed)) {
    D <- .explanatory_candidate(out, it, "location")
    if (.explanatory_addable(out$est$B, D) > 0L)
      out <- relax_explanatory(out, it, "location")
    if (out$m[match(it, colnames(out$X))] > 1L) {
      D <- .explanatory_candidate(out, it, "thresholds")
      if (.explanatory_addable(out$est$B, D) > 0L)
        out <- relax_explanatory(out, it, "thresholds")
    }
  }
  out$mc <- .explanatory_inherit_mc(
    fit, source, inherit, exclude = unique(fully_relaxed))
  out
}

#' @export
print.rasch_explanatory <- function(x, ...) {
  cat(sprintf("rasch %s analysis: %d items, %d persons\n",
              x$explanatory_model, ncol(x$X), nrow(x$X)))
  cat("Formula: ", x$explanatory$formula_text, "\n", sep = "")
  cat(sprintf("Conditional calibration: %d explanatory parameter(s), %d fixed departure(s)\n",
              nrow(x$est$coefficients),
              nrow(x$explanatory$relaxations)))
  tst <- explanatory_test(x)
  if (tst$df > 0L)
    cat(sprintf("Free calibration comparison: adjusted chi-square %.3f on %d df, p = %s\n",
                tst$chisq_kent, tst$df, .fmt_p(tst$p_kent)))
  invisible(x)
}
