# rasch :: data simulation
#
# Generate data from the model family with dial-in departures from it, so a
# known pathology can be planted and the matching diagnostic watched as it
# fires. Every simulator returns data ready for its fit function, with the
# true generating parameters attached as attr(x, "truth") and a class that
# prints a summary of what was planted.
#
# Shared truth schema (all simulators populate what applies):
#   list(layout, n_*, theta/locations, difficulty/thresholds, discrimination,
#        guessing, groups, planted = <character, human-readable pathologies>)

# null-coalescing helper (package-internal; base R gained %||% only in 4.4)
`%||%` <- function(a, b) if (is.null(a)) b else a

# A seeded simulator call should be reproducible without commandeering the
# caller's random number stream: capture the stream, seed, and restore on
# exit, so simulate_*(seed = s) twice gives the same data while code after
# the call draws exactly what it would have drawn anyway. Restoring an
# absent stream means removing the one set.seed() created.
.sim_seed_capture <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  else NULL
}

.sim_seed_restore <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", old, envir = globalenv())
  }
  invisible(NULL)
}

.sim_count <- function(x, name, min = 1L) {
  if (length(x) != 1L || !is.numeric(x) || is.complex(x) ||
      !is.null(dim(x)) || !is.null(oldClass(x)) || !is.finite(x) ||
      x != floor(x) || x < min ||
      x > .Machine$integer.max)
    stop(name, " must be one whole number >= ", min)
  as.integer(x)
}

.sim_scalar <- function(x, name, lower = -Inf, upper = Inf,
                        lower_open = FALSE, upper_open = FALSE) {
  ok <- length(x) == 1L && is.numeric(x) && !is.complex(x) &&
    is.null(dim(x)) && is.null(oldClass(x)) && is.finite(x) &&
    if (lower_open) x > lower else x >= lower
  ok <- ok && if (upper_open) x < upper else x <= upper
  if (!ok) {
    left <- if (lower_open) "(" else "["
    right <- if (upper_open) ")" else "]"
    stop(name, " must be one finite value in ", left, lower, ", ", upper, right)
  }
  as.numeric(x)
}

.sim_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  .sim_count(seed, "seed", 0L)
}

.sim_vector <- function(x, name, lengths, lower = -Inf, upper = Inf,
                        lower_open = FALSE, upper_open = FALSE) {
  ok <- is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
    is.null(oldClass(x)) && length(x) %in% lengths && all(is.finite(x))
  if (ok) {
    ok <- all(if (lower_open) x > lower else x >= lower) &&
      all(if (upper_open) x < upper else x <= upper)
  }
  if (!ok)
    stop(name, " must contain plain finite numeric values with length ",
         paste(lengths, collapse = " or "))
  as.numeric(x)
}

.sim_structure <- function(x, name, allowed, required = character()) {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || is.data.frame(x))
    stop(name, " must be NULL or a named list")
  if (!length(x)) {
    if (length(required))
      stop(name, " must contain ", paste(required, collapse = ", "))
    return(x)
  }
  nm <- names(x)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(trimws(nm))) || anyDuplicated(nm))
    stop(name, " must be a named list with unique, non-missing names")
  unknown <- setdiff(nm, allowed)
  if (length(unknown))
    stop(name, " has unknown component(s): ", paste(unknown, collapse = ", "))
  missing <- setdiff(required, nm)
  if (length(missing))
    stop(name, " must contain ", paste(missing, collapse = ", "))
  x
}

# Construct a fixed departure that the explanatory design cannot absorb. The
# residual of a unit vector is orthogonal to every fitted design column; scale
# it so the best-supported focal row carries the requested magnitude. A
# saturated design has no such departure and must be refused rather than
# advertised as planted misfit.
.sim_explanatory_departure <- function(B, magnitude) {
  B <- as.matrix(B)
  magnitude <- .sim_scalar(magnitude, "explanatory departure", lower = 0,
                           lower_open = TRUE)
  if (!is.numeric(B) || !nrow(B) || !ncol(B) || any(!is.finite(B)))
    stop("the explanatory simulation design must be a finite numeric matrix")
  q <- qr(B, tol = 1e-10)
  if (q$rank >= nrow(B))
    stop("the explanatory design is saturated; increase the number of items ",
         "or objects before planting a fixed departure")
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  residual_leverage <- pmax(1 - rowSums(Q^2), 0)
  focal <- which.max(residual_leverage)
  if (!length(focal) || residual_leverage[focal] <= 1e-10)
    stop("the explanatory design has no stable row on which to plant a fixed departure")
  direction <- -drop(Q %*% Q[focal, ])
  direction[focal] <- direction[focal] + 1
  direction <- direction * magnitude / direction[focal]
  list(values = direction, index = focal)
}

# Construct a finite object vector with the requested sample correlation to x.
# The random component is centred and projected off x before being rescaled,
# rather than merely combined with x and expected to be orthogonal in a small
# realised object set.
.sim_correlated_fixed <- function(x, rho, target_sd,
                                  draw = function(n) stats::rnorm(n)) {
  if (abs(rho) == 1) return(rho * x)
  xx <- sum(x^2)
  for (attempt in seq_len(20L)) {
    z <- draw(length(x))
    z <- z - mean(z)
    z <- z - x * sum(z * x) / xx
    zsd <- stats::sd(z)
    if (is.finite(zsd) && zsd > sqrt(.Machine$double.eps)) {
      z <- z / zsd * target_sd
      return(rho * x + sqrt(1 - rho^2) * z)
    }
  }
  stop("could not construct the requested finite object correlation")
}

# Allocate a finite set of observations as evenly as possible across the
# requested IDs. Unlike sampling with replacement, this guarantees that every
# declared judge contributes data when the design contains enough rows.
.sim_balanced_ids <- function(ids, n, what = "levels") {
  if (n < length(ids))
    stop("the generated design has ", n, " observations, fewer than the ",
         length(ids), " requested ", what)
  out <- rep(ids, n %/% length(ids))
  rem <- n %% length(ids)
  if (rem) out <- c(out, sample(ids, rem, replace = FALSE))
  sample(out, length(out), replace = FALSE)
}

.sim_planted_count <- function(prop, n) {
  if (prop <= 0) return(0L)
  min(n, max(1L, as.integer(round(prop * n))))
}

# person locations from one of a few distributions
.sim_theta <- function(n, mean, sd, dist = "normal") {
  n <- .sim_count(n, "n")
  mean <- .sim_scalar(mean, "mean")
  sd <- .sim_scalar(sd, "sd", lower = 0)
  if (sd == 0) return(rep(mean, n))
  z <- switch(dist,
    normal  = stats::rnorm(n),
    uniform = stats::runif(n, -sqrt(3), sqrt(3)),
    skew    = { u <- stats::rgamma(n, 2, 1); (u - 2) / sqrt(2) },
    bimodal = { s <- sample(c(-1, 1), n, TRUE); s * 1.1 + stats::rnorm(n, 0, 0.5) },
    stats::rnorm(n))
  if (n == 1L) return(mean)
  mean + sd * as.numeric(scale(z))
}

# draw one item's responses for a vector of person locations. tau are the
# item's thresholds (length m; dichotomous m = 1). disc scales the whole
# exponent (a departure when != 1); guess is a lower asymptote (dichotomous).
.sim_item <- function(theta, tau, disc = 1, guess = 0) {
  if (!length(theta) || any(!is.finite(theta)) ||
      !length(tau) || any(!is.finite(tau)))
    stop("theta and tau must contain finite values")
  disc <- .sim_scalar(disc, "disc", lower = 0, lower_open = TRUE)
  guess <- .sim_scalar(guess, "guess", lower = 0, upper = 1,
                       upper_open = TRUE)
  m <- length(tau); xs <- 0:m
  cum <- c(0, cumsum(tau))
  eta <- disc * (outer(theta, xs) - matrix(cum, length(theta), m + 1L, byrow = TRUE))
  eta <- eta - apply(eta, 1, max)
  P <- exp(eta); P <- P / rowSums(P)
  if (guess > 0 && m == 1L) {
    P[, 2] <- guess + (1 - guess) * P[, 2]; P[, 1] <- 1 - P[, 2]
  }
  cs <- t(apply(P, 1, cumsum))
  as.integer(rowSums(stats::runif(length(theta)) > cs))     # category 0..m
}

# thresholds around an item location: evenly spread (rating-scale pattern),
# or with a supplied relative pattern (partial-credit: pattern varies by
# item); optionally deliberately disordered (one threshold dropped below its
# predecessor, re-centred so the item's mean location is preserved)
.sim_thresholds <- function(delta, m, spread, disordered = FALSE,
                            pattern = NULL) {
  if (m == 1L) return(delta)
  step <- if (is.null(pattern)) seq(-1, 1, length.out = m) * spread
          else pattern
  tau <- delta + step - mean(step)
  if (disordered && m >= 2L) {
    i <- max(2L, ceiling(m / 2))           # always has a predecessor to undercut
    # Set the adjacent reversal directly. Subtracting a fixed amount from the
    # original threshold did not guarantee disorder for PCM patterns whose
    # randomly generated gap happened to be wider than that amount.
    tau[i] <- tau[i - 1L] - 0.2 * spread
    tau <- tau - mean(tau) + delta         # keep the item's location honest
  }
  tau
}

#' Simulate person-by-item Rasch data
#'
#' Generates dichotomous, partial credit, or rating scale data. Optional
#' arguments introduce item misfit, guessing, multidimensionality, local
#' dependence, DIF, response styles, or missingness. Generating values are
#' stored in \code{attr(x, "truth")}. A positive planted proportion selects
#' at least one person or response cell when the requested departures can
#' coexist.
#'
#' @param n_persons,n_items Sample size and test length.
#' @param model \code{"dichotomous"}, \code{"PCM"}, or \code{"RSM"}. Under
#'   \code{"RSM"} every item shares one category-threshold pattern (items
#'   differ by location only); under \code{"PCM"} each item's threshold
#'   spacings and span are drawn afresh, as the partial credit model allows.
#' @param n_categories Response categories for polytomous models (>= 3).
#' @param theta_mean,theta_sd Mean and standard deviation of the person
#'   distribution.
#' @param theta_dist Shape of the person distribution: \code{"normal"},
#'   \code{"uniform"}, \code{"skew"}, or \code{"bimodal"}.
#' @param difficulty Either the two endpoints of an evenly spaced location
#'   range, or one location per item.
#' @param threshold_spread Half-range of the category thresholds about each
#'   item location (polytomous).
#' @param discrimination The item slope, supplied as one value or one per item.
#'   Values above 1 produce steeper responses and negative fit residuals.
#'   Values below 1 produce flatter responses and positive fit residuals.
#' @param guessing Scalar or length-\code{n_items} lower asymptote
#'   (dichotomous): low-location persons answer correctly by chance.
#' @param second_dim \code{NULL}, or \code{list(items=, rho=)}: the named items
#'   load on a second trait whose realised sample correlation with the first
#'   is \code{rho}. At least three persons are needed unless \code{rho} is
#'   -1 or 1. Each item is named once.
#' @param dependence \code{NULL}, or \code{list(pairs=, strength=)}: each pair's
#'   second item responds partly to the first. This departure feeds the
#'   residual-dependence diagnostics. Each directed pair is listed once.
#' @param dif \code{NULL}, or \code{list(items=, uniform=, nonuniform=)}: the
#'   named items function differently for the last person group: a location
#'   shift (\code{uniform}) and/or a slope change (\code{nonuniform}). Needs
#'   \code{n_groups >= 2}; each item is named once.
#' @param careless Proportion of persons who answer at random. Careless and
#'   response-style assignments are disjoint; their requested counts must fit.
#' @param response_style \code{NULL}, or \code{list(type=, prop=, strength=)}
#'   with \code{type} \code{"extreme"} or \code{"middle"}: a proportion
#'   \code{prop} of persons favour the end (or middle) categories regardless
#'   of the trait, with distortion \code{strength} (default 1.6) on the
#'   log-probability scale (polytomous).
#' @param speeded Proportion not reached at the last item: a growing tail of
#'   missing responses over the final items. These cells are kept distinct
#'   from any completely-at-random missing cells.
#' @param disordered \code{NULL} or item names/indices given disordered
#'   thresholds (polytomous; feeds the threshold diagnostics).
#' @param n_groups Number of equal person groups (a \code{group} factor column
#'   is added when > 1, for DIF).
#' @param missing Proportion of responses set missing completely at random,
#'   drawn from cells not already missing through speededness. The requested
#'   count must fit among those cells and leave at least one observed response.
#' @param seed Optional non-negative whole-number RNG seed.
#' @return A data frame of class \code{"rasch_sim"} (item columns
#'   \code{I01}..., an \code{id} column, and a \code{group} column when
#'   grouped), with \code{attr(x, "truth")} holding the generating parameters
#'   and the planted departures.
#' @examples
#' # a clean scale with one over-discriminating item and one DIF item
#' d <- simulate_rasch(400, 12, discrimination = c(3, rep(1, 11)),
#'                     dif = list(items = "I06", uniform = 1), n_groups = 2,
#'                     seed = 1)
#' fit <- rasch(d, id = "id", factors = "group")
#' fit$items[c("item", "infit_ms", "outfit_ms")]   # item 1 misfits
#' dif_anova(fit)$summary                           # item 6 flags
#' @export
simulate_rasch <- function(n_persons = 500, n_items = 20,
                           model = c("dichotomous", "PCM", "RSM"),
                           n_categories = 3, theta_mean = 0, theta_sd = 1,
                           theta_dist = "normal", difficulty = c(-2.5, 2.5),
                           threshold_spread = 1.2, discrimination = 1,
                           guessing = 0, second_dim = NULL, dependence = NULL,
                           dif = NULL, careless = 0, response_style = NULL,
                           speeded = 0, disordered = NULL,
                           n_groups = 1, missing = 0, seed = NULL) {
  seed <- .sim_seed(seed)
  if (!is.null(seed)) {
    .old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(.old_stream), add = TRUE)
    set.seed(seed)
  }
  model <- match.arg(model)
  N <- .sim_count(n_persons, "n_persons", 2L)
  I <- .sim_count(n_items, "n_items", 2L)
  n_groups <- .sim_count(n_groups, "n_groups")
  if (n_groups > N) stop("n_groups cannot exceed n_persons")
  theta_mean <- .sim_scalar(theta_mean, "theta_mean")
  theta_sd <- .sim_scalar(theta_sd, "theta_sd", lower = 0)
  threshold_spread <- .sim_scalar(threshold_spread, "threshold_spread",
                                  lower = 0, lower_open = TRUE)
  careless <- .sim_scalar(careless, "careless", lower = 0, upper = 1)
  speeded <- .sim_scalar(speeded, "speeded", lower = 0, upper = 1)
  missing <- .sim_scalar(missing, "missing", lower = 0, upper = 1,
                         upper_open = TRUE)
  theta_dist <- match.arg(theta_dist, c("normal", "uniform", "skew", "bimodal"))
  if (model != "dichotomous")
    n_categories <- .sim_count(n_categories, "n_categories", 3L)
  m <- if (model == "dichotomous") 1L else as.integer(n_categories) - 1L
  second_dim <- .sim_structure(second_dim, "second_dim",
                               c("items", "rho"), "items")
  dependence <- .sim_structure(dependence, "dependence",
                               c("pairs", "strength"), "pairs")
  dif <- .sim_structure(dif, "dif",
                        c("items", "uniform", "nonuniform"), "items")
  response_style <- .sim_structure(response_style, "response_style",
                                   c("type", "prop", "strength"))
  if (!is.null(response_style) && m == 1L) {
    warning("response_style applies to polytomous items only; ignored for dichotomous data")
    response_style <- NULL
  }
  inm <- sprintf("I%02d", seq_len(I))
  as_idx <- function(x, allow_duplicates = FALSE) {
    if (is.null(x) || !length(x)) return(integer(0))
    if (!is.null(dim(x)))
      stop("item selectors must be plain vectors, not matrices or arrays")
    if (is.character(x)) {
      i <- match(x, inm)
    } else if (is.numeric(x) && !is.complex(x) && is.null(oldClass(x))) {
      if (any(!is.finite(x)) || any(x != round(x)))
        stop("item index(es) must be finite whole numbers")
      i <- as.integer(x)
    } else {
      stop("item selectors must be item names or plain numeric indices")
    }
    if (anyNA(i) || any(i < 1L | i > I))
      stop("unknown item name(s)/index(es): ",
           paste(x[is.na(i) | i < 1L | i > I], collapse = ", "))
    if (!allow_duplicates && anyDuplicated(i))
      stop("item selectors must name each generated item at most once")
    i
  }

  # item locations, thresholds, slopes, guessing (with per-item overrides)
  difficulty <- .sim_vector(difficulty, "difficulty", unique(c(2L, I)))
  delta <- setNames(if (length(difficulty) == I) difficulty
                    else seq(difficulty[1], difficulty[2], length.out = I), inm)
  discrimination <- .sim_vector(discrimination, "discrimination",
                                unique(c(1L, I)), lower = 0,
                                lower_open = TRUE)
  guessing <- .sim_vector(guessing, "guessing", unique(c(1L, I)),
                          lower = 0, upper = 1, upper_open = TRUE)
  disc <- if (length(discrimination) == I) discrimination else
    rep(discrimination, I)
  guess <- if (length(guessing) == I) guessing else rep(guessing, I)
  if (m > 1L && any(guess > 0)) {
    warning("guessing applies to dichotomous items only; ignored for ", model)
    guess[] <- 0
  }
  dis_items <- as_idx(disordered)
  # a dichotomous item has one threshold and cannot be disordered: planting
  # it would record a generating feature the data do not carry
  if (m == 1L && length(dis_items)) {
    warning("disordered thresholds apply to polytomous items only; ",
            "ignored for dichotomous data")
    dis_items <- integer(0)
  }
  # RSM: one common threshold pattern (per-item location only). PCM: the
  # pattern itself varies across items, as the model allows -- each item's
  # spacings are drawn afresh (gated so the RNG stream of the other models
  # is untouched)
  patterns <- vector("list", I)
  if (model == "PCM" && m > 1L) patterns <- lapply(seq_len(I), function(i) {
    # ordered spacings with varied gaps and a varied overall span per item
    p <- cumsum(stats::runif(m, 0.5, 1.5))
    p <- (p - mean(p)) / (max(p) - min(p)) * 2 * threshold_spread *
      stats::runif(1, 0.85, 1.15)
    p
  })
  tau <- lapply(seq_len(I), function(i)
    .sim_thresholds(delta[i], m, threshold_spread, i %in% dis_items,
                    pattern = patterns[[i]]))
  names(tau) <- inm

  # person locations (primary) and groups
  theta <- .sim_theta(N, theta_mean, theta_sd, theta_dist)
  group <- if (n_groups > 1L)
    factor(sprintf("g%d", (seq_len(N) - 1L) %% n_groups + 1L)) else NULL

  # a second dimension for the nominated items: a correlated latent trait
  theta2 <- NULL; dim_items <- integer(0)
  if (!is.null(second_dim)) {
    dim_items <- as_idx(second_dim$items)
    if (!length(dim_items))
      stop("second_dim$items must name at least one item")
    rho <- second_dim$rho %||% 0.5
    if (length(rho) != 1L || !is.numeric(rho) || is.complex(rho) ||
        !is.null(dim(rho)) || !is.null(oldClass(rho)) || !is.finite(rho) ||
        abs(rho) > 1)
      stop("second_dim$rho must be a single correlation in [-1, 1]")
    if (theta_sd == 0)
      stop("second_dim requires theta_sd > 0 because a latent correlation is otherwise undefined")
    if (N < 3L && abs(rho) < 1)
      stop("second_dim needs at least three persons to realise a finite ",
           "correlation strictly between -1 and 1")
    # Construct the orthogonal component in the realised sample. This makes
    # the requested mean, spread and correlation exact rather than relying on
    # an independent random draw to be orthogonal in a small sample.
    theta2 <- theta_mean + .sim_correlated_fixed(
      theta - theta_mean, rho, theta_sd,
      draw = function(n) .sim_theta(n, 0, theta_sd, theta_dist))
  }

  X <- matrix(NA_integer_, N, I, dimnames = list(NULL, inm))
  dif_items <- as_idx(if (is.null(dif)) NULL else dif$items)
  if (!is.null(dif) && !length(dif_items))
    stop("dif$items must name at least one generated item")
  if (!is.null(dif)) {
    du <- .sim_scalar(dif$uniform %||% 0, "dif$uniform")
    dn <- .sim_scalar(dif$nonuniform %||% 0, "dif$nonuniform")
    if (length(dif_items) && any(disc[dif_items] + dn <= 0))
      stop("dif$nonuniform makes a planted item discrimination non-positive")
    dif$uniform <- du; dif$nonuniform <- dn
    if (du == 0 && dn == 0) {
      dif <- NULL
      dif_items <- integer(0)
    }
  }
  if (length(dif_items) && n_groups < 2L)
    stop("dif needs n_groups >= 2 (the last group carries the DIF)")
  dif_grp <- if (n_groups > 1L) levels(group)[n_groups] else NA

  # every regeneration of an item must honour that item's OWN generating
  # structure -- its trait (second dimension), its group-shifted thresholds
  # and slope (DIF), its guessing -- or a later misfit layer would silently
  # erase an earlier planted one
  item_pars <- function(i, p) {
    dif_here <- i %in% dif_items && !is.na(dif_grp) && group[p] == dif_grp
    list(tau = tau[[i]] + if (dif_here) dif$uniform %||% 0 else 0,
         disc = disc[i] + if (dif_here) dif$nonuniform %||% 0 else 0)
  }
  gen_item <- function(i, shift = rep(0, N)) {
    th <- (if (i %in% dim_items) theta2 else theta) + shift
    if (i %in% dif_items && !is.na(dif_grp)) {
      g2 <- group == dif_grp
      out <- integer(N)
      out[!g2] <- .sim_item(th[!g2], tau[[i]], disc[i], guess[i])
      out[g2]  <- .sim_item(th[g2], tau[[i]] + (dif$uniform %||% 0),
                            disc[i] + (dif$nonuniform %||% 0), guess[i])
      out
    } else .sim_item(th, tau[[i]], disc[i], guess[i])
  }
  # the model expectation of item i for every person, under the same
  # generating structure (trait, DIF shift, guessing) the responses used --
  # including any dependence carry-over already applied to that item, or in
  # a chain the residual of the middle item would carry the first pair's
  # shift as a systematic mean and leak it into the second pair
  exp_item <- function(i, shift = NULL) {
    th <- (if (i %in% dim_items) theta2 else theta) + (shift %||% 0)
    E <- vapply(seq_len(N), function(p) {
      pp <- item_pars(i, p)
      sum((0:m) * .p_item(th[p], pp$tau, pp$disc))
    }, 0)
    if (m == 1L && guess[i] > 0) E <- guess[i] + (1 - guess[i]) * E
    E
  }

  for (i in seq_len(I)) X[, i] <- gen_item(i)

  # response dependence: the second item of each pair partly follows the
  # first (adds d*(x1 - E1) to its exponent, inducing residual correlation);
  # the regeneration keeps i2's own DIF / second-dimension structure
  dep_pairs <- list()
  dep_pair_idx <- list()
  dep_shift <- vector("list", I)
  dep_sources <- integer(0)
  dep_keys <- character(0)
  if (!is.null(dependence)) {
    d_str <- .sim_scalar(dependence$strength %||% 1,
                         "dependence$strength")
    if (d_str == 0) dependence <- NULL
    for (pp in if (is.null(dependence)) list() else dependence$pairs) {
      ij <- as_idx(pp, allow_duplicates = TRUE)
      if (length(ij) != 2L || ij[1L] == ij[2L])
        stop("each dependence pair must name two different items")
      i1 <- ij[1]; i2 <- ij[2]
      pair_key <- paste(i1, i2, sep = "->")
      if (pair_key %in% dep_keys)
        stop("dependence$pairs repeats the directed pair ", inm[i1], " -> ",
             inm[i2], "; list each directed pair once")
      dep_keys <- c(dep_keys, pair_key)
      # regenerating i2 replaces the responses any earlier pair drew its
      # carry-over from, so a source may not become a later target
      if (i2 %in% dep_sources)
        stop("item ", inm[i2], " is the source of an earlier dependence ",
             "pair and the target of a later one: regenerating it would ",
             "erase the earlier dependence -- order the pairs so that each ",
             "item is a target before it is a source")
      dep_sources <- c(dep_sources, i1)
      # the expectation must match X1's actual generating structure
      # (guessing, DIF, second dimension), or the "residual" x1 - E1 has a
      # systematic mean and leaks an unplanted shift into the second item
      # the source's expectation is taken under the structure its responses
      # were actually drawn from, its own carry-over included
      shift <- d_str * (X[, i1] - exp_item(i1, dep_shift[[i1]])) / m
      # an item named as the second of several pairs carries every one of
      # them: regenerating from the new pair alone would erase the earlier
      # dependence and plant only the last
      dep_shift[[i2]] <- (dep_shift[[i2]] %||% 0) + shift
      X[, i2] <- gen_item(i2, shift = dep_shift[[i2]])
      dep_pairs[[length(dep_pairs) + 1L]] <- inm[ij]
      dep_pair_idx[[length(dep_pair_idx) + 1L]] <- ij
    }
  }

  # response styles (polytomous): a proportion of persons distort toward the
  # end or the middle categories regardless of the trait; the base
  # probabilities keep each item's own structure (trait, DIF) per person
  style_idx <- integer(0)
  if (!is.null(response_style) && m >= 2L) {
    stype <- response_style$type %||% "extreme"
    if (!is.character(stype) || length(stype) != 1L || is.na(stype) ||
        !is.null(dim(stype)) || !is.null(oldClass(stype)))
      stop("response_style$type must be one of 'extreme' or 'middle'")
    stype <- match.arg(stype, c("extreme", "middle"))
    response_style$type <- stype
    sprop <- .sim_scalar(response_style$prop %||% 0.15,
                         "response_style$prop", lower = 0, upper = 1)
    ss <- .sim_scalar(response_style$strength %||% 1.6,
                      "response_style$strength", lower = 0)
    # a style of zero strength distorts nothing and no one: drawing the
    # persons anyway would record an effect the data do not carry
    if (ss == 0 || sprop == 0) {
      response_style <- NULL
    } else {
    style_idx <- sample(N, .sim_planted_count(sprop, N))
    mid <- m / 2
    dev2 <- ((0:m - mid) / mid)^2
    # Apply the style on the log-probability scale. Constructing exp(ss *
    # dev2) first overflows at perfectly finite strengths, even though the
    # normalised categorical probabilities remain well defined.
    log_w <- if (stype == "extreme") ss * dev2 else -ss * dev2
    style_prob <- function(i, p, shift = 0) {
      th_p <- if (i %in% dim_items) theta2[p] else theta[p]
      pp <- item_pars(i, p)
      eta <- pp$disc * ((0:m) * (th_p + shift) -
        c(0, cumsum(pp$tau))) + log_w
      eta <- eta - max(eta)
      pr <- exp(eta)
      pr / sum(pr)
    }
    for (p in style_idx) {
      # Draw the styled marginal responses first. Dependence is then rebuilt
      # in its declared order from these responses. Reusing dep_shift here
      # would condition a target on the source response drawn before the
      # source itself acquired its response style.
      for (i in seq_len(I)) {
        if (is.na(X[p, i])) next
        X[p, i] <- sample.int(m + 1L, 1L,
                              prob = style_prob(i, p)) - 1L
      }
      styled_shift <- numeric(I)
      for (ij in dep_pair_idx) {
        i1 <- ij[1L]; i2 <- ij[2L]
        if (is.na(X[p, i1]) || is.na(X[p, i2])) next
        source_prob <- style_prob(i1, p, styled_shift[i1])
        source_mean <- sum((0:m) * source_prob)
        styled_shift[i2] <- styled_shift[i2] +
          d_str * (X[p, i1] - source_mean) / m
        X[p, i2] <- sample.int(
          m + 1L, 1L, prob = style_prob(i2, p, styled_shift[i2])) - 1L
      }
    }
    }
  }

  # careless responders: answer uniformly at random
  careless_idx <- integer(0)
  if (careless > 0) {
    n_careless <- .sim_planted_count(careless, N)
    # Careless responses replace the whole response vector, so overlap would
    # erase a requested response style and make its realised proportion
    # smaller than the argument states.
    ordinary <- setdiff(seq_len(N), style_idx)
    if (n_careless > length(ordinary))
      stop("the requested careless and response-style proportions cannot ",
           "coexist without overlap; reduce one of them")
    careless_idx <- if (n_careless)
      ordinary[sample.int(length(ordinary), n_careless)] else integer(0)
    X[careless_idx, ] <- matrix(sample(0:m, length(careless_idx) * I, TRUE),
                                length(careless_idx), I)
  }

  # speededness: a contiguous not-reached tail over the last items. `speeded`
  # persons drop out somewhere in the final zone, so the missing rate grows
  # linearly to `speeded` at the last item
  # the not-reached tail spans the last 40% of the test: with fewer than
  # three items there is no tail, so the request cannot be planted and must
  # not be recorded as though it were
  if (speeded > 0 && I < 3L) {
    warning("speededness needs at least three items to place a ",
            "not-reached tail; ignored")
    speeded <- 0
  }
  speeded_idx <- integer(0)
  if (speeded > 0 && I >= 3L) {
    k <- max(1L, round(0.4 * I)); z0 <- I - k
    speeded_idx <- sample(N, .sim_planted_count(speeded, N))
    for (p in speeded_idx) {
      stop_at <- z0 + sample.int(k, 1L) - 1L              # drop point in zone
      if (stop_at < I) X[p, (stop_at + 1L):I] <- NA
    }
  }
  missing_cells <- integer(0)
  if (missing > 0) {
    n_missing <- .sim_planted_count(missing, length(X))
    available <- which(!is.na(X))
    if (n_missing >= length(available))
      stop("the requested speededness and completely-at-random missingness ",
           "would remove every remaining response; reduce one of them")
    missing_cells <- available[sample.int(length(available), n_missing)]
    X[missing_cells] <- NA
  }

  out <- data.frame(id = sprintf("P%04d", seq_len(N)), X,
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(group)) out$group <- group

  planted <- character(0)
  if (any(disc != 1)) planted <- c(planted, sprintf("discrimination != 1 on %s",
    paste(inm[disc != 1], collapse = ", ")))
  if (any(guess > 0)) planted <- c(planted, sprintf("guessing on %s",
    paste(inm[guess > 0], collapse = ", ")))
  if (length(dim_items)) planted <- c(planted, sprintf(
    "second dimension (rho %.2f) on %s", second_dim$rho %||% 0.5,
    paste(inm[dim_items], collapse = ", ")))
  if (length(dep_pairs)) planted <- c(planted, sprintf(
    "response dependence: %s", paste(vapply(dep_pairs, paste, "",
                                            collapse = "-"), collapse = "; ")))
  if (length(dif_items)) planted <- c(planted, sprintf(
    "DIF (group %s) on %s: uniform %.2f, non-uniform %.2f", dif_grp,
    paste(inm[dif_items], collapse = ", "), dif$uniform %||% 0,
    dif$nonuniform %||% 0))
  if (length(careless_idx)) planted <- c(planted, sprintf(
    "%d careless responder(s)", length(careless_idx)))
  if (length(style_idx)) planted <- c(planted, sprintf(
    "%s response style on %d person(s)",
    response_style$type %||% "extreme", length(style_idx)))
  if (speeded > 0) planted <- c(planted, sprintf(
    "speededness (%.0f%% not-reached at the last item)", 100 * speeded))
  if (length(dis_items)) planted <- c(planted, sprintf(
    "disordered thresholds on %s", paste(inm[dis_items], collapse = ", ")))
  if (length(missing_cells)) planted <- c(planted, sprintf(
    "%d response%s missing completely at random (%.1f%%)",
    length(missing_cells), if (length(missing_cells) == 1L) "" else "s",
    100 * length(missing_cells) / length(X)))

  attr(out, "truth") <- list(
    layout = "rasch",
    description = sprintf("%s, %d persons x %d items%s", model, N, I,
      if (!is.null(group)) sprintf(", %d groups", nlevels(group)) else ""),
    model = model, n_persons = N, n_items = I,
    person_id = out$id, theta = theta, theta2 = theta2,
    difficulty = delta, thresholds = tau,
    discrimination = disc, guessing = guess,
    groups = group, dim_items = inm[dim_items], dif_items = inm[dif_items],
    careless_idx = careless_idx, style_idx = style_idx,
    speeded_idx = speeded_idx, missing_cells = missing_cells,
    planted = planted)
  class(out) <- c("rasch_sim", "data.frame")
  out
}

# category probabilities for one location (used by the dependence term)
.p_item <- function(theta, tau, disc = 1) {
  m <- length(tau); cum <- c(0, cumsum(tau))
  e <- disc * ((0:m) * theta - cum); e <- e - max(e)
  p <- exp(e); p / sum(p)
}

#' @export
print.rasch_sim <- function(x, ...) {
  tr <- attr(x, "truth")
  cat(sprintf("Simulated %s data: %s\n", tr$layout,
              tr$description %||% sprintf("%d rows", nrow(x))))
  if (length(tr$planted)) {
    cat("Planted departures:\n")
    for (p in tr$planted) cat(paste0("  - ", p, "\n"))
  } else cat("Model-conforming (no departures planted).\n")
  invisible(x)
}

#' Simulate paired-comparison data
#'
#' Generates dichotomous or ordered paired comparisons from the
#' Bradley--Terry--Luce model. Optional arguments introduce a second object
#' attribute, erratic judges, or within-judge dependence. Generating values are
#' stored in \code{attr(x, "truth")}.
#'
#' @param n_objects,n_judges Objects to scale and judges comparing them. Every
#'   judge is allocated at least one comparison; the simulator refuses a
#'   design with fewer comparisons than judges.
#' @param reps_per_pair Comparisons made of each object pair.
#' @param model \code{"dichotomous"} (a winner) or \code{"polytomous"} (a rated
#'   margin in \code{n_categories} categories; an earlier development-era value
#'   \code{"graded"} is accepted as an alias).
#' @param n_categories Categories for the polytomous model.
#' @param object_sd Spread of the object locations (evenly spaced, sum-zero).
#' @param object_locations Optional numeric vector of generated object
#'   locations. It must have length \code{n_objects}; names, when supplied,
#'   must identify the generated objects. Values are centred to identify the
#'   origin and take precedence over \code{object_sd}.
#' @param second_attribute \code{NULL}, or \code{list(rho=)}: half the judges
#'   rank by a second object attribute whose realised correlation with the
#'   first is \code{rho}.
#'   This introduces residual dimensionality and possible intransitivity.
#' @param erratic_judges Proportion of judges who choose at random.
#' @param dependence \code{NULL}, or \code{list(exposure=, carry_over=)}:
#'   within-judge order effects (a seen-before advantage and a pull from the
#'   judge's own earlier verdicts). Adds an \code{order} column. Feeds the
#'   dependence effects fitted by \code{\link{btl}}.
#' @param seed Optional non-negative whole-number RNG seed.
#' @return A data frame of class \code{"rasch_sim"}: \code{object_a},
#'   \code{object_b}, \code{winner} (or \code{response} when polytomous),
#'   \code{judge}, and \code{order} when dependence is planted; with
#'   \code{attr(x, "truth")}.
#' @examples
#' d <- simulate_btl(8, 12, erratic_judges = 0.15, seed = 1)
#' bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
#' bt$judges          # the erratic judges carry large fit residuals
#' @export
simulate_btl <- function(n_objects = 8, n_judges = 12, reps_per_pair = 25,
                         model = c("dichotomous", "polytomous", "graded"),
                         n_categories = 4,
                         object_sd = 1, second_attribute = NULL,
                         erratic_judges = 0, dependence = NULL, seed = NULL,
                         object_locations = NULL) {
  seed <- .sim_seed(seed)
  if (!is.null(seed)) {
    .old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(.old_stream), add = TRUE)
    set.seed(seed)
  }
  model <- match.arg(model)
  # "graded" is an earlier development name for the polytomous comparison model,
  # kept as a working alias for released user code
  if (model == "graded") model <- "polytomous"
  K <- .sim_count(n_objects, "n_objects", 3L)
  J <- .sim_count(n_judges, "n_judges", 2L)
  reps_per_pair <- .sim_count(reps_per_pair, "reps_per_pair")
  object_sd <- .sim_scalar(object_sd, "object_sd", lower = 0)
  erratic_judges <- .sim_scalar(erratic_judges, "erratic_judges",
                                lower = 0, upper = 1)
  second_attribute <- .sim_structure(second_attribute, "second_attribute",
                                     "rho")
  dependence <- .sim_structure(dependence, "dependence",
                               c("exposure", "carry_over"))
  if (model == "polytomous")
    n_categories <- .sim_count(n_categories, "n_categories", 3L)
  m <- if (model == "polytomous") as.integer(n_categories) - 1L else 1L
  objs <- sprintf("O%d", seq_len(K)); jids <- sprintf("J%d", seq_len(J))
  if (is.null(object_locations)) {
    beta <- setNames(as.numeric(scale(seq_len(K))) * object_sd, objs)
  } else {
    location_names <- names(object_locations)
    object_locations <- .sim_vector(object_locations, "object_locations", K)
    if (!is.null(location_names) &&
        (!setequal(location_names, objs) || anyDuplicated(location_names)))
      stop("named object_locations must identify every generated object exactly")
    if (!is.null(location_names))
      object_locations <- stats::setNames(object_locations, location_names)[objs]
    beta <- setNames(object_locations - mean(object_locations), objs)
  }
  tau <- if (m > 1L) .sim_thresholds(0, m, 1.2) else NULL

  # a second object attribute (orthogonal part) for the two-camp design
  beta2 <- NULL; camp <- NULL
  if (!is.null(second_attribute)) {
    beta_sd <- stats::sd(beta)
    if (!is.finite(beta_sd) || beta_sd <= sqrt(.Machine$double.eps))
      stop("second_attribute requires object locations with positive spread because an object correlation is otherwise undefined")
    rho <- second_attribute$rho %||% 0.3
    rho <- .sim_scalar(rho, "second_attribute$rho", lower = -1, upper = 1)
    beta2 <- setNames(.sim_correlated_fixed(beta, rho, beta_sd), objs)
    camp <- setNames(rep(c("a", "b"), length.out = J), jids)
  }
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(object_a = rep(pr[, 1], each = reps_per_pair),
                  object_b = rep(pr[, 2], each = reps_per_pair),
                  stringsAsFactors = FALSE)
  d$judge <- .sim_balanced_ids(jids, nrow(d), "judges")
  n_erratic <- .sim_planted_count(erratic_judges, J)
  erratic <- if (n_erratic) sample(jids, n_erratic, replace = FALSE) else
    character(0)

  win_prob <- function(a, b, jd) {
    ba <- if (!is.null(camp) && camp[jd] == "b") beta2[a] else beta[a]
    bb <- if (!is.null(camp) && camp[jd] == "b") beta2[b] else beta[b]
    ba - bb
  }
  # dependence needs a per-judge judgment order and running history
  if (!is.null(dependence)) {
    exq <- .sim_scalar(dependence$exposure %||% 0, "dependence$exposure")
    cry <- .sim_scalar(dependence$carry_over %||% 0,
                       "dependence$carry_over")
    if (exq == 0 && cry == 0) dependence <- NULL
  }
  if (!is.null(dependence)) {
    d <- d[order(d$judge), ]
    d$order <- stats::ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
    seen <- new.env(parent = emptyenv()); hs <- new.env(parent = emptyenv())
    hc <- new.env(parent = emptyenv())
    g0 <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
    resp <- integer(nrow(d))
    for (r in seq_len(nrow(d))) {
      j <- d$judge[r]; a <- d$object_a[r]; b <- d$object_b[r]
      ka <- paste(j, a); kb <- paste(j, b)
      lp <- win_prob(a, b, j) +
        exq * (as.numeric(g0(seen, ka) > 0) - as.numeric(g0(seen, kb) > 0)) +
        cry * ((if (g0(hc, ka) > 0) g0(hs, ka) / g0(hc, ka) else 0) -
               (if (g0(hc, kb) > 0) g0(hs, kb) / g0(hc, kb) else 0))
      x <- if (j %in% erratic) sample(0:m, 1)
           else if (m == 1L) as.integer(stats::runif(1) < stats::plogis(lp))
           else sample(0:m, 1, prob = .p_item(lp, tau))
      resp[r] <- x
      assign(ka, g0(seen, ka) + 1, seen); assign(kb, g0(seen, kb) + 1, seen)
      assign(ka, g0(hc, ka) + 1, hc);     assign(kb, g0(hc, kb) + 1, hc)
      assign(ka, g0(hs, ka) + (2 * x / m - 1), hs)
      assign(kb, g0(hs, kb) + (2 * (m - x) / m - 1), hs)
    }
  } else {
    lp <- vapply(seq_len(nrow(d)), function(r)
      win_prob(d$object_a[r], d$object_b[r], d$judge[r]), 0)
    resp <- integer(nrow(d))
    reg <- !(d$judge %in% erratic)
    resp[reg] <- if (m == 1L) as.integer(stats::runif(sum(reg)) < stats::plogis(lp[reg]))
                 else vapply(which(reg), function(r) sample(0:m, 1, prob = .p_item(lp[r], tau)), 0L)
    if (any(!reg)) resp[!reg] <- sample(0:m, sum(!reg), TRUE)
  }

  if (m == 1L) d$winner <- ifelse(resp == 1L, d$object_a, d$object_b)
  else d$response <- resp
  rownames(d) <- NULL

  planted <- character(0)
  if (length(erratic)) planted <- c(planted,
    sprintf("%d erratic judge(s): %s", length(erratic), paste(erratic, collapse = ", ")))
  if (!is.null(second_attribute)) planted <- c(planted,
    sprintf("second object attribute (rho %.2f), two judge camps",
            second_attribute$rho %||% 0.3))
  if (!is.null(dependence)) planted <- c(planted, sprintf(
    "within-judge dependence: exposure %.2f, carry-over %.2f",
    dependence$exposure %||% 0, dependence$carry_over %||% 0))

  attr(d, "truth") <- list(
    layout = "btl",
    description = sprintf("%s, %d objects, %d judges, %d comparisons",
                          model, K, J, nrow(d)),
    model = model, location = beta, location2 = beta2,
    attribute_correlation = if (is.null(beta2)) NULL else
      stats::cor(beta, beta2), camp = camp,
    erratic = erratic, planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}

#' Simulate many-facet Rasch data
#'
#' Generates fully crossed ratings from a many-facet Rasch model (Linacre
#' 1989), with optional erratic raters, item-by-rater interaction, or halo. A
#' positive rater proportion selects at least one rater when the requested
#' departures can coexist.
#'
#' @param n_persons,n_items,n_raters Facet sizes (fully crossed).
#' @param n_categories Rating categories.
#' @param theta_sd,item_sd Spread of person ability and item difficulty.
#' @param rater_severity_sd Spread of rater severities (the core facet;
#'   recovered in \code{facet_effects}).
#' @param erratic_raters Proportion of raters who rate at random (feeds the
#'   rater fit residual). Erratic and halo raters are disjoint.
#' @param interaction \code{NULL}, or \code{list(rater=, item=, bias=)}: one
#'   rater is unusually harsh (positive) or lenient (negative) on one item.
#'   Feeds the item-by-rater interaction (fit with \code{interaction = }).
#' @param halo Proportion of raters showing a halo effect: they rate by the
#'   person's overall level and barely differentiate items (feeds the rater
#'   fit residual and the item-by-rater interaction). Its requested count
#'   must fit among the non-erratic raters.
#' @param seed Optional non-negative whole-number RNG seed.
#' @return A long data frame of class \code{"rasch_sim"} (\code{person},
#'   \code{item}, \code{rater}, \code{score}) ready for
#'   \code{\link{rasch_mfrm}}, with the truth attached.
#' @examples
#' d <- simulate_mfrm(60, 5, 6, rater_severity_sd = 0.8, seed = 1)
#' mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
#'                  facets = "rater")
#' cor(mf$facet_effects$rater$severity, attr(d, "truth")$severity)  # recovered
#' @export
simulate_mfrm <- function(n_persons = 80, n_items = 5, n_raters = 6,
                          n_categories = 4, theta_sd = 1.2, item_sd = 1,
                          rater_severity_sd = 0.6, erratic_raters = 0,
                          interaction = NULL, halo = 0, seed = NULL) {
  seed <- .sim_seed(seed)
  if (!is.null(seed)) {
    .old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(.old_stream), add = TRUE)
    set.seed(seed)
  }
  N <- .sim_count(n_persons, "n_persons", 2L)
  I <- .sim_count(n_items, "n_items", 2L)
  R <- .sim_count(n_raters, "n_raters", 2L)
  n_categories <- .sim_count(n_categories, "n_categories", 2L)
  theta_sd <- .sim_scalar(theta_sd, "theta_sd", lower = 0)
  item_sd <- .sim_scalar(item_sd, "item_sd", lower = 0)
  rater_severity_sd <- .sim_scalar(rater_severity_sd, "rater_severity_sd",
                                   lower = 0)
  erratic_raters <- .sim_scalar(erratic_raters, "erratic_raters",
                                lower = 0, upper = 1)
  halo <- .sim_scalar(halo, "halo", lower = 0, upper = 1)
  interaction <- .sim_structure(interaction, "interaction",
                                c("rater", "item", "bias"),
                                c("rater", "item", "bias"))
  m <- n_categories - 1L
  pids <- sprintf("P%03d", seq_len(N)); iids <- sprintf("I%d", seq_len(I))
  rids <- sprintf("R%d", seq_len(R))
  theta <- .sim_theta(N, 0, theta_sd)
  delta <- setNames(seq(-item_sd, item_sd, length.out = I), iids)
  lambda <- setNames(as.numeric(scale(stats::rnorm(R))) * rater_severity_sd, rids)
  base_tau <- .sim_thresholds(0, m, 1.2)
  erratic <- if (erratic_raters > 0)
    rids[seq_len(.sim_planted_count(erratic_raters, R))] else character(0)
  # Halo raters are drawn from the end and remain disjoint from erratic
  # raters, whose random scores would erase the halo mechanism.
  halo_r <- if (halo > 0) {
    pool <- setdiff(rev(rids), erratic)
    n_halo <- .sim_planted_count(halo, R)
    if (n_halo > length(pool))
      stop("the requested erratic-rater and halo proportions cannot coexist ",
           "without overlap; reduce one of them")
    pool[seq_len(n_halo)]
  } else character(0)
  int_bias <- matrix(0, I, R, dimnames = list(iids, rids))
  if (!is.null(interaction)) {
    plain_level <- function(x)
      is.atomic(x) && is.null(dim(x)) && is.null(oldClass(x)) &&
        length(x) == 1L && !is.na(x)
    if (!plain_level(interaction$item) || !(interaction$item %in% iids) ||
        !plain_level(interaction$rater) || !(interaction$rater %in% rids))
      stop("interaction$item and interaction$rater must each name one generated level")
    interaction$bias <- .sim_scalar(interaction$bias, "interaction$bias")
    if (interaction$bias == 0) {
      interaction <- NULL
    } else {
      # an erratic rater answers at random, discarding the whole rating model
      # for that rater: a bias planted on one would not be in the data, while
      # the recorded truth would still claim it
      if (interaction$rater %in% erratic)
        stop("interaction$rater '", interaction$rater, "' is one of the ",
             "erratic raters, who answer at random: the planted bias would ",
             "not appear in the data. Nominate another rater, or lower ",
             "erratic_raters")
      int_bias[interaction$item, interaction$rater] <- interaction$bias
    }
  }

  grid <- expand.grid(p = seq_len(N), i = seq_len(I), r = seq_len(R))
  score <- integer(nrow(grid))
  for (i in seq_len(I)) for (r in seq_len(R)) {
    rows <- grid$i == i & grid$r == r
    # item difficulty and rater severity shift the person's thresholds; a halo
    # rater ignores the item's own difficulty (uses the mean instead)
    di <- if (rids[r] %in% halo_r) mean(delta) else delta[i]
    tau_ir <- base_tau + di + lambda[r] + int_bias[i, r]
    score[rows] <- if (rids[r] %in% erratic) sample(0:m, sum(rows), TRUE)
                   else .sim_item(theta[grid$p[rows]], tau_ir)
  }
  d <- data.frame(person = pids[grid$p], item = iids[grid$i],
                  rater = rids[grid$r], score = score, stringsAsFactors = FALSE)

  planted <- character(0)
  if (length(erratic)) planted <- c(planted,
    sprintf("%d erratic rater(s): %s", length(erratic), paste(erratic, collapse = ", ")))
  if (length(halo_r)) planted <- c(planted,
    sprintf("%d halo rater(s): %s", length(halo_r), paste(halo_r, collapse = ", ")))
  if (!is.null(interaction)) planted <- c(planted, sprintf(
    "rater-by-item bias: %s on %s (%.2f)", interaction$rater,
    interaction$item, interaction$bias))
  planted <- c(planted, sprintf("rater severities SD %.2f", stats::sd(lambda)))

  attr(d, "truth") <- list(
    layout = "mfrm",
    description = sprintf("%d persons x %d items x %d raters (%d ratings)",
                          N, I, R, nrow(d)),
    person_id = pids, theta = theta, difficulty = delta, severity = lambda,
    erratic = erratic, halo = halo_r, planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}

#' Simulate extended frame-of-reference data with differing units
#'
#' Generates data whose latent unit differs across item-set by person-group
#' frames (Humphry 2005): a person in group g responding to an item in set s
#' does so at the frame unit rho = alpha_set * phi_group scaling the whole
#' exponent. The planted set- and group-unit ratios are recovered by
#' \code{\link{rasch_efrm}}. A positive careless-response or missingness
#' proportion selects at least one person or response cell.
#'
#' @param n_per_group Persons in each group.
#' @param items_per_set Items in each set.
#' @param n_sets,n_groups Numbers of item sets and person groups.
#' @param set_unit_ratio,group_unit_ratio Geometric span of the set and group
#'   units across their levels (1 = equal units, i.e. an ordinary Rasch fit).
#' @param n_categories Response categories per item: 2 (the default) gives
#'   dichotomous items; larger values give partial credit items whose
#'   evenly spaced thresholds are centred on the item locations, with the
#'   frame unit scaling the whole exponent as in the dichotomous case.
#' @param theta_sd Spread of person ability.
#' @param item_drift Optional \code{list(items=, group=, shift=)}. The named
#'   item or items move by \code{shift} logits in one generated person group,
#'   violating item invariance across frames.
#' @param careless Proportion of persons whose complete response vectors are
#'   replaced by random category choices.
#' @param missing Proportion of response cells set missing completely at
#'   random after the responses are generated. It must leave at least one
#'   observed response.
#' @param seed Optional non-negative whole-number RNG seed.
#' @return A wide data frame of class \code{"rasch_sim"}, containing an ID,
#'   item columns, and group. Its truth attribute contains the item-set map
#'   required by \code{\link{rasch_efrm}}.
#' @examples
#' d <- simulate_efrm(200, 6, set_unit_ratio = 1.3, seed = 1)
#' tr <- attr(d, "truth")
#' ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
#'                  boot_reps = 0)    # point estimates only
#' ef$alpha_table   # planted ratio 1.3, recovered within small-sample noise
#' @export
simulate_efrm <- function(n_per_group = 300, items_per_set = 8, n_sets = 2,
                          n_groups = 2, set_unit_ratio = 1.3,
                          group_unit_ratio = 1, n_categories = 2,
                          theta_sd = 1.3, seed = NULL, item_drift = NULL,
                          careless = 0, missing = 0) {
  # recorded before validation assigns to the formals, which would make
  # missing() report every argument as supplied
  set_ratio_given <- !missing(set_unit_ratio)
  group_ratio_given <- !missing(group_unit_ratio)
  seed <- .sim_seed(seed)
  if (!is.null(seed)) {
    .old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(.old_stream), add = TRUE)
    set.seed(seed)
  }
  S <- .sim_count(n_sets, "n_sets")
  G <- .sim_count(n_groups, "n_groups")
  K <- .sim_count(items_per_set, "items_per_set", 2L)
  npg <- .sim_count(n_per_group, "n_per_group", 2L)
  n_categories <- .sim_count(n_categories, "n_categories", 2L)
  m <- as.integer(n_categories) - 1L
  set_unit_ratio <- .sim_scalar(set_unit_ratio, "set_unit_ratio",
                                lower = 0, lower_open = TRUE)
  group_unit_ratio <- .sim_scalar(group_unit_ratio, "group_unit_ratio",
                                  lower = 0, lower_open = TRUE)
  theta_sd <- .sim_scalar(theta_sd, "theta_sd", lower = 0)
  careless <- .sim_scalar(careless, "careless", lower = 0, upper = 1)
  missing <- .sim_scalar(missing, "missing", lower = 0, upper = 1,
                         upper_open = TRUE)
  item_drift <- .sim_structure(item_drift, "item_drift",
                               c("items", "group", "shift"),
                               c("items", "group", "shift"))
  # a ratio is a comparison between frames: with one frame there is nothing
  # to span, and the span would silently return a unit while the recorded
  # truth claimed the requested ratio
  # only an EXPLICIT request is refused: the defaults describe the ordinary
  # multi-frame design and resolve to the reference when there is one frame
  if (S == 1L) {
    if (set_ratio_given && set_unit_ratio != 1)
      stop("`set_unit_ratio` must be 1 when `n_sets` is 1: a unit ratio is ",
           "a comparison between sets, so none can be planted")
    set_unit_ratio <- 1
  }
  if (G == 1L) {
    if (group_ratio_given && group_unit_ratio != 1)
      stop("`group_unit_ratio` must be 1 when `n_groups` is 1: a unit ratio ",
           "is a comparison between groups, so none can be planted")
    group_unit_ratio <- 1
  }
  # set and group units span the ratio geometrically, normalised to mean 1
  gspan <- function(ratio, n) { u <- exp(seq(0, log(ratio), length.out = n)); u / exp(mean(log(u))) }
  alpha <- gspan(set_unit_ratio, S)
  phi <- gspan(group_unit_ratio, G)
  set_items <- lapply(seq_len(S), function(s) sprintf("S%dI%02d", s, seq_len(K)))
  inm <- unlist(set_items)
  delta <- setNames(rep(seq(-1.5, 1.5, length.out = K), S), inm)
  set_of <- setNames(rep(seq_len(S), each = K), inm)

  tau_list <- lapply(inm, function(nm) .sim_thresholds(delta[nm], m, 0.8))
  names(tau_list) <- inm

  grp <- factor(rep(sprintf("g%d", seq_len(G)), each = npg))
  gidx <- rep(seq_len(G), each = npg)                  # group NUMBER per person
  names(phi) <- sprintf("g%d", seq_len(G))             # unit by group label
  drift <- rep(0, length(inm)); names(drift) <- inm
  drift_group <- NULL
  if (!is.null(item_drift)) {
    if (!is.atomic(item_drift$items) || !is.null(dim(item_drift$items)) ||
        !is.null(oldClass(item_drift$items)))
      stop("item_drift$items must be a plain vector of generated item names")
    drift_items <- as.character(item_drift$items)
    if (!length(drift_items) || anyNA(drift_items) ||
        any(!nzchar(drift_items)) || anyDuplicated(drift_items) ||
        any(!drift_items %in% inm))
      stop("item_drift$items must name one or more generated items exactly once")
    if (!is.atomic(item_drift$group) || !is.null(dim(item_drift$group)) ||
        !is.null(oldClass(item_drift$group)) ||
        length(item_drift$group) != 1L || is.na(item_drift$group) ||
        !as.character(item_drift$group) %in% levels(grp))
      stop("item_drift$group must name one generated group")
    item_drift$shift <- .sim_scalar(item_drift$shift, "item_drift$shift")
    if (item_drift$shift == 0) {
      item_drift <- NULL
    } else {
      drift[drift_items] <- item_drift$shift
      drift_group <- as.character(item_drift$group)
    }
  }
  N <- length(grp); theta <- .sim_theta(N, 0, theta_sd)
  X <- matrix(NA_integer_, N, length(inm), dimnames = list(NULL, inm))
  for (col in seq_along(inm)) {
    # phi is indexed by group NUMBER: as.integer() on the factor would order
    # the levels as strings (g1, g10, g2, ...) and attach each unit to the
    # wrong group as soon as there are ten of them
    s <- set_of[inm[col]]
    rho <- alpha[s] * phi[gidx]                        # per-person unit
    shift <- if (is.null(drift_group)) rep(0, N) else
      ifelse(as.character(grp) == drift_group, drift[inm[col]], 0)
    if (m == 1L) {
      X[, col] <- as.integer(stats::runif(N) <
                    stats::plogis(rho * (theta - delta[inm[col]] - shift)))
    } else {
      ct <- cumsum(tau_list[[inm[col]]])
      E <- cbind(0, sapply(seq_len(m), function(k)
        rho * (k * theta - ct[k] - k * shift)))
      P <- exp(E - apply(E, 1, max)); P <- P / rowSums(P)
      cum <- P %*% upper.tri(diag(m + 1L), diag = TRUE)
      X[, col] <- as.integer(rowSums(stats::runif(N) > cum))
    }
  }
  careless_idx <- integer(0)
  if (careless > 0) {
    careless_idx <- sample.int(N, .sim_planted_count(careless, N))
    X[careless_idx, ] <- matrix(
      sample(0:m, length(careless_idx) * ncol(X), replace = TRUE),
      nrow = length(careless_idx), ncol = ncol(X))
  }
  missing_cells <- integer(0)
  if (missing > 0) {
    n_missing <- .sim_planted_count(missing, length(X))
    if (n_missing >= length(X))
      stop("the requested missingness would remove every response; reduce `missing`")
    missing_cells <- sample.int(
      length(X), n_missing)
    X[missing_cells] <- NA_integer_
  }
  out <- data.frame(id = sprintf("P%04d", seq_len(N)), X, group = grp,
                    check.names = FALSE, stringsAsFactors = FALSE)

  planted <- sprintf("set-unit ratio %.2f across %d sets", set_unit_ratio, S)
  if (group_unit_ratio != 1)
    planted <- c(planted, sprintf("group-unit ratio %.2f across %d groups",
                                  group_unit_ratio, G))
  if (!is.null(item_drift) && item_drift$shift != 0)
    planted <- c(planted, sprintf("item drift in %s: %s shifted %.2f logits",
      drift_group, paste(names(drift)[drift != 0], collapse = ", "),
      item_drift$shift))
  if (length(careless_idx)) planted <- c(planted,
    sprintf("%d careless responder(s)", length(careless_idx)))
  if (length(missing_cells)) planted <- c(planted,
    sprintf("%.1f%% responses missing completely at random",
            100 * length(missing_cells) / length(X)))
  attr(out, "truth") <- list(
    layout = "efrm",
    description = sprintf("%d persons, %d sets x %d groups, %d items (%d categories)",
                          N, S, G, length(inm), m + 1L),
    person_id = out$id, theta = theta, difficulty = delta,
    thresholds = if (m > 1L) tau_list else NULL,
    alpha = alpha, phi = phi,
    item_sets = setNames(set_items, sprintf("set%d", seq_len(S))),
    groups = grp, item_drift = if (is.null(item_drift)) NULL else
      list(items = names(drift)[drift != 0], group = drift_group,
           shift = item_drift$shift),
    careless_idx = careless_idx, missing_cells = missing_cells,
    planted = planted)
  class(out) <- c("rasch_sim", "data.frame")
  out
}

#' Replicate a simulation for Monte Carlo studies
#'
#' Calls one of the \code{simulate_*} functions \code{n} times with successive
#' seeds, returning the datasets as a list -- for power, Type-I, or
#' parameter-recovery studies.
#'
#' @param FUN A simulator, e.g. \code{\link{simulate_rasch}}.
#' @param n Number of datasets.
#' @param ... Arguments passed to \code{FUN} (the same each replicate).
#' @param seed Seed of the first replicate (each subsequent one increments it).
#' @return A list of class \code{"rasch_sim_batch"}, one simulated dataset per
#'   element.
#' @examples
#' # 8 datasets with a planted DIF item; how often is it flagged?
#' batch <- sim_replicate(simulate_rasch, 4, n_persons = 300, n_items = 8,
#'                        dif = list(items = "I05", uniform = 0.8), n_groups = 2,
#'                        seed = 1)
#' # sim_apply() is resilient: a replicate the estimator refuses (e.g. a
#' # small or disconnected draw) contributes NA instead of aborting the run
#' flagged <- sim_apply(batch, function(d)
#'   dif_anova(rasch(d, id = "id", factors = "group"))$summary$uniform_DIF[5])
#' mean(flagged, na.rm = TRUE)
#' @export
sim_replicate <- function(FUN, n, ..., seed = NULL) {
  if (!is.function(FUN)) stop("FUN must be a simulation function")
  n <- .sim_count(n, "n")
  base <- if (is.null(seed)) sample.int(1e6, 1L)
          else .sim_seed(seed)
  if (base > .Machine$integer.max - n + 1L)
    stop("seed plus the requested replicate count exceeds the integer range")
  # Subtract before adding: at the upper integer boundary, `base + k - 1L`
  # overflows at the intermediate addition even when k is one and the final
  # seed is valid.
  reps <- lapply(seq_len(n), function(k) FUN(..., seed = base + (k - 1L)))
  structure(reps, class = "rasch_sim_batch", n = n,
            layout = attr(reps[[1]], "truth")$layout)
}

#' Apply a statistic across a simulation batch
#'
#' Applies \code{FUN} to each replicate of a \code{\link{sim_replicate}}
#' batch, catching replicates on which \code{FUN} errors -- for example a
#' small or disconnected draw the estimator refuses as unidentified -- so a
#' single failure does not abort the whole Monte-Carlo run. Failed
#' replicates contribute \code{NA}; the number of failures and the distinct
#' error messages are attached as attributes.
#'
#' @param batch A \code{"rasch_sim_batch"} from \code{\link{sim_replicate}}
#'   (or any list of datasets).
#' @param FUN A function of one dataset returning a scalar statistic.
#' @param ... Further arguments passed to \code{FUN}.
#' @return A vector of per-replicate statistics, with \code{NA} where the
#'   function failed. Attribute \code{n_failed} gives the failure count;
#'   \code{failure_messages} contains the distinct messages.
#' @examples
#' batch <- sim_replicate(simulate_rasch, 10, n_persons = 300, n_items = 8,
#'                        seed = 1)
#' psi <- sim_apply(batch, function(d) rasch(d)$psi$PSI)
#' mean(psi, na.rm = TRUE)
#' @export
sim_apply <- function(batch, FUN, ...) {
  if (!is.list(batch) || !length(batch)) stop("batch must be a non-empty list")
  if (!is.function(FUN)) stop("FUN must be a function")
  res <- lapply(batch, function(d)
    tryCatch(list(ok = TRUE, v = FUN(d, ...)),
             error = function(e) list(ok = FALSE, v = NA, msg = conditionMessage(e))))
  # a one-column data frame and a list holding a vector both have length 1:
  # unlisting them later would expand one dataset into several results while
  # the failure count stayed at zero
  scalar <- function(v) !is.null(v) && is.atomic(v) && !is.list(v) &&
    length(v) == 1L
  valid <- vapply(res, function(r) isTRUE(r$ok) && scalar(r$v), TRUE)
  for (i in which(!valid & vapply(res, `[[`, logical(1), "ok"))) {
    res[[i]]$ok <- FALSE
    res[[i]]$msg <- "FUN must return one atomic scalar value"
  }
  ok <- vapply(res, `[[`, logical(1), "ok")
  vals <- lapply(res, function(r) {
    v <- r$v; if (!scalar(v)) NA else v[[1]]
  })
  out <- tryCatch(unlist(vals, use.names = FALSE),
                  error = function(e) vals)
  msgs <- unique(vapply(res[!ok], function(r) r$msg, ""))
  structure(out, n_failed = sum(!ok), failure_messages = msgs)
}

#' @export
print.rasch_sim_batch <- function(x, ...) {
  cat(sprintf("%d simulated %s datasets (Monte Carlo batch)\n",
              attr(x, "n"), attr(x, "layout")))
  cat("Each element is a simulated dataset; fit and summarise across them, e.g.\n")
  cat("  vapply(batch, function(d) <statistic of fit>, 0)\n")
  invisible(x)
}

#' Compare fitted and generating parameters
#'
#' Compares fitted parameters with the generating values from a
#' \code{simulate_*} function (carried on the data as
#' \code{attr(sim, "truth")}): item difficulties and person abilities for a
#' Rasch fit, object locations for a paired-comparison fit, rater severities
#' (with item and person measures) for a many-facet fit, set and group units
#' for a Rasch frames fit, and common object locations, panel and set units,
#' and set origins for a paired-comparison frames fit. Locations are
#' mean-centred where the model identifies them only up to an origin. An
#' externally anchored Rasch or paired-comparison fit retains its identified
#' origin, so recovery and bias are reported on the anchored scale.
#' The fit must be from the simulated model family and must have converged.
#' For a many-facet simulation, the planted rater facet must be identifiable
#' uniquely by its name or level labels.
#' EFRM set parameters are matched by their item or object membership, not by
#' the spelling of the set labels; a different fitted partition is refused.
#'
#' @param fit A fit of the simulated data (\code{\link{rasch}},
#'   \code{\link{btl}}, \code{\link{rasch_mfrm}}, \code{\link{rasch_efrm}},
#'   or \code{\link{btl_efrm}}).
#' @param sim The simulated data (from a \code{simulate_*} function).
#' @return A list of class \code{"rasch_recovery"}: \code{summary} (per
#'   parameter type: n, correlation, RMSE, bias) and \code{pieces} (the true
#'   and estimated values behind each).
#' @examples
#' d <- simulate_rasch(500, 12, seed = 1)
#' sim_recovery(rasch(d, id = "id"), d)$summary
#' @export
sim_recovery <- function(fit, sim) {
  tr <- attr(sim, "truth")
  if (is.null(tr)) stop("`sim` carries no simulation truth")
  lay <- tr$layout
  if (!is.character(lay) || length(lay) != 1L || is.na(lay) ||
      !lay %in% c("rasch", "btl", "mfrm", "efrm", "btl_efrm"))
    stop("`sim` carries an unsupported or malformed simulation layout")
  family_ok <- switch(lay,
    rasch = inherits(fit, "rasch") &&
      !inherits(fit, c("rasch_btl", "rasch_mfrm", "rasch_efrm")),
    btl = inherits(fit, "rasch_btl") &&
      !inherits(fit, "rasch_btl_efrm"),
    mfrm = inherits(fit, "rasch_mfrm"),
    efrm = inherits(fit, "rasch_efrm") &&
      !inherits(fit, "rasch_btl_efrm"),
    btl_efrm = inherits(fit, "rasch_btl_efrm"))
  if (!isTRUE(family_ok))
    stop("`fit` does not match the simulated ", lay, " model family")
  converged <- if (inherits(fit, "rasch_btl")) fit$converged else
    fit$est$converged
  if (!isTRUE(converged))
    stop("`fit` did not converge; recovery summaries are unavailable")
  pieces <- list(); centred <- list()
  add <- function(name, true, est, label = NULL, centre = FALSE) {
    true <- as.numeric(true); est <- as.numeric(est)
    keep <- is.finite(true) & is.finite(est)
    if (!any(keep)) return(invisible())
    if (centre) {                       # location-type: identified up to origin
      true <- true - mean(true[keep]); est <- est - mean(est[keep])
    }
    centred[[name]] <<- isTRUE(centre)
    pieces[[name]] <<- data.frame(
      parameter = name,
      label = if (is.null(label)) NA_character_ else as.character(label)[keep],
      true = true[keep], estimated = est[keep], stringsAsFactors = FALSE)
  }
  add_person <- function(centre = TRUE) {
    true_id <- tr$person_id %||% names(tr$theta)
    est_id <- if (!is.null(fit$person) && "id" %in% names(fit$person))
      as.character(fit$person$id) else NULL
    if (!is.null(true_id) && !is.null(est_id)) {
      at <- match(est_id, as.character(true_id))
      keep <- !is.na(at)
      if (any(keep))
        add("person ability", tr$theta[at[keep]], fit$person$theta[keep],
            est_id[keep], centre = centre)
    } else if (!is.null(fit$person) &&
               length(tr$theta) == length(fit$person$theta)) {
      # Compatibility with simulation objects created before person IDs were
      # recorded in their truth attributes.
      add("person ability", tr$theta, fit$person$theta, centre = centre)
    }
  }
  # Set labels are presentation metadata. A caller may give the same item or
  # object partition different names when fitting the simulated data. Recover
  # the fitted label from membership rather than assuming the simulator's
  # conventional set1, set2, ... names survived unchanged.
  fitted_set_labels <- function(truth_sets) {
    fallback <- names(truth_sets)
    if (is.null(fallback)) fallback <- sprintf("set%d", seq_along(truth_sets))
    fitted_map <- fit$set_of
    if ((is.null(fitted_map) || is.null(names(fitted_map))) &&
        is.data.frame(fit$objects) &&
        all(c("object", "set") %in% names(fit$objects)))
      fitted_map <- stats::setNames(as.character(fit$objects$set),
                                    as.character(fit$objects$object))
    if (is.null(fitted_map) || is.null(names(fitted_map))) return(fallback)
    mapped <- vapply(truth_sets, function(members) {
      z <- fitted_map[as.character(members)]
      if (length(z) != length(members) || anyNA(z)) return(NA_character_)
      z <- unique(as.character(z))
      if (length(z) == 1L) z else NA_character_
    }, character(1))
    if (anyNA(mapped) || anyDuplicated(mapped))
      stop("the fitted set partition does not match the simulated set ",
           "membership; set-parameter recovery would compare different ",
           "estimands")
    unname(mapped)
  }
  if (lay == "rasch") {
    anchored <- !is.null(fit$refit_spec$anchors) &&
      nrow(fit$refit_spec$anchors) > 0L
    ei <- setNames(fit$items$location, fit$items$item)
    cm <- intersect(names(tr$difficulty), names(ei))
    add("item difficulty", tr$difficulty[cm], ei[cm], cm,
        centre = !anchored)
    add_person(centre = !anchored)
  } else if (lay == "btl") {
    ot <- fit$objects
    # recovery is judged on calibrated locations; an extrapolated boundary
    # row is a reporting value, not an estimate of the planted location
    if ("extreme" %in% names(ot)) ot <- ot[!(ot$extreme %in% TRUE), ]
    eo <- setNames(ot$location, ot$object)
    cm <- intersect(names(tr$location), names(eo))
    add("object location", tr$location[cm], eo[cm], cm,
        centre = is.null(fit$anchors))
  } else if (lay == "mfrm") {
    # The simulated severity belongs to the rater facet.  Do not assume that
    # this is the first fitted facet: callers may add or reorder facets before
    # passing the fit here.  Prefer the simulator's conventional facet name;
    # otherwise use the single facet whose levels best match the planted
    # rater labels.  An ambiguous match is not a valid recovery comparison.
    fes <- fit$facet_effects
    overlap <- if (length(fes) && length(names(tr$severity)))
      vapply(fes, function(z) {
        if (!is.data.frame(z) ||
            !all(c("level", "severity") %in% names(z))) return(0L)
        sum(names(tr$severity) %in% as.character(z$level))
      }, integer(1)) else integer(0)
    named_rater <- match("rater", names(fes))
    best <- if (length(overlap) && max(overlap) > 0L)
      which(overlap == max(overlap)) else integer(0)
    fi <- if (!is.na(named_rater) && named_rater %in% best) {
      named_rater
    } else if (length(best) == 1L) best else integer(0)
    if (length(fi) != 1L)
      stop("the planted rater severity cannot be matched to one fitted ",
           "facet; recovery would be ambiguous")
    fe <- fes[[fi]]
    es <- setNames(fe$severity, fe$level)
    cm <- intersect(names(tr$severity), names(es))
    add("rater severity", tr$severity[cm], es[cm], cm, centre = TRUE)
    # per-item margins live in item_effects (fit$items holds the virtual
    # item-by-facet combinations, whose names never match)
    ie <- fit$item_effects
    if (!is.null(ie)) {
      ei <- setNames(ie$location, ie$item)
      ci <- intersect(names(tr$difficulty), names(ei))
      add("item difficulty", tr$difficulty[ci], ei[ci], ci, centre = TRUE)
    }
    add_person()
  } else if (lay == "efrm") {
    at <- fit$alpha_table
    set_labels <- fitted_set_labels(tr$item_sets)
    # units are identified up to a common scale, so compare on the centred
    # log scale (a ratio); the planted alpha is normalised the same way
    add("set unit (log)", log(tr$alpha), log(at$alpha[match(
      set_labels, at$set)]), set_labels, centre = TRUE)
    # the person-group units phi are a fitted, reported quantity too --
    # recover them, not only the set units
    if (!is.null(tr$phi) && !is.null(fit$phi_table)) {
      glab <- if (!is.null(names(tr$phi))) names(tr$phi) else
        sprintf("g%d", seq_along(tr$phi))
      ephi <- fit$phi_table$phi[match(glab, fit$phi_table$group)]
      add("group unit (log)", log(tr$phi), log(ephi), glab, centre = TRUE)
    }
  } else if (lay == "btl_efrm") {
    # The simulator and estimator use the same identifying conventions:
    # geometric-mean-one panel units, alpha_1 = 1 and kappa_1 = 0. These
    # parameters can therefore be compared directly rather than re-centred.
    ot <- fit$objects
    ev <- setNames(ot$v %||% ot$location, ot$object)
    cm <- intersect(names(tr$v), names(ev))
    add("object location", tr$v[cm], ev[cm], cm)

    ep <- setNames(fit$phi_table$phi, fit$phi_table$panel)
    cm <- intersect(names(tr$phi), names(ep))
    add("panel unit (log)", log(tr$phi[cm]), log(ep[cm]), cm)

    set_labels <- fitted_set_labels(tr$object_sets)
    ea <- setNames(fit$alpha_table$alpha, fit$alpha_table$set)
    add("set unit (log)", log(unname(tr$alpha)), log(ea[set_labels]),
        set_labels)

    ek <- setNames(fit$kappa_table$kappa, fit$kappa_table$set)
    add("set origin", unname(tr$kappa), ek[set_labels], set_labels)
  } else stop("unsupported layout: ", lay)

  # bias after centring is structurally zero for any parameter identified
  # only up to an origin/scale convention: it is not identifiable, so report
  # NA rather than a misleading ~0. Correlation and RMSE (scatter about the
  # aligned scale) remain meaningful.
  if (!length(pieces))
    stop("the fit and simulation truth have no comparable parameters")
  summ <- do.call(rbind, lapply(pieces, function(d) {
    nm <- d$parameter[1]
    data.frame(
    parameter = nm, n = nrow(d),
    correlation = if (nrow(d) > 2) .safe_cor(d$true, d$estimated) else NA_real_,
    rmse = sqrt(mean((d$estimated - d$true)^2)),
    bias = if (isTRUE(centred[[nm]])) NA_real_ else mean(d$estimated - d$true),
    stringsAsFactors = FALSE)}))
  rownames(summ) <- NULL
  structure(list(summary = summ, pieces = pieces, layout = lay),
            class = "rasch_recovery")
}

#' @export
print.rasch_recovery <- function(x, ...) {
  cat("Parameter recovery (planted vs recovered):\n")
  s <- x$summary
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-16s n=%-4d r=%.3f  RMSE=%.3f  bias=%+.3f\n",
                s$parameter[i], s$n[i], s$correlation[i], s$rmse[i], s$bias[i]))
  invisible(x)
}

#' Plot fitted against generating parameters
#'
#' One true-versus-estimated panel per parameter type, with the identity line
#' and the correlation and RMSE.
#'
#' @param x A \code{"rasch_recovery"} object.
#' @param ... Unused.
#' @return Called for its plotting side effect.
#' @examples
#' \donttest{
#' d <- simulate_rasch(300, 8, seed = 1)
#' fit <- rasch(d, id = "id")
#' plot_recovery(sim_recovery(fit, d))
#' }
#' @export
plot_recovery <- function(x, ...) {
  if (!inherits(x, "rasch_recovery"))
    stop("`x` must be a recovery result from sim_recovery()", call. = FALSE)
  np <- length(x$pieces)
  op <- par(mfrow = c(1, np), mar = c(4.2, 4.2, 2.4, 1), mgp = c(2.4, 0.7, 0),
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            cex.main = 1.0, font.main = 2)
  on.exit(par(op))
  for (i in seq_len(np)) {
    d <- x$pieces[[i]]; s <- x$summary[i, ]
    rng <- range(c(d$true, d$estimated))
    plot(d$true, d$estimated, xlim = rng, ylim = rng, xlab = "planted",
         ylab = "recovered", main = d$parameter[1], axes = FALSE,
         pch = 21, bg = .rr$blue, col = "white", cex = 1.1)
    abline(0, 1, col = .rr$red, lty = 2, lwd = 1.5)
    .rr_axis(1)
    .rr_axis(2)
    mtext(sprintf("r = %.3f   RMSE = %.2f", s$correlation, s$rmse), 3,
          line = 0.2, cex = 0.8, col = .rr$soft)
  }
}

#' Simulate paired-comparison EFRM data with differing frame units
#'
#' Generates dichotomous paired comparisons whose latent unit differs across
#' judge-panel by object-set frames -- the paired-comparison extension of the
#' extended frame of reference model (Humphry 2005) fitted by
#' \code{\link{btl_efrm}}. Objects in set \code{s} have a within-set
#' calibration location \code{beta}; their common-scale value is
#' \code{v = alpha_s beta + kappa_s}. A comparison judged in panel \code{g}
#' carries the panel unit \code{phi_g}: within a set the comparison logit is
#' \code{phi_g (beta_a - beta_b)}, across sets it is \code{phi_g (v_a - v_b)}.
#' The planted panel units, set units and origins are recovered by
#' \code{\link{btl_efrm}}.
#'
#' @param n_objects_per_set,n_sets Objects in each set and number of sets.
#' @param n_judges_per_panel,n_panels Judges in each panel and number of panels.
#'   Comparisons are balanced across panels and judges; a design with fewer
#'   comparisons than judges is refused.
#' @param reps_within Replications of each within-set object pair.
#' @param reps_cross Replications of each cross-set object pair.
#' @param panel_units Panel units \code{phi} (length \code{n_panels}); the
#'   default is all one, and any supplied vector is rescaled to geometric
#'   mean one.
#' @param set_units Set units \code{alpha} (length \code{n_sets}); the default
#'   is all one, and \code{alpha_1} is forced to one (the reference set).
#' @param set_origins Set origins \code{kappa} (length \code{n_sets}); the
#'   default is all zero, and \code{kappa_1} is forced to zero.
#' @param object_sd Spread of the within-set calibration locations.
#' @param erratic_judges Proportion of judges who choose between the two
#'   objects at random.
#' @param seed Optional non-negative whole-number RNG seed.
#' @return A data frame of class \code{"rasch_sim"} with columns
#'   \code{object_a}, \code{object_b}, \code{winner}, \code{judge} and
#'   \code{panel}, and \code{attr(x, "truth")} holding the common-scale values
#'   \code{v}, the per-set \code{beta}, the units \code{phi}, \code{alpha},
#'   \code{kappa}, and the \code{object_sets} map to pass to
#'   \code{\link{btl_efrm}}.
#' @examples
#' d <- simulate_btl_efrm(6, 2, set_units = c(1, 1.4), seed = 1)
#' bt <- btl_efrm(d, "object_a", "object_b", winner = "winner",
#'                judge = "judge", panels = "panel",
#'                object_sets = attr(d, "truth")$object_sets,
#'                se_method = "conditional")
#' bt$alpha_table   # recovers the ~1.4 set unit
#' @export
simulate_btl_efrm <- function(n_objects_per_set = 8, n_sets = 2,
                              n_judges_per_panel = 6, n_panels = 2,
                              reps_within = 20, reps_cross = 20,
                              panel_units = NULL, set_units = NULL,
                              set_origins = NULL, object_sd = 1,
                              seed = NULL, erratic_judges = 0) {
  seed <- .sim_seed(seed)
  if (!is.null(seed)) {
    .old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(.old_stream), add = TRUE)
    set.seed(seed)
  }
  S <- .sim_count(n_sets, "n_sets")
  G <- .sim_count(n_panels, "n_panels")
  Kp <- .sim_count(n_objects_per_set, "n_objects_per_set", 2L)
  Jp <- .sim_count(n_judges_per_panel, "n_judges_per_panel")
  reps_within <- .sim_count(reps_within, "reps_within")
  reps_cross <- .sim_count(reps_cross, "reps_cross")
  object_sd <- .sim_scalar(object_sd, "object_sd", lower = 0,
                           lower_open = TRUE)
  erratic_judges <- .sim_scalar(erratic_judges, "erratic_judges",
                                lower = 0, upper = 1)

  phi <- if (is.null(panel_units)) rep(1, G) else panel_units
  if (!is.numeric(phi) || is.complex(phi) || !is.null(dim(phi)) ||
      !is.null(oldClass(phi)) ||
      length(phi) != G || any(!is.finite(phi) | phi <= 0))
    stop("panel_units must contain n_panels positive finite values")
  phi <- as.numeric(phi)
  alpha <- if (is.null(set_units)) rep(1, S) else set_units
  if (!is.numeric(alpha) || is.complex(alpha) || !is.null(dim(alpha)) ||
      !is.null(oldClass(alpha)) ||
      length(alpha) != S || any(!is.finite(alpha) | alpha <= 0))
    stop("set_units must contain n_sets positive finite values")
  alpha <- as.numeric(alpha)
  kappa <- if (is.null(set_origins)) rep(0, S) else set_origins
  if (!is.numeric(kappa) || is.complex(kappa) || !is.null(dim(kappa)) ||
      !is.null(oldClass(kappa)) ||
      length(kappa) != S || any(!is.finite(kappa)))
    stop("set_origins must contain n_sets finite values")
  kappa <- as.numeric(kappa)
  # a unit and an origin are relative to the other frames: with one frame
  # the normalisation below would silently return them to the reference
  # while the recorded truth claimed the requested values
  if (G == 1L && !isTRUE(all.equal(phi, 1)))
    stop("`panel_units` must be 1 when `n_panels` is 1: a panel unit is ",
         "relative to the other panels, so none can be planted")
  if (S == 1L && !isTRUE(all.equal(alpha, 1)))
    stop("`set_units` must be 1 when `n_sets` is 1: a set unit is relative ",
         "to the other sets, so none can be planted")
  if (S == 1L && !isTRUE(all.equal(kappa, 0)))
    stop("`set_origins` must be 0 when `n_sets` is 1: an origin is relative ",
         "to the other sets, so none can be planted")
  phi <- phi / exp(mean(log(phi)))                    # geometric mean one
  alpha <- alpha / alpha[1]                            # alpha_1 = 1
  kappa <- kappa - kappa[1]                            # kappa_1 = 0

  set_nm <- sprintf("set%d", seq_len(S))
  panel_nm <- sprintf("panel%d", seq_len(G))
  objs_by_set <- lapply(seq_len(S), function(s) sprintf("S%dO%02d", s, seq_len(Kp)))
  names(objs_by_set) <- set_nm
  beta <- numeric(0)
  for (s in seq_len(S)) {
    bs <- as.numeric(scale(seq_len(Kp))) * object_sd   # sum-zero, spread object_sd
    beta <- c(beta, setNames(bs, objs_by_set[[s]]))
  }
  set_of <- setNames(rep(seq_len(S), each = Kp), unlist(objs_by_set))
  v <- alpha[set_of] * beta + kappa[set_of]

  judges <- sprintf("J%03d", seq_len(G * Jp))
  panel_of <- setNames(panel_nm[rep(seq_len(G), each = Jp)], judges)

  # assemble the object pairs (within each set, then across every set pair)
  aa <- bb <- character(0)
  for (s in seq_len(S)) {
    pr <- t(utils::combn(objs_by_set[[s]], 2))
    aa <- c(aa, rep(pr[, 1], reps_within)); bb <- c(bb, rep(pr[, 2], reps_within))
  }
  if (S > 1L) for (i in seq_len(S - 1L)) for (j in (i + 1L):S) {
    grid <- expand.grid(oa = objs_by_set[[i]], ob = objs_by_set[[j]],
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    aa <- c(aa, rep(grid$oa, reps_cross)); bb <- c(bb, rep(grid$ob, reps_cross))
  }
  R <- length(aa)
  # Balance panels within every set-pair stratum and then balance judges
  # within their panels. Requested panels and judges therefore cannot vanish
  # from a finite generated sample merely by chance.
  stratum <- paste(pmin(set_of[aa], set_of[bb]),
                   pmax(set_of[aa], set_of[bb]), sep = ":")
  if (R < length(judges))
    stop("the generated design has ", R, " comparisons, fewer than the ",
         length(judges), " requested judges")
  pan <- character(R)
  # Continue one randomly ordered panel cycle across strata. Counts differ by
  # at most one both overall and, when a stratum is large enough, within it.
  panel_cycle <- rep(sample(panel_nm, G, replace = FALSE), length.out = R)
  cursor <- 0L
  for (ss in unique(stratum)) {
    rows <- which(stratum == ss)
    take <- cursor + seq_along(rows)
    pan[rows] <- sample(panel_cycle[take], length(rows), replace = FALSE)
    cursor <- cursor + length(rows)
  }
  jd <- character(R)
  for (g in panel_nm) {
    rows <- which(pan == g)
    pool <- names(panel_of)[panel_of == g]
    jd[rows] <- .sim_balanced_ids(pool, length(rows),
                                  paste0("judges in ", g))
  }

  n_erratic <- .sim_planted_count(erratic_judges, length(judges))
  if (n_erratic) {
    # Distribute erratic judges across panels as evenly as their number
    # permits, then select them at random within each panel.
    take <- rep(n_erratic %/% G, G)
    extra <- n_erratic %% G
    if (extra) {
      extra_panels <- sample(seq_len(G), extra, replace = FALSE)
      take[extra_panels] <- take[extra_panels] + 1L
    }
    erratic <- unlist(lapply(seq_len(G), function(g) {
      pool <- names(panel_of)[panel_of == panel_nm[g]]
      if (take[g]) sample(pool, take[g], replace = FALSE) else character(0)
    }), use.names = FALSE)
  } else erratic <- character(0)
  # frame-dependent logit: within-set uses beta, cross-set uses the common v
  same <- set_of[aa] == set_of[bb]
  lp <- ifelse(same, phi[match(pan, panel_nm)] * (beta[aa] - beta[bb]),
               phi[match(pan, panel_nm)] * (v[aa] - v[bb]))
  win_a <- stats::runif(R) < stats::plogis(lp)
  bad <- jd %in% erratic
  if (any(bad)) win_a[bad] <- sample(c(FALSE, TRUE), sum(bad), replace = TRUE)
  d <- data.frame(object_a = aa, object_b = bb,
                  winner = ifelse(win_a, aa, bb),
                  judge = jd, panel = unname(pan),
                  stringsAsFactors = FALSE)
  rownames(d) <- NULL

  planted <- character(0)
  if (any(phi != 1)) planted <- c(planted, sprintf(
    "panel units phi = (%s)", paste(sprintf("%.2f", phi), collapse = ", ")))
  if (any(alpha != 1)) planted <- c(planted, sprintf(
    "set units alpha = (%s)", paste(sprintf("%.2f", alpha), collapse = ", ")))
  if (any(kappa != 0)) planted <- c(planted, sprintf(
    "set origins kappa = (%s)", paste(sprintf("%.2f", kappa), collapse = ", ")))
  if (length(erratic)) planted <- c(planted, sprintf(
    "%d erratic judge(s): %s", length(erratic), paste(erratic, collapse = ", ")))
  if (!length(planted)) planted <- "equal units (phi = alpha = 1)"

  attr(d, "truth") <- list(
    layout = "btl_efrm",
    description = sprintf("%d objects (%d sets) x %d panels, %d comparisons",
                          S * Kp, S, G, R),
    v = v, beta = beta, phi = setNames(phi, panel_nm),
    alpha = setNames(alpha, set_nm), kappa = setNames(kappa, set_nm),
    set_of = set_of, object_sets = objs_by_set, erratic = erratic,
    planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}
