#' Draw a Wright map with WrightMap
#'
#' Prepares person estimates and item locations from a fitted model and passes
#' them to \code{WrightMap::wrightMap()}. Person panels may be formed from
#' variables retained in the fit or supplied as a matrix. Item panels use the
#' \code{item.groups} facility in WrightMap 1.5.
#'
#' For partial credit and rating scale models, \code{type = "thresholds"}
#' displays each item's estimated category thresholds. \code{type =
#' "locations"} displays one location per item. In MFRM and EFRM fits, the
#' rows are the calibrated item-by-facet or item-by-frame response columns.
#' For EFRM fits, \code{person_panels = "groups"} and \code{item_panels =
#' "sets"} use the fitted frame design; both may be specified together.
#'
#' @param fit A fitted object from \code{\link{rasch}}, including explanatory,
#'   MFRM and EFRM fits. Comparative-judgement models do not estimate person
#'   locations and are not supported.
#' @param type Plot category \code{"thresholds"} or one \code{"locations"}
#'   estimate per item.
#' @param person_panels Optional person-panel specification. Supply the name
#'   of one or more variables retained in \code{fit$person}, a vector or factor
#'   of panel memberships with one value per person, or a numeric matrix or
#'   data frame whose columns contain the estimates for the panels. Missing
#'   estimates are permitted. For EFRM fits, \code{"groups"} uses the fitted
#'   person groups.
#' @param item_panels Optional vector or factor assigning each item row to a
#'   panel. A named vector is matched to item names; an unnamed vector is used
#'   in item order. A named list may instead map panel names to item names.
#'   For EFRM fits, use \code{"sets"}, \code{"groups"}, or
#'   \code{c("sets", "groups")} to arrange the calibrated response columns by
#'   item set, person group, or frame. This option requires WrightMap 1.5 or
#'   later.
#' @param ... Further arguments passed to \code{WrightMap::wrightMap()}, such
#'   as \code{person.side}, \code{item.side}, \code{main.title}, or graphical
#'   settings.
#' @return Invisibly, the threshold matrix returned by
#'   \code{WrightMap::wrightMap()} when its \code{return.thresholds} argument
#'   is true; otherwise \code{NULL}.
#' @seealso \code{\link{plot_wright}}
#' @references Torres Irribarra, D., and Freund, R. (2025). WrightMap: IRT
#'   item-person map with ConQuest integration. R package version 1.5.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' fit <- rasch(X)
#' if (requireNamespace("WrightMap", quietly = TRUE)) {
#'   wright_map(fit)
#' }
#' @export
wright_map <- function(fit, type = c("thresholds", "locations"),
                       person_panels = NULL, item_panels = NULL, ...) {
  if (!requireNamespace("WrightMap", quietly = TRUE)) {
    stop("Package 'WrightMap' is required. Install it with ",
         "install.packages(\"WrightMap\") before using wright_map().",
         call. = FALSE)
  }
  type <- match.arg(type)
  dat <- .wright_map_data(fit, type = type, person_panels = person_panels,
                          item_panels = item_panels)

  dots <- list(...)
  reserved <- intersect(names(dots), c("thetas", "thresholds", "item.groups"))
  if (length(reserved)) {
    stop("Use fit, person_panels, and item_panels rather than passing ",
         paste(reserved, collapse = ", "), " through '...'.", call. = FALSE)
  }
  if (!is.null(dat$item_panels) &&
      !("item.groups" %in% names(formals(WrightMap::wrightMap)))) {
    stop("Item panels require WrightMap 1.5 or later. Update WrightMap before ",
         "using item_panels.", call. = FALSE)
  }

  args <- c(list(thetas = dat$persons, thresholds = dat$items), dots)
  if (ncol(dat$persons) > 1L && !("dim.names" %in% names(dots)))
    args$dim.names <- colnames(dat$persons)
  if (!is.null(dat$item_panels)) args$item.groups <- dat$item_panels
  invisible(do.call(WrightMap::wrightMap, args))
}

.wright_map_data <- function(fit, type = c("thresholds", "locations"),
                             person_panels = NULL, item_panels = NULL) {
  if (inherits(fit, "rasch_btl"))
    stop("Wright maps require person estimates; comparative-judgement models ",
         "estimate object locations but not judge locations.", call. = FALSE)
  if (!inherits(fit, "rasch"))
    stop("fit must be a fitted rasch object.", call. = FALSE)
  type <- match.arg(type)

  theta <- fit$person$theta
  if (is.null(theta) || !length(theta))
    stop("The fitted object does not contain person estimates.", call. = FALSE)

  persons <- .wright_person_panels(fit, person_panels)
  items <- .wright_item_matrix(fit, type)
  panels <- .wright_item_panels(fit, item_panels, rownames(items))
  list(persons = persons, items = items, item_panels = panels)
}

.wright_person_panels <- function(fit, panels) {
  theta <- as.numeric(fit$person$theta)
  n <- length(theta)
  if (is.null(panels)) return(matrix(theta, ncol = 1L))

  if (is.matrix(panels) || is.data.frame(panels)) {
    ans <- as.matrix(panels)
    if (!is.numeric(ans))
      stop("A person-panel matrix must contain numeric estimates.", call. = FALSE)
    if (!nrow(ans)) stop("person_panels has no rows.", call. = FALSE)
    return(ans)
  }

  if (inherits(fit, "rasch_efrm") && is.character(panels) &&
      length(panels) == 1L && panels == "groups") {
    panels <- fit$frame_group[1L]
  }
  if (is.character(panels) && length(panels) <= ncol(fit$person) &&
      all(panels %in% names(fit$person))) {
    gdat <- fit$person[panels]
    group <- if (ncol(gdat) == 1L) gdat[[1L]] else
      interaction(gdat, drop = TRUE, sep = " x ")
  } else {
    if (length(panels) != n)
      stop("person_panels must have one value per person.", call. = FALSE)
    group <- panels
  }
  if (anyNA(group))
    stop("Person-panel memberships cannot be missing.", call. = FALSE)
  group <- if (is.factor(group)) droplevels(group) else
    factor(group, levels = unique(group))
  if (nlevels(group) < 1L)
    stop("person_panels defines no panels.", call. = FALSE)
  ans <- vapply(levels(group), function(g) {
    out <- rep(NA_real_, n)
    out[group == g] <- theta[group == g]
    out
  }, numeric(n))
  if (is.null(dim(ans))) ans <- matrix(ans, ncol = 1L)
  colnames(ans) <- levels(group)
  ans
}

.wright_item_matrix <- function(fit, type) {
  item_names <- as.character(fit$items$item)
  if (!length(item_names))
    stop("The fitted object does not contain item estimates.", call. = FALSE)

  if (type == "locations") {
    loc <- fit$items$location
    if (is.null(loc) || length(loc) != length(item_names))
      stop("The fitted object does not contain item locations.", call. = FALSE)
    return(matrix(as.numeric(loc), ncol = 1L,
                  dimnames = list(item_names, "Location")))
  }

  thr <- fit$thresholds
  needed <- c("item", "k", "tau")
  if (is.null(thr) || !all(needed %in% names(thr)))
    stop("The fitted object does not contain item thresholds.", call. = FALSE)
  kmax <- max(as.integer(thr$k), na.rm = TRUE)
  ans <- matrix(NA_real_, nrow = length(item_names), ncol = kmax,
                dimnames = list(item_names, paste0("Threshold ", seq_len(kmax))))
  ii <- as.integer(thr$item)
  kk <- as.integer(thr$k)
  keep <- is.finite(ii) & is.finite(kk) & ii >= 1L &
    ii <= nrow(ans) & kk >= 1L & kk <= ncol(ans)
  ans[cbind(ii[keep], kk[keep])] <- as.numeric(thr$tau[keep])
  ans
}

.wright_item_panels <- function(fit, panels, item_names) {
  if (is.null(panels)) return(NULL)

  if (is.list(panels)) {
    if (is.null(names(panels)) || any(!nzchar(names(panels))))
      stop("A list supplied to item_panels must name every panel.", call. = FALSE)
    members <- as.character(unlist(panels, use.names = FALSE))
    if (anyDuplicated(members) || !setequal(members, item_names))
      stop("A list supplied to item_panels must contain each item exactly once.",
           call. = FALSE)
    lookup <- stats::setNames(rep(names(panels), lengths(panels)), members)
    panels <- unname(lookup[item_names])
  }

  if (inherits(fit, c("rasch_mfrm", "rasch_efrm")) &&
      is.character(panels)) {
    map <- fit$virtual_map
    fields <- panels
    fields[fields == "sets"] <- "set"
    fields[fields == "groups"] <- "group"
    if (!is.null(map) && length(fields) <= 2L &&
        all(fields %in% names(map))) {
      if (nrow(map) != length(item_names) ||
          !identical(as.character(map$vkey), item_names))
        stop("The fitted design map is not aligned with the response columns.",
             call. = FALSE)
      values <- lapply(fields, function(x) as.character(map[[x]]))
      panels <- if (length(values) == 1L) values[[1L]] else
        interaction(values, drop = TRUE, sep = " x ")
    }
  }

  if (!is.atomic(panels) || is.matrix(panels))
    stop("item_panels must be a vector or factor.", call. = FALSE)
  if (!is.null(names(panels))) {
    if (anyDuplicated(names(panels)) || !all(item_names %in% names(panels)))
      stop("Named item_panels must contain each item name exactly once.",
           call. = FALSE)
    panels <- panels[item_names]
  } else if (length(panels) != length(item_names)) {
    stop("item_panels must have one value per item.", call. = FALSE)
  }
  if (anyNA(panels))
    stop("Item-panel memberships cannot be missing.", call. = FALSE)
  if (is.factor(panels)) droplevels(panels) else
    factor(panels, levels = unique(panels))
}
