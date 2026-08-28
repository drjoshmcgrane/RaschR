# rasch :: Guttman display
# ===========================================================================
# The Guttman (1944) scalogram: persons and items are sorted by location, so
# under a perfectly deterministic (Guttman) response pattern the matrix splits
# into a block of high scores at the top-left and low scores at the
# bottom-right. The coefficient of reproducibility measures how close the
# observed data is to that deterministic ideal. The probabilistic Rasch model
# relaxes the deterministic Guttman pattern; the display remains a useful
# qualitative check on ordering and aberrant responses.
# ===========================================================================

#' Guttman-ordered response matrix and reproducibility
#'
#' Orders persons by descending location and dichotomous items by ascending
#' location (the Guttman scalogram arrangement), computing the
#' reproducibility coefficient against the deterministic pattern implied by each
#' person's total score. PCM category steps can interleave across items, so
#' polytomous fits are rejected rather than forced into an invalid whole-item
#' deterministic order.
#'
#' @param fit A fitted object from \code{\link{rasch}} whose columns form one
#'   administered item set. Expanded EFRM and MFRM response-cell matrices are
#'   not accepted; one-cell-per-item reductions are.
#' @return A list with the ordered score matrix \code{matrix} (persons by
#'   items, row and column names carrying the ID and item labels), the person
#'   and item orderings, and the coefficient of reproducibility \code{CR}.
#' @references
#' Guttman, L. (1944). A basis for scaling qualitative data. American
#' Sociological Review, 9(2), 139--150.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(200 * 6, 1, plogis(outer(rnorm(200), d, "-"))), 200, 6)
#' colnames(X) <- paste0("I", 1:6)
#' guttman_table(rasch(X))$CR
#' @export
guttman_table <- function(fit) {
  if (!inherits(fit, "rasch")) stop("guttman_table needs a rasch fit")
  structural <- inherits(fit, c("rasch_efrm", "rasch_mfrm"))
  if (!.classical_design_applicable(fit))
    stop("the whole-item scalogram is not defined when an item is ",
         "represented by several frame or facet response cells; construct ",
         "it for one observable item design instead")
  X <- fit$X; m <- fit$m
  if (structural) {
    item_names <- fit$virtual_map$item[
      match(colnames(X), fit$virtual_map$vkey)]
    if (anyNA(item_names) || anyDuplicated(item_names))
      stop("the one-cell-per-item reduction could not be matched to its ",
           "item names")
    colnames(X) <- item_names
  }
  if (any(m != 1L))
    stop("guttman_table currently supports dichotomous items only; ",
         "polytomous category steps can interleave and require a step-level ",
         "scalogram")
  th <- fit$person$theta
  porder <- order(th, fit$person$raw, decreasing = TRUE)
  iorder <- order(fit$items$location)
  G <- X[porder, iorder, drop = FALSE]
  rownames(G) <- as.character(fit$person$id)[porder]
  colnames(G) <- colnames(X)[iorder]

  # coefficient of reproducibility: 1 - errors / responses, where the
  # error count compares each observed score with the deterministic pattern
  # that reproduces the person's total from the easiest categories up.
  errors <- 0L; total <- 0L
  mi <- m[iorder]
  for (r in seq_len(nrow(G))) {
    obs <- G[r, ]; ok <- !is.na(obs)
    if (!any(ok)) next
    R <- sum(obs[ok])
    # deterministic fill: allocate the raw score to the easiest items first
    pred <- numeric(sum(ok)); cap <- mi[ok]; left <- R
    for (j in seq_along(pred)) { take <- min(cap[j], left); pred[j] <- take; left <- left - take }
    errors <- errors + sum(abs(obs[ok] - pred))
    total  <- total + sum(cap)
  }
  structure(class = "rasch_guttman", list(matrix = G, person_order = porder, item_order = iorder,
       CR = 1 - errors / total, errors = errors, responses = total))
}

#' Plot the Guttman scalogram
#'
#' Displays the dichotomous Guttman-ordered response matrix as a heatmap (dark
#' for a score of one), with persons sorted by location down the rows and items by
#' location across the columns. The coefficient of reproducibility is shown in
#' the subtitle.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param max_persons Persons are thinned to at most this many evenly spaced
#'   rows for legibility on large samples.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(200 * 6, 1, plogis(outer(rnorm(200), d, "-"))), 200, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_guttman(rasch(X))
#' @export
plot_guttman <- function(fit, max_persons = 80) {
  max_persons <- .check_whole(max_persons, "max_persons", 1)
  g <- guttman_table(fit); G <- g$matrix; m <- max(fit$m)
  N <- nrow(G)
  if (N > max_persons) G <- G[round(seq(1, N, length.out = max_persons)), , drop = FALSE]
  L <- ncol(G); rN <- nrow(G)
  pal <- colorRampPalette(c("#f1f5f9", "#93c5fd", "#1e3a8a"))(m + 1)
  op <- par(mar = c(6, 1.5, 5, 1.5), col.main = .rr$ink, font.main = 2,
            cex.main = 1.15)
  on.exit(par(op))
  image(seq_len(L), seq_len(rN), t(G[rN:1, , drop = FALSE]), col = pal,
        zlim = c(0, m), axes = FALSE, xlab = "", ylab = "", main = "")
  mtext(sprintf("persons by location (high at top), items easy to hard; coefficient of reproducibility %.3f",
                g$CR), side = 3, line = 2.5, adj = 0, cex = 0.8, col = .rr$soft)
  axis(1, seq_len(L), colnames(G), las = 2, cex.axis = 0.62, col = NA, col.ticks = NA)
  na_cells <- which(is.na(t(G[rN:1, , drop = FALSE])), arr.ind = TRUE)
  if (nrow(na_cells)) points(na_cells[, 1], na_cells[, 2], pch = 4,
                             cex = 0.5, col = .rr$red)
  box(col = .rr$grid)
  legend("topright", inset = c(0, -0.1), xpd = TRUE, horiz = TRUE,
         legend = 0:m, fill = pal, border = NA, bty = "n", cex = 0.75,
         title = "score")
  invisible(NULL)
}


#' @export
print.rasch_guttman <- function(x, ...) {
  cat(sprintf("Guttman scalogram: %d persons x %d items (difficulty-ordered)\n",
              nrow(x$matrix), ncol(x$matrix)))
  cat(sprintf("Coefficient of reproducibility %.3f (%d error(s) in %d responses)\n",
              x$CR, x$errors, x$responses))
  cat("(the ordered response matrix is on $matrix; plot with plot_guttman())\n")
  invisible(x)
}
