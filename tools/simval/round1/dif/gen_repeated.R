# Custom generator for repeated-measures (within-person) DIF designs.
# Follows the same construction as the dif_anova() roxygen mixed-design
# example (dichotomous Rasch, logistic response function), extended to
# allow: a between-person factor with g_levels levels, a within-person
# occasion factor with occ_levels levels, an optional uniform-DIF shift on
# one item for one (non-reference) GROUP level, and an optional uniform
# shift on one item for one (non-reference) OCCASION level (applied to
# every person at that occasion -- a break of invariance over time, i.e.
# within-person / occasion DIF).
gen_repeated <- function(N = 300, n_items = 8, g_levels = 2, occ_levels = 2,
                          group_dif_item = NULL, group_dif_shift = 0,
                          occ_dif_item = NULL, occ_dif_shift = 0,
                          theta_sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  d <- seq(-2, 2, length.out = n_items)
  theta <- rnorm(N, 0, theta_sd)
  group <- factor(rep(paste0("g", seq_len(g_levels)), length.out = N))
  occ_lv <- paste0("T", seq_len(occ_levels))

  make_wave <- function(occ) {
    shift <- matrix(0, N, n_items)
    if (!is.null(group_dif_item) && g_levels >= 2)
      shift[group == paste0("g", g_levels), group_dif_item] <- group_dif_shift
    if (!is.null(occ_dif_item) && occ == occ_lv[occ_levels])
      shift[, occ_dif_item] <- occ_dif_shift
    matrix(rbinom(N * n_items, 1,
                   plogis(outer(theta, d, "-") - shift)), N, n_items)
  }
  Xl <- lapply(occ_lv, make_wave)
  X <- do.call(rbind, Xl)
  colnames(X) <- sprintf("I%02d", seq_len(n_items))
  data.frame(X, person = rep(seq_len(N), occ_levels),
             group = rep(as.character(group), occ_levels),
             occasion = rep(occ_lv, each = N),
             stringsAsFactors = FALSE)
}
