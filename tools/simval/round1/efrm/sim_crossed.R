# custom generator: single item-set, 2x2 crossed person design (region x cohort)
# phi varies by region only (region_effect); cohort has no effect (null factor)
sim_crossed_phi <- function(n_per_cell = 150, n_items = 10, region_effect = log(1.3),
                             cohort_effect = 0, theta_sd = 1.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cells <- expand.grid(region = c("N", "S"), cohort = c("Y", "O"),
                        stringsAsFactors = FALSE)
  delta <- seq(-1.5, 1.5, length.out = n_items)
  names(delta) <- sprintf("I%02d", seq_len(n_items))
  rows <- list(); phi_cell <- numeric(nrow(cells))
  for (c in seq_len(nrow(cells))) {
    phi <- exp(region_effect * (cells$region[c] == "S") +
               cohort_effect * (cells$cohort[c] == "O"))
    phi_cell[c] <- phi
    theta <- rnorm(n_per_cell, 0, theta_sd)
    X <- sapply(delta, function(dd)
      rbinom(n_per_cell, 1, plogis(phi * (theta - dd))))
    rows[[c]] <- data.frame(X, region = cells$region[c], cohort = cells$cohort[c],
                             theta_true = theta, stringsAsFactors = FALSE)
  }
  # centre phi_cell geometrically for reporting truth (matches package convention)
  phi_cell <- phi_cell / exp(mean(log(phi_cell)))
  d <- do.call(rbind, rows)
  attr(d, "phi_cell") <- phi_cell
  attr(d, "cells") <- cells
  attr(d, "item_names") <- names(delta)
  d
}
