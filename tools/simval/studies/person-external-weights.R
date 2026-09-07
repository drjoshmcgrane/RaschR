# STUDY: person-external-weights
#
# Conditional bias, sandwich-SE calibration and interval coverage for person
# estimates formed with externally imposed item weights. Item parameters are
# held at their generating values so the study isolates the weighted person
# estimator from calibration error. Run from the package root:
#
#   Rscript tools/simval/studies/person-external-weights.R

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

STUDY <- "person-external-weights"
N <- 5000L

make_bank <- function(model, n_items) {
  location <- seq(-2, 2, length.out = n_items)
  tau <- if (model == "PCM")
    lapply(location, function(b) b + c(-0.9, 0, 0.9)) else
    lapply(location, identity)
  names(tau) <- sprintf("I%02d", seq_len(n_items))
  tau
}

run_one <- function(theta_true, model, weight_shape, varying_units = FALSE,
                    seed) {
  L <- if (model == "PCM") 15L else 30L
  tau <- make_bank(model, L)
  disc <- if (varying_units) rep(c(0.75, 1.35), length.out = L) else rep(1, L)
  q <- switch(weight_shape,
    equal = rep(1, L),
    moderate = c(rep(2, ceiling(L / 2)), rep(0.5, floor(L / 2))),
    strong = c(rep(3, ceiling(L / 2)), rep(0.2, floor(L / 2))),
    zero = c(rep(1, ceiling(2 * L / 3)), rep(0, floor(L / 3))))

  set.seed(seed)
  X <- vapply(seq_len(L), function(j)
    rasch:::.sim_item(rep(theta_true, N), tau[[j]], disc[j]), integer(N))
  colnames(X) <- names(tau)
  fit <- list(X = X, tau_list = tau, disc = disc,
              person = data.frame(id = sprintf("P%05d", seq_len(N))),
              factors = NULL, est = list(converged = TRUE))
  class(fit) <- "rasch"
  z <- weighted_person_estimates(fit, stats::setNames(q, names(tau)))
  ok <- is.finite(z$theta) & is.finite(z$se) & z$se > 0
  est <- z$theta[ok]
  se <- z$se[ok]
  scenario <- sprintf("%s; weights=%s; units=%s; theta=%+.1f", model,
                      weight_shape, if (varying_units) "varying" else "equal",
                      theta_true)
  sv_row(
    STUDY, scenario, "person location", n_reps = sum(ok),
    n_attempted = N, n_refused = N - sum(ok), n_nonconv = 0L,
    effect = theta_true, bias = mean(est - theta_true), emp_sd = stats::sd(est),
    mean_se = mean(se), coverage95 = mean(abs(est - theta_true) <= 1.96 * se),
    notes = paste("fixed generating calibration; externally weighted Warm",
                  "score with sandwich standard error"))
}

design <- rbind(
  expand.grid(theta = c(-1, 0, 1), model = "dichotomous",
              weights = c("equal", "moderate", "strong", "zero"),
              varying = FALSE, stringsAsFactors = FALSE),
  data.frame(theta = c(-1, 0, 1), model = "PCM", weights = "strong",
             varying = FALSE, stringsAsFactors = FALSE),
  data.frame(theta = c(-1, 0, 1), model = "dichotomous", weights = "strong",
             varying = TRUE, stringsAsFactors = FALSE))

rows <- lapply(seq_len(nrow(design)), function(i)
  run_one(design$theta[i], design$model[i], design$weights[i],
          design$varying[i], seed = 880000L + i))

sv_write(do.call(rbind, rows), "person-external-weights")
