# STUDY: alpha-npml-coverage
#
# Sampling calibration of the EFRM semiparametric item-set link under normal
# and wide bimodal person distributions. Each replicate re-estimates the grid
# masses, and the hybrid bootstrap re-estimates them after resampling persons.
# Run from the package root. Set SV_CORES on Unix-like systems to parallelise.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
STUDY <- "alpha-npml-coverage"

gen <- function(n, ratio, population, seed, ips = 8L) {
  set.seed(seed)
  alpha <- c(set1 = ratio^(-0.5), set2 = ratio^0.5)
  if (population == "contrasting groups") {
    n1 <- n %/% 2L
    group <- c(rep("lower", n1), rep("upper", n - n1))
    theta <- c(rnorm(n1, -1.4, 0.8),
               sample(c(0.7, 2.4), n - n1, replace = TRUE) +
                 rnorm(n - n1, 0, 0.45))
    phi <- c(lower = 1.5^(-0.5), upper = 1.5^0.5)
  } else {
    group <- rep("g1", n); phi <- c(g1 = 1)
    theta <- if (population == "normal") rnorm(n, 0, 1.3) else
      sample(c(-2.2, 2.2), n, replace = TRUE) + rnorm(n, 0, 0.6)
  }
  sets <- rep(names(alpha), each = ips)
  delta <- rep(seq(-1.5, 1.5, length.out = ips), 2L)
  X <- vapply(seq_along(sets), function(i)
    rbinom(n, 1, plogis(alpha[sets[i]] * phi[group] *
                         (theta - delta[i]))), numeric(n))
  colnames(X) <- sprintf("%sI%02d", sets, seq_along(sets))
  d <- data.frame(id = sprintf("P%05d", seq_len(n)), X, group = group,
                  check.names = FALSE)
  target <- ips + ceiling(ips / 2)
  list(data = d, sets = split(colnames(X), sets),
       item = colnames(X)[target], item_truth = delta[target])
}

one_cell <- function(population, ratio, seed0, R = 150L, n = 500L,
                     boot_reps = 80L, se_method = "hybrid") {
  one <- function(r) {
    d <- gen(n, ratio, population, seed0 + r)
    f <- tryCatch(rasch_efrm(d$data, item_sets = d$sets, groups = "group",
                             id = "id", boot_reps = boot_reps,
                             se_method = se_method),
                  error = function(e) NULL)
    if (is.null(f)) return(c(est = NA, se = NA, p = NA,
                             item_est = NA, item_se = NA, item_truth = NA,
                             phi_est = NA, phi_se = NA, phi_truth = NA,
                             refused = 1, nonconv = 0))
    link_converged <- nrow(f$linking$alpha_edges) == 0L ||
      all(f$linking$alpha_edges$converged %in% TRUE)
    if (!isTRUE(f$est$converged) || !link_converged)
      return(c(est = NA, se = NA, p = NA,
               item_est = NA, item_se = NA, item_truth = NA,
               phi_est = NA, phi_se = NA, phi_truth = NA,
               refused = 0, nonconv = 1))
    if (!identical(f$se_method, se_method))
      return(c(est = NA, se = NA, p = NA,
               item_est = NA, item_se = NA, item_truth = NA,
               phi_est = NA, phi_se = NA, phi_truth = NA,
               refused = 1, nonconv = 0))
    at <- f$alpha_table
    it <- f$item_arbitrary[f$item_arbitrary$item == d$item, ]
    pt <- f$phi_table[f$phi_table$group == "upper", ]
    c(est = log(at$alpha[at$set == "set2"]),
      se = at$se_log_alpha[at$set == "set2"],
      p = f$efrm_vs_rasch$unit_omnibus$p[
        f$efrm_vs_rasch$unit_omnibus$term == "set units (alpha)"],
      item_est = it$location, item_se = it$se, item_truth = d$item_truth,
      phi_est = if (nrow(pt)) log(pt$phi) else NA_real_,
      phi_se = if (nrow(pt)) pt$se_log_phi else NA_real_,
      phi_truth = if (nrow(pt)) log(1.5) / 2 else NA_real_,
      refused = 0, nonconv = 0)
  }
  cores <- suppressWarnings(as.integer(Sys.getenv("SV_CORES", "1")))
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  z <- if (cores > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(R), one, mc.cores = cores, mc.set.seed = FALSE)
  else lapply(seq_len(R), one)
  z <- do.call(rbind, z)
  ok <- is.finite(z[, "est"]) & is.finite(z[, "se"]) & z[, "se"] > 0
  truth <- log(ratio) / 2
  out_alpha <- sv_row(STUDY,
    sprintf("%s persons; ratio %.1f; N=%d; %s SE; %d bootstrap replicates",
            population, ratio, n, se_method, boot_reps),
    "log alpha[set2] bias, SE calibration, coverage and omnibus rejection",
    sum(ok), effect = ratio, bias = mean(z[ok, "est"]) - truth,
    emp_sd = sd(z[ok, "est"]), mean_se = mean(z[ok, "se"]),
    coverage95 = mean(abs(z[ok, "est"] - truth) <= 1.96 * z[ok, "se"]),
    type1 = if (ratio == 1) mean(z[ok, "p"] < 0.05) else NA_real_,
    power = if (ratio != 1) mean(z[ok, "p"] < 0.05) else NA_real_,
    n_attempted = R, n_refused = sum(z[, "refused"]),
    n_nonconv = sum(z[, "nonconv"]))
  ok_item <- is.finite(z[, "item_est"]) & is.finite(z[, "item_se"]) &
    z[, "item_se"] > 0
  out_item <- sv_row(STUDY,
    sprintf("%s persons; ratio %.1f; N=%d; %s SE; %d bootstrap replicates",
            population, ratio, n, se_method, boot_reps),
    "common-scale item location bias, SE calibration and coverage",
    sum(ok_item), effect = if (any(ok_item))
      z[which(ok_item)[1L], "item_truth"] else NA_real_,
    bias = mean(z[ok_item, "item_est"] - z[ok_item, "item_truth"]),
    emp_sd = sd(z[ok_item, "item_est"]), mean_se = mean(z[ok_item, "item_se"]),
    coverage95 = mean(abs(z[ok_item, "item_est"] - z[ok_item, "item_truth"]) <=
                        1.96 * z[ok_item, "item_se"]),
    n_attempted = R, n_refused = sum(z[, "refused"]),
    n_nonconv = sum(z[, "nonconv"]))
  cat(sprintf("[%s] %s ratio %.1f %s: n=%d bias %+.4f SD/SE %.3f coverage %.3f\n",
              format(Sys.time(), "%H:%M"), population, ratio, se_method,
              sum(ok), out_alpha$bias, out_alpha$se_ratio,
              out_alpha$coverage95))
  ok_phi <- is.finite(z[, "phi_est"]) & is.finite(z[, "phi_se"]) &
    z[, "phi_se"] > 0
  out_phi <- if (any(ok_phi)) sv_row(STUDY,
    sprintf("%s persons; ratio %.1f; N=%d; %s SE; %d bootstrap replicates",
            population, ratio, n, se_method, boot_reps),
    "group-unit log phi bias, SE calibration and coverage",
    sum(ok_phi), effect = z[which(ok_phi)[1L], "phi_truth"],
    bias = mean(z[ok_phi, "phi_est"] - z[ok_phi, "phi_truth"]),
    emp_sd = sd(z[ok_phi, "phi_est"]), mean_se = mean(z[ok_phi, "phi_se"]),
    coverage95 = mean(abs(z[ok_phi, "phi_est"] - z[ok_phi, "phi_truth"]) <=
                        1.96 * z[ok_phi, "phi_se"]),
    n_attempted = R, n_refused = sum(z[, "refused"]),
    n_nonconv = sum(z[, "nonconv"])) else NULL
  do.call(rbind, Filter(Negate(is.null), list(out_alpha, out_item, out_phi)))
}

rows <- list(
  one_cell("normal", 1.0, 210000L, R = 200L),
  one_cell("normal", 1.4, 220000L),
  one_cell("bimodal", 1.0, 230000L, R = 200L),
  one_cell("bimodal", 1.4, 240000L),
  one_cell("contrasting groups", 1.0, 245000L, R = 150L, n = 600L),
  one_cell("contrasting groups", 1.4, 247000L, R = 150L, n = 600L),
  one_cell("normal", 1.0, 250000L, R = 80L, se_method = "bootstrap"),
  one_cell("normal", 1.4, 260000L, R = 80L, se_method = "bootstrap")
)
sv_write(do.call(rbind, rows), STUDY)
