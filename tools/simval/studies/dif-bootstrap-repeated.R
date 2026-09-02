# Operating characteristics of the public conditional DIF bootstrap for a
# genuinely repeated design. The conformance study checks orchestration but its
# two-level occasion factor has no non-trivial Greenhouse--Geisser correction.
# Here four occasions exercise that correction under balanced, nonspherical and
# differentially incomplete null designs, followed by uniform and non-uniform
# within-person DIF alternatives.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NULL_REPS <- as.integer(Sys.getenv("SV_NULL_REPS", "350"))
POWER_REPS <- as.integer(Sys.getenv("SV_POWER_REPS", "125"))
NBOOT <- as.integer(Sys.getenv("SV_BOOT", "99"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "4")))
ONLY <- Sys.getenv("SV_ONLY", "")

scenarios <- list(
  balanced_null = list(kind = "null", covariance = "balanced",
                       missing = FALSE, effect = 0),
  nonspherical_null = list(kind = "null", covariance = "nonspherical",
                           missing = FALSE, effect = 0),
  incomplete_mixed_null = list(kind = "null", covariance = "balanced",
                               missing = TRUE, effect = 0),
  uniform_occasion_0.70 = list(kind = "uniform", covariance = "balanced",
                               missing = FALSE, effect = .70),
  uniform_occasion_1.20 = list(kind = "uniform", covariance = "balanced",
                               missing = FALSE, effect = 1.20),
  nonuniform_occasion_0.80 = list(kind = "nonuniform",
                                  covariance = "balanced",
                                  missing = FALSE, effect = .80),
  nonuniform_occasion_1.50 = list(kind = "nonuniform",
                                  covariance = "balanced",
                                  missing = FALSE, effect = 1.50))
if (nzchar(ONLY)) {
  wanted <- trimws(strsplit(ONLY, ",", fixed = TRUE)[[1]])
  scenarios <- scenarios[intersect(wanted, names(scenarios))]
}
if (!length(scenarios)) stop("SV_ONLY selected no repeated-DIF scenario")

failure <- function(status, stage, condition = NULL) list(
  status = status, stage = stage,
  reason = if (is.null(condition)) "" else conditionMessage(condition),
  B_attempted = condition$B %||% 0L,
  B_used = condition$B_used %||% 0L,
  B_nonconverged = condition$B_nonconverged %||% 0L,
  B_errors = condition$B_errors %||% 0L)

make_data <- function(seed, scenario) {
  set.seed(seed)
  np <- if (identical(scenario$kind, "null")) 120L else 180L
  K <- 4L; L <- 7L
  person <- sprintf("P%03d", seq_len(np))
  id <- rep(person, each = K)
  occasion <- factor(rep(paste0("T", seq_len(K)), times = np),
                     levels = paste0("T", seq_len(K)))
  group0 <- if (isTRUE(scenario$missing))
    factor(c(rep("A", 36L), rep("B", np - 36L))) else
      factor(rep(c("A", "B"), length.out = np))
  group <- factor(rep(group0, each = K), levels = levels(group0))

  base <- rnorm(np)
  e <- matrix(rnorm(np * K), np, K)
  # Unequal variances and correlations produce a non-spherical trajectory.
  # The row-specific person locations are still absorbed by the Rasch score;
  # no item changes its response function under either null.
  if (identical(scenario$covariance, "nonspherical")) {
    A <- matrix(c(
      1.00, .72, .30, .05,
       .72, 1.00, .45, .12,
       .30, .45, 1.00, .68,
       .05, .12, .68, 1.00), K, K, byrow = TRUE)
    e <- e %*% chol(A)
    e <- sweep(e, 2L, c(.15, .35, .70, 1.10), "*")
  } else e <- e * .30
  time_shift <- c(-.35, -.05, .25, .55)
  theta <- rep(base, each = K) + as.vector(t(e)) +
    rep(time_shift, times = np)

  difficulty <- seq(-1.5, 1.5, length.out = L)
  eta <- outer(theta, difficulty, "-")
  target <- "I4"
  target_col <- 4L
  last <- occasion == "T4"
  if (identical(scenario$kind, "uniform"))
    eta[last, target_col] <- eta[last, target_col] + scenario$effect
  if (identical(scenario$kind, "nonuniform"))
    eta[last, target_col] <- (1 + scenario$effect) * theta[last]
  prob <- stats::plogis(eta)
  if (identical(scenario$covariance, "nonspherical")) {
    # A Gaussian copula preserves every marginal Rasch probability while
    # making repeated responses to the same item non-spherical. It is a
    # deliberate local-dependence stress condition, separate from the exact
    # fitted Rasch null in balanced_null.
    R <- matrix(c(
      1.00, .85, 0, 0,
       .85, 1.00, 0, 0,
       0, 0, 1.00, -.45,
       0, 0, -.45, 1.00), K, K, byrow = TRUE)
    root <- chol(R)
    X <- matrix(NA_integer_, nrow(eta), ncol(eta))
    for (n in seq_len(np)) {
      rows <- (n - 1L) * K + seq_len(K)
      for (j in seq_len(L)) {
        z <- drop(stats::rnorm(K) %*% root)
        X[rows, j] <- as.integer(z < stats::qnorm(prob[rows, j]))
      }
    }
  } else {
    X <- matrix(rbinom(length(eta), 1L, prob), nrow(eta), ncol(eta))
  }
  colnames(X) <- paste0("I", seq_len(L))

  # A common occasion shift plus group-dependent panel loss is a difficult
  # null for the between-person terms. No item parameter depends on group.
  if (isTRUE(scenario$missing)) {
    lose <- (group == "B" & occasion == "T4" & runif(length(id)) < .45) |
      (group == "A" & occasion == "T2" & runif(length(id)) < .15)
    keep <- !lose
    X <- X[keep, , drop = FALSE]
    id <- id[keep]; group <- droplevels(group[keep])
    occasion <- droplevels(occasion[keep])
  }
  # Preserve a modest item-level missingness burden in every design.
  miss <- matrix(runif(length(X)) < .02, nrow(X), ncol(X))
  X[miss] <- NA_integer_
  list(X = X, id = id, factors = data.frame(group, occasion),
       target = target)
}

one <- function(seed, scenario) {
  d <- make_data(seed, scenario)
  fit <- tryCatch(rasch(
    d$X, id = d$id, factors = d$factors, n_groups = 3), error = identity)
  if (inherits(fit, "error"))
    return(failure(if (inherits(fit, "rasch_refusal")) "refused" else
      "error", "observed fit", fit))
  if (!isTRUE(fit$est$converged))
    return(failure("nonconverged", "observed fit"))

  da <- tryCatch(dif_anova(
    fit, effects = "factorial", within = "occasion", n_groups = 3,
    p_adjust = "holm"), error = identity)
  if (inherits(da, "error"))
    return(failure(if (inherits(da, "rasch_refusal")) "refused" else
      "error", "observed DIF", da))
  target_row <- which(as.character(da$summary$item) == d$target &
                      as.character(da$summary$term) == "occasion")
  if (length(target_row) != 1L)
    return(failure("error", "target term"))

  db <- tryCatch(suppressWarnings(dif_bootstrap(
    fit, da, B = NBOOT, workers = 1L, seed = seed + 1000000L)),
    error = identity)
  if (inherits(db, "error"))
    return(failure(if (inherits(db, "rasch_refusal")) "refused" else
      "error", "DIF bootstrap", db))

  primary <- is.finite(da$terms$p_adj) & da$terms$p_adj < .05
  booted <- is.finite(db$terms$p_boot_adj) & db$terms$p_boot_adj < .05
  invariant <- as.character(da$terms$item) != d$target
  eps <- da$terms$gg_epsilon[is.finite(da$terms$gg_epsilon)]
  list(
    status = "ok", stage = "analysed", reason = "",
    primary_fwer = any(primary), bootstrap_fwer = any(booted),
    primary_invariant_fwer = any(primary[invariant]),
    bootstrap_invariant_fwer = any(booted[invariant]),
    primary_uniform = isTRUE(da$summary$uniform_DIF[target_row]),
    bootstrap_uniform = isTRUE(db$summary$uniform_DIF_boot[target_row]),
    primary_nonuniform = isTRUE(da$summary$nonuniform_DIF[target_row]),
    bootstrap_nonuniform = isTRUE(db$summary$nonuniform_DIF_boot[target_row]),
    gg_epsilon = if (length(eps)) mean(eps) else NA_real_,
    B_attempted = db$B, B_used = db$B_used,
    B_nonconverged = db$B_nonconverged, B_errors = db$B_errors)
}

set.seed(9.04e7)
rows <- list(); row_id <- 0L
for (nm in names(scenarios)) {
  scenario <- scenarios[[nm]]
  n_attempted <- if (identical(scenario$kind, "null"))
    NULL_REPS else POWER_REPS
  message(nm, " (", n_attempted, " datasets; B = ", NBOOT, ")")
  seeds <- sample.int(1e8, n_attempted)
  z <- parallel::mclapply(
    seeds, one, scenario = scenario, mc.cores = NCORE,
    mc.preschedule = TRUE)
  status <- vapply(z, `[[`, "", "status")
  ok <- status == "ok"; zz <- z[ok]; n <- length(zz)
  if (!n) stop("no analysed replicate in ", nm)
  val <- function(name) vapply(zz, `[[`, FALSE, name)
  num <- function(name) vapply(zz, `[[`, 0, name)
  pf <- val("primary_fwer"); bf <- val("bootstrap_fwer")
  pif <- val("primary_invariant_fwer")
  bif <- val("bootstrap_invariant_fwer")
  pu <- val("primary_uniform"); bu <- val("bootstrap_uniform")
  pn <- val("primary_nonuniform"); bn <- val("bootstrap_nonuniform")
  eps <- num("gg_epsilon")
  b_attempted <- vapply(z, `[[`, 0, "B_attempted")
  b_used <- vapply(z, `[[`, 0, "B_used")
  b_nonconv <- vapply(z, `[[`, 0, "B_nonconverged")
  b_errors <- vapply(z, `[[`, 0, "B_errors")
  common <- list(
    study = "dif-bootstrap-repeated", scenario = nm, n_reps = n,
    n_attempted = n_attempted, n_refused = sum(status == "refused"),
    n_nonconv = sum(status == "nonconverged"),
    n_error = sum(status == "error"),
    n_boot_attempted = sum(b_attempted), n_boot_used = sum(b_used),
    n_boot_nonconv = sum(b_nonconv), n_boot_errors = sum(b_errors),
    effect = scenario$effect,
    notes = paste0(
      "four-occasion mixed residual ANOVA; Holm primary and complete-family ",
      "minimum-p bootstrap; performance conditional on an analysable observed ",
      "fit and the public usable-refit rule"))
  mk <- function(quantity, ...) do.call(
    sv_row, c(common, list(quantity = quantity), list(...)))
  block <- list(
    mk("mean Greenhouse-Geisser epsilon", bias = mean(eps, na.rm = TRUE),
       emp_sd = stats::sd(eps, na.rm = TRUE)),
    mk("mean usable bootstrap proportion", bias = mean(b_used[ok] / NBOOT),
       emp_sd = stats::sd(b_used[ok] / NBOOT)))
  if (identical(scenario$kind, "null")) {
    block <- c(list(
      mk("primary adjusted FWER", familywise = mean(pf)),
      mk("bootstrap minimum-p FWER", familywise = mean(bf)),
      mk("paired bootstrap-minus-primary FWER", bias = mean(bf - pf),
         mean_se = stats::sd(bf - pf) / sqrt(n), se_ratio = NA_real_)),
      block)
  } else {
    target_primary <- if (identical(scenario$kind, "uniform")) pu else pn
    target_boot <- if (identical(scenario$kind, "uniform")) bu else bn
    block <- c(block, list(
      mk("primary affected-item power", power = mean(target_primary)),
      mk("bootstrap affected-item power", power = mean(target_boot)),
      mk("primary invariant-item FWER", familywise = mean(pif)),
      mk("bootstrap invariant-item FWER", familywise = mean(bif)),
      mk("paired bootstrap-minus-primary power",
         bias = mean(target_boot - target_primary),
         mean_se = stats::sd(target_boot - target_primary) / sqrt(n),
         se_ratio = NA_real_)))
  }
  row_id <- row_id + 1L
  rows[[row_id]] <- do.call(rbind, block)
  sv_write(do.call(rbind, rows), "dif-bootstrap-repeated")
}
