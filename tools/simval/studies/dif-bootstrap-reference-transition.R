# Narrow validation of the reference-probability DIF bootstrap.
#
# Earlier releases compared raw F statistics for the marginal bootstrap
# probability and then used that value as a floor on the minimum-p result.
# That comparison is equivalent only while a term's F reference is fixed.
# This study retains both F and reference-p matrices from the same public
# dif_bootstrap() run, so the old and current calculations are compared on
# identical observed datasets and identical conditional draws.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

dif_source_md5 <- unname(tools::md5sum("R/dif-bootstrap.R"))

NREP <- as.integer(Sys.getenv("SV_REPS", "100"))
NBOOT <- as.integer(Sys.getenv("SV_BOOT", "99"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "4")))
ONLY <- Sys.getenv("SV_ONLY", "")

scenarios <- c("balanced", "nonspherical", "incomplete")
if (nzchar(ONLY)) {
  wanted <- trimws(strsplit(ONLY, ",", fixed = TRUE)[[1L]])
  scenarios <- intersect(wanted, scenarios)
}
if (!length(scenarios)) stop("SV_ONLY selected no transition scenario")

make_data <- function(seed, scenario) {
  set.seed(seed)
  np <- 120L; occasions <- 4L; items <- 7L
  id <- rep(sprintf("P%03d", seq_len(np)), each = occasions)
  occasion <- factor(rep(paste0("T", seq_len(occasions)), times = np),
                     levels = paste0("T", seq_len(occasions)))
  group0 <- if (identical(scenario, "incomplete"))
    factor(c(rep("A", 36L), rep("B", np - 36L))) else
      factor(rep(c("A", "B"), length.out = np))
  group <- factor(rep(group0, each = occasions), levels = levels(group0))

  base <- rnorm(np)
  e <- matrix(rnorm(np * occasions), np, occasions)
  if (identical(scenario, "nonspherical")) {
    cor_mat <- matrix(c(
      1.00, .72, .30, .05,
       .72, 1.00, .45, .12,
       .30, .45, 1.00, .68,
       .05, .12, .68, 1.00), occasions, occasions, byrow = TRUE)
    e <- e %*% chol(cor_mat)
    e <- sweep(e, 2L, c(.15, .35, .70, 1.10), "*")
  } else {
    e <- e * .30
  }
  theta <- rep(base, each = occasions) + as.vector(t(e)) +
    rep(c(-.35, -.05, .25, .55), times = np)
  difficulty <- seq(-1.5, 1.5, length.out = items)
  X <- matrix(rbinom(length(theta) * items, 1L,
    stats::plogis(outer(theta, difficulty, "-"))), length(theta), items)
  colnames(X) <- paste0("I", seq_len(items))

  if (identical(scenario, "incomplete")) {
    lose <- (group == "B" & occasion == "T4" & runif(length(id)) < .45) |
      (group == "A" & occasion == "T2" & runif(length(id)) < .15)
    keep <- !lose
    X <- X[keep, , drop = FALSE]
    id <- id[keep]
    group <- droplevels(group[keep])
    occasion <- droplevels(occasion[keep])
  }
  miss <- matrix(runif(length(X)) < .02, nrow(X), ncol(X))
  X[miss] <- NA_integer_
  list(X = X, id = id, factors = data.frame(group, occasion))
}

failure <- function(status) list(status = status, B = 0L, used = 0L,
  nonconv = 0L, errors = 0L)

one <- function(seed, scenario) {
  d <- make_data(seed, scenario)
  fit <- tryCatch(rasch(d$X, id = d$id, factors = d$factors, n_groups = 3),
                  error = identity)
  if (inherits(fit, "error"))
    return(failure(if (inherits(fit, "rasch_refusal")) "refused" else
      "error"))
  if (!isTRUE(fit$est$converged)) return(failure("nonconverged"))
  dif <- tryCatch(dif_anova(fit, effects = "factorial", within = "occasion",
                            n_groups = 3, p_adjust = "holm"),
                  error = identity)
  if (inherits(dif, "error"))
    return(failure(if (inherits(dif, "rasch_refusal")) "refused" else
      "error"))
  boot <- tryCatch(suppressWarnings(dif_bootstrap(
    fit, dif, B = NBOOT, workers = 1L, seed = seed + 1000000L)),
    error = identity)
  if (inherits(boot, "error"))
    return(failure(if (inherits(boot, "rasch_refusal")) "refused" else
      "error"))

  family <- !dif$term_ids %in% c("Residuals", "ci")
  observed <- dif$terms[family, , drop = FALSE]
  reference_marginal <- boot$terms$p_boot[family]
  reference_adjusted <- boot$terms$p_boot_adj[family]
  raw_f_marginal <- vapply(seq_len(ncol(boot$replicates$F)), function(j)
    (1 + sum(boot$replicates$F[, j] >= observed$F_value[j])) /
      (boot$B_used + 1), numeric(1))
  legacy_adjusted <- pmax(raw_f_marginal, reference_adjusted)
  list(
    status = "ok", B = boot$B, used = boot$B_used,
    nonconv = boot$B_nonconverged, errors = boot$B_errors,
    reference_fwer = any(reference_adjusted < .05),
    legacy_fwer = any(legacy_adjusted < .05),
    reference_marginal_rate = mean(reference_marginal < .05),
    legacy_marginal_rate = mean(raw_f_marginal < .05),
    marginal_max_change = max(abs(reference_marginal - raw_f_marginal)),
    adjusted_max_change = max(abs(reference_adjusted - legacy_adjusted)),
    marginal_decision_changes = sum(
      (reference_marginal < .05) != (raw_f_marginal < .05)),
    adjusted_decision_changes = sum(
      (reference_adjusted < .05) != (legacy_adjusted < .05)))
}

set.seed(9.06e7)
rows <- vector("list", length(scenarios))
for (j in seq_along(scenarios)) {
  scenario <- scenarios[j]
  message(scenario, " (", NREP, " datasets; B = ", NBOOT, ")")
  seeds <- sample.int(1e8, NREP)
  z <- parallel::mclapply(seeds, one, scenario = scenario, mc.cores = NCORE,
                          mc.preschedule = TRUE)
  status <- vapply(z, `[[`, "", "status")
  ok <- status == "ok"
  zz <- z[ok]
  n <- length(zz)
  if (!n) stop("no analysed replicate in ", scenario)
  val <- function(name) vapply(zz, `[[`, 0, name)
  reference_fwer <- as.logical(val("reference_fwer"))
  legacy_fwer <- as.logical(val("legacy_fwer"))
  reference_marginal_rate <- val("reference_marginal_rate")
  legacy_marginal_rate <- val("legacy_marginal_rate")
  marginal_max <- val("marginal_max_change")
  adjusted_max <- val("adjusted_max_change")
  marginal_changed <- val("marginal_decision_changes")
  adjusted_changed <- val("adjusted_decision_changes")
  used <- val("used")
  common <- list(
    study = "dif-bootstrap-reference-transition", scenario = scenario,
    n_reps = n, n_attempted = NREP,
    n_refused = sum(status == "refused"),
    n_nonconv = sum(status == "nonconverged"),
    n_error = sum(status == "error"),
    n_boot_attempted = sum(vapply(z, `[[`, 0, "B")),
    n_boot_used = sum(used),
    n_boot_nonconv = sum(vapply(z, `[[`, 0, "nonconv")),
    n_boot_errors = sum(vapply(z, `[[`, 0, "errors")),
    notes = paste(
      "same observed datasets and conditional draws under both calculations;",
      "current algorithm uses reference probabilities for marginal and",
      "single-step minimum-p inference; R/dif-bootstrap.R md5",
      dif_source_md5))
  mk <- function(quantity, ...) do.call(
    sv_row, c(common, list(quantity = quantity), list(...)))
  rows[[j]] <- rbind(
    mk("current reference-p FWER", familywise = mean(reference_fwer)),
    mk("legacy raw-F-floor FWER", familywise = mean(legacy_fwer)),
    mk("paired current-minus-legacy FWER",
       bias = mean(reference_fwer - legacy_fwer),
       mean_se = stats::sd(reference_fwer - legacy_fwer) / sqrt(n),
       se_ratio = NA_real_),
    mk("current reference-p marginal term Type I",
       type1 = mean(reference_marginal_rate),
       mc_override = list(
         type1 = stats::sd(reference_marginal_rate) / sqrt(n))),
    mk("legacy raw-F marginal term Type I",
       type1 = mean(legacy_marginal_rate),
       mc_override = list(
         type1 = stats::sd(legacy_marginal_rate) / sqrt(n))),
    mk("datasets with changed marginal decision",
       type1 = mean(marginal_changed > 0L)),
    mk("datasets with changed adjusted decision",
       familywise = mean(adjusted_changed > 0L)),
    mk("maximum marginal probability change",
       bias = max(marginal_max), emp_sd = stats::sd(marginal_max)),
    mk("maximum adjusted probability change",
       bias = max(adjusted_max), emp_sd = stats::sd(adjusted_max)),
    mk("mean usable bootstrap proportion",
       bias = mean(used / NBOOT), emp_sd = stats::sd(used / NBOOT)))
  sv_write(do.call(rbind, rows), "dif-bootstrap-reference-transition")
}
